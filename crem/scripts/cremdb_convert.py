#!/usr/bin/env python3
"""
Script to convert existing CReM databases to the new deduplicated format.

Old format (per radius table):
    CREATE TABLE radiusN(env TEXT, core_smi TEXT, core_sma TEXT, dist2 INTEGER, ...)

New format:
    CREATE TABLE radiusN(env_id INTEGER NOT NULL, core_smi_id INTEGER NOT NULL[, set_name ...])
    CREATE TABLE envs(env_id INTEGER PRIMARY KEY AUTOINCREMENT, env TEXT NOT NULL UNIQUE)
    CREATE TABLE frags(core_smi_id INTEGER PRIMARY KEY AUTOINCREMENT,
                      core_smi TEXT NOT NULL UNIQUE,
                      dist2 INTEGER,
                      core_smi_h_id INTEGER NOT NULL)
    CREATE TABLE frags_h(core_smi_h_id INTEGER PRIMARY KEY AUTOINCREMENT,
                        smi TEXT NOT NULL UNIQUE)

Usage:
    python convert_crem_db.py old_database.db new_database.db [--radii 1 2 3] [--set-name NAME]
"""

import sqlite3
import argparse
import traceback
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import sys
import re
from tqdm import tqdm
from rdkit import Chem, RDLogger
from crem.mol_context import combine_core_env_to_rxn_smarts
from crem.scripts.cremdb_create import create_indices, _replace_attachment_points_with_h


def replace_attachment_points_with_h(smiles: str) -> str:
    """
    Replace attachment point markers (*) with hydrogen.

    Args:
        smiles: SMILES string with attachment points marked as *

    Returns:
        SMILES string with hydrogens instead of attachment points
    """
    return _replace_attachment_points_with_h(smiles)


def _validate_set_name(set_name: str) -> str:
    if not re.match(r'^[A-Za-z_][A-Za-z0-9_]*$', set_name):
        raise ValueError(
            "set name must be a valid SQLite identifier (letters, numbers, underscores; cannot start with a number)"
        )
    if set_name in ('env_id', 'core_smi_id'):
        raise ValueError("set name cannot be env_id or core_smi_id")
    return set_name


def _get_freq_column_type(old_conn: sqlite3.Connection, radius: int) -> str:
    table_name = f"radius{radius}"
    table_check = old_conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (table_name,)
    ).fetchone()
    if not table_check:
        return "INTEGER"
    schema = old_conn.execute(f"PRAGMA table_info({table_name})").fetchall()
    for col in schema:
        if col[1] == 'freq':
            return col[2] or "INTEGER"
    raise ValueError(f"Missing required column 'freq' in {table_name}")


def create_new_schema(
    new_conn: sqlite3.Connection,
    old_conn: sqlite3.Connection,
    radii: List[int],
    set_name: str = "undefined",
):
    """
    Create the new database schema.

    Args:
        new_conn: SQLite connection to new database
        old_conn: SQLite connection to old database
        radii: List of radius values to create tables for
        set_name: Optional column name to add to radius tables
    """
    cur = new_conn.cursor()

    cur.execute("PRAGMA user_version = 1;")

    # Create envs table
    cur.execute("""
        CREATE TABLE IF NOT EXISTS envs(
            env_id INTEGER PRIMARY KEY AUTOINCREMENT,
            env TEXT NOT NULL UNIQUE
        )
    """)

    # Create frags_h table (must be created before frags due to foreign key)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS frags_h(
            core_smi_h_id INTEGER PRIMARY KEY AUTOINCREMENT,
            smi TEXT NOT NULL UNIQUE
        )
    """)

    # Create frags table
    cur.execute("""
        CREATE TABLE IF NOT EXISTS frags(
            core_smi_id INTEGER PRIMARY KEY AUTOINCREMENT,
            core_smi TEXT NOT NULL UNIQUE,
            core_smi_h_id INTEGER NOT NULL,
            FOREIGN KEY (core_smi_h_id) REFERENCES frags_h(core_smi_h_id)
        )
    """)

    # Create radius tables
    # core_num_atoms and dist2 live on radius{N} so the hot-path filter
    # queries do not have to join frags.
    for radius in radii:
        column_defs = [
            "env_id INTEGER NOT NULL",
            "core_smi_id INTEGER NOT NULL",
            "core_num_atoms INTEGER NOT NULL",
            "dist2 INTEGER NOT NULL",
            "is_ring_closure INTEGER NOT NULL DEFAULT 0",
        ]
        if set_name:
            freq_type = _get_freq_column_type(old_conn, radius)
            column_defs.append(f"{set_name} {freq_type} DEFAULT 0")
        cur.execute(f"""
            CREATE TABLE IF NOT EXISTS radius{radius}(
                {", ".join(column_defs)},
                FOREIGN KEY (env_id) REFERENCES envs(env_id),
                FOREIGN KEY (core_smi_id) REFERENCES frags(core_smi_id),
                UNIQUE (env_id, core_smi_id, is_ring_closure)
            )
        """)

    new_conn.commit()

    # Add columns, which are not pre-defined, to transfer them
    predefined_colnames = ['env', 'core_smi', 'core_num_atoms', 'core_sma', 'dist2', 'freq']
    old_tables = old_conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'radius%'").fetchall()
    for row in old_tables:
        col_data_old = old_conn.execute(f"PRAGMA table_info({row[0]})").fetchall()
        col_data_new = new_conn.execute(f"PRAGMA table_info(frags)").fetchall()
        col_names_new = [row[1] for row in col_data_new]
        for row in col_data_old:
            if row[1] not in predefined_colnames and row[1] not in col_names_new:
                sql = f'ALTER TABLE frags ADD {row[1]} {row[2]} {"NOT NULL" if row[3] == 1 else ""} {f"DEFAULT {row[4]}" if row[4] is not None else ""}'
                print(sql)
                new_conn.execute(sql)


def convert_database(
    old_db_path: str,
    new_db_path: str,
    radii: List[int],
    batch_size: int = 10000,
    set_name: str = "undefined",
    verbose: bool = True,
):
    """
    Convert old database format to new deduplicated format.

    Args:
        old_db_path: Path to existing database
        new_db_path: Path to new database (will be created)
        radii: List of radius values to convert
        batch_size: Number of rows to process in each batch
        set_name: column name to add to radius tables and fill from freq column from old db
        verbose: Print progress information
    """

    RDLogger.DisableLog('rdApp.warning')
    set_name = _validate_set_name(set_name)

    if verbose:
        print(f"Converting database from {old_db_path} to {new_db_path}")
        print(f"Processing radii: {radii}")

    # Connect to databases
    old_conn = sqlite3.connect(old_db_path)
    new_conn = sqlite3.connect(new_db_path)

    old_cur = old_conn.cursor()
    new_cur = new_conn.cursor()

    try:
        # Create new schema
        create_new_schema(new_conn, old_conn, radii, set_name)

        # Dictionaries to map unique values to IDs
        env_to_id: Dict[str, int] = {}
        core_smi_to_id: Dict[str, int] = {}
        core_smi_h_to_id: Dict[str, int] = {}

        env_counter = 0
        frag_counter = 0
        frag_h_counter = 0

        # Process each radius table
        for radius in radii:
            if verbose:
                print(f"\nProcessing radius{radius}...")

            # Check if table exists in old database
            table_check = old_cur.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                (f'radius{radius}',)
            ).fetchone()

            if not table_check:
                print(f"Warning: radius{radius} table not found in old database, skipping")
                continue

            # Get total row count for progress bar
            total_rows = old_cur.execute(f"SELECT COUNT(rowid) FROM radius{radius}").fetchone()[0]

            if verbose:
                print(f"Total rows: {total_rows}")

            # Get schema of old table
            schema = old_cur.execute(f"PRAGMA table_info(radius{radius})").fetchall()
            column_names = [col[1] for col in schema]

            # Verify required columns exist
            required_cols = ['env', 'core_smi', 'core_num_atoms', 'core_sma', 'dist2']
            missing_cols = [col for col in required_cols if col not in column_names]
            if missing_cols:
                raise ValueError(f"Missing required columns in radius{radius}: {missing_cols}")
            if set_name and 'freq' not in column_names:
                raise ValueError(f"Missing required column 'freq' in radius{radius} for set-name")

            # Process in batches
            offset = 0
            pbar = tqdm(total=total_rows, disable=not verbose, desc=f"radius{radius}")

            # Columns we SELECT from v0 radius{N}: everything pivoted out of the
            # row plus dist2 / core_num_atoms (which now live only on the v1
            # radius table).
            v0_extra_columns = [col for col in column_names if col not in ('env', 'core_smi', 'core_sma', 'freq')]
            # Columns we INSERT into v1 frags: same minus dist2 and
            # core_num_atoms (both denormalized into radius{N}).
            frags_columns = [col for col in v0_extra_columns if col not in ('dist2', 'core_num_atoms')]
            column_names_transfer = ['env', 'core_smi'] + v0_extra_columns
            if set_name:
                column_names_transfer.append('freq')

            # Position of v0_extra_columns inside `row` (offset by 2 for env, core_smi).
            v0_col_pos = {col: 2 + i for i, col in enumerate(v0_extra_columns)}
            core_num_atoms_pos = v0_col_pos['core_num_atoms']
            dist2_pos = v0_col_pos['dist2']

            while offset < total_rows:
                # Fetch batch
                rows = old_cur.execute(
                    f'SELECT {",".join(column_names_transfer)} FROM radius{radius} LIMIT ? OFFSET ?',
                    (batch_size, offset)
                ).fetchall()

                if not rows:
                    break

                # Process each row
                new_rows_radius = []
                new_rows_frags = []
                new_rows_frags_h = []
                new_rows_envs = []

                for row in rows:
                    env = row[0]
                    core_smi = row[1]
                    core_num_atoms = row[core_num_atoms_pos]
                    dist2 = row[dist2_pos]

                    # Get or create env_id
                    if env not in env_to_id:
                        env_counter += 1
                        env_to_id[env] = env_counter
                        new_rows_envs.append((env_counter, env))
                    env_id = env_to_id[env]

                    # Get or create core_smi_h and core_smi_h_id (natural key: core_smi_h).
                    core_smi_h = replace_attachment_points_with_h(core_smi)
                    if core_smi_h not in core_smi_h_to_id:
                        frag_h_counter += 1
                        core_smi_h_to_id[core_smi_h] = frag_h_counter
                        new_rows_frags_h.append((frag_h_counter, core_smi_h))
                    core_smi_h_id = core_smi_h_to_id[core_smi_h]

                    # Get or create core_smi_id
                    if core_smi not in core_smi_to_id:
                        frag_counter += 1
                        core_smi_to_id[core_smi] = frag_counter
                        new_rows_frags.append(
                            tuple(row[v0_col_pos[c]] for c in frags_columns)
                            + (core_smi_h_id, frag_counter, core_smi)
                        )
                    core_smi_id = core_smi_to_id[core_smi]

                    # Add to radius table (denormalized core_num_atoms / dist2).
                    if set_name:
                        freq_value = row[2 + len(v0_extra_columns)]
                        new_rows_radius.append((env_id, core_smi_id, core_num_atoms, dist2, freq_value))
                    else:
                        new_rows_radius.append((env_id, core_smi_id, core_num_atoms, dist2))

                # Batch insert into radius table
                new_cur.executemany("INSERT INTO envs (env_id, env) VALUES (?, ?)",
                                    new_rows_envs)

                new_cur.executemany("INSERT INTO frags_h (core_smi_h_id, smi) VALUES (?, ?)",
                                    new_rows_frags_h)

                all_frag_cols = frags_columns + ['core_smi_h_id', 'core_smi_id', 'core_smi']
                sql = (f"INSERT INTO frags ({','.join(all_frag_cols)}) "
                       f"VALUES ({','.join('?' * len(all_frag_cols))})")
                new_cur.executemany(sql, new_rows_frags)

                if set_name:
                    new_cur.executemany(
                        f"INSERT INTO radius{radius} "
                        f"(env_id, core_smi_id, core_num_atoms, dist2, {set_name}) "
                        f"VALUES (?, ?, ?, ?, ?)",
                        new_rows_radius
                    )
                else:
                    new_cur.executemany(
                        f"INSERT INTO radius{radius} "
                        f"(env_id, core_smi_id, core_num_atoms, dist2) "
                        f"VALUES (?, ?, ?, ?)",
                        new_rows_radius
                    )

                # Commit batch
                new_conn.commit()

                offset += batch_size
                pbar.update(len(rows))

            pbar.close()

            if verbose:
                # Get statistics
                new_count = new_cur.execute(f"SELECT COUNT(*) FROM radius{radius}").fetchone()[0]
                print(f"Converted {new_count} rows")

        # Create indices for optimized queries
        if verbose:
            print("\nCreating indices...")

        create_indices(new_conn, radii, verbose)

        # Print final statistics
        if verbose:
            print("\n" + "="*60)
            print("Conversion Summary:")
            print(f"  Unique environments: {len(env_to_id)}")
            print(f"  Unique fragments (with attachment points): {len(core_smi_to_id)}")
            print(f"  Unique fragments (H-replaced): {len(core_smi_h_to_id)}")

            # Calculate size reduction
            import os
            if os.path.exists(old_db_path) and os.path.exists(new_db_path):
                old_size = os.path.getsize(old_db_path) / (1024**2)  # MB
                new_size = os.path.getsize(new_db_path) / (1024**2)  # MB
                reduction = ((old_size - new_size) / old_size) * 100
                print(f"\nDatabase sizes:")
                print(f"  Old: {old_size:.2f} MB")
                print(f"  New: {new_size:.2f} MB")
                print(f"  Reduction: {reduction:.1f}%")
            print("="*60)

    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        raise

    finally:
        old_conn.close()
        new_conn.close()
        RDLogger.EnableLog('rdApp.warning')


def verify_conversion(old_db_path: str, new_db_path: str, radius: int = 3,
                     sample_size: int = 100, verbose: bool = True):
    """
    Verify that conversion was successful by comparing sample data.

    Args:
        old_db_path: Path to old database
        new_db_path: Path to new database
        radius: Radius to check
        sample_size: Number of rows to sample for verification
        verbose: Print details
    """
    if verbose:
        print(f"\nVerifying conversion for radius{radius}...")

    old_conn = sqlite3.connect(old_db_path)
    new_conn = sqlite3.connect(new_db_path)

    old_cur = old_conn.cursor()
    new_cur = new_conn.cursor()

    try:
        table_check = old_cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (f'radius{radius}',)
        ).fetchone()

        if not table_check:
            print(f"Warning: radius{radius} table not found in old database, skipping verification")
            return False

        # Get total counts
        old_count = old_cur.execute(f"SELECT COUNT(*) FROM radius{radius}").fetchone()[0]
        new_count = new_cur.execute(f"SELECT COUNT(*) FROM radius{radius}").fetchone()[0]

        if verbose:
            print(f"Row counts - Old: {old_count}, New: {new_count}")

        if old_count != new_count:
            print(f"WARNING: Row count mismatch! Old: {old_count}, New: {new_count}")
            return False

        # Sample and compare data
        sample_rows = old_cur.execute(
            f"SELECT env, core_smi, core_sma, dist2 FROM radius{radius} "
            f"ORDER BY RANDOM() LIMIT ?",
            (sample_size,)
        ).fetchall()

        mismatches = 0
        for env, core_smi, core_sma, dist2 in sample_rows:
            # Query new database
            result = new_cur.execute(f"""
                SELECT f.core_smi, r.dist2
                FROM radius{radius} r
                JOIN frags f ON r.core_smi_id = f.core_smi_id
                JOIN envs e ON r.env_id = e.env_id
                WHERE e.env = ? AND f.core_smi = ?
            """, (env, core_smi)).fetchone()

            if result is None:
                mismatches += 1
                if verbose:
                    print(f"Missing entry: env={env[:50]}..., core_smi={core_smi}")
            else:
                new_core_sma = combine_core_env_to_rxn_smarts(core_smi, env, False)
                if result != (core_smi, dist2) or new_core_sma != core_sma:
                    mismatches += 1
                    if verbose:
                        print(f"Data mismatch: expected {(core_smi, core_sma, dist2)}, got {result}")

        if mismatches == 0:
            if verbose:
                print(f"✓ Verification passed! Checked {sample_size} random samples.")
            return True
        else:
            print(f"✗ Verification failed! {mismatches}/{sample_size} samples had issues.")
            return False

    finally:
        old_conn.close()
        new_conn.close()


def main():
    parser = argparse.ArgumentParser(
        description='Convert CReM database to deduplicated format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  cremdb_convert -i old.db -o new.db
  cremdb_convert -i old.db -o new.db --radii 1 2 3 4 5
  cremdb_convert -i old.db -o new.db --batch-size 5000 --verify
  cremdb_convert -i old.db -o new.db --set-name my_set
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Path to existing database')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to new database (will be created)')
    parser.add_argument(
        '--radii',
        nargs='+',
        type=int,
        default=[1, 2, 3, 4, 5],
        help='Space-separated list of radii to convert (default: 1 2 3 4 5)',
    )
    parser.add_argument('--batch-size', type=int, default=10000,
                       help='Number of rows to process per batch (default: 10000)')
    parser.add_argument('--set-name', default="undefined",
                        help='Name of the new column to create in radius tables and populate from the freq column of old_db (default: undefined')
    parser.add_argument('--verify', action='store_true',
                       help='Verify conversion after completion')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Parse radii
    radii = sorted(set(args.radii))

    # Check if output file exists
    import os
    if os.path.exists(args.output):
        response = input(f"{args.output} already exists. Overwrite? (y/N): ")
        if response.lower() != 'y':
            print("Conversion cancelled.")
            return
        os.remove(args.output)

    # Perform conversion
    verbose = not args.quiet

    try:
        convert_database(old_db_path=args.input,
                         new_db_path=args.output,
                         radii=radii,
                         batch_size=args.batch_size,
                         set_name=args.set_name,
                         verbose=verbose)

        # Verify if requested
        if args.verify:
            for radius in radii:
                verify_conversion(args.input, args.output, radius, verbose=verbose)

        print("\n✓ Conversion completed successfully!")

    except Exception as e:
        print(f"\n✗ Conversion failed: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
