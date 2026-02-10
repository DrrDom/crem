#!/usr/bin/env python3
"""
Script to convert existing CReM databases to the new deduplicated format.

Old format (per radius table):
    CREATE TABLE radiusN(env TEXT, core_smi TEXT, core_sma TEXT, dist2 INTEGER, ...)

New format:
    CREATE TABLE radiusN(env_id INTEGER NOT NULL, core_smi_id INTEGER NOT NULL)
    CREATE TABLE envs(env_id INTEGER PRIMARY KEY AUTOINCREMENT, env TEXT NOT NULL UNIQUE)
    CREATE TABLE frags(core_smi_id INTEGER PRIMARY KEY AUTOINCREMENT,
                      core_smi TEXT NOT NULL UNIQUE,
                      core_sma TEXT NOT NULL UNIQUE,
                      dist2 INTEGER,
                      core_smi_h_id INTEGER NOT NULL)
    CREATE TABLE frags_h(core_smi_h_id INTEGER PRIMARY KEY AUTOINCREMENT,
                        smi TEXT NOT NULL UNIQUE)

Usage:
    python convert_crem_db.py old_database.db new_database.db [--radii 1,2,3]
"""

import sqlite3
import argparse
import traceback
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import sys
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi


def replace_attachment_points_with_h(smiles: str) -> str:
    """
    Replace attachment point markers (*) with hydrogen.

    Args:
        smiles: SMILES string with attachment points marked as *

    Returns:
        SMILES string with hydrogens instead of attachment points
    """
    # Simple replacement - in production you might want to use RDKit for proper canonicalization
    return Chem.CanonSmiles(smiles.replace('*', 'H'))


def create_new_schema(new_conn: sqlite3.Connection, old_conn: sqlite3.Connection, radii: List[int]):
    """
    Create the new database schema.

    Args:
        new_conn: SQLite connection to new database
        radii: List of radius values to create tables for
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
            smi TEXT NOT NULL,
            inchi TEXT NOT NULL UNIQUE
        )
    """)

    # Create frags table
    cur.execute("""
        CREATE TABLE IF NOT EXISTS frags(
            core_smi_id INTEGER PRIMARY KEY AUTOINCREMENT,
            core_smi TEXT NOT NULL UNIQUE,
            core_num_atoms INTEGER NOT NULL,
            core_sma TEXT NOT NULL,
            dist2 INTEGER NOT NULL,
            core_smi_h_id INTEGER NOT NULL,
            FOREIGN KEY (core_smi_h_id) REFERENCES frags_h(core_smi_h_id)
        )
    """)

    # Create radius tables
    for radius in radii:
        cur.execute(f"""
            CREATE TABLE IF NOT EXISTS radius{radius}(
                env_id INTEGER NOT NULL,
                core_smi_id INTEGER NOT NULL,
                FOREIGN KEY (env_id) REFERENCES envs(env_id),
                FOREIGN KEY (core_smi_id) REFERENCES frags(core_smi_id),
                UNIQUE (env_id, core_smi_id)
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


def convert_database(old_db_path: str, new_db_path: str, radii: List[int],
                     batch_size: int = 10000, verbose: bool = True):
    """
    Convert old database format to new deduplicated format.

    Args:
        old_db_path: Path to existing database
        new_db_path: Path to new database (will be created)
        radii: List of radius values to convert
        batch_size: Number of rows to process in each batch
        verbose: Print progress information
    """

    RDLogger.DisableLog('rdApp.warning')

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
        create_new_schema(new_conn, old_conn, radii)

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

            # Process in batches
            offset = 0
            pbar = tqdm(total=total_rows, disable=not verbose, desc=f"radius{radius}")

            # make a specific order: env, core_smi, other columns ...
            column_names_tranfer = column_names
            column_names_tranfer.remove('freq')
            column_names_tranfer.remove('env')
            column_names_tranfer.insert(0, 'env')
            column_names_tranfer.remove('core_smi')
            column_names_tranfer.insert(1, 'core_smi')

            while offset < total_rows:
                # Fetch batch
                rows = old_cur.execute(
                    f'SELECT {",".join(column_names_tranfer)} FROM radius{radius} LIMIT ? OFFSET ?',
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

                    # Get or create env_id
                    if env not in env_to_id:
                        env_counter += 1
                        env_to_id[env] = env_counter
                        new_rows_envs.append((env_counter, env))
                    env_id = env_to_id[env]

                    # Get or create core_smi_h and core_smi_h_id
                    core_smi_h = replace_attachment_points_with_h(core_smi)
                    inchi_value = inchi.MolToInchi(Chem.MolFromSmiles(core_smi_h))
                    if inchi_value not in core_smi_h_to_id:
                        frag_h_counter += 1
                        core_smi_h_to_id[inchi_value] = frag_h_counter
                        new_rows_frags_h.append((frag_h_counter, core_smi_h, inchi_value))
                    core_smi_h_id = core_smi_h_to_id[inchi_value]

                    # Get or create core_smi_id
                    if core_smi not in core_smi_to_id:
                        frag_counter += 1
                        core_smi_to_id[core_smi] = frag_counter
                        new_rows_frags.append(row[2:] + (core_smi_h_id, frag_counter, core_smi))
                    core_smi_id = core_smi_to_id[core_smi]

                    # Add to radius table
                    new_rows_radius.append((env_id, core_smi_id))

                # Batch insert into radius table
                new_cur.executemany("INSERT INTO envs (env_id, env) VALUES (?, ?)",
                                    new_rows_envs)

                new_cur.executemany("INSERT INTO frags_h (core_smi_h_id, smi, inchi) VALUES (?, ?, ?)",
                                    new_rows_frags_h)

                sql = (f"INSERT INTO frags ({','.join(column_names_tranfer[2:])}, core_smi_h_id, core_smi_id, core_smi) "
                       f"VALUES ({','.join('?' * (len(column_names_tranfer[2:]) + 3))})")
                new_cur.executemany(sql, new_rows_frags)

                new_cur.executemany(
                    f"INSERT INTO radius{radius} (env_id, core_smi_id) VALUES (?, ?)",
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


def create_indices(conn: sqlite3.Connection, radii: List[int], verbose: bool = True):
    """
    Create optimized indices on the new database.

    Args:
        conn: SQLite connection
        radii: List of radius values
        verbose: Print progress
    """
    cur = conn.cursor()

    indices = [
        ("idx_envs_env", "CREATE INDEX IF NOT EXISTS idx_envs_env ON envs(env)"),
        ("idx_frags_core_smi", "CREATE INDEX IF NOT EXISTS idx_frags_core_smi ON frags(core_smi)"),
        # ("idx_frags_core_smi_h_id", "CREATE INDEX IF NOT EXISTS idx_frags_core_smi_h_id ON frags(core_smi_h_id)"),
        # ("idx_frags_dist2", "CREATE INDEX IF NOT EXISTS idx_frags_dist2 ON frags(dist2)"),
        # ("idx_frags_h_id", "CREATE INDEX IF NOT EXISTS idx_frags_h_id ON frags_h(core_smi_h_id)"),
        ("idx_frags_h_smi", "CREATE INDEX IF NOT EXISTS idx_frags_h_smi ON frags_h(smi)"),
    ]

    # Add indices for each radius table
    for radius in radii:
        indices.extend([
            (f"idx_radius{radius}_env_id",
             f"CREATE INDEX IF NOT EXISTS idx_radius{radius}_env_id ON radius{radius}(env_id)"),
            (f"idx_radius{radius}_core_smi_id",
             f"CREATE INDEX IF NOT EXISTS idx_radius{radius}_core_smi_id ON radius{radius}(core_smi_id)"),
            (f"idx_radius{radius}_both",
             f"CREATE INDEX IF NOT EXISTS idx_radius{radius}_both ON radius{radius}(env_id, core_smi_id)"),
        ])

    for idx_name, sql in tqdm(indices, desc="Creating indices", disable=not verbose):
        cur.execute(sql)

    conn.commit()

    # Analyze database for query optimization
    if verbose:
        print("Analyzing database...")
    cur.execute("ANALYZE")
    conn.commit()


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
                SELECT f.core_smi, f.core_sma, f.dist2
                FROM radius{radius} r
                JOIN frags f ON r.core_smi_id = f.core_smi_id
                JOIN envs e ON r.env_id = e.env_id
                WHERE e.env = ? AND f.core_smi = ?
            """, (env, core_smi)).fetchone()

            if result is None:
                mismatches += 1
                if verbose:
                    print(f"Missing entry: env={env[:50]}..., core_smi={core_smi}")
            elif result != (core_smi, core_sma, dist2):
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
  python convert_crem_db.py old.db new.db
  python convert_crem_db.py old.db new.db --radii 1,2,3,4,5
  python convert_crem_db.py old.db new.db --batch-size 5000 --verify
        """
    )

    parser.add_argument('old_db', help='Path to existing database')
    parser.add_argument('new_db', help='Path to new database (will be created)')
    parser.add_argument('--radii', default='1,2,3,4,5',
                       help='Comma-separated list of radii to convert (default: 1,2,3,4,5)')
    parser.add_argument('--batch-size', type=int, default=10000,
                       help='Number of rows to process per batch (default: 10000)')
    parser.add_argument('--verify', action='store_true',
                       help='Verify conversion after completion')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Parse radii
    radii = [int(r.strip()) for r in args.radii.split(',')]

    # Check if output file exists
    import os
    if os.path.exists(args.new_db):
        response = input(f"{args.new_db} already exists. Overwrite? (y/N): ")
        if response.lower() != 'y':
            print("Conversion cancelled.")
            return
        os.remove(args.new_db)

    # Perform conversion
    verbose = not args.quiet

    try:
        convert_database(args.old_db, args.new_db, radii, args.batch_size, verbose)

        # Verify if requested
        if args.verify:
            for radius in radii:
                verify_conversion(args.old_db, args.new_db, radius, verbose=verbose)

        print("\n✓ Conversion completed successfully!")

    except Exception as e:
        print(f"\n✗ Conversion failed: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
