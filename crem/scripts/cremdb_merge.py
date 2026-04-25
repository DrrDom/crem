#!/usr/bin/env python3
"""
Merge multiple CReM fragment databases (v1 schema) into a single database.

Typical use: combine shard databases created by cremdb_create.py --shard-size.

Usage:
    cremdb_merge.py -t base.db -i shard_001.db shard_002.db ...

The target (-t) must already exist and contain the base data (e.g. shard_000.db).
Each source database is merged into it in order.  Indices are recreated on the
target after all sources are merged unless --no-index is given.
"""

import argparse
import sqlite3
import sys
from typing import List

from tqdm import tqdm


def _get_radii(conn: sqlite3.Connection) -> List[int]:
    """Return sorted list of radius values found in the database."""
    radii = []
    for (name,) in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'radius%'"
    ):
        try:
            radii.append(int(name[6:]))
        except ValueError:
            pass
    return sorted(radii)


def merge_into(
    target_conn: sqlite3.Connection,
    source_paths: List[str],
    verbose: bool = True,
) -> None:
    """Merge source shard databases into target_conn.

    Drops radius query indices on the target before merging to speed up bulk
    inserts.  Does NOT recreate indices — the caller must call create_indices()
    afterwards.

    Args:
        target_conn: Open SQLite connection to the target (base) database.
        source_paths: Paths to source shard databases to merge in.
        verbose: Print per-shard progress to stderr.
    """
    radii = _get_radii(target_conn)

    # Drop radius query indices on target to avoid per-row index maintenance.
    # The UNIQUE autoindex (sqlite_autoindex_radiusN_1) is kept — it drives
    # ON CONFLICT DO UPDATE.
    for radius in radii:
        for suffix in ("env_id", "core_smi_id", "both", "lookup"):
            target_conn.execute(f"DROP INDEX IF EXISTS idx_radius{radius}_{suffix}")
    target_conn.execute("DROP INDEX IF EXISTS idx_frags_core_num_atoms")
    target_conn.execute("DROP INDEX IF EXISTS idx_frags_h_smi")
    target_conn.commit()

    for source_path in tqdm(source_paths, desc="Merging shards", disable=not verbose):
        if verbose:
            sys.stderr.write(f"  Merging {source_path}...\n")
            sys.stderr.flush()

        target_conn.execute("ATTACH DATABASE ? AS src", (source_path,))
        try:
            target_conn.execute("PRAGMA foreign_keys = OFF")
            target_conn.execute("BEGIN")

            # 1. Merge envs — natural key is env text, no ID translation needed.
            target_conn.execute(
                "INSERT OR IGNORE INTO main.envs(env) SELECT env FROM src.envs"
            )

            # 2. Merge frags_h — natural key is the H-canonical SMILES (smi).
            target_conn.execute(
                "INSERT OR IGNORE INTO main.frags_h(smi) "
                "SELECT smi FROM src.frags_h"
            )

            # 3. Merge frags — core_smi_h_id must be translated via smi.
            target_conn.execute("""
                INSERT OR IGNORE INTO main.frags(core_smi, core_num_atoms, dist2, core_smi_h_id)
                SELECT sf.core_smi, sf.core_num_atoms, sf.dist2, mh.core_smi_h_id
                FROM src.frags sf
                JOIN src.frags_h sh ON sf.core_smi_h_id = sh.core_smi_h_id
                JOIN main.frags_h mh ON sh.smi = mh.smi
            """)

            # Build integer translation tables once per shard (one text JOIN each)
            # so that the per-radius merges use fast integer-to-integer lookups.
            target_conn.execute("""
                CREATE TEMP TABLE _env_map AS
                SELECT se.env_id AS src_id, me.env_id AS dst_id
                FROM src.envs se
                JOIN main.envs me ON se.env = me.env
            """)
            target_conn.execute(
                "CREATE INDEX temp._env_map_src ON _env_map(src_id)"
            )
            target_conn.execute("""
                CREATE TEMP TABLE _frag_map AS
                SELECT sf.core_smi_id AS src_id, mf.core_smi_id AS dst_id
                FROM src.frags sf
                JOIN main.frags mf ON sf.core_smi = mf.core_smi
            """)
            target_conn.execute(
                "CREATE INDEX temp._frag_map_src ON _frag_map(src_id)"
            )

            # 4. Merge each radius table using integer translation maps (no text JOINs).
            for radius in radii:
                src_cols = {
                    row[1]
                    for row in target_conn.execute(
                        f"PRAGMA src.table_info(radius{radius})"
                    ).fetchall()
                }
                src_set_cols = sorted(src_cols - {"env_id", "core_smi_id"})
                if not src_set_cols:
                    continue

                # Ensure target has every set column present in the source.
                dst_cols = {
                    row[1]
                    for row in target_conn.execute(
                        f"PRAGMA main.table_info(radius{radius})"
                    ).fetchall()
                }
                for col in src_set_cols:
                    if col not in dst_cols:
                        target_conn.execute(
                            f"ALTER TABLE main.radius{radius} "
                            f"ADD COLUMN {col} INTEGER NOT NULL DEFAULT 0"
                        )

                col_list = ", ".join(src_set_cols)
                src_vals = ", ".join(f"r.{c}" for c in src_set_cols)
                conflict_upd = ", ".join(
                    f"{c} = {c} + excluded.{c}" for c in src_set_cols
                )

                target_conn.execute(f"""
                    INSERT INTO main.radius{radius}(env_id, core_smi_id, core_num_atoms, dist2, {col_list})
                    SELECT em.dst_id, fm.dst_id, r.core_num_atoms, r.dist2, {src_vals}
                    FROM src.radius{radius} r
                    JOIN temp._env_map em  ON r.env_id      = em.src_id
                    JOIN temp._frag_map fm ON r.core_smi_id = fm.src_id
                    ORDER BY em.dst_id, fm.dst_id
                    ON CONFLICT(env_id, core_smi_id) DO UPDATE SET {conflict_upd}
                """)

            target_conn.commit()
            target_conn.execute("PRAGMA foreign_keys = ON")

        except Exception:
            target_conn.rollback()
            raise

        finally:
            target_conn.execute("DROP TABLE IF EXISTS temp._env_map")
            target_conn.execute("DROP TABLE IF EXISTS temp._frag_map")
            target_conn.execute("DETACH DATABASE src")


def main():
    parser = argparse.ArgumentParser(
        description="Merge individual CReM fragment databases into a single database"
    )
    parser.add_argument(
        "-t", "--target",
        required=True,
        help="Target database to merge into. Must already exist with schema and data.",
    )
    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="Source (shard) database files to merge into the target.",
    )
    parser.add_argument(
        "--no-index",
        action="store_true",
        help="Skip index creation after merge (useful when more shards will be merged later).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=True,
    )
    args = parser.parse_args()

    with sqlite3.connect(args.target) as conn:
        conn.execute("PRAGMA cache_size = -262144")   # 256 MB
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA temp_store = MEMORY")

        radii = _get_radii(conn)
        merge_into(conn, args.input, verbose=args.verbose)

        if not args.no_index:
            # Lazy import to avoid circular dependency with cremdb_create.
            from crem.scripts.cremdb_create import create_indices
            create_indices(conn, radii, verbose=args.verbose)

    print(f"Merged {len(args.input)} shard(s) into {args.target}.")


if __name__ == "__main__":
    main()