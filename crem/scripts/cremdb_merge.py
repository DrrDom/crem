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
import os
import shutil
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor
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
    target_conn.execute("DROP INDEX IF EXISTS idx_frags_h_smi")

    # _merge_history makes repeated merges of the same source idempotent.
    # Resume after a crash between commit and source-unlink would otherwise
    # double-count radius rows on the second invocation.
    target_conn.execute(
        "CREATE TABLE IF NOT EXISTS _merge_history("
        "source_id TEXT PRIMARY KEY, merged_at TEXT)"
    )
    target_conn.commit()

    for source_path in tqdm(source_paths, desc="Merging shards", disable=not verbose):
        source_id = os.path.basename(source_path)

        # Cheap pre-check: skip already-absorbed sources before opening the
        # source DB. The authoritative check is inside the transaction below
        # so that a concurrent absorber and we don't both insert + merge.
        already = target_conn.execute(
            "SELECT 1 FROM _merge_history WHERE source_id = ?",
            (source_id,),
        ).fetchone()
        if already:
            if verbose:
                sys.stderr.write(f"  Skipping {source_path} (already absorbed)\n")
                sys.stderr.flush()
            continue

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
                INSERT OR IGNORE INTO main.frags(core_smi, core_smi_h_id)
                SELECT sf.core_smi, mh.core_smi_h_id
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
            # Metadata columns (env_id, core_smi_id, core_num_atoms, dist2,
            # is_ring_closure) are handled explicitly; everything else is a
            # per-set count column to be summed on conflict.
            metadata_cols = {
                "env_id", "core_smi_id", "core_num_atoms", "dist2", "is_ring_closure",
            }
            for radius in radii:
                src_cols = {
                    row[1]
                    for row in target_conn.execute(
                        f"PRAGMA src.table_info(radius{radius})"
                    ).fetchall()
                }
                src_set_cols = sorted(src_cols - metadata_cols)
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

                # Older v1 shards may lack is_ring_closure; default to 0 for them.
                src_has_ring = "is_ring_closure" in src_cols
                ring_select = "r.is_ring_closure" if src_has_ring else "0"

                col_list = ", ".join(src_set_cols)
                src_vals = ", ".join(f"r.{c}" for c in src_set_cols)
                conflict_upd = ", ".join(
                    f"{c} = {c} + excluded.{c}" for c in src_set_cols
                )

                target_conn.execute(f"""
                    INSERT INTO main.radius{radius}(
                        env_id, core_smi_id, core_num_atoms, dist2,
                        is_ring_closure, {col_list}
                    )
                    SELECT em.dst_id, fm.dst_id, r.core_num_atoms, r.dist2,
                           {ring_select}, {src_vals}
                    FROM src.radius{radius} r
                    JOIN temp._env_map em  ON r.env_id      = em.src_id
                    JOIN temp._frag_map fm ON r.core_smi_id = fm.src_id
                    ORDER BY em.dst_id, fm.dst_id
                    ON CONFLICT(env_id, core_smi_id, is_ring_closure)
                    DO UPDATE SET {conflict_upd}
                """)

            # Record absorption in the same transaction as the merge so that
            # commit-or-nothing semantics make merge_into idempotent.
            target_conn.execute(
                "INSERT OR IGNORE INTO _merge_history(source_id, merged_at) "
                "VALUES (?, strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))",
                (source_id,),
            )

            target_conn.commit()
            target_conn.execute("PRAGMA foreign_keys = ON")

        except Exception:
            target_conn.rollback()
            raise

        finally:
            target_conn.execute("DROP TABLE IF EXISTS temp._env_map")
            target_conn.execute("DROP TABLE IF EXISTS temp._frag_map")
            target_conn.execute("DETACH DATABASE src")


def _pair_merge(target_path: str, source_path: str) -> None:
    """Worker entry point: absorb source_path into target_path.

    Runs in a child process under ProcessPoolExecutor. Each call opens its
    own SQLite connection on a different file pair so the parallel workers
    do not contend. Indices are not rebuilt — the orchestrator calls
    create_indices once on the final survivor.
    """
    with sqlite3.connect(target_path) as conn:
        conn.execute("PRAGMA cache_size = -262144")
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA temp_store = MEMORY")
        merge_into(conn, [source_path], verbose=False)


def binary_tree_merge_into(
    final_target_path: str,
    shard_paths: List[str],
    *,
    max_parallel: int,
    rebuild_index: bool = True,
    verbose: bool = True,
) -> None:
    """Merge shard_paths into final_target_path via parallel pair reduction.

    Pairs of shards are merged concurrently in each round; the surviving
    shards advance to the next round; an odd one out is carried forward.
    After the final round the single survivor is moved to final_target_path
    and create_indices is called once.

    Args:
        final_target_path: Destination DB. Overwritten on success.
        shard_paths: Shard DB paths. Sources are unlinked as they're absorbed.
        max_parallel: Max concurrent pair-merges per round.
        rebuild_index: Rebuild covering indices on the final survivor.
        verbose: Print per-round progress to stderr.
    """
    if not shard_paths:
        raise ValueError("shard_paths is empty")

    paths = sorted(str(p) for p in shard_paths)
    round_idx = 0
    while len(paths) > 1:
        round_idx += 1
        pairs = []
        for i in range(0, len(paths) - 1, 2):
            pairs.append((paths[i], paths[i + 1]))
        carry = paths[-1] if len(paths) % 2 else None

        if verbose:
            sys.stderr.write(
                f"  Tree merge round {round_idx}: {len(pairs)} pair(s)"
                f"{', 1 carry' if carry else ''}\n"
            )
            sys.stderr.flush()

        with ProcessPoolExecutor(max_workers=max(1, max_parallel)) as ex:
            futures = [ex.submit(_pair_merge, t, s) for t, s in pairs]
            for f in futures:
                f.result()  # raises on worker failure

        # Sources have been absorbed (and recorded in _merge_history on each
        # target). Unlink them; if a unlink fails because the file is already
        # gone (resumed run), that's fine.
        for _, s in pairs:
            try:
                os.unlink(s)
            except FileNotFoundError:
                pass

        paths = [t for t, _ in pairs]
        if carry:
            paths.append(carry)

    survivor = paths[0]
    if os.path.abspath(survivor) != os.path.abspath(final_target_path):
        # On overwrite, remove an existing destination first to make shutil.move
        # behave consistently across same-FS and cross-FS cases.
        if os.path.exists(final_target_path):
            os.unlink(final_target_path)
        shutil.move(survivor, final_target_path)

    if rebuild_index:
        with sqlite3.connect(final_target_path) as conn:
            conn.execute("PRAGMA cache_size = -262144")
            conn.execute("PRAGMA synchronous = NORMAL")
            conn.execute("PRAGMA journal_mode = WAL")
            conn.execute("PRAGMA temp_store = MEMORY")
            radii = _get_radii(conn)
            from crem.scripts.cremdb_create import create_indices  # lazy: avoids circular import
            create_indices(conn, radii, verbose=verbose)


def run(
    target_path: str,
    source_paths: List[str],
    rebuild_index: bool = True,
    verbose: bool = True,
    parallel: int = 1,
) -> None:
    """Merge source shard databases into target and optionally rebuild indices.

    Args:
        target_path: Path to the target (base) database. Must already exist.
        source_paths: Paths to source shard databases to merge in.
        rebuild_index: Recreate covering indices on the target after merge.
        verbose: Print per-shard progress to stderr.
        parallel: When > 1, merge with binary-tree reduction using up to
            this many concurrent pair-merges per round. The target is treated
            as one of the contributors and the final survivor is moved back
            to target_path.
    """
    if parallel < 1:
        raise ValueError("parallel must be >= 1")

    if parallel > 1:
        binary_tree_merge_into(
            target_path,
            [str(target_path), *(str(s) for s in source_paths)],
            max_parallel=parallel,
            rebuild_index=rebuild_index,
            verbose=verbose,
        )
        if verbose:
            print(
                f"Tree-merged {len(source_paths)} shard(s) into {target_path}."
            )
        return

    with sqlite3.connect(target_path) as conn:
        conn.execute("PRAGMA cache_size = -262144")
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA temp_store = MEMORY")

        radii = _get_radii(conn)
        merge_into(conn, [str(s) for s in source_paths], verbose=verbose)

        if rebuild_index:
            from crem.scripts.cremdb_create import create_indices  # lazy: avoids circular import
            create_indices(conn, radii, verbose=verbose)

    if verbose:
        print(f"Merged {len(source_paths)} shard(s) into {target_path}.")


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
        "--parallel",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Merge with binary-tree reduction using up to N concurrent "
            "pair-merges per round. Default 1 (serial). With N>1 the target "
            "is one of the contributors and the final survivor is moved back "
            "to the target path."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=True,
    )
    args = parser.parse_args()
    run(
        target_path=args.target,
        source_paths=args.input,
        rebuild_index=not args.no_index,
        verbose=args.verbose,
        parallel=args.parallel,
    )


if __name__ == "__main__":
    main()
