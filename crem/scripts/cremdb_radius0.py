#!/usr/bin/env python3
"""
Add a radius0 table to an existing v2 CReM database, derived from `frags`.

Radius 0 imposes no environment constraint: the only thing left of an environment is how
many attachment points a fragment has and which of them close a ring, encoded as the
`RxAy` env class. Everything a radius-0 row needs other than occurrence counts is a
function of `frags.core_smi` alone, which is what makes this script possible at all - under
the v2 convention the isotope label on ring cuts carries the provenance, so `frags` is a
self-sufficient description of the fragment universe.

Occurrence counts are NOT recoverable this way. A radius table's counts are occurrences
multiplied by the env symmetry orbit size, and that factor depends on the radius, so no
sum over them is a radius-independent occurrence count. What IS exactly recoverable is
per-set *membership*, because membership is a boolean OR over rows rather than a sum. Rows
written here therefore carry 1 where the fragment occurs in a set and 0 where it does not -
1 being the true lower bound, since anything in `frags` was observed at least once.

    min_freq=0 and min_freq=1 behave as "no filter" on such a table; min_freq >= 2 matches
    nothing, because the counts needed to answer it were never recorded.

For real counts, build the database with `cremdb_create --radii 0 ...`, which counts one
increment per fragmentation event.

Usage:
    cremdb_radius0 -i fragments.db
"""

import argparse
import sqlite3
import sys
from itertools import islice
from multiprocessing import Pool

from rdkit import RDLogger

from crem.db import _RESERVED_RADIUS_COLUMNS
from crem.mol_context import RADIUS0_ENV_CLASSES, get_radius0_rows
from crem.scripts.cremdb_create import (_core_dist2, _count_heavy_atoms,
                                        _replace_attachment_points_with_h, create_indices)


def _radius_tables(con):
    """Existing radius tables with radius > 0, sorted."""
    out = []
    for (name,) in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'radius%'"):
        try:
            radius = int(name[6:])
        except ValueError:
            continue
        if radius > 0:
            out.append(radius)
    return sorted(out)


def _set_columns(con, radius):
    cols = {row[1] for row in con.execute(f"PRAGMA table_info(radius{radius})")}
    return sorted(cols - _RESERVED_RADIUS_COLUMNS)


def _check_database(con, force):
    version = con.execute("PRAGMA user_version").fetchone()[0]
    if version < 2:
        raise ValueError(
            f"this database is v{version}; deriving radius 0 needs v2, where the ring-cut "
            f"isotope label makes frags self-sufficient. On v1 the two-attachment ring arcs "
            f"are unlabelled and their provenance cannot be read from the fragment string."
        )
    radii = _radius_tables(con)
    if not radii:
        raise ValueError("no radius tables found; is this a CReM fragment database?")
    existing = con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='radius0'").fetchone()
    if existing:
        n = con.execute("SELECT COUNT(*) FROM radius0").fetchone()[0]
        if n and not force:
            raise ValueError(
                f"radius0 already exists with {n} rows. It may have been built by "
                f"cremdb_create with real occurrence counts, which this script would "
                f"replace with membership flags. Pass --force to overwrite it."
            )
    return radii


def _create_radius0(con, set_cols, freq_types):
    con.execute("DROP TABLE IF EXISTS radius0")
    defs = ["env_id INTEGER NOT NULL", "core_smi_id INTEGER NOT NULL",
            "core_num_atoms INTEGER NOT NULL", "dist2 INTEGER NOT NULL"]
    defs += [f"{c} {freq_types.get(c, 'INTEGER')} NOT NULL DEFAULT 0" for c in set_cols]
    con.execute(f"""CREATE TABLE radius0({', '.join(defs)},
                    FOREIGN KEY (env_id) REFERENCES envs(env_id),
                    FOREIGN KEY (core_smi_id) REFERENCES frags(core_smi_id),
                    UNIQUE (env_id, core_smi_id))""")
    con.executemany("INSERT OR IGNORE INTO envs(env) VALUES (?)",
                    [(e,) for e in RADIUS0_ENV_CLASSES])


def _worker(core_smi):
    """Fragment identity and its stored orientations. Runs in a worker process."""
    env, rows = get_radius0_rows(core_smi)
    if not env or not rows:
        return core_smi, None, None, ()
    return core_smi, env, min(rows), rows


def _init_worker():
    RDLogger.DisableLog('rdApp.*')


def _chunks(iterable, size):
    it = iter(iterable)
    while True:
        block = list(islice(it, size))
        if not block:
            return
        yield block


def build_radius0(db_path, ncpu=1, batch_size=50000, verbose=True, force=False):
    """Derive a radius0 table from `frags` in an existing v2 database."""
    RDLogger.DisableLog('rdApp.*')
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA synchronous = NORMAL")
    con.execute("PRAGMA cache_size = -262144")
    con.execute("PRAGMA temp_store = MEMORY")

    radii = _check_database(con, force)
    set_cols = _set_columns(con, radii[0])
    if not set_cols:
        raise ValueError(f"radius{radii[0]} has no per-set count columns; database is malformed")
    freq_types = {row[1]: row[2] for row in con.execute(f"PRAGMA table_info(radius{radii[0]})")}
    if verbose:
        sys.stderr.write(f"radius tables: {radii}; fragment sets: {set_cols}\n")

    _create_radius0(con, set_cols, freq_types)
    con.execute("""CREATE TEMP TABLE _frag_map(
                       core_smi_id INTEGER PRIMARY KEY, frag_key TEXT NOT NULL)""")
    con.execute("""CREATE TEMP TABLE _rows(
                       frag_key TEXT NOT NULL, env TEXT NOT NULL, core_smi TEXT NOT NULL,
                       UNIQUE (frag_key, core_smi))""")
    con.commit()

    # 1. map every frags row to its fragment, and collect that fragment's orientations once.
    #    get_radius0_rows is labelling-invariant (see tests), so all labellings of one
    #    fragment agree on both the key and the orientation set.
    total = con.execute("SELECT COUNT(*) FROM frags").fetchone()[0]
    cur = con.execute("SELECT core_smi_id, core_smi FROM frags")
    seen_fragments = set()
    done = 0
    pool = Pool(ncpu, initializer=_init_worker) if ncpu > 1 else None
    try:
        for block in _chunks(cur, batch_size):
            smis = [smi for _cid, smi in block]
            results = (pool.map(_worker, smis, chunksize=200) if pool
                       else [_worker(s) for s in smis])
            by_smi = {smi: (env, key, rows) for smi, env, key, rows in results}
            mapping, new_rows = [], []
            for cid, smi in block:
                env, key, rows = by_smi[smi]
                if key is None:
                    continue
                mapping.append((cid, key))
                if key not in seen_fragments:
                    seen_fragments.add(key)
                    new_rows.extend((key, env, o) for o in rows)
            con.executemany("INSERT OR IGNORE INTO _frag_map VALUES (?,?)", mapping)
            con.executemany("INSERT OR IGNORE INTO _rows VALUES (?,?,?)", new_rows)
            con.commit()
            done += len(block)
            if verbose:
                sys.stderr.write(f"\r  mapped {done}/{total} frags rows, "
                                 f"{len(seen_fragments)} fragments")
                sys.stderr.flush()
    finally:
        if pool:
            pool.close()
            pool.join()
    if verbose:
        sys.stderr.write("\n")

    # 2. per-set membership per fragment: a boolean OR over every labelling of the fragment
    #    and every radius table. Unlike a count this is orbit-independent, so it is exact.
    con.execute(f"""CREATE TEMP TABLE _member(
                        frag_key TEXT PRIMARY KEY,
                        {', '.join(f'{c} INTEGER NOT NULL DEFAULT 0' for c in set_cols)})""")
    for radius in radii:
        cols = _set_columns(con, radius)
        present = [c for c in cols if c in set_cols]
        if not present:
            continue
        sel = ", ".join(f"MAX(CASE WHEN r.{c} > 0 THEN 1 ELSE 0 END)" for c in present)
        upd = ", ".join(f"{c} = MAX({c}, excluded.{c})" for c in present)
        con.execute(f"""INSERT INTO _member(frag_key, {', '.join(present)})
                        SELECT m.frag_key, {sel}
                        FROM radius{radius} r JOIN _frag_map m ON r.core_smi_id = m.core_smi_id
                        GROUP BY m.frag_key
                        ON CONFLICT(frag_key) DO UPDATE SET {upd}""")
        con.commit()
        if verbose:
            sys.stderr.write(f"  membership from radius{radius}\n")

    # 3. every orientation needs a frags row. All orientations of a fragment share the
    #    H-collapsed form, so they can reuse the existing frags_h id.
    con.execute("""CREATE TEMP TABLE _need AS
                   SELECT DISTINCT w.core_smi FROM _rows w
                   WHERE NOT EXISTS (SELECT 1 FROM frags f WHERE f.core_smi = w.core_smi)""")
    missing = [r[0] for r in con.execute("SELECT core_smi FROM _need")]
    if verbose:
        sys.stderr.write(f"  {len(missing)} orientations absent from frags; inserting\n")
    for block in _chunks(missing, batch_size):
        rows = []
        for smi in block:
            h = _replace_attachment_points_with_h(smi)
            if h is None:
                continue
            rows.append((smi, h))
        con.executemany("INSERT OR IGNORE INTO frags_h(smi) VALUES (?)",
                        [(h,) for _s, h in rows])
        con.executemany("""INSERT OR IGNORE INTO frags(core_smi, core_smi_h_id)
                           SELECT ?, core_smi_h_id FROM frags_h WHERE smi = ?""", rows)
        con.commit()

    # 4. write the rows: count 1 where the fragment belongs to the set, 0 where it does not
    member_sel = ", ".join(f"COALESCE(mb.{c}, 0)" for c in set_cols)
    con.execute(f"""INSERT OR IGNORE INTO radius0(env_id, core_smi_id, core_num_atoms,
                                                 dist2, {', '.join(set_cols)})
                    SELECT e.env_id, f.core_smi_id, 0, 0, {member_sel}
                    FROM _rows w
                    JOIN envs e ON e.env = w.env
                    JOIN frags f ON f.core_smi = w.core_smi
                    LEFT JOIN _member mb ON mb.frag_key = w.frag_key""")
    con.commit()

    # core_num_atoms and dist2 are pure functions of core_smi; fill them in one pass.
    # dist2 follows the existing convention - only cores with exactly two attachment
    # points carry it; the ring-cut pair distance for 3-/4-point cores is a separate change.
    todo = con.execute("""SELECT r.rowid, f.core_smi FROM radius0 r
                          JOIN frags f ON r.core_smi_id = f.core_smi_id""").fetchall()
    for block in _chunks(todo, batch_size):
        con.executemany("UPDATE radius0 SET core_num_atoms = ?, dist2 = ? WHERE rowid = ?",
                        [(_count_heavy_atoms(smi), _core_dist2(smi), rid)
                         for rid, smi in block])
        con.commit()

    n_rows = con.execute("SELECT COUNT(*) FROM radius0").fetchone()[0]
    create_indices(con, [0], verbose=verbose)
    con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    con.close()
    if verbose:
        sys.stderr.write(f"radius0: {n_rows} rows from {len(seen_fragments)} fragments\n")
    return n_rows


def main():
    parser = argparse.ArgumentParser(
        description="Derive a radius0 table from frags in an existing v2 CReM database",
        epilog="Counts are membership flags (1/0), not occurrence counts - see the module "
               "docstring. Build with cremdb_create --radii 0 for real counts.")
    parser.add_argument('-i', '--input', required=True, help='v2 database to extend')
    parser.add_argument('-c', '--ncpu', type=int, default=1, help='worker processes')
    parser.add_argument('--batch-size', type=int, default=50000, help='rows per batch')
    parser.add_argument('--force', action='store_true',
                        help='overwrite an existing radius0 table')
    parser.add_argument('--quiet', action='store_true', help='suppress progress output')
    args = parser.parse_args()
    build_radius0(args.input, ncpu=args.ncpu, batch_size=args.batch_size,
                  verbose=not args.quiet, force=args.force)


if __name__ == '__main__':
    main()
