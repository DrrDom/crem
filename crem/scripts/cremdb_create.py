#!/usr/bin/env python3
"""
Stream-like fragment DB creation for the new schema.

Processes input SMILES in chunks, aggregates counts per chunk, and updates
the output DB on the fly. A processed-chunks file can be used to resume.
The set name can be a single column name or one or more files with SMILES ids.
"""

import argparse
import io
import os
import re
import sqlite3
import sys
import time
from collections import defaultdict
from itertools import permutations
from multiprocessing import Pool, cpu_count

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi, rdMMPA

from crem.mol_context import get_std_context_core_permutations
from crem.scripts.cremdb_convert import create_indices


def _validate_set_name(set_name):
    if not re.match(r'^[A-Za-z_][A-Za-z0-9_]*$', set_name):
        raise ValueError(
            "set_name must be a valid SQLite identifier (letters, numbers, underscores; cannot start with a number)"
        )
    if set_name in ('env_id', 'core_smi_id'):
        raise ValueError("set_name cannot be env_id or core_smi_id")
    return set_name


def _read_ids(path):
    ids = []
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            ids.append(line)
    return set(ids)


def _resolve_set_names(values):
    if not values:
        raise ValueError("At least one set name or set file must be specified")

    existing_files = [v for v in values if os.path.exists(v)]
    missing_files = [v for v in values if not os.path.exists(v)]

    if existing_files:
        if len(missing_files) > 1:
            raise ValueError(
                "Mixed set names and files are not allowed. "
                f"Missing files: {', '.join(missing_files)}"
            )

        set_names = []
        set_ids = {}
        for path in existing_files:
            base = os.path.basename(path)
            name = os.path.splitext(base)[0]
            name = _validate_set_name(name)
            if name in set_ids:
                raise ValueError(f"Duplicate set name derived from files: {name}")
            set_names.append(name)
            set_ids[name] = _read_ids(path)

        if len(missing_files) == 1:
            base = os.path.basename(missing_files[0])
            default_name = _validate_set_name(os.path.splitext(base)[0])
            if default_name in set_ids:
                raise ValueError(f"Duplicate set name: {default_name}")
            set_names.append(default_name)
            set_ids[default_name] = None

        return set_names, set_ids

    if len(values) != 1:
        raise ValueError("Provide a single set name or one or more set files")
    return [_validate_set_name(values[0])], None


def _replace_attachment_points_with_h(smiles):
    return Chem.CanonSmiles(smiles.replace('*', 'H'))


def _iter_chunks(handle, chunk_size):
    buf = []
    for line in handle:
        buf.append(line)
        if len(buf) >= chunk_size:
            yield buf
            buf = []
    if buf:
        yield buf


def _parse_line(line, sep):
    line = line.strip()
    if not line:
        return None, None
    if sep is None:
        parts = line.split()
    else:
        parts = line.split(sep)
    if not parts:
        return None, None
    smi = parts[0]
    smi_id = parts[1] if len(parts) > 1 else ''
    return smi, smi_id


def _fragment_mol(smi, smi_id, mode):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return set()

    outlines = set()

    if mode == 0 or mode == 1:
        frags = rdMMPA.FragmentMol(
            mol,
            pattern="[!#1]!@!=!#[!#1]",
            maxCuts=4,
            resultsAsMols=False,
            maxCutBonds=30,
        )
        frags += rdMMPA.FragmentMol(
            mol,
            pattern="[!#1]!@!=!#[!#1]",
            maxCuts=3,
            resultsAsMols=False,
            maxCutBonds=30,
        )
        frags = set(frags)
        for core, chains in frags:
            outlines.add((smi, smi_id, core, chains))

    if mode == 0 or mode == 2:
        mol = Chem.AddHs(mol)
        n = mol.GetNumAtoms() - mol.GetNumHeavyAtoms()
        if n < 60:  # TODO: remove this limit, it is not very reasonable
            frags = rdMMPA.FragmentMol(
                mol,
                pattern="[#1]!@!=!#[!#1]",
                maxCuts=1,
                resultsAsMols=False,
                maxCutBonds=100,   # TODO: why we need this?
            )
            for core, chains in frags:
                outlines.add((smi, smi_id, core, chains))

    return outlines


def _core_dist2(core_smi):
    if core_smi.count('*') == 2:
        mol = Chem.MolFromSmiles(core_smi, sanitize=False)
        if mol is None:
            return 0
        mat = Chem.GetDistanceMatrix(mol)
        ids = [a.GetIdx() for a in mol.GetAtoms() if not a.GetAtomicNum()]
        if len(ids) == 2:
            return int(mat[ids[0], ids[1]])
    return 0


def _env_core_from_fragment(core, context, radius, keep_stereo, max_heavy_atoms):
    output = []

    if not core:
        residues = context.split('.')
        if len(residues) == 2:
            for ctx, c in permutations(residues, 2):
                if ctx == '[H][*:1]':
                    continue
                mm = Chem.MolFromSmiles(c, sanitize=False)
                num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float('inf')
                if num_heavy_atoms <= max_heavy_atoms:
                    env, cores = get_std_context_core_permutations(ctx, c, radius, keep_stereo)
                    if env and cores:
                        output.append((env, cores[0], num_heavy_atoms))  # only one item in cores
    else:
        mm = Chem.MolFromSmiles(core, sanitize=False)
        num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float('inf')
        if num_heavy_atoms <= max_heavy_atoms:
            env, cores = get_std_context_core_permutations(context, core, radius, keep_stereo)
            if env and cores:
                for core_smi in cores:
                    output.append((env, core_smi, num_heavy_atoms))

    return output


def _init_worker(radii, keep_stereo, max_heavy_atoms, mode, sep, set_names, set_ids):
    global _RADII
    global _KEEP_STEREO
    global _MAX_HEAVY_ATOMS
    global _MODE
    global _SEP
    global _SET_NAMES
    global _SET_IDS
    _RADII = radii
    _KEEP_STEREO = keep_stereo
    _MAX_HEAVY_ATOMS = max_heavy_atoms
    _MODE = mode
    _SEP = sep
    _SET_NAMES = set_names
    _SET_IDS = set_ids
    RDLogger.DisableLog('rdApp.warning')


def _process_chunk(task):
    chunk_id, lines = task
    envs = set()
    core_info = {}
    counts = {name: {r: defaultdict(int) for r in _RADII} for name in _SET_NAMES}
    stats = {
        "lines": 0,
        "fragments": 0,
        "pairs": 0,
    }

    for line in lines:
        smi, smi_id = _parse_line(line, _SEP)
        if not smi:
            continue
        stats["lines"] += 1
        if _SET_IDS is None:
            member_sets = _SET_NAMES
        else:
            member_sets = []
            for name, ids in _SET_IDS.items():
                if ids is None:
                    member_sets.append(name)
                elif smi_id in ids:
                    member_sets.append(name)
            if not member_sets:
                continue

        frags = _fragment_mol(smi, smi_id, _MODE)
        stats["fragments"] += len(frags)
        for _, _, core, context in frags:
            for radius in _RADII:
                for env, core_smi, core_num_atoms in _env_core_from_fragment(
                    core,
                    context,
                    radius,
                    _KEEP_STEREO,
                    _MAX_HEAVY_ATOMS,
                ):
                    for set_name in member_sets:
                        counts[set_name][radius][(env, core_smi)] += 1
                    envs.add(env)
                    if core_smi not in core_info:
                        dist2 = _core_dist2(core_smi)
                        core_smi_h = _replace_attachment_points_with_h(core_smi)
                        mol_h = Chem.MolFromSmiles(core_smi_h)
                        inchi_val = inchi.MolToInchi(mol_h) if mol_h else None
                        core_info[core_smi] = (core_num_atoms, dist2, core_smi_h, inchi_val)
                    stats["pairs"] += 1

    return chunk_id, envs, core_info, counts, stats


def _read_chunk_ids(path):
    if not path or not os.path.exists(path):
        return set()
    ids = set()
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                ids.add(int(line))
            except ValueError:
                continue
    return ids


def _ensure_schema(conn, radii, set_names):
    conn.execute("PRAGMA user_version = 1")
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA temp_store = MEMORY")

    conn.execute("""
        CREATE TABLE IF NOT EXISTS envs(
            env_id INTEGER PRIMARY KEY AUTOINCREMENT,
            env TEXT NOT NULL UNIQUE
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS frags_h(
            core_smi_h_id INTEGER PRIMARY KEY AUTOINCREMENT,
            smi TEXT NOT NULL,
            inchi TEXT NOT NULL UNIQUE
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS frags(
            core_smi_id INTEGER PRIMARY KEY AUTOINCREMENT,
            core_smi TEXT NOT NULL UNIQUE,
            core_num_atoms INTEGER NOT NULL,
            dist2 INTEGER NOT NULL,
            core_smi_h_id INTEGER NOT NULL,
            FOREIGN KEY (core_smi_h_id) REFERENCES frags_h(core_smi_h_id)
        )
    """)

    for radius in radii:
        conn.execute(f"""
            CREATE TABLE IF NOT EXISTS radius{radius}(
                env_id INTEGER NOT NULL,
                core_smi_id INTEGER NOT NULL,
                FOREIGN KEY (env_id) REFERENCES envs(env_id),
                FOREIGN KEY (core_smi_id) REFERENCES frags(core_smi_id),
                UNIQUE (env_id, core_smi_id)
            )
        """)
        cols = {row[1] for row in conn.execute(f"PRAGMA table_info(radius{radius})")}
        for set_name in set_names:
            if set_name not in cols:
                conn.execute(f"ALTER TABLE radius{radius} ADD COLUMN {set_name} INTEGER NOT NULL DEFAULT 0")

    conn.commit()


def _init_temp_tables(conn):
    conn.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_envs(env TEXT PRIMARY KEY)")
    conn.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_cores(
            core_smi TEXT PRIMARY KEY,
            core_num_atoms INTEGER NOT NULL,
            dist2 INTEGER NOT NULL,
            core_smi_h TEXT NOT NULL,
            inchi TEXT NOT NULL
        )
    """)
    conn.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_counts(
            env TEXT NOT NULL,
            core_smi TEXT NOT NULL,
            cnt INTEGER NOT NULL
        )
    """)


def _clear_temp_tables(conn):
    conn.execute("DELETE FROM tmp_envs")
    conn.execute("DELETE FROM tmp_cores")
    conn.execute("DELETE FROM tmp_counts")


def _load_tmp_envs(conn, envs):
    if envs:
        conn.executemany("INSERT OR IGNORE INTO tmp_envs (env) VALUES (?)", [(e,) for e in envs])


def _load_tmp_cores(conn, core_info):
    rows = []
    for core_smi, (core_num_atoms, dist2, core_smi_h, inchi_val) in core_info.items():
        if not inchi_val:
            continue
        rows.append((core_smi, core_num_atoms, dist2, core_smi_h, inchi_val))
    if rows:
        conn.executemany(
            "INSERT OR IGNORE INTO tmp_cores (core_smi, core_num_atoms, dist2, core_smi_h, inchi) "
            "VALUES (?, ?, ?, ?, ?)",
            rows,
        )


def _ensure_env_ids(conn):
    conn.execute("INSERT OR IGNORE INTO envs (env) SELECT env FROM tmp_envs")


def _ensure_inchi_ids(conn):
    conn.execute("INSERT OR IGNORE INTO frags_h (smi, inchi) SELECT core_smi_h, inchi FROM tmp_cores")


def _ensure_core_ids(conn):
    conn.execute("""
        INSERT OR IGNORE INTO frags (core_smi, core_num_atoms, dist2, core_smi_h_id)
        SELECT t.core_smi, t.core_num_atoms, t.dist2, h.core_smi_h_id
        FROM tmp_cores t
        JOIN frags_h h ON t.inchi = h.inchi
    """)


def _update_radius_counts(conn, counts, set_names):
    for set_name in set_names:
        per_set = counts.get(set_name, {})
        for radius, mapping in per_set.items():
            if not mapping:
                continue
            conn.execute("DELETE FROM tmp_counts")
            conn.executemany(
                "INSERT INTO tmp_counts (env, core_smi, cnt) VALUES (?, ?, ?)",
                [(env, core_smi, cnt) for (env, core_smi), cnt in mapping.items()],
            )
            conn.execute(f"""
                INSERT INTO radius{radius} (env_id, core_smi_id, {set_name})
                SELECT e.env_id, f.core_smi_id, SUM(tc.cnt)
                FROM tmp_counts tc
                JOIN envs e ON e.env = tc.env
                JOIN frags f ON f.core_smi = tc.core_smi
                GROUP BY e.env_id, f.core_smi_id
                ON CONFLICT(env_id, core_smi_id)
                DO UPDATE SET {set_name} = {set_name} + excluded.{set_name}
            """)


def _open_input(path, force_zstd):
    if force_zstd or path.endswith(".zst"):
        try:
            import zstandard as zstd
        except ImportError as exc:
            raise RuntimeError("zstandard is required to read .zst files") from exc
        fh = open(path, "rb")
        dctx = zstd.ZstdDecompressor()
        stream = dctx.stream_reader(fh)
        return io.TextIOWrapper(stream), fh
    return open(path, "rt", encoding="utf-8", errors="replace"), None


def main():
    parser = argparse.ArgumentParser(description="Stream-like fragment DB creation for new schema")
    parser.add_argument("input", help="Input SMILES file (text or .zst)")
    parser.add_argument("output_db", help="Output SQLite DB file")
    parser.add_argument(
        "set_name",
        nargs="+",
        help=(
            "Set name (single) and/or one or more files with SMILES ids; file basenames "
            "become set names in the DB, so they must follow SQLite column naming rules "
            "(letters, numbers, underscores; cannot start with a number)"
        ),
    )
    parser.add_argument(
        "--radii",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4, 5],
        help="Space-separated radii (default: 1 2 3 4 5)",
    )
    parser.add_argument("--chunk-size", type=int, default=100, help="Lines per chunk (default: 100)")
    parser.add_argument("--max-heavy-atoms", type=int, default=15, help="Max heavy atoms in core (default: 15)")
    parser.add_argument("--keep-stereo", action="store_true", help="Keep stereo in env/core")
    parser.add_argument("--mode", type=int, choices=[0, 1, 2], default=0, help="Fragmentation mode: 0 - all atoms, 1 - only heavy atoms, 2 - only hydrogen atoms (default: 0)")
    parser.add_argument("--sep", default=None, help="Input delimiter (default: whitespace)")
    parser.add_argument("--processed-chunks", default=None, help="Path to processed chunks file (append)")
    parser.add_argument("--exclude-chunks", default=None, help="Path to chunk list to skip")
    parser.add_argument("--zstd", action="store_true", help="Force zstd input")
    parser.add_argument("--log-every", type=int, default=None, help="Log progress every N chunks (default: None)")
    parser.add_argument("--ncpu", type=int, default=1, help="Number of worker processes (default: 1)")
    args = parser.parse_args()

    set_names, set_ids = _resolve_set_names(args.set_name)
    radii = sorted(set(args.radii))
    if not radii:
        raise ValueError("At least one radius must be specified")

    exclude_chunks = _read_chunk_ids(args.exclude_chunks)
    processed_chunks = _read_chunk_ids(args.processed_chunks)
    skip_chunks = exclude_chunks | processed_chunks

    RDLogger.DisableLog("rdApp.warning")

    with sqlite3.connect(args.output_db) as conn:
        _ensure_schema(conn, radii, set_names)
        _init_temp_tables(conn)
        cur = conn.cursor()

        processed_handle = None
        if args.processed_chunks:
            processed_handle = open(args.processed_chunks, "a", encoding="utf-8")

        input_handle, zstd_handle = _open_input(args.input, args.zstd)

        total_stats = {"lines": 0, "fragments": 0, "pairs": 0}
        chunks_processed = 0
        chunks_skipped = 0
        start_time = time.time()

        ncpu = min(cpu_count(), max(args.ncpu, 1))
        pool = Pool(
            ncpu,
            initializer=_init_worker,
            initargs=(
                radii,
                args.keep_stereo,
                args.max_heavy_atoms,
                args.mode,
                args.sep,
                set_names,
                set_ids,
            ),
        )

        try:
            def task_iter():
                nonlocal chunks_skipped
                for chunk_id, lines in enumerate(_iter_chunks(input_handle, args.chunk_size)):
                    if chunk_id in skip_chunks:
                        chunks_skipped += 1
                        continue
                    yield (chunk_id, lines)

            for chunk_id, envs, core_info, counts, stats in pool.imap_unordered(
                _process_chunk, task_iter(), chunksize=1
            ):
                chunk_start_time = time.time()

                conn.execute("BEGIN")
                _clear_temp_tables(conn)
                _load_tmp_envs(conn, envs)
                _load_tmp_cores(conn, core_info)
                _ensure_env_ids(conn)
                _ensure_inchi_ids(conn)
                _ensure_core_ids(conn)
                _update_radius_counts(conn, counts, set_names)
                conn.commit()

                chunks_processed += 1
                total_stats["lines"] += stats["lines"]
                total_stats["fragments"] += stats["fragments"]
                total_stats["pairs"] += stats["pairs"]

                if processed_handle:
                    processed_handle.write(f"{chunk_id}\n")
                    processed_handle.flush()

                if args.log_every and chunks_processed % args.log_every == 0:
                    elapsed = time.time() - start_time
                    rate = total_stats["lines"] / elapsed if elapsed > 0 else 0
                    sys.stderr.write(
                        f"\rChunks: {chunks_processed} processed, {chunks_skipped} skipped | "
                        f"mols: {total_stats['lines']} | frags: {total_stats['fragments']} | "
                        f"pairs: {total_stats['pairs']} | {rate:.1f} mol/s"
                    )
                    sys.stderr.flush()

            create_indices(conn, radii, True)

        finally:
            pool.close()
            pool.join()
            input_handle.close()
            if zstd_handle:
                zstd_handle.close()
            if processed_handle:
                processed_handle.close()

        sys.stderr.write("\n")

        env_count = cur.execute("SELECT COUNT(*) FROM envs").fetchone()[0]
        frag_count = cur.execute("SELECT COUNT(*) FROM frags").fetchone()[0]
        frag_h_count = cur.execute("SELECT COUNT(*) FROM frags_h").fetchone()[0]
        radius_counts = {}
        for radius in radii:
            radius_counts[radius] = cur.execute(f"SELECT COUNT(rowid) FROM radius{radius}").fetchone()[0]

        elapsed = time.time() - start_time
        print("Done.")
        print(f"Chunks processed: {chunks_processed}")
        print(f"Chunks skipped: {chunks_skipped}")
        print(f"Molecules processed: {total_stats['lines']}")
        print(f"Fragments generated: {total_stats['fragments']}")
        print(f"Env/core pairs: {total_stats['pairs']}")
        print(f"Unique envs: {env_count}")
        print(f"Unique frags: {frag_count}")
        print(f"Unique frags_h: {frag_h_count}")
        for radius, count in radius_counts.items():
            print(f"radius{radius} rows: {count}")
        print(f"Elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
