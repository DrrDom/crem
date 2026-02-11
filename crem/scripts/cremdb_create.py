#!/usr/bin/env python3
"""
Stream-like fragment DB creation for the new schema.

Processes input SMILES in chunks, aggregates counts per chunk, and updates
the output DB on the fly. A processed-chunks file can be used to resume.
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


def _validate_set_name(set_name):
    if not re.match(r'^[A-Za-z_][A-Za-z0-9_]*$', set_name):
        raise ValueError(
            "set_name must be a valid SQLite identifier (letters, numbers, underscores; cannot start with a number)"
        )
    if set_name in ('env_id', 'core_smi_id'):
        raise ValueError("set_name cannot be env_id or core_smi_id")
    return set_name


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

    if mode == 1 or mode == 2:
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


def _init_worker(radii, keep_stereo, max_heavy_atoms, mode, sep):
    global _RADII
    global _KEEP_STEREO
    global _MAX_HEAVY_ATOMS
    global _MODE
    global _SEP
    _RADII = radii
    _KEEP_STEREO = keep_stereo
    _MAX_HEAVY_ATOMS = max_heavy_atoms
    _MODE = mode
    _SEP = sep
    RDLogger.DisableLog('rdApp.warning')


def _process_chunk(task):
    chunk_id, lines = task
    envs = set()
    core_info = {}
    counts = {r: defaultdict(int) for r in _RADII}
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
                    counts[radius][(env, core_smi)] += 1
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


def _ensure_schema(conn, radii, set_name):
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
                {set_name} INTEGER NOT NULL DEFAULT 0,
                FOREIGN KEY (env_id) REFERENCES envs(env_id),
                FOREIGN KEY (core_smi_id) REFERENCES frags(core_smi_id),
                UNIQUE (env_id, core_smi_id)
            )
        """)
        cols = {row[1] for row in conn.execute(f"PRAGMA table_info(radius{radius})")}
        if set_name not in cols:
            conn.execute(f"ALTER TABLE radius{radius} ADD COLUMN {set_name} INTEGER NOT NULL DEFAULT 0")

    conn.commit()


def _batched(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def _ensure_env_ids(conn, envs, cache, batch_size):
    new_envs = [env for env in envs if env not in cache]
    if new_envs:
        conn.executemany("INSERT OR IGNORE INTO envs (env) VALUES (?)", [(e,) for e in new_envs])
        for batch in _batched(new_envs, batch_size):
            placeholders = ",".join("?" * len(batch))
            for env_id, env in conn.execute(
                f"SELECT env_id, env FROM envs WHERE env IN ({placeholders})",
                batch,
            ):
                cache[env] = env_id


def _ensure_inchi_ids(conn, core_info, cache, batch_size):
    new_inchis = []
    for core_smi, (_, _, core_smi_h, inchi_val) in core_info.items():
        if inchi_val and inchi_val not in cache:
            new_inchis.append((core_smi_h, inchi_val))

    if new_inchis:
        conn.executemany(
            "INSERT OR IGNORE INTO frags_h (smi, inchi) VALUES (?, ?)",
            new_inchis,
        )
        inchi_vals = [i for _, i in new_inchis]
        for batch in _batched(inchi_vals, batch_size):
            placeholders = ",".join("?" * len(batch))
            for core_smi_h_id, inchi_val in conn.execute(
                f"SELECT core_smi_h_id, inchi FROM frags_h WHERE inchi IN ({placeholders})",
                batch,
            ):
                cache[inchi_val] = core_smi_h_id


def _ensure_core_ids(conn, core_info, core_cache, inchi_cache, batch_size):
    new_cores = []
    for core_smi, (core_num_atoms, dist2, _, inchi_val) in core_info.items():
        if core_smi in core_cache:
            continue
        if not inchi_val:
            continue
        core_smi_h_id = inchi_cache.get(inchi_val)
        if core_smi_h_id is None:
            continue
        new_cores.append((core_smi, core_num_atoms, dist2, core_smi_h_id))

    if new_cores:
        conn.executemany(
            "INSERT OR IGNORE INTO frags (core_smi, core_num_atoms, dist2, core_smi_h_id) VALUES (?, ?, ?, ?)",
            new_cores,
        )
        core_vals = [c for c, _, _, _ in new_cores]
        for batch in _batched(core_vals, batch_size):
            placeholders = ",".join("?" * len(batch))
            for core_smi_id, core_smi in conn.execute(
                f"SELECT core_smi_id, core_smi FROM frags WHERE core_smi IN ({placeholders})",
                batch,
            ):
                core_cache[core_smi] = core_smi_id


def _update_radius_counts(conn, counts, env_cache, core_cache, set_name):
    for radius, mapping in counts.items():
        rows = []
        for (env, core_smi), cnt in mapping.items():
            env_id = env_cache.get(env)
            core_id = core_cache.get(core_smi)
            if env_id is None or core_id is None:
                continue
            rows.append((env_id, core_id, cnt))
        if rows:
            conn.executemany(
                f"""
                INSERT INTO radius{radius} (env_id, core_smi_id, {set_name})
                VALUES (?, ?, ?)
                ON CONFLICT(env_id, core_smi_id)
                DO UPDATE SET {set_name} = {set_name} + excluded.{set_name}
                """,
                rows,
            )


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
    parser.add_argument("set_name", help="Set name (column) to accumulate frequencies")
    parser.add_argument("--radii", default="1,2,3,4,5", help="Comma-separated radii (default: 1,2,3,4,5)")
    parser.add_argument("--chunk-size", type=int, default=10000, help="Lines per chunk (default: 10000)")
    parser.add_argument("--ncpu", type=int, default=1, help="Number of worker processes (default: 1)")
    parser.add_argument("--max-heavy-atoms", type=int, default=15, help="Max heavy atoms in core (default: 15)")
    parser.add_argument("--keep-stereo", action="store_true", help="Keep stereo in env/core")
    parser.add_argument("--mode", type=int, choices=[0, 1, 2], default=1, help="Fragmentation mode (default: 1)")
    parser.add_argument("--sep", default=None, help="Input delimiter (default: whitespace)")
    parser.add_argument("--processed-chunks", default=None, help="Path to processed chunks file (append)")
    parser.add_argument("--exclude-chunks", default=None, help="Path to chunk list to skip")
    parser.add_argument("--zstd", action="store_true", help="Force zstd input")
    parser.add_argument("--log-every", type=int, default=1, help="Log progress every N chunks (default: 1)")
    parser.add_argument("--batch-size", type=int, default=500, help="DB lookup batch size (default: 500)")
    args = parser.parse_args()

    set_name = _validate_set_name(args.set_name)
    radii = [int(r.strip()) for r in args.radii.split(",") if r.strip()]
    if not radii:
        raise ValueError("At least one radius must be specified")

    exclude_chunks = _read_chunk_ids(args.exclude_chunks)
    processed_chunks = _read_chunk_ids(args.processed_chunks)
    skip_chunks = exclude_chunks | processed_chunks

    RDLogger.DisableLog("rdApp.warning")

    conn = sqlite3.connect(args.output_db)
    _ensure_schema(conn, radii, set_name)
    cur = conn.cursor()

    env_cache = {}
    inchi_cache = {}
    core_cache = {}

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
        initargs=(radii, args.keep_stereo, args.max_heavy_atoms, args.mode, args.sep),
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
            conn.execute("BEGIN")
            _ensure_env_ids(conn, envs, env_cache, args.batch_size)
            _ensure_inchi_ids(conn, core_info, inchi_cache, args.batch_size)
            _ensure_core_ids(conn, core_info, core_cache, inchi_cache, args.batch_size)
            _update_radius_counts(conn, counts, env_cache, core_cache, set_name)
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
