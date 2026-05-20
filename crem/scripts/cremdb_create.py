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
import signal
import sqlite3
import subprocess
import sys
import time
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, islice, permutations
from multiprocessing import Pool, cpu_count
from operator import itemgetter
from pathlib import Path
from typing import List

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMMPA
from tqdm import tqdm

from crem.mol_context import get_std_context_core_permutations


_SQLITE_BATCH = 32000

# Magic value written into PRAGMA application_id at the end of a stride-mode
# shard build. The parallel-shards orchestrator reads this to decide which
# children to skip on resume.
_STRIDE_SHARD_SENTINEL = 0xC0DECAFE

# How many MMPA failures each worker is allowed to log to stderr before going
# silent. Bad-molecule edge cases in rdMMPA.FragmentMol can come in clusters
# (millions of structures), so we print the first few smi_ids for traceability
# then rely on the aggregate counter in the end-of-run summary.
_MMPA_LOG_BUDGET = 10
_mmpa_log_remaining = _MMPA_LOG_BUDGET


def _safe_mmpa_fragment(mol, smi_id, **kwargs):
    """Call rdMMPA.FragmentMol; swallow the two RDKit edge-case exceptions.

    Returns ``(frags, n_failures)``. The fragmenter is observed to occasionally
    raise ``IndexError: map::at`` or ``RuntimeError: Pre-condition Violation
    (endIdx not connected to end atom of bond)`` deep in
    ``MMPA::detail::addBondsFromTemplate``; both kill the whole worker if
    propagated. We log the first few offending smi_ids per worker, count the
    failure, and let the build continue with whatever fragments were produced
    by the other passes.
    """
    global _mmpa_log_remaining
    try:
        return rdMMPA.FragmentMol(mol, **kwargs), 0
    except (RuntimeError, IndexError, KeyError) as exc:
        if _mmpa_log_remaining > 0:
            sys.stderr.write(
                f"[mmpa-fail] {type(exc).__name__}: smi_id={smi_id or '<no_id>'}\n"
            )
            sys.stderr.flush()
            _mmpa_log_remaining -= 1
        return [], 1

def create_indices(conn: sqlite3.Connection, radii: List[int], verbose: bool = True):
    """
    Create optimized indices on the new database.

    Args:
        conn: SQLite connection
        radii: List of radius values
        verbose: Print progress

    Notes:
        Redundant indices are NOT created — the UNIQUE constraints on
        envs(env), frags(core_smi), frags_h(smi) and
        radius{N}(env_id, core_smi_id, is_ring_closure) already produce
        autoindices. One explicit covering index per radius table is added:
        (env_id, is_ring_closure, core_num_atoms, dist2). This serves the
        hot query path: every call filters by env_id and provenance, and most
        also by size/dist.
    """
    cur = conn.cursor()

    indices = []
    for radius in radii:
        indices.append((
            f"idx_radius{radius}_lookup",
            f"CREATE INDEX IF NOT EXISTS idx_radius{radius}_lookup "
            f"ON radius{radius}(env_id, is_ring_closure, core_num_atoms, dist2)",
        ))

    for idx_name, sql in tqdm(indices, desc="Creating indices", disable=not verbose):
        cur.execute(sql)

    conn.commit()

    # Analyze database for query optimization
    if verbose:
        print("Analyzing database...")
    cur.execute("ANALYZE")
    conn.commit()


def _shard_db_path(output_db: str, shard_idx: int) -> str:
    """Return file path for the given shard index.

    Shard 0 is the base output file; shards 1+ get a zero-padded numeric suffix.
    Example: output.db → output.db (0), output_001.db (1), output_002.db (2), ...
    """
    if shard_idx == 0:
        return output_db
    p = Path(output_db)
    return str(p.parent / f"{p.stem}_{shard_idx:03d}{p.suffix}")


def _find_last_shard_idx(output_db: str) -> int:
    """Return the index of the last existing shard file, or -1 if none exist."""
    if not os.path.exists(output_db):
        return -1
    last = 0
    idx = 1
    while os.path.exists(_shard_db_path(output_db, idx)):
        last = idx
        idx += 1
    return last


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

    existing_files = [v for v in values if os.path.isfile(v)]
    missing_files = [v for v in values if not os.path.isfile(v)]

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


def _build_membership(set_names, set_ids):
    """Fold ``set_ids`` into a single per-molecule lookup structure.

    Returns ``(membership, all_mols_mask)`` where:

      * ``membership`` is ``dict[smi_id, int]`` — bit ``i`` is set when
        ``smi_id`` belongs to ``set_names[i]``. Sets that have no ID file
        (membership is universal) are not represented in ``membership``;
        their bits live in ``all_mols_mask``. Returns ``None`` if there is
        no per-ID filtering at all (single set name without files).
      * ``all_mols_mask`` is the bitmask of sets that include every input
        molecule (``set_ids[name] is None``).

    Folding into one dict deduplicates IDs that appear in multiple sets,
    which materially shrinks the in-memory index when set files are nested
    subsets of each other.
    """
    if set_ids is None:
        return None, (1 << len(set_names)) - 1
    membership: dict = {}
    all_mols_mask = 0
    for i, name in enumerate(set_names):
        ids = set_ids.get(name)
        bit = 1 << i
        if ids is None:
            all_mols_mask |= bit
            continue
        for smi_id in ids:
            membership[smi_id] = membership.get(smi_id, 0) | bit
    return membership, all_mols_mask


# Workers are forked Python processes; each owns an independent cache. The
# same core_smi typically appears many times across chunks within a worker
# (it's the SMILES of a fragment core, of which there are far fewer than
# (env, core) pairs), so caching skips the parse+canonicalize round-trip on
# repeats. The cap bounds peak memory: ~32-byte keys/values × 200k ≈ a few MB
# per worker.
@lru_cache(maxsize=200_000)
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


def _fragment_mol_acyclic(smi, smi_id, mode):
    """Acyclic-bond MMPA fragmentation (the historical CReM fragmenter).

    Returns ``(outlines_set, n_mmpa_failures)``. Each MMPA call is wrapped in
    ``_safe_mmpa_fragment`` so that an RDKit edge-case exception only drops
    the offending pass for the offending molecule rather than killing the
    whole worker.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return set(), 0

    outlines = set()
    n_failures = 0

    if mode == 0 or mode == 1:
        frags4, f1 = _safe_mmpa_fragment(
            mol, smi_id,
            pattern="[!#1]!@!=!#[!#1]",
            maxCuts=4,
            resultsAsMols=False,
            maxCutBonds=30,
        )
        frags3, f2 = _safe_mmpa_fragment(
            mol, smi_id,
            pattern="[!#1]!@!=!#[!#1]",
            maxCuts=3,
            resultsAsMols=False,
            maxCutBonds=30,
        )
        n_failures += f1 + f2
        for core, chains in set(list(frags4) + list(frags3)):
            outlines.add((core, chains))

    if mode == 0 or mode == 2:
        mol = Chem.AddHs(mol)
        n = mol.GetNumAtoms() - mol.GetNumHeavyAtoms()
        if n < 60:  # TODO: remove this limit, it is not very reasonable
            frags, f3 = _safe_mmpa_fragment(
                mol, smi_id,
                pattern="[#1]!@!=!#[!#1]",
                maxCuts=1,
                resultsAsMols=False,
                maxCutBonds=100,   # TODO: why we need this?
            )
            n_failures += f3
            for core, chains in frags:
                outlines.add((core, chains))

    return outlines, n_failures


def _fragment_mol_ring(smi, smi_id):
    """Ring-bond fragmentation: cut every pair of SINGLE bonds that lies
    inside the same ring, yielding 2-AP arc fragments.

    Aromatic / double / triple ring bonds are NOT cut. Cuts that don't
    disconnect the graph (e.g. spanning two rings of a fused system without
    severing both paths) are skipped — only pairs producing exactly two
    fragments, each with exactly two * atoms, are kept.

    Both arcs are emitted as candidate cores: (arc_a, arc_b) and (arc_b, arc_a).
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return set(), 0

    try:
        bond_rings = mol.GetRingInfo().BondRings()
    except (RuntimeError, IndexError):
        # Defensive: RDKit corner cases have been observed to throw here.
        # Treat as "no rings to cut" and continue.
        return set(), 0
    if not bond_rings:
        return set(), 0

    outlines = set()
    seen_pairs = set()
    for ring_bond_ids in bond_rings:
        for b1, b2 in combinations(sorted(ring_bond_ids), 2):
            if (b1, b2) in seen_pairs:
                continue
            seen_pairs.add((b1, b2))
            if mol.GetBondWithIdx(b1).GetBondType() != Chem.BondType.SINGLE:
                continue
            if mol.GetBondWithIdx(b2).GetBondType() != Chem.BondType.SINGLE:
                continue
            try:
                cut = Chem.FragmentOnBonds(mol, [b1, b2], dummyLabels=[(1, 1), (2, 2)])
                pieces = Chem.GetMolFrags(cut, asMols=True)
            except Exception:
                continue
            if len(pieces) != 2:
                continue  # didn't disconnect (fused ring system)
            valid = True
            for piece in pieces:
                stars = [a for a in piece.GetAtoms() if a.GetAtomicNum() == 0]
                if len(stars) != 2:
                    valid = False
                    break
                # FragmentOnBonds stores cut labels as isotopes; convert to
                # atom map numbers so SMILES carry [*:N] as elsewhere in CReM.
                for a in stars:
                    iso = a.GetIsotope()
                    if iso:
                        a.SetIsotope(0)
                        a.SetAtomMapNum(iso)
            if not valid:
                continue
            try:
                smi_a = Chem.MolToSmiles(pieces[0])
                smi_b = Chem.MolToSmiles(pieces[1])
            except Exception:
                continue
            outlines.add((smi_a, smi_b))
            outlines.add((smi_b, smi_a))

    return outlines, 0


def _fragment_mol(smi, smi_id, mode, frag_mode):
    """Drive acyclic and/or ring-bond fragmentation per `frag_mode`.

    Returns ``(out, n_mmpa_failures)`` where ``out`` is a set of
    ``(smi, smi_id, core, chains, is_ring_closure)`` tuples and
    ``n_mmpa_failures`` is how many RDKit fragmentation calls were skipped
    due to edge-case exceptions on this molecule.
    """
    out = set()
    n_failures = 0
    if frag_mode in ('acyclic', 'both'):
        outlines, f = _fragment_mol_acyclic(smi, smi_id, mode)
        n_failures += f
        for core, chains in outlines:
            out.add((smi, smi_id, core, chains, 0))
    if frag_mode in ('ring', 'both'):
        outlines, f = _fragment_mol_ring(smi, smi_id)
        n_failures += f
        for core, chains in outlines:
            out.add((smi, smi_id, core, chains, 1))
    return out, n_failures


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


# `_env_core_from_fragment` is called once per (core, context) MMPA fragment
# at every radius — typically with heavy repetition of the SMILES strings
# across chunks within one worker. The Chem.MolFromSmiles(..., sanitize=False)
# parse used to be done solely to read GetNumHeavyAtoms(); caching collapses
# the repeated parses to a single dict lookup.
@lru_cache(maxsize=200_000)
def _count_heavy_atoms(smi):
    mm = Chem.MolFromSmiles(smi, sanitize=False)
    return mm.GetNumHeavyAtoms() if mm else float('inf')


def _env_core_from_fragment(core, context, radius, keep_stereo, max_heavy_atoms):
    output = []

    if not core:
        residues = context.split('.')
        if len(residues) == 2:
            for ctx, c in permutations(residues, 2):
                if ctx == '[H][*:1]':
                    continue
                num_heavy_atoms = _count_heavy_atoms(c)
                if num_heavy_atoms <= max_heavy_atoms:
                    env, cores = get_std_context_core_permutations(ctx, c, radius, keep_stereo)
                    if env and cores:
                        output.append((env, cores[0], num_heavy_atoms))  # only one item in cores
    else:
        num_heavy_atoms = _count_heavy_atoms(core)
        if num_heavy_atoms <= max_heavy_atoms:
            env, cores = get_std_context_core_permutations(context, core, radius, keep_stereo)
            if env and cores:
                for core_smi in cores:
                    output.append((env, core_smi, num_heavy_atoms))

    return output


def _init_worker(radii, keep_stereo, max_heavy_atoms, mode, set_names, frag_mode):
    """Initialise fragmentation workers.

    The set-membership index is intentionally NOT shipped here — set decisions
    are made in the orchestrating main process so workers don't fork-inherit
    huge per-file ID sets that Python's refcount-on-access would COW-duplicate
    into every worker.
    """
    global _RADII
    global _KEEP_STEREO
    global _MAX_HEAVY_ATOMS
    global _MODE
    global _SET_NAMES
    global _FRAG_MODE
    _RADII = radii
    _KEEP_STEREO = keep_stereo
    _MAX_HEAVY_ATOMS = max_heavy_atoms
    _MODE = mode
    _SET_NAMES = set_names
    _FRAG_MODE = frag_mode
    RDLogger.DisableLog('rdApp.warning')


def _process_chunk(task):
    """Fragment one chunk of pre-parsed, pre-filtered molecules.

    The main process performs line parsing and set-membership filtering and
    hands workers ``(smi, smi_id, member_sets)`` triples. Workers therefore
    never see the (potentially huge) per-set ID index.
    """
    chunk_id, items = task
    envs = set()
    core_info = {}
    counts = {name: {r: defaultdict(int) for r in _RADII} for name in _SET_NAMES}
    stats = {
        "lines": 0,
        "fragments": 0,
        "pairs": 0,
        "skipped_mmpa": 0,
    }

    for smi, smi_id, member_sets in items:
        stats["lines"] += 1
        frags, n_mmpa_fail = _fragment_mol(smi, smi_id, _MODE, _FRAG_MODE)
        stats["fragments"] += len(frags)
        stats["skipped_mmpa"] += n_mmpa_fail
        for _, _, core, context, is_ring_closure in frags:
            for radius in _RADII:
                for env, core_smi, core_num_atoms in _env_core_from_fragment(
                    core,
                    context,
                    radius,
                    _KEEP_STEREO,
                    _MAX_HEAVY_ATOMS,
                ):
                    for set_name in member_sets:
                        counts[set_name][radius][(env, core_smi, is_ring_closure)] += 1
                    envs.add(env)
                    if core_smi not in core_info:
                        dist2 = _core_dist2(core_smi)
                        core_smi_h = _replace_attachment_points_with_h(core_smi)
                        core_info[core_smi] = (core_num_atoms, dist2, core_smi_h)
                    stats["pairs"] += 1

    return chunk_id, envs, core_info, counts, stats


def _read_chunk_ids(path):
    if not path or not os.path.exists(path):
        return set()
    ids = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                try:
                    ids.append(int(line))
                except ValueError:
                    continue
    return set(ids)


def _ensure_schema(conn, radii, set_names):
    conn.execute("PRAGMA page_size = 16384")     # larger pages for better B-tree packing (new DBs only)
    conn.execute("PRAGMA user_version = 1")
    conn.execute("PRAGMA foreign_keys = ON")
    # Bulk-load journaling: OFF/OFF is materially faster than WAL/NORMAL when
    # there is a single writer and recovery is acceptable (a crash mid-build
    # is recovered by re-running with --processed-chunks, or by re-creating
    # the shard from scratch). _finalize_shard restores WAL/NORMAL before
    # the shard is released.
    # When the DB is already in WAL mode (e.g. updating an existing DB) and
    # another connection is alive, switching to OFF requires exclusive access
    # and will raise OperationalError — fall back to WAL/NORMAL gracefully.
    try:
        conn.execute("PRAGMA journal_mode = OFF")
        conn.execute("PRAGMA synchronous = OFF")
    except Exception:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("PRAGMA cache_size = -65536")   # 64 MB page cache

    conn.execute("""
        CREATE TABLE IF NOT EXISTS envs(
            env_id INTEGER PRIMARY KEY AUTOINCREMENT,
            env TEXT NOT NULL UNIQUE
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS frags_h(
            core_smi_h_id INTEGER PRIMARY KEY AUTOINCREMENT,
            smi TEXT NOT NULL UNIQUE
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS frags(
            core_smi_id INTEGER PRIMARY KEY AUTOINCREMENT,
            core_smi TEXT NOT NULL UNIQUE,
            core_smi_h_id INTEGER NOT NULL,
            FOREIGN KEY (core_smi_h_id) REFERENCES frags_h(core_smi_h_id)
        )
    """)

    for radius in radii:
        conn.execute(f"""
            CREATE TABLE IF NOT EXISTS radius{radius}(
                env_id INTEGER NOT NULL,
                core_smi_id INTEGER NOT NULL,
                core_num_atoms INTEGER NOT NULL,
                dist2 INTEGER NOT NULL,
                is_ring_closure INTEGER NOT NULL DEFAULT 0,
                FOREIGN KEY (env_id) REFERENCES envs(env_id),
                FOREIGN KEY (core_smi_id) REFERENCES frags(core_smi_id),
                UNIQUE (env_id, core_smi_id, is_ring_closure)
            )
        """)
        cols = {row[1] for row in conn.execute(f"PRAGMA table_info(radius{radius})")}
        for set_name in set_names:
            if set_name not in cols:
                conn.execute(f"ALTER TABLE radius{radius} ADD COLUMN {set_name} INTEGER NOT NULL DEFAULT 0")
        # Drop query indices on radius tables before bulk loading; they will be
        # recreated by create_indices() at the end. The UNIQUE autoindex is kept
        # because it is required for ON CONFLICT DO UPDATE upserts.
        for suffix in ("env_id", "core_smi_id", "both", "lookup"):
            conn.execute(f"DROP INDEX IF EXISTS idx_radius{radius}_{suffix}")

    conn.commit()


def _merge_core_info(acc, new):
    for k, v in new.items():
        if k not in acc:
            acc[k] = v


def _merge_counts(acc, new, set_names, radii):
    for set_name in set_names:
        new_per_set = new.get(set_name, {})
        acc_per_set = acc[set_name]
        for radius in radii:
            src = new_per_set.get(radius, {})
            dst = acc_per_set[radius]
            for key, cnt in src.items():
                dst[key] += cnt


def _flush_to_db(conn, envs, core_info, counts, set_names, radii, timings=None):
    """Flush accumulated data to DB using per-flush local ID maps.

    For every flush we resolve (env, core_smi_h, core_smi) → row-id mappings
    just for the rows we are about to touch. ``INSERT OR IGNORE … RETURNING``
    returns IDs for newly-inserted rows; a fallback ``SELECT … WHERE x IN
    (...)`` covers rows that were IGNORED because they already exist on disk.

    The maps go out of scope at the end of the call, so there is no lifetime
    state mirroring the DB in Python.
    """

    env_ids: dict = {}        # env_str        -> env_id
    smi_h_ids: dict = {}      # H-canonical SMI -> core_smi_h_id
    core_smi_ids: dict = {}   # core_smi       -> core_smi_id

    # Step 1: resolve env_ids for every env touched in this flush.
    t0 = time.perf_counter()
    env_list = list(envs)
    for i in range(0, len(env_list), _SQLITE_BATCH):
        batch = env_list[i:i + _SQLITE_BATCH]
        ph = ",".join(["(?)"] * len(batch))
        inserted = set()
        for env_str, env_id in conn.execute(
            f"INSERT OR IGNORE INTO envs (env) VALUES {ph} "
            f"RETURNING env, env_id",
            batch,
        ):
            env_ids[env_str] = env_id
            inserted.add(env_str)
        missing = [e for e in batch if e not in inserted]
        if missing:
            ph_sel = ",".join(["?"] * len(missing))
            for env_str, env_id in conn.execute(
                f"SELECT env, env_id FROM envs WHERE env IN ({ph_sel})",
                missing,
            ):
                env_ids[env_str] = env_id
    if timings is not None:
        timings["envs"] = time.perf_counter() - t0

    # Step 2: resolve core_smi_h_ids for every H-canonical SMILES touched.
    t0 = time.perf_counter()
    smi_h_list: list = []
    seen_smi_h: set = set()
    for core_smi, (_core_num_atoms, _dist2, core_smi_h) in core_info.items():
        if core_smi_h and core_smi_h not in seen_smi_h:
            smi_h_list.append(core_smi_h)
            seen_smi_h.add(core_smi_h)

    for i in range(0, len(smi_h_list), _SQLITE_BATCH):
        batch = smi_h_list[i:i + _SQLITE_BATCH]
        ph = ",".join(["(?)"] * len(batch))
        inserted = set()
        for smi_val, h_id in conn.execute(
            f"INSERT OR IGNORE INTO frags_h (smi) VALUES {ph} "
            f"RETURNING smi, core_smi_h_id",
            batch,
        ):
            smi_h_ids[smi_val] = h_id
            inserted.add(smi_val)
        missing = [s for s in batch if s not in inserted]
        if missing:
            ph_sel = ",".join(["?"] * len(missing))
            for smi_val, h_id in conn.execute(
                f"SELECT smi, core_smi_h_id FROM frags_h WHERE smi IN ({ph_sel})",
                missing,
            ):
                smi_h_ids[smi_val] = h_id
    if timings is not None:
        timings["smi_h"] = time.perf_counter() - t0

    # Step 3: resolve core_smi_ids for every core touched.
    t0 = time.perf_counter()
    core_inserts: list = []  # (core_smi, core_smi_h_id) — INSERT payload
    core_smis_touched: list = []  # core_smis we need to resolve, regardless of new/old
    for core_smi, (_core_num_atoms, _dist2, core_smi_h) in core_info.items():
        if not core_smi_h or core_smi_h not in smi_h_ids:
            continue
        core_smis_touched.append(core_smi)
        core_inserts.append((core_smi, smi_h_ids[core_smi_h]))

    # 2 columns per row → at most _SQLITE_BATCH/2 rows per chunk to stay
    # under SQLITE_MAX_VARIABLE_NUMBER.
    rows_per_chunk = max(1, _SQLITE_BATCH // 2)
    for i in range(0, len(core_inserts), rows_per_chunk):
        batch = core_inserts[i:i + rows_per_chunk]
        ph = ",".join(["(?,?)"] * len(batch))
        flat = [v for row in batch for v in row]
        inserted = set()
        for core_smi, core_smi_id in conn.execute(
            f"INSERT OR IGNORE INTO frags "
            f"(core_smi, core_smi_h_id) "
            f"VALUES {ph} "
            f"RETURNING core_smi, core_smi_id",
            flat,
        ):
            core_smi_ids[core_smi] = core_smi_id
            inserted.add(core_smi)
        missing = [r[0] for r in batch if r[0] not in inserted]
        if missing:
            ph_sel = ",".join(["?"] * len(missing))
            for core_smi, core_smi_id in conn.execute(
                f"SELECT core_smi, core_smi_id FROM frags WHERE core_smi IN ({ph_sel})",
                missing,
            ):
                core_smi_ids[core_smi] = core_smi_id
    if timings is not None:
        timings["cores"] = time.perf_counter() - t0

    # Step 4: upsert radius counts using resolved IDs (no JOINs).
    # core_num_atoms and dist2 are denormalized into radius{N} so the hot
    # query path can filter directly on the radius table without joining
    # frags. They are written on first INSERT and not touched on conflict.
    # The (env, core, is_ring_closure) triple is the conflict key — the same
    # (env, core) pair can carry both an acyclic-cut row (is_ring_closure=0)
    # and a ring-cut row (is_ring_closure=1) with independent per-set counts.
    t0 = time.perf_counter()
    for set_name in set_names:
        per_set = counts.get(set_name, {})
        for radius in radii:
            mapping = per_set.get(radius, {})
            if not mapping:
                continue
            rows = []
            for (env, core_smi, is_ring_closure), cnt in mapping.items():
                env_id = env_ids.get(env)
                core_smi_id = core_smi_ids.get(core_smi)
                if env_id is None or core_smi_id is None:
                    continue
                core_entry = core_info.get(core_smi)
                if core_entry is None:
                    continue
                core_num_atoms, dist2, _ = core_entry
                rows.append((env_id, core_smi_id, core_num_atoms, dist2,
                             is_ring_closure, cnt))
            # Sort only by (env_id, core_smi_id, is_ring_closure) — the
            # conflict key — to drive sequential B-tree writes; the trailing
            # columns are irrelevant for write order and full-tuple comparison
            # would cost noticeably more per call.
            rows.sort(key=itemgetter(0, 1, 4))
            if rows:
                conn.executemany(
                    f"INSERT INTO radius{radius} "
                    f"(env_id, core_smi_id, core_num_atoms, dist2, "
                    f"is_ring_closure, {set_name}) "
                    f"VALUES (?, ?, ?, ?, ?, ?) "
                    f"ON CONFLICT(env_id, core_smi_id, is_ring_closure) "
                    f"DO UPDATE SET {set_name} = {set_name} + excluded.{set_name}",
                    rows,
                )
    if timings is not None:
        timings["radius"] = time.perf_counter() - t0


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


def run(
    input_path: str,
    output_db: str,
    set_name: list,
    radii=(1, 2, 3, 4, 5),
    chunk_size: int = 100,
    max_heavy_atoms: int = 15,
    keep_stereo: bool = False,
    mode: int = 0,
    sep=None,
    processed_chunks=None,
    force_zstd: bool = False,
    log_every=None,
    ncpu: int = 1,
    flush_every: int = 100,
    prefetch: int = 4,
    timings: bool = False,
    shard_size=None,
    verbose: bool = True,
    frag_mode: str = 'both',
    stride_mod: int = 1,
    stride_idx: int = 0,
) -> None:
    """Build or extend a v1 CReM fragment database.

    Args:
        input_path: Path to input SMILES file (plain text or .zst).
        output_db: Path to output SQLite DB (created if absent, extended if present).
        set_name: CLI-style list: ``["name"]`` for a single set, or
            ``["name1", "ids1.txt", "name2", "ids2.txt"]`` for multiple sets.
        radii: Fragment radii to compute (default 1–5).
        chunk_size: Input lines per worker task.
        max_heavy_atoms: Maximum heavy atoms in a core fragment.
        keep_stereo: Preserve stereocentres in env/core SMILES.
        mode: Fragmentation mode — 0 all atoms, 1 heavy only, 2 H only.
        sep: Column separator in the input file (None = whitespace).
        processed_chunks: Path to a file tracking completed chunk IDs for resumable runs.
        force_zstd: Force zstd decompression regardless of file extension.
        log_every: Print a progress line every N chunks (None = silent).
        ncpu: Worker processes (capped at available CPUs).
        flush_every: Chunks to accumulate before each DB flush.
        prefetch: In-flight task batches per worker.
        timings: Print per-flush timing breakdown to stderr.
        shard_size: Max input structures per shard DB (None = single DB).
        frag_mode: Fragmentation source — 'acyclic' (today's MMPA cuts of
            acyclic bonds), 'ring' (cuts of pairs of SINGLE bonds inside the
            same ring; used by make_cycle(ring_closures=True)), or 'both'.
        stride_mod: When > 1, only process chunks whose ``chunk_id %
            stride_mod == stride_idx``. Used by ``run_parallel_shards`` to
            split the input across N concurrent shard builders.
        stride_idx: Stride index for this build; meaningful only when
            ``stride_mod > 1``.
    """
    if frag_mode not in ('acyclic', 'ring', 'both'):
        raise ValueError(
            f"frag_mode must be one of 'acyclic', 'ring', 'both' (got {frag_mode!r})"
        )
    if stride_mod < 1:
        raise ValueError("stride_mod must be >= 1")
    if not (0 <= stride_idx < stride_mod):
        raise ValueError("stride_idx must be in [0, stride_mod)")
    set_names, set_ids = _resolve_set_names(set_name)
    radii = sorted(set(radii))
    if not radii:
        raise ValueError("At least one radius must be specified")

    # Build the per-molecule membership index once, in the main process.
    # Workers never see set_ids — see _init_worker's docstring for the why.
    membership, all_mols_mask = _build_membership(set_names, set_ids)
    # The full-membership case (all bits set, no per-ID filtering) lets us
    # ship a single shared tuple as `member_sets` for every molecule.
    full_member_sets = tuple(set_names) if membership is None else None
    set_names_tuple = tuple(set_names)

    skip_chunks = _read_chunk_ids(processed_chunks)

    RDLogger.DisableLog("rdApp.warning")

    # Determine which shard to start from.
    # In shard mode on resume: skip the last (possibly partial) shard and open
    # a fresh one.  In non-shard mode always use index 0 (the output file).
    if shard_size is not None:
        last_shard_idx = _find_last_shard_idx(output_db)
        current_shard_idx = last_shard_idx + 1
    else:
        current_shard_idx = 0

    def _open_shard(idx: int) -> sqlite3.Connection:
        path = _shard_db_path(output_db, idx)
        c = sqlite3.connect(path)
        _ensure_schema(c, radii, set_names)
        return c

    def _finalize_shard(c: sqlite3.Connection) -> None:
        """Restore WAL/NORMAL and build query indices, making the shard a
        valid standalone DB. During the build it ran with journal_mode=OFF
        (see _ensure_schema)."""
        c.execute("PRAGMA synchronous = NORMAL")
        c.execute("PRAGMA journal_mode = WAL")
        create_indices(c, radii, verbose)
        # In stride mode, mark the shard as fully built so the orchestrator
        # can skip re-launching this child on resume.
        if stride_mod > 1:
            c.execute(f"PRAGMA application_id = {_STRIDE_SHARD_SENTINEL}")
            c.commit()

    conn = _open_shard(current_shard_idx)

    # No lifetime Python-side ID caches: _flush_to_db builds per-flush local
    # maps via INSERT OR IGNORE RETURNING with a SELECT fallback for rows that
    # already exist on disk. This keeps memory bounded by per-flush size
    # instead of growing with the cumulative DB size.

    processed_handle = None
    if processed_chunks:
        processed_handle = open(processed_chunks, "a", encoding="utf-8")

    input_handle, zstd_handle = _open_input(input_path, force_zstd)

    total_stats = {"lines": 0, "fragments": 0, "pairs": 0, "skipped_mmpa": 0}
    chunks_processed = 0
    chunks_skipped = 0
    start_time = time.time()

    _ncpu = min(cpu_count(), max(ncpu, 1))
    batch_size = _ncpu * max(prefetch, 1)

    pool = Pool(
        _ncpu,
        initializer=_init_worker,
        initargs=(
            radii,
            keep_stereo,
            max_heavy_atoms,
            mode,
            set_names,
            frag_mode,
        ),
    )

    # Accumulators for flush batching
    acc_envs = set()
    acc_core_info = {}
    acc_counts = {name: {r: defaultdict(int) for r in radii} for name in set_names}
    acc_chunk_ids = []
    flush_counter = 0
    structures_in_current_shard = 0  # reset on each shard rotation
    sources_to_merge: list = []       # populated at merge time; kept for stats

    # _do_flush and _rotate_shard access `conn` via late binding: rebinding
    # `conn` in this scope is picked up by subsequent _do_flush calls.

    def _do_flush():
        nonlocal chunks_processed, flush_counter
        _timings = {} if timings else None

        conn.execute("PRAGMA foreign_keys = OFF")  # IDs are guaranteed to exist
        t0 = time.perf_counter()
        conn.execute("BEGIN")
        if _timings is not None:
            _timings["begin"] = time.perf_counter() - t0

        _flush_to_db(
            conn, acc_envs, acc_core_info, acc_counts,
            set_names, radii,
            timings=_timings,
        )

        t0 = time.perf_counter()
        conn.commit()
        conn.execute("PRAGMA foreign_keys = ON")
        if _timings is not None:
            _timings["commit"] = time.perf_counter() - t0

        chunks_processed += len(acc_chunk_ids)
        if processed_handle:
            for cid in acc_chunk_ids:
                processed_handle.write(f"{cid}\n")
            processed_handle.flush()

        flush_counter += 1
        # No WAL checkpoint here: build runs with journal_mode=OFF.

        if _timings is not None:
            parts = "  ".join(f"{k}={v*1000:.1f}ms" for k, v in _timings.items())
            total = sum(_timings.values())
            sys.stderr.write(f"[flush #{chunks_processed}] {parts}  total={total*1000:.1f}ms\n")
            sys.stderr.flush()

    def _rotate_shard():
        nonlocal conn, current_shard_idx, flush_counter, structures_in_current_shard
        _finalize_shard(conn)
        conn.close()
        current_shard_idx += 1
        conn = _open_shard(current_shard_idx)
        flush_counter = 0
        structures_in_current_shard = 0
        sys.stderr.write(
            f"\n  Shard full → opened {_shard_db_path(output_db, current_shard_idx)}\n"
        )
        sys.stderr.flush()

    try:
        def task_iter():
            """Yield (chunk_id, items) where items are pre-parsed and
            pre-filtered (smi, smi_id, member_sets) triples.

            All line parsing and set-membership decisions happen here in the
            main process, so workers don't need the (potentially huge) ID
            index. Chunks whose every molecule is filtered out (no member
            sets) are dropped entirely — the pool never sees them.
            """
            nonlocal chunks_skipped
            for chunk_id, lines in enumerate(_iter_chunks(input_handle, chunk_size)):
                if stride_mod > 1 and chunk_id % stride_mod != stride_idx:
                    continue
                if chunk_id in skip_chunks:
                    chunks_skipped += 1
                    continue
                items = []
                for line in lines:
                    smi, smi_id = _parse_line(line, sep)
                    if not smi:
                        continue
                    if membership is None:
                        # No per-ID filtering — every molecule is in every set.
                        items.append((smi, smi_id, full_member_sets))
                        continue
                    mask = membership.get(smi_id, 0) | all_mols_mask
                    if mask == 0:
                        continue
                    member_sets = tuple(
                        name for i, name in enumerate(set_names_tuple)
                        if mask & (1 << i)
                    )
                    items.append((smi, smi_id, member_sets))
                # Yield even all-filtered chunks so chunk_id still flows
                # through the pool and gets persisted in processed_chunks;
                # without that, resume re-parses these chunks every time.
                yield (chunk_id, items)

        tasks = task_iter()

        while True:
            # Back-pressure: only submit batch_size tasks at a time
            current_batch = list(islice(tasks, batch_size))
            if not current_batch:
                break

            for chunk_id, envs, core_info, counts, stats in pool.imap_unordered(
                _process_chunk, current_batch, chunksize=1
            ):
                t0 = time.perf_counter()
                acc_envs.update(envs)
                _merge_core_info(acc_core_info, core_info)
                _merge_counts(acc_counts, counts, set_names, radii)
                acc_chunk_ids.append(chunk_id)
                if timings:
                    sys.stderr.write(f"[merge chunk {chunk_id}] {(time.perf_counter()-t0)*1000:.1f}ms\n")
                    sys.stderr.flush()

                total_stats["lines"] += stats["lines"]
                total_stats["fragments"] += stats["fragments"]
                total_stats["pairs"] += stats["pairs"]
                total_stats["skipped_mmpa"] += stats.get("skipped_mmpa", 0)
                structures_in_current_shard += stats["lines"]

                if len(acc_chunk_ids) >= flush_every:
                    _do_flush()

                    if log_every and chunks_processed % log_every == 0:
                        elapsed = time.time() - start_time
                        rate = total_stats["lines"] / elapsed if elapsed > 0 else 0
                        sys.stderr.write(
                            f"\rChunks: {chunks_processed} processed, {chunks_skipped} skipped | "
                            f"mols: {total_stats['lines']} | frags: {total_stats['fragments']} | "
                            f"pairs: {total_stats['pairs']} | {rate:.1f} mol/s"
                        )
                        sys.stderr.flush()

                    # Reset accumulators
                    acc_envs = set()
                    acc_core_info = {}
                    acc_counts = {name: {r: defaultdict(int) for r in radii} for name in set_names}
                    acc_chunk_ids = []

                    # Rotate shard if size threshold reached
                    if shard_size is not None and structures_in_current_shard >= shard_size:
                        _rotate_shard()

        # Final flush for any remaining accumulated chunks
        if acc_chunk_ids:
            _do_flush()

        # Finalize the last (or only) shard
        _finalize_shard(conn)
        conn.close()
        conn = None

        # In shard mode: collect ALL shard files from index 1 upward (including
        # those created in previous sessions) and merge them into shard 0.
        if shard_size is not None:
            sources_to_merge = []
            idx = 1
            while os.path.exists(_shard_db_path(output_db, idx)):
                sources_to_merge.append(_shard_db_path(output_db, idx))
                idx += 1

            if sources_to_merge:
                sys.stderr.write(
                    f"\nMerging {len(sources_to_merge) + 1} shards into {output_db}...\n"
                )
                sys.stderr.flush()
                from crem.scripts.cremdb_merge import merge_into  # lazy: avoids circular import
                with sqlite3.connect(output_db) as merge_conn:
                    merge_conn.execute("PRAGMA cache_size = -262144")  # 256 MB
                    merge_conn.execute("PRAGMA synchronous = NORMAL")
                    merge_conn.execute("PRAGMA journal_mode = WAL")
                    merge_conn.execute("PRAGMA temp_store = MEMORY")
                    merge_into(merge_conn, sources_to_merge, verbose=True)
                    create_indices(merge_conn, radii, True)

    finally:
        pool.close()
        pool.join()
        input_handle.close()
        if zstd_handle:
            zstd_handle.close()
        if processed_handle:
            processed_handle.close()
        if conn is not None:
            conn.close()

    if verbose:
        sys.stderr.write("\n")

    elapsed = time.time() - start_time

    if verbose:
        # Print statistics from the final merged (or single) database
        with sqlite3.connect(output_db) as stats_conn:
            cur = stats_conn.cursor()
            env_count = cur.execute("SELECT COUNT(*) FROM envs").fetchone()[0]
            frag_count = cur.execute("SELECT COUNT(*) FROM frags").fetchone()[0]
            frag_h_count = cur.execute("SELECT COUNT(*) FROM frags_h").fetchone()[0]
            radius_counts = {}
            for radius in radii:
                radius_counts[radius] = cur.execute(f"SELECT COUNT(rowid) FROM radius{radius}").fetchone()[0]

        print("Done.")
        print(f"Chunks processed: {chunks_processed}")
        print(f"Chunks skipped: {chunks_skipped}")
        print(f"Molecules processed: {total_stats['lines']}")
        print(f"Fragments generated: {total_stats['fragments']}")
        print(f"Env/core pairs: {total_stats['pairs']}")
        print(f"MMPA failures: {total_stats['skipped_mmpa']}")
        print(f"Unique envs: {env_count}")
        print(f"Unique frags: {frag_count}")
        print(f"Unique frags_h: {frag_h_count}")
        for radius, count in radius_counts.items():
            print(f"radius{radius} rows: {count}")
        print(f"Elapsed: {elapsed:.1f}s")
        if shard_size is not None:
            print(f"Shards on disk: {len(sources_to_merge) + 1}")


def _shard_is_finalised(shard_path: str) -> bool:
    """Return True if shard_path has the stride-mode finalisation sentinel."""
    if not os.path.exists(shard_path):
        return False
    try:
        with sqlite3.connect(shard_path) as c:
            app_id = c.execute("PRAGMA application_id").fetchone()[0]
    except sqlite3.DatabaseError:
        return False
    return int(app_id) == _STRIDE_SHARD_SENTINEL


def run_parallel_shards(
    input_path: str,
    output_db: str,
    set_name: list,
    parallel_shards: int,
    ncpu: int,
    radii=(1, 2, 3, 4, 5),
    chunk_size: int = 100,
    max_heavy_atoms: int = 15,
    keep_stereo: bool = False,
    mode: int = 0,
    sep=None,
    force_zstd: bool = False,
    log_every=None,
    flush_every: int = 100,
    prefetch: int = 4,
    timings: bool = False,
    verbose: bool = True,
    frag_mode: str = 'both',
    merge_parallel: int = None,
) -> None:
    """Build a v1 CReM fragment DB using N concurrent shard builders.

    Each child runs ``cremdb_create`` itself with ``--stride-mod N --stride-idx i``
    so it processes only ``chunk_id % N == i`` of the input. Available CPUs
    are split evenly across children. After all children finish, the parts
    are merged into ``output_db`` via binary-tree reduction.

    Resume: re-running the same command re-launches only those children
    whose shard DB does not yet carry the finalisation sentinel.

    Args:
        See ``run`` for shared kwargs.
        parallel_shards: Number of concurrent shard builders.
        merge_parallel: Max concurrent pair-merges during the final tree
            merge. Defaults to ``parallel_shards``.
    """
    if parallel_shards < 2:
        raise ValueError("parallel_shards must be >= 2 (use run() for serial)")

    start_time = time.time()

    parts_dir = Path(f"{output_db}.parts")
    parts_dir.mkdir(parents=True, exist_ok=True)

    ncpu_per_shard = max(1, ncpu // parallel_shards)
    if merge_parallel is None:
        merge_parallel = parallel_shards

    # Build the child argv from the current command-line, removing the
    # orchestration-only flags and adding per-child stride/output flags.
    # The child runs the same Python entrypoint as the parent.
    base_argv = [sys.executable, "-m", "crem.scripts.cremdb_create"]
    base_argv.extend(["--input", input_path])
    base_argv.extend(["--set-name", *set_name])
    base_argv.extend(["--radii", *(str(r) for r in radii)])
    base_argv.extend(["--chunk-size", str(chunk_size)])
    base_argv.extend(["--max-heavy-atoms", str(max_heavy_atoms)])
    if keep_stereo:
        base_argv.append("--keep-stereo")
    base_argv.extend(["--mode", str(mode)])
    base_argv.extend(["--frag-mode", frag_mode])
    if sep is not None:
        base_argv.extend(["--sep", sep])
    if force_zstd:
        base_argv.append("--zstd")
    if log_every is not None:
        base_argv.extend(["--log-every", str(log_every)])
    base_argv.extend(["--flush-every", str(flush_every)])
    base_argv.extend(["--prefetch", str(prefetch)])
    if timings:
        base_argv.append("--timings")

    procs: List[subprocess.Popen] = []
    log_handles: List = []
    try:
        for i in range(parallel_shards):
            shard_db = str(parts_dir / f"shard_{i:03d}.db")
            chunk_file = str(parts_dir / f"processed_chunks_{i:03d}.txt")
            log_file = str(parts_dir / f"shard_{i:03d}.log")

            if _shard_is_finalised(shard_db):
                if verbose:
                    sys.stderr.write(f"  Shard {i:03d}: already finalised, skipping\n")
                    sys.stderr.flush()
                continue

            child_argv = list(base_argv)
            child_argv.extend(["--output", shard_db])
            child_argv.extend(["--processed-chunks", chunk_file])
            child_argv.extend(["--ncpu", str(ncpu_per_shard)])
            child_argv.extend(["--stride-mod", str(parallel_shards)])
            child_argv.extend(["--stride-idx", str(i)])

            log_fh = open(log_file, "ab")
            log_handles.append(log_fh)
            if verbose:
                sys.stderr.write(
                    f"  Launching shard {i:03d}: {ncpu_per_shard} CPUs, "
                    f"log {log_file}\n"
                )
                sys.stderr.flush()
            procs.append(subprocess.Popen(
                child_argv,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
            ))

        # Wait on all children. On Ctrl-C, propagate SIGTERM and re-raise.
        failures = []
        try:
            for p in procs:
                rc = p.wait()
                if rc != 0:
                    failures.append((p.pid, rc))
        except KeyboardInterrupt:
            for p in procs:
                if p.poll() is None:
                    try:
                        p.send_signal(signal.SIGTERM)
                    except ProcessLookupError:
                        pass
            for p in procs:
                try:
                    p.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    p.kill()
            raise

        if failures:
            details = ", ".join(f"pid={pid} rc={rc}" for pid, rc in failures)
            raise RuntimeError(
                f"{len(failures)} shard builder(s) failed: {details}. "
                f"Parts left in {parts_dir} for inspection/resume."
            )
    finally:
        for fh in log_handles:
            try:
                fh.close()
            except Exception:
                pass

    # All children finished successfully. Collect every shard (including any
    # that we skipped because already finalised) for the final merge.
    shard_paths = [
        str(parts_dir / f"shard_{i:03d}.db") for i in range(parallel_shards)
    ]
    missing = [p for p in shard_paths if not os.path.exists(p)]
    if missing:
        raise RuntimeError(
            f"Shard files missing after build: {missing}. Re-run to resume."
        )

    if verbose:
        sys.stderr.write(
            f"\nMerging {len(shard_paths)} shards into {output_db} "
            f"(parallel={merge_parallel})...\n"
        )
        sys.stderr.flush()

    from crem.scripts.cremdb_merge import binary_tree_merge_into  # lazy import
    binary_tree_merge_into(
        final_target_path=output_db,
        shard_paths=shard_paths,
        max_parallel=merge_parallel,
        rebuild_index=True,
        verbose=verbose,
    )

    elapsed = time.time() - start_time

    if verbose:
        # Orchestrator-side summary. Build-time counters (chunks processed,
        # molecules processed, fragments generated, env/core pairs) only exist
        # inside each child process and are deliberately omitted here.
        with sqlite3.connect(output_db) as stats_conn:
            cur = stats_conn.cursor()
            env_count = cur.execute("SELECT COUNT(*) FROM envs").fetchone()[0]
            frag_count = cur.execute("SELECT COUNT(*) FROM frags").fetchone()[0]
            frag_h_count = cur.execute("SELECT COUNT(*) FROM frags_h").fetchone()[0]
            radius_counts = {}
            for radius in radii:
                radius_counts[radius] = cur.execute(
                    f"SELECT COUNT(rowid) FROM radius{radius}"
                ).fetchone()[0]

        print("Done.")
        print(f"Unique envs: {env_count}")
        print(f"Unique frags: {frag_count}")
        print(f"Unique frags_h: {frag_h_count}")
        for radius, count in radius_counts.items():
            print(f"radius{radius} rows: {count}")
        print(f"Elapsed: {elapsed:.1f}s")


def main():
    parser = argparse.ArgumentParser(description="Stream-like fragment DB creation for new schema")
    parser.add_argument("-i", "--input", required=True, help="Input SMILES file (text or .zst)")
    parser.add_argument("-o", "--output", dest="output_db", required=True, help="Output SQLite DB file")
    parser.add_argument(
        "-s", "--set-name",
        dest="set_name",
        nargs="+",
        required=True,
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
    parser.add_argument(
        "--frag-mode",
        choices=["acyclic", "ring", "both"],
        default="both",
        help=(
            "Fragmentation source: 'acyclic' (cuts of acyclic bonds, today's "
            "CReM behaviour); 'ring' (cuts of pairs of SINGLE bonds inside "
            "the same ring, used by make_cycle(ring_closures=True)); "
            "'both' includes both kinds, tagging each row with is_ring_closure "
            "(default: both)"
        ),
    )
    parser.add_argument("--sep", default=None, help="Input delimiter (default: whitespace)")
    parser.add_argument("--processed-chunks", default=None, help="Path to processed chunks file (append)")
    parser.add_argument("--zstd", action="store_true", help="Force zstd input")
    parser.add_argument("--log-every", type=int, default=None, help="Log progress every N chunks (default: None)")
    parser.add_argument("--ncpu", type=int, default=1, help="Number of worker processes (default: 1)")
    parser.add_argument(
        "--flush-every",
        type=int,
        default=100,
        help="Accumulate this many chunks in memory before flushing to DB (default: 100)",
    )
    parser.add_argument(
        "--prefetch",
        type=int,
        default=4,
        help="In-flight task batches per worker (batch_size = ncpu * prefetch). "
             "Lower values reduce peak memory at the cost of potentially underutilising workers (default: 4)",
    )
    parser.add_argument("--timings", action="store_true",
        help="Print per-flush timing breakdown to stderr")
    parser.add_argument(
        "--shard-size",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Max number of input structures per shard database. "
            "If omitted, all data goes to one DB (default behaviour). "
            "Shard 0 is the output file itself; subsequent shards get a "
            "numeric suffix (e.g. output_001.db, output_002.db, ...). "
            "On resume, the last existing shard is detected automatically "
            "and a fresh new shard is started — the partial last shard is "
            "left as-is and will be merged at the end."
        ),
    )
    parser.add_argument(
        "--parallel-shards",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Build N shards concurrently, each on a stride of the input. "
            "CPUs from --ncpu are split evenly across them. Each shard DB "
            "lives in <output>.parts/shard_NNN.db; the parts are merged into "
            "<output> via a parallel binary-tree reduction at the end. "
            "Default 1 (current single-process behaviour)."
        ),
    )
    # Internal flags used by the orchestrator to drive its children. Not
    # intended for direct invocation.
    parser.add_argument("--stride-mod", type=int, default=1, help=argparse.SUPPRESS)
    parser.add_argument("--stride-idx", type=int, default=0, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.parallel_shards < 1:
        parser.error("--parallel-shards must be >= 1")
    if args.parallel_shards > 1:
        if args.shard_size is not None:
            parser.error("--parallel-shards is incompatible with --shard-size")
        if args.stride_mod > 1:
            parser.error("--parallel-shards is incompatible with --stride-mod (internal)")
        run_parallel_shards(
            input_path=args.input,
            output_db=args.output_db,
            set_name=args.set_name,
            parallel_shards=args.parallel_shards,
            ncpu=args.ncpu,
            radii=args.radii,
            chunk_size=args.chunk_size,
            max_heavy_atoms=args.max_heavy_atoms,
            keep_stereo=args.keep_stereo,
            mode=args.mode,
            sep=args.sep,
            force_zstd=args.zstd,
            log_every=args.log_every,
            flush_every=args.flush_every,
            prefetch=args.prefetch,
            timings=args.timings,
            frag_mode=args.frag_mode,
        )
        return

    run(
        input_path=args.input,
        output_db=args.output_db,
        set_name=args.set_name,
        radii=args.radii,
        chunk_size=args.chunk_size,
        max_heavy_atoms=args.max_heavy_atoms,
        keep_stereo=args.keep_stereo,
        mode=args.mode,
        sep=args.sep,
        processed_chunks=args.processed_chunks,
        force_zstd=args.zstd,
        log_every=args.log_every,
        ncpu=args.ncpu,
        flush_every=args.flush_every,
        prefetch=args.prefetch,
        timings=args.timings,
        shard_size=args.shard_size,
        frag_mode=args.frag_mode,
        stride_mod=args.stride_mod,
        stride_idx=args.stride_idx,
    )


if __name__ == "__main__":
    main()
