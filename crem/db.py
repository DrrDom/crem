"""Python API for CReM fragment database management.

All database operations (creation, merging, property annotation, inspection) are
available as plain Python functions importable from this module:

    from crem.db import create_db, merge_dbs, add_fragment_props, get_db_info

This module also owns the schema column constants shared by the generation code
and the command-line tools; it is deliberately free of RDKit imports at module
level so that importing a constant costs nothing.
"""

import os
import pickle
import sqlite3
import tempfile
from multiprocessing import Pool
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Union

PathLike = Union[str, Path]

_BUILTIN_PROPS = ('mw', 'logp', 'rtb', 'tpsa', 'fcsp3')

# Table metadata: (id_column, smiles_column)
_TABLE_COLS = {
    'frags':   ('core_smi_id', 'core_smi'),
    'frags_h': ('core_smi_h_id', 'smi'),
}

#: Columns of radius{N} that are schema metadata rather than per-set occurrence counts.
#: `is_ring_closure` is absent from v2 tables; naming it here is harmless because set
#: discovery is a set difference.
_RESERVED_RADIUS_COLUMNS = frozenset(
    {'env_id', 'core_smi_id', 'core_num_atoms', 'dist2', 'is_ring_closure'}
)

#: Schema columns of the v1/v2 fragment tables. Everything else in them is a property
#: column added after the build (by cremdb_add_prop or add_fragment_props), so property
#: discovery is a set difference against these.
#: `core_num_atoms` is listed for `frags` because older v1 builds carry it there while
#: current builds denormalize it into radius{N} only - naming a column that may be
#: absent is harmless in a set difference, and omitting it would report it as a property.
BASE_TABLE_COLUMNS = {
    'frags':   frozenset({'core_smi_id', 'core_smi', 'core_smi_h_id', 'core_num_atoms'}),
    'frags_h': frozenset({'core_smi_h_id', 'smi'}),
}

#: Schema columns of a v0 radius{N} table (`freq` only when built with counts). v0 has no
#: separate fragment table, so property columns live on the radius tables themselves.
V0_BASE_COLUMNS = frozenset({'env', 'core_smi', 'core_num_atoms', 'core_sma', 'dist2', 'freq'})

_CUSTOM_WRITE_BATCH = 10_000

# Sentinel used by add_fragment_props to distinguish "argument not supplied"
# from explicit None / [] — see that function's docstring.
_PROPS_DEFAULT = object()


def create_db(
    input: Union[PathLike, Iterable[str]],
    output: PathLike,
    set_name: Union[str, Dict[str, Optional[set]]],
    radii=(1, 2, 3, 4, 5),
    *,
    ncpu: int = 1,
    max_heavy_atoms: int = 15,
    keep_stereo: bool = False,
    mode: int = 0,
    chunk_size: int = 100,
    flush_every: int = 100,
    shard_size: Optional[int] = None,
    parallel_shards: int = 1,
    frag_mode: str = 'both_optimal',
    verbose: bool = True,
    sep: Optional[str] = None,
    processed_chunks: Optional[PathLike] = None,
    force_zstd: bool = False,
    log_every: Optional[int] = None,
    prefetch: int = 4,
    timings: bool = False,
    merge_parallel: Optional[int] = None,
    fragment_error_log: bool = False,
) -> None:
    """Create or extend a v2 CReM fragment database.

    Calling on an existing database is safe and additive: ``_ensure_schema``
    uses ``CREATE TABLE IF NOT EXISTS`` and incremental ``ALTER TABLE``, so
    any new set names or radii are added and existing data is preserved.

    :param input: path to a SMILES file (``str`` / ``Path``) **or** an iterable
        of ``"SMILES [ID]"`` strings (one molecule per item).
    :param output: path to the output SQLite database.
    :param set_name: a single set name (``str``), or a ``dict`` mapping each set
        name to either ``None`` (all molecules) or a ``set`` of molecule IDs
        that belong to that set.
    :param radii: fragment radii to build (default 1–5).
    :param ncpu: worker processes.
    :param max_heavy_atoms: maximum heavy atoms in a core fragment.
    :param keep_stereo: preserve stereocentres in env/core SMILES.
    :param mode: fragmentation mode — 0 all atoms, 1 heavy only, 2 H only.
    :param chunk_size: input lines per worker task.
    :param flush_every: chunks to accumulate before each DB flush.
    :param shard_size: max input structures per shard DB (``None`` = single DB).
        Incompatible with ``parallel_shards > 1``.
    :param parallel_shards: when > 1, run N shard builders concurrently, each
        fragmenting a stride of the input. CPUs from ``ncpu`` are split evenly
        across them. Shard DBs live in ``<output>.parts/`` and are merged into
        ``output`` via a parallel binary-tree reduction. Default 1
        (single-process build).
    :param frag_mode: fragmentation source: ``'acyclic'``, ``'ring'``,
        ``'both'``, ``'ring_optimal'``, or ``'both_optimal'``. Default
        ``'both_optimal'`` matches ``cremdb_create``.
    :param verbose: print progress and statistics to stdout/stderr.
    :param sep: input delimiter (``None`` = whitespace).
    :param processed_chunks: path to a processed-chunks file for resumable
        non-parallel builds from file input. Ignored when ``input`` is an
        iterable. Also ignored for ``parallel_shards > 1``; parallel builds
        manage per-shard processed-chunk files internally.
    :param force_zstd: force zstd input decompression regardless of file suffix.
    :param log_every: print a progress line every N chunks (``None`` = silent).
    :param prefetch: in-flight task batches per worker.
    :param timings: print per-flush timing breakdown to stderr.
    :param merge_parallel: max concurrent pair-merges for ``parallel_shards > 1``.
    :param fragment_error_log: write defensive fragment validation issues to
        ``<output>.errors``. If false, issues are written to stderr.
    """
    if parallel_shards < 1:
        raise ValueError("parallel_shards must be >= 1")
    if parallel_shards > 1 and shard_size is not None:
        raise ValueError("parallel_shards > 1 is incompatible with shard_size")
    from crem.scripts.cremdb_create import run as _run, run_parallel_shards as _run_parallel

    tmp_input: Optional[str] = None
    tmp_ids: List[str] = []

    try:
        # --- resolve input ---------------------------------------------------
        if isinstance(input, (str, Path)):
            input_path = str(input)
            processed_chunks_arg = (
                str(processed_chunks) if processed_chunks is not None else None
            )
        else:
            processed_chunks_arg = None
            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.smi', delete=False, encoding='utf-8'
            ) as fh:
                tmp_input = fh.name
                for line in input:
                    fh.write(line.rstrip('\n') + '\n')
            input_path = tmp_input

        # --- resolve set_name ------------------------------------------------
        if isinstance(set_name, str):
            set_name_arg = [set_name]
        elif isinstance(set_name, dict):
            set_name_arg = []
            for name, ids in set_name.items():
                set_name_arg.append(name)
                if ids is not None:
                    with tempfile.NamedTemporaryFile(
                        mode='w', suffix='.txt', delete=False, encoding='utf-8'
                    ) as fh:
                        tmp_ids.append(fh.name)
                        for mol_id in ids:
                            fh.write(str(mol_id) + '\n')
                    set_name_arg.append(tmp_ids[-1])
        else:
            raise TypeError("set_name must be a str or dict")

        if parallel_shards > 1:
            _run_parallel(
                input_path=input_path,
                output_db=str(output),
                set_name=set_name_arg,
                parallel_shards=parallel_shards,
                ncpu=ncpu,
                radii=list(radii),
                chunk_size=chunk_size,
                max_heavy_atoms=max_heavy_atoms,
                keep_stereo=keep_stereo,
                mode=mode,
                flush_every=flush_every,
                verbose=verbose,
                frag_mode=frag_mode,
                sep=sep,
                force_zstd=force_zstd,
                log_every=log_every,
                prefetch=prefetch,
                timings=timings,
                merge_parallel=merge_parallel,
                fragment_error_log=fragment_error_log,
            )
        else:
            _run(
                input_path=input_path,
                output_db=str(output),
                set_name=set_name_arg,
                radii=list(radii),
                chunk_size=chunk_size,
                max_heavy_atoms=max_heavy_atoms,
                keep_stereo=keep_stereo,
                mode=mode,
                flush_every=flush_every,
                shard_size=shard_size,
                ncpu=ncpu,
                verbose=verbose,
                frag_mode=frag_mode,
                sep=sep,
                processed_chunks=processed_chunks_arg,
                force_zstd=force_zstd,
                log_every=log_every,
                prefetch=prefetch,
                timings=timings,
                fragment_error_log=fragment_error_log,
            )

    finally:
        if tmp_input and os.path.exists(tmp_input):
            os.unlink(tmp_input)
        for p in tmp_ids:
            if os.path.exists(p):
                os.unlink(p)


def merge_dbs(
    target: PathLike,
    sources: List[PathLike],
    *,
    rebuild_index: bool = True,
    parallel: int = 1,
    verbose: bool = True,
) -> None:
    """Merge source shard databases into ``target``.

    :param target: path to the target (base) database. Must already exist.
    :param sources: list of source shard database paths to merge in.
    :param rebuild_index: recreate covering indices on the target after merge.
    :param parallel: when > 1, merge with binary-tree reduction using up to this
        many concurrent pair-merges per round. The target is treated as one of
        the contributors; the final survivor is moved back to ``target``.
        Default 1 (serial).
    :param verbose: print per-shard progress.
    """
    if parallel < 1:
        raise ValueError("parallel must be >= 1")
    from crem.scripts.cremdb_merge import run as _run
    _run(
        target_path=str(target),
        source_paths=[str(s) for s in sources],
        rebuild_index=rebuild_index,
        verbose=verbose,
        parallel=parallel,
    )


def add_fragment_props(
    db: PathLike,
    properties=_PROPS_DEFAULT,
    *,
    custom_props: Optional[Dict[str, Callable[[str], float]]] = None,
    table: str = 'frags',
    ncpu: int = 1,
    verbose: bool = False,
) -> None:
    """Add molecular properties to a CReM fragment database.

    Only rows with ``NULL`` property values are processed, so calling this
    function after adding new fragments fills only the newly added rows.

    Built-in properties are computed on the ``frags`` table (``core_smi``
    column) using RDKit descriptors.  Custom properties can target either
    ``'frags'`` (``core_smi``) or ``'frags_h'`` (H-replaced SMILES ``smi``).

    :param db: path to the fragment database.
    :param properties: built-in property names to compute (``'mw'``, ``'logp'``,
        ``'rtb'``, ``'tpsa'``, ``'fcsp3'``). Accepted values: if **omitted**, all
        built-ins are computed when ``custom_props`` is not given and **no**
        built-ins when ``custom_props`` is given (so
        ``add_fragment_props(db, custom_props={...})`` adds only the custom
        columns, while the usual ``add_fragment_props(db)`` is unchanged);
        ``None`` or ``'all'`` forces all built-ins (combine with ``custom_props``
        to add both at once); a list/tuple computes that subset; and ``[]``
        skips built-ins entirely.
    :param custom_props: mapping of ``{column_name: func(smi) -> value}``.
        Picklable functions (named functions, ``functools.partial``) use
        ``ncpu`` workers; non-picklable ones (lambdas, closures) are processed
        serially.
    :param table: target table for ``custom_props`` — ``'frags'`` or
        ``'frags_h'``.
    :param ncpu: workers for built-in and picklable custom properties.
    :param verbose: print progress to stderr.
    """
    if table not in _TABLE_COLS:
        raise ValueError(f"table must be one of {list(_TABLE_COLS)}, got {table!r}")

    if properties is _PROPS_DEFAULT:
        compute_builtins = custom_props is None
        builtins_arg: Optional[List[str]] = None
    elif properties is None or properties == 'all':
        compute_builtins = True
        builtins_arg = None
    elif isinstance(properties, (list, tuple)) and properties:
        compute_builtins = True
        builtins_arg = list(properties)
    else:  # [] / () / explicit empty → skip built-ins
        compute_builtins = False
        builtins_arg = None

    if compute_builtins:
        from crem.scripts.cremdb_add_prop import run as _run
        _run(db_path=str(db), properties=builtins_arg, ncpu=ncpu, verbose=verbose)

    if custom_props:
        _add_custom_props(str(db), custom_props, table=table, ncpu=ncpu, verbose=verbose)


def get_db_info(db: PathLike) -> Dict:
    """Summarize the schema of a CReM fragment database.

    Reads only PRAGMAs and `sqlite_master`, so it is cheap on databases of any size,
    and opens the file read-only: a database is never created or modified by this call.

    :param db: path to the fragment database.
    :return: a dict with the keys

        * ``path`` - the path as given;
        * ``version`` - ``PRAGMA user_version`` (0 legacy, 1 deprecated, 2 current);
        * ``has_sets`` - whether the schema supports fragment sets (False for v0);
        * ``radius_tables`` - ``{'radius1': [set names], ...}``, ordered by radius; the
          lists are empty for v0, which has no sets;
        * ``properties`` - ``{table: [property columns]}`` in column order, covering
          ``frags`` and ``frags_h`` for v1/v2 and each radius table for v0. These are
          the columns that may be used as property filters at generation time;
        * ``is_shard`` - True for an unmerged stride-mode shard written by
          ``cremdb_create --parallel-shards``.

    :raises FileNotFoundError: if the path does not exist.
    :raises ValueError: if the file is not a recognizable CReM database.
    :raises sqlite3.OperationalError: if the file cannot be opened read-only (for
        instance a database with a hot write-ahead log needing recovery).
    """
    from crem.scripts.cremdb_create import _STRIDE_SHARD_SENTINEL

    path = str(db)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{path}: file not found")

    con = _connect_ro(path)
    try:
        version = con.execute("PRAGMA user_version").fetchone()[0]
        # application_id is a 32-bit *signed* field, so a sentinel with the high bit set
        # reads back as its negative twin (and SQLite ignores the write altogether when
        # the value does not fit). Compare the low 32 bits, which is true either way.
        app_id = con.execute("PRAGMA application_id").fetchone()[0]
        is_shard = (app_id & 0xFFFFFFFF) == (_STRIDE_SHARD_SENTINEL & 0xFFFFFFFF)
        tables = {row[0] for row in
                  con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        radius_names = _radius_table_names(tables)

        info = {'path': path, 'version': version, 'is_shard': is_shard}

        if version >= 1:
            # v1/v2 (and any future member of that family): sets are the radius{N}
            # columns that are not schema, properties live on the shared fragment tables.
            if 'frags' not in tables:
                raise ValueError(
                    f"{path}: user_version={version} implies the normalized schema, "
                    f"but the frags table is missing - this is not a CReM database")
            info['has_sets'] = True
            info['radius_tables'] = {
                name: sorted(_columns(con, name) - _RESERVED_RADIUS_COLUMNS)
                for name in radius_names
            }
            info['properties'] = {
                name: [c for c in _columns(con, name, ordered=True)
                       if c not in BASE_TABLE_COLUMNS[name]]
                for name in ('frags', 'frags_h') if name in tables
            }
        else:
            # v0, or a SQLite file that simply is not a CReM database: user_version
            # defaults to 0, so the radius tables are the only evidence either way.
            if not radius_names:
                raise ValueError(
                    f"{path}: no radius tables and no schema version - "
                    f"this is not a CReM database")
            info['has_sets'] = False
            info['radius_tables'] = {name: [] for name in radius_names}
            info['properties'] = {
                name: [c for c in _columns(con, name, ordered=True)
                       if c not in V0_BASE_COLUMNS]
                for name in radius_names
            }
        return info
    finally:
        con.close()


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _connect_ro(path: str) -> sqlite3.Connection:
    """Open a database read-only, so that inspecting it can neither create the file
    (which a plain connect() to a mistyped path does) nor write to it."""
    from urllib.request import pathname2url
    return sqlite3.connect(f"file:{pathname2url(os.path.abspath(path))}?mode=ro", uri=True)


def _columns(con: sqlite3.Connection, table: str, ordered: bool = False):
    """Column names of `table`, in declaration order when `ordered`, else as a set."""
    names = [row[1] for row in con.execute(f"PRAGMA table_info({table})")]
    return names if ordered else set(names)


def _radius_table_names(tables: Iterable[str]) -> List[str]:
    """`radiusN` table names sorted by N. Names with a non-numeric suffix are ignored
    rather than raising: nothing stops a user table from starting with "radius"."""
    found = []
    for name in tables:
        if name.startswith('radius'):
            try:
                found.append((int(name[6:]), name))
            except ValueError:
                continue
    return [name for _, name in sorted(found)]

def _is_picklable(obj) -> bool:
    try:
        pickle.dumps(obj)
        return True
    except Exception:
        return False


def _custom_prop_worker_init(funcs):
    global _CUSTOM_FUNCS
    _CUSTOM_FUNCS = funcs


def _custom_prop_worker(item):
    row_id, smi = item
    return row_id, [f(smi) for f in _CUSTOM_FUNCS]


def _add_custom_props(
    db_path: str,
    custom_props: Dict[str, Callable],
    table: str,
    ncpu: int,
    verbose: bool,
) -> None:
    id_col, smi_col = _TABLE_COLS[table]
    col_names = list(custom_props.keys())
    funcs = list(custom_props.values())

    picklable = all(_is_picklable(f) for f in funcs)

    with sqlite3.connect(db_path) as conn:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA cache_size = -65536")

        # Add columns (silently skip if already present).
        for col in col_names:
            try:
                conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} NUMERIC DEFAULT NULL")
            except sqlite3.OperationalError:
                pass
        conn.commit()

        null_filter = " OR ".join(f"{c} IS NULL" for c in col_names)
        rows = conn.execute(
            f"SELECT {id_col}, {smi_col} FROM {table} WHERE {null_filter}"
        ).fetchall()

        if not rows:
            if verbose:
                import sys
                sys.stderr.write(f"  No unlabeled rows in {table} for custom props\n")
            return

        if verbose:
            import sys
            sys.stderr.write(f"  Computing {col_names} for {len(rows)} rows in {table}\n")

        update_sql = (
            f"UPDATE {table} SET "
            + ", ".join(f"{c} = ?" for c in col_names)
            + f" WHERE {id_col} = ?"
        )

        if picklable and ncpu > 1:
            pool = Pool(ncpu, initializer=_custom_prop_worker_init, initargs=(funcs,))
            try:
                results = pool.map(_custom_prop_worker, rows)
            finally:
                pool.close()
                pool.join()
        else:
            results = [(row_id, [f(smi) for f in funcs]) for row_id, smi in rows]

        write_buf = []
        for row_id, vals in results:
            write_buf.append((*vals, row_id))
            if len(write_buf) >= _CUSTOM_WRITE_BATCH:
                conn.executemany(update_sql, write_buf)
                conn.commit()
                write_buf = []

        if write_buf:
            conn.executemany(update_sql, write_buf)
            conn.commit()

        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
