"""Python API for CReM fragment database management.

All three database operations (creation, merging, property annotation) are
available as plain Python functions importable from this module:

    from crem.db import create_db, merge_dbs, add_fragment_props
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

_CUSTOM_WRITE_BATCH = 10_000

# Sentinel used by add_fragment_props to distinguish "argument not supplied"
# from explicit None / [] — see that function's docstring.
_PROPS_DEFAULT = object()


def create_db(output: PathLike, input: Union[PathLike, Iterable[str]], set_name: Union[str, Dict[str, Optional[set]]],
              radii=(1, 2, 3, 4, 5), *, ncpu: int = 1, max_heavy_atoms: int = 15, keep_stereo: bool = False,
              mode: int = 0, chunk_size: int = 100, flush_every: int = 100, shard_size: Optional[int] = None,
              verbose: bool = True) -> None:
    """Create or extend a v1 CReM fragment database.

    Calling on an existing database is safe and additive: ``_ensure_schema``
    uses ``CREATE TABLE IF NOT EXISTS`` and incremental ``ALTER TABLE``, so
    any new set names or radii are added and existing data is preserved.

    Args:
        input: Path to a SMILES file (``str`` / ``Path``) **or** an iterable
            of ``"SMILES [ID]"`` strings (one molecule per item).
        output: Path to the output SQLite database.
        set_name: A single set name (``str``), or a ``dict`` mapping each set
            name to either ``None`` (all molecules) or a ``set`` of molecule
            IDs that belong to that set.
        radii: Fragment radii to build (default 1–5).
        ncpu: Worker processes.
        max_heavy_atoms: Maximum heavy atoms in a core fragment.
        keep_stereo: Preserve stereocentres in env/core SMILES.
        mode: Fragmentation mode — 0 all atoms, 1 heavy only, 2 H only.
        chunk_size: Input lines per worker task.
        flush_every: Chunks to accumulate before each DB flush.
        shard_size: Max input structures per shard DB (``None`` = single DB).
        verbose: Print progress and statistics to stdout/stderr.
    """
    from crem.scripts.cremdb_create import run as _run

    tmp_input: Optional[str] = None
    tmp_ids: List[str] = []

    try:
        # --- resolve input ---------------------------------------------------
        if isinstance(input, (str, Path)):
            input_path = str(input)
        else:
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
            # pass-through defaults for rarely needed options
            sep=None,
            processed_chunks=None,
            force_zstd=False,
            log_every=None,
            prefetch=4,
            timings=False,
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
    verbose: bool = True,
) -> None:
    """Merge source shard databases into *target*.

    Args:
        target: Path to the target (base) database. Must already exist.
        sources: List of source shard database paths to merge in.
        rebuild_index: Recreate covering indices on the target after merge.
        verbose: Print per-shard progress.
    """
    from crem.scripts.cremdb_merge import run as _run
    _run(
        target_path=str(target),
        source_paths=[str(s) for s in sources],
        rebuild_index=rebuild_index,
        verbose=verbose,
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

    Args:
        db: Path to the fragment database.
        properties: Built-in property names to compute
            (``'mw'``, ``'logp'``, ``'rtb'``, ``'tpsa'``, ``'fcsp3'``).
            Behaviour by value:

            * **omitted** — defaults to all built-ins when *custom_props* is
              not given, and to **no** built-ins when *custom_props* is given.
              This makes ``add_fragment_props(db, custom_props={...})`` add
              only the requested custom columns; the typical built-in case
              ``add_fragment_props(db)`` is unchanged.
            * ``None`` or ``'all'`` — explicitly compute all built-ins
              (use this together with *custom_props* to add both at once).
            * a list/tuple of names — compute that subset.
            * ``[]`` — explicitly skip built-ins.
        custom_props: Mapping of ``{column_name: func(smi) -> value}``.
            Picklable functions (named functions, ``functools.partial``) use
            *ncpu* workers; non-picklable ones (lambdas, closures) are
            processed serially.
        table: Target table for *custom_props* — ``'frags'`` or ``'frags_h'``.
        ncpu: Workers for built-in and picklable custom properties.
        verbose: Print progress to stderr.
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


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

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
