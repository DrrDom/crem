#!/usr/bin/env python3

import argparse
import sqlite3
import sys
from multiprocessing import Pool

from crem.arg_types import filepath_type, cpu_type
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds, CalcTPSA, CalcFractionCSP3


_ALL_PROPS = ['mw', 'logp', 'rtb', 'tpsa', 'fcsp3']

_FETCH_BATCH = 50000
_WRITE_BATCH = 10000
_IMAP_CHUNK = 500

_SELECTED_PROPS = []   # set in worker initializer


def _init_worker(selected_props):
    global _SELECTED_PROPS
    _SELECTED_PROPS = selected_props


def _calc(item):
    row_id, smi = item
    mol = Chem.MolFromSmiles(smi)
    vals = []
    for prop in _SELECTED_PROPS:
        if mol is None:
            vals.append(None)
        elif prop == 'mw':
            vals.append(round(MolWt(mol), 2))
        elif prop == 'logp':
            vals.append(round(MolLogP(mol), 2))
        elif prop == 'rtb':
            vals.append(CalcNumRotatableBonds(Chem.RemoveHs(mol)))
        elif prop == 'tpsa':
            vals.append(CalcTPSA(mol))
        elif prop == 'fcsp3':
            vals.append(round(CalcFractionCSP3(mol), 3))
    return row_id, vals


def _add_columns(conn, table, selected_props):
    for prop in selected_props:
        try:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {prop} NUMERIC DEFAULT NULL")
        except sqlite3.OperationalError as e:
            sys.stderr.write(str(e) + '\n')
    conn.commit()


def _process_table(conn, pool, table, id_col, smi_col, selected_props,
                   verbose, fetch_batch, write_batch, imap_chunk):
    null_filter = " OR ".join(f"{p} IS NULL" for p in selected_props)

    total = conn.execute(f"SELECT COUNT(*) FROM {table} WHERE {null_filter}").fetchone()[0]
    if total == 0:
        if verbose:
            sys.stderr.write(f"  No unlabeled rows in {table}\n")
        return

    if verbose:
        sys.stderr.write(f"  {total} rows to process in {table}\n")

    cur = conn.cursor()
    cur.execute(f"SELECT {id_col}, {smi_col} FROM {table} WHERE {null_filter}")

    update_sql = (
        f"UPDATE {table} SET "
        + ", ".join(f"{p} = ?" for p in selected_props)
        + f" WHERE {id_col} = ?"
    )

    write_buf = []
    processed = 0

    while True:
        rows = cur.fetchmany(fetch_batch)
        if not rows:
            break
        for row_id, vals in pool.imap_unordered(_calc, rows, chunksize=imap_chunk):
            write_buf.append((*vals, row_id))
            if len(write_buf) >= write_batch:
                conn.executemany(update_sql, write_buf)
                conn.commit()
                processed += len(write_buf)
                write_buf = []
                if verbose:
                    sys.stderr.write(f"\r  {processed}/{total}")
                    sys.stderr.flush()

    if write_buf:
        conn.executemany(update_sql, write_buf)
        conn.commit()
        processed += len(write_buf)

    if verbose:
        sys.stderr.write(f"\r  {processed}/{total}\n")


def run(
    db_path: str,
    properties=None,
    ncpu: int = 1,
    verbose: bool = False,
    fetch_batch: int = _FETCH_BATCH,
    write_batch: int = _WRITE_BATCH,
    imap_chunk: int = _IMAP_CHUNK,
) -> None:
    """Add built-in molecular properties to a CReM fragment database.

    Only rows with NULL property values are processed, so re-running after
    adding new fragments fills only the new rows.

    Args:
        db_path: Path to the SQLite fragment database.
        properties: List of property names to compute (``'mw'``, ``'logp'``,
            ``'rtb'``, ``'tpsa'``, ``'fcsp3'``). ``None`` computes all.
        ncpu: Worker processes.
        verbose: Print progress to stderr.
        fetch_batch: Rows fetched per batch.
        write_batch: Rows per executemany write.
        imap_chunk: Multiprocessing imap chunksize.
    """
    selected_props = list(properties) if properties is not None else list(_ALL_PROPS)
    if not selected_props:
        sys.stderr.write('No valid property names supplied.\n')
        return

    pool = Pool(ncpu, initializer=_init_worker, initargs=(selected_props,))

    with sqlite3.connect(db_path) as conn:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA cache_size = -65536")

        version = conn.execute("PRAGMA user_version").fetchone()[0]

        if version == 1:
            all_tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if 'frags' not in all_tables:
                raise RuntimeError("Detected new schema (user_version=1) but frags table is missing")
            if verbose:
                sys.stderr.write("Schema version 1: adding properties to frags table\n")
            _add_columns(conn, 'frags', selected_props)
            _process_table(
                conn, pool, 'frags', 'core_smi_id', 'core_smi', selected_props,
                verbose, fetch_batch, write_batch, imap_chunk,
            )

        elif version == 0:
            radius_tables = [
                row[0] for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'radius%'"
                )
            ]
            if verbose:
                sys.stderr.write(f"Old schema: adding properties to radius tables: {radius_tables}\n")
            for table in radius_tables:
                if verbose:
                    sys.stderr.write(f"\nTable {table}\n")
                _add_columns(conn, table, selected_props)
                _process_table(
                    conn, pool, table, 'rowid', 'core_smi', selected_props,
                    verbose, fetch_batch, write_batch, imap_chunk,
                )

        else:
            raise ValueError(f"Unsupported database version: {version}")

        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")

    pool.close()
    pool.join()

    if verbose:
        sys.stderr.write(f"\nFinished: {db_path}\n")


def entry_point():
    parser = argparse.ArgumentParser(
        description='Add columns with values of chosen properties to CReM fragment database.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=filepath_type,
                        help='SQLite DB with CReM fragments.')
    parser.add_argument('-p', '--properties', metavar='NAMES', required=False, nargs='*',
                        default=_ALL_PROPS, choices=_ALL_PROPS,
                        help='properties to compute.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')
    parser.add_argument('--fetch-batch', metavar='INTEGER', default=_FETCH_BATCH, type=int,
                        help='rows fetched from DB per batch.')
    parser.add_argument('--write-batch', metavar='INTEGER', default=_WRITE_BATCH, type=int,
                        help='rows per executemany write.')
    parser.add_argument('--imap-chunk', metavar='INTEGER', default=_IMAP_CHUNK, type=int,
                        help='multiprocessing imap chunksize.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    if not args.properties:
        sys.stderr.write('No valid property names supplied.\n')
        sys.exit(1)

    run(
        db_path=args.input,
        properties=args.properties,
        ncpu=args.ncpu,
        verbose=args.verbose,
        fetch_batch=args.fetch_batch,
        write_batch=args.write_batch,
        imap_chunk=args.imap_chunk,
    )


if __name__ == '__main__':
    entry_point()
