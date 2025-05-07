#!/usr/bin/env python3

import argparse
import sqlite3
import sys
from functools import partial
from multiprocessing import Pool

from crem.arg_types import filepath_type, cpu_type
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds, CalcTPSA, CalcFractionCSP3


props = ['mw', 'logp', 'rtb', 'tpsa', 'fcsp3']

# SQLite DB will be updated by chunks
if sqlite3.sqlite_version_info[:2] <= (3, 32):
    CHUNK_SIZE = 999
else:
    CHUNK_SIZE = 32766


def property_type(x):
    return [item.lower() for item in x if item.lower() in props]


def calc(items, mw=False, logp=False, rtb=False, tpsa=False, fcsp3=False):
    rowid, smi = items
    res = dict()
    mol = Chem.MolFromSmiles(smi)
    if mol:
        if mw:
            res['mw'] = round(MolWt(mol), 2)
        if logp:
            res['logp'] = round(MolLogP(mol), 2)
        if rtb:
            res['rtb'] = CalcNumRotatableBonds(Chem.RemoveHs(mol))
        if tpsa:
            res['tpsa'] = CalcTPSA(mol)
        if fcsp3:
            res['fcsp3'] = round(CalcFractionCSP3(mol), 3)
    upd_str = ','.join(f'{k} = {v}' for k, v in res.items())
    return rowid, upd_str


def entry_point():
    parser = argparse.ArgumentParser(description='Add columns with values of chosen properties to CReM fragment '
                                                 'database.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=filepath_type,
                        help='SQLite DB with CReM fragments.')
    parser.add_argument('-p', '--properties', metavar='NAMES', required=False, nargs='*', default=props, choices=props,
                        help='properties to compute.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    if not args.properties:
        sys.stderr.write(f'No valid names of properties were supplied. Check them please: {", ".join(args.properties)}\n')
        exit()

    pool = Pool(args.ncpu)

    mw = 'mw' in args.properties
    logp = 'logp' in args.properties
    rtb = 'rtb' in args.properties
    tpsa = 'tpsa' in args.properties
    fcsp3 = 'fcsp3' in args.properties

    with sqlite3.connect(args.input) as conn:
        cur = conn.cursor()
        tables = cur.execute("SELECT name FROM sqlite_master WHERE type = 'table' AND name LIKE 'radius%'")
        tables = [i[0] for i in tables]

        for table in tables:
            for prop in args.properties:
                try:
                    cur.execute(f"ALTER TABLE {table} ADD COLUMN {prop} NUMERIC DEFAULT NULL")
                    conn.commit()
                except sqlite3.OperationalError as e:
                    sys.stderr.write(str(e) + '\n')
            sql = f"SELECT rowid, core_smi FROM {table} WHERE " + \
                  " OR ".join([f"{prop} IS NULL" for prop in args.properties])
            cur.execute(sql)
            res = cur.fetchall()

            for i, (rowid, upd_str) in enumerate(pool.imap_unordered(partial(calc, mw=mw, logp=logp, rtb=rtb, tpsa=tpsa, fcsp3=fcsp3), res), 1):
                cur.execute(f"UPDATE {table} SET {upd_str} WHERE rowid = '{rowid}'")
                if i % 10000 == 0:
                    conn.commit()
                    if args.verbose:
                        sys.stderr.write(f'\r{i} fragments processed')
            conn.commit()

            sys.stderr.write(f'\nProperties were successfully added to {args.input}\n')


if __name__ == '__main__':
    entry_point()
