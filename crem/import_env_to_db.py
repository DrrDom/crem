import argparse
import sqlite3
import sys
import re
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from .mol_context import combine_core_env_to_rxn_smarts

__author__ = 'pavel'


def __calc(env, core):
    sma = combine_core_env_to_rxn_smarts(core, env, False)
    if core.count('*') == 2:
        mol = Chem.MolFromSmiles(core, sanitize=False)
        mat = Chem.GetDistanceMatrix(mol)
        ids = []
        for a in mol.GetAtoms():
            if not a.GetAtomicNum():
                ids.append(a.GetIdx())
        dist2 = mat[ids[0], ids[1]]
    else:
        dist2 = 0
    return sma, dist2


def __calc_mp(items):
    return __calc(*items)


def __get_additional_data(data, pool):
    # data is a list of tuples (env, core)
    if pool:
        res = [items for items in pool.imap(__calc_mp, data, chunksize=100)]
    else:
        res = [__calc(*items) for items in data]
    return res


def main(input_fname, output_fname, radius, counts, ncpu, verbose):

    pool = Pool(min(ncpu, cpu_count())) if ncpu > 1 else None

    table_name = 'radius%i' % radius

    with sqlite3.connect(output_fname) as conn:
        cur = conn.cursor()

        cur.execute("DROP TABLE IF EXISTS %s" % table_name)
        if counts:
            cur.execute("CREATE TABLE %s("
                        "env TEXT NOT NULL, "
                        "core_smi TEXT NOT NULL, "
                        "core_num_atoms INTEGER NOT NULL, "
                        "core_sma TEXT NOT NULL, "
                        "dist2 INTEGER NOT NULL, "
                        "freq INTEGER NOT NULL)" % table_name)
        else:
            cur.execute("CREATE TABLE %s("
                        "env TEXT NOT NULL, "
                        "core_smi TEXT NOT NULL, "
                        "core_num_atoms INTEGER NOT NULL, "
                        "core_sma TEXT NOT NULL,"
                        "dist2 INTEGER NOT NULL)" % table_name)
        conn.commit()

        buf = []
        with open(input_fname) as f:
            for i, line in enumerate(f):
                if counts:
                    tmp = re.split(',| ', line.strip())
                    tmp.append(tmp.pop(0))  # move the first item to the end
                    buf.append(tuple(tmp))
                else:
                    buf.append(tuple(line.strip().split(",")))
                if (i + 1) % 100000 == 0:
                    adata = __get_additional_data((items[:2] for items in buf), pool)
                    if counts:
                        buf = [a[:-1] + b + (a[-1],) for a, b in zip(buf, adata)]
                        cur.executemany("INSERT INTO %s VALUES (?, ?, ?, ?, ?, ?)" % table_name, buf)
                    else:
                        buf = [a + b for a, b in zip(buf, adata)]
                        cur.executemany("INSERT INTO %s VALUES (?, ?, ?, ?, ?)" % table_name, buf)
                    conn.commit()
                    buf = []
                    if verbose:
                        sys.stderr.write("\r%i lines proceed" % (i + 1))

        if buf:
            adata = __get_additional_data((items[:2] for items in buf), pool)
            if counts:
                buf = [a[:-1] + b + (a[-1],) for a, b in zip(buf, adata)]
                cur.executemany("INSERT INTO %s VALUES (?, ?, ?, ?, ?, ?)" % table_name, buf)
            else:
                buf = [a + b for a, b in zip(buf, adata)]
                cur.executemany("INSERT INTO %s VALUES (?, ?, ?, ?, ?)" % table_name, buf)
            conn.commit()

        idx_name = "%s_env_idx" % table_name
        cur.execute("DROP INDEX IF EXISTS %s" % idx_name)
        cur.execute("CREATE INDEX %s ON %s (env)" % (idx_name, table_name))
        conn.commit()

    if pool is not None:
        pool.close()


def entry_point():
    parser = argparse.ArgumentParser(description='Create SQLite DB from a text file containing env_smi, core_smi, '
                                                 'core_atom_num and core_sma.')
    parser.add_argument('-i', '--input', metavar='env_frags.txt', required=True,
                        help='a comma-separated  text file with env_smi, core_smi, core_atom_num and core_sma.')
    parser.add_argument('-o', '--out', metavar='output.db', required=True,
                        help='output SQLite DB file.')
    parser.add_argument('-r', '--radius', metavar='RADIUS', required=True,
                        help='radius of environment. If table for this radius value exists in output DB '
                             'it will be dropped.')
    parser.add_argument('-c', '--counts', action='store_true', default=False,
                        help='set if the input file contains number of occurrences as a first column '
                             '(output of sort | uniq -c). This will add a column freq to the output DB.')
    parser.add_argument('-n', '--ncpu', default=1,
                        help='number of cpus. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "out": output_fname = v
        if o == "verbose": verbose = v
        if o == "radius": radius = int(v)
        if o == "counts": counts = v
        if o == "ncpu": ncpu = int(v)

    main(input_fname=input_fname,
         output_fname=output_fname,
         radius=radius,
         counts=counts,
         ncpu=ncpu,
         verbose=verbose)


if __name__ == '__main__':
    entry_point()
