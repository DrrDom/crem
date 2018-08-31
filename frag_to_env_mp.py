__author__ = 'pavel'

import argparse
import sqlite3
from pprint import pprint
import sys
import re

from itertools import permutations
from multiprocessing import Pool, cpu_count
from rdkit import Chem

from mol_context import get_std_context_core_permutations, combine_core_env_to_rxn_smarts
from functions import smiles_to_smarts


def create_db(conn):
    cur = conn.cursor()
    cur.execute("""CREATE TABLE IF NOT EXISTS core_table
                  (
                    std_context_smi TEXT NOT NULL,
                    std_core_smi TEXT NOT NULL,
                    smart_core_smi TEXT NOT NULL
                  )""")
    cur.execute("DELETE FROM core_table")
    cur.execute("DROP INDEX IF EXISTS std_context_smi_idx")
    cur.execute("DROP INDEX IF EXISTS std_core_smi_idx")
    conn.commit()


def insert_db(cur, items):
    cur.executemany("INSERT INTO core_table VALUES(?,?,?)", items)


def compress_db(conn):
    cur = conn.cursor()
    cur.execute("DELETE FROM core_table WHERE rowid NOT IN (SELECT MIN(rowid) "
                "FROM core_table GROUP BY std_context_smi, std_core_smi)")
    conn.commit()


def index_db(conn):
    cur = conn.cursor()
    cur.execute("CREATE INDEX std_context_smi_idx ON core_table (std_context_smi)")
    cur.execute("CREATE INDEX std_core_smi_idx ON core_table (std_core_smi)")
    conn.commit()


def mol_to_smarts(mol):

    # change the isotope to 42
    for atom in mol.GetAtoms():
        atom.SetIsotope(42)

    # print out the smiles - all the atom attributes will be fully specified
    smarts = Chem.MolToSmiles(mol, isomericSmiles=True)
    # remove the 42 isotope labels
    smarts = re.sub(r'\[42', "[", smarts)
    # now have a fully specified SMARTS - simples!

    return smarts


def process_line(line):
    # returns env_smi, core_smi, heavy_atoms_num, core_smarts

    output = []
    smi, id, core, context = line.strip().split(',')

    if (not core and not context) or (keep_mols_set and id not in keep_mols_set):
        return output
    else:
        # one split
        if not core:
            residues = context.split('.')
            if len(residues) == 2:
                for context, core in permutations(residues, 2):
                    if context == '[H][*:1]':   # ignore such cases
                        continue
                    mm = Chem.MolFromSmiles(core, sanitize=False)
                    num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float('inf')
                    if num_heavy_atoms <= max_heavy_atoms:
                        env, cores = get_std_context_core_permutations(context, core, radius, keep_stereo)
                        if env and cores:
                            # for 1 cut cores will always contain 1 item
                            output.append((env, cores[0], num_heavy_atoms, combine_core_env_to_rxn_smarts(cores[0], env)))
            else:
                sys.stderr.write('more than two fragments in context (%s) where core is empty' % context)
                sys.stderr.flush()
        # two or more splits
        else:
            mm = Chem.MolFromSmiles(core, sanitize=False)
            num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float('inf')
            if num_heavy_atoms <= max_heavy_atoms:
                env, cores = get_std_context_core_permutations(context, core, radius, keep_stereo)
                if env and cores:
                    for c in cores:
                        output.append((env, c, num_heavy_atoms, combine_core_env_to_rxn_smarts(c, env)))
        return output


def init(keep_mols):
    global keep_mols_set
    keep_mols_set = set([line.strip() for line in open(keep_mols).readlines()]) if keep_mols else set()


def main(input_fname, output_fname, keep_mols, radius, keep_stereo, max_heavy_atoms, ncpu, verbose):

    # radius and remove_stereo are supplied to process_context_core via global environment (ugly but working solution)

    # if keep_mols:
    #     keep_mols = set([line.strip() for line in open(keep_mols).readlines()])

    ncpu = min(cpu_count(), max(ncpu, 1))
    p = Pool(ncpu, initializer=init, initargs=(keep_mols,))

    # conn = sqlite3.connect(output_fname)
    # create_db(conn)
    # cur = conn.cursor()
    #
    # cache = []
    #
    # try:
    #
    #     with open(input_fname) as f:
    #
    #         for i, res in enumerate(p.imap_unordered(process_line, f, chunksize=1000), 1):
    #
    #             cache.extend((i for i in res if i))
    #             # insert_db(cur, tuple(i for i in res if i))
    #
    #             if verbose and i % 1000 == 0:
    #                 sys.stderr.write('\r%i lines passed' % i)
    #                 sys.stderr.flush()
    #
    #             if i % 10000:
    #                 insert_db(cur, set(cache))
    #                 cache = []
    #                 conn.commit()
    #
    #             if i > 10000:
    #                 break
    #
    # finally:
    #     insert_db(cur, set(cache))
    #     p.close()
    #     conn.commit()
    #     compress_db(conn)
    #     index_db(conn)

    try:
        with open(output_fname, 'wt') as out:

            with open(input_fname) as f:

                for i, res in enumerate(p.imap_unordered(process_line, f, chunksize=1000), 1):

                    for item in res:
                        if item:
                            out.write('%s,%s,%i,%s\n' % item)

                    if verbose and i % 1000 == 0:
                        sys.stderr.write('\r%i lines passed' % i)
                        sys.stderr.flush()

    finally:
        p.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create text file for fragment replacement from fragmented molecules '
                                                 'obtained with fragmentation.py. '
                                                 'The output may contain duplicated lines which should be filtered out '
                                                 'externally.')
    parser.add_argument('-i', '--input', metavar='frags.txt', required=True,
                        help='fragmented molecules.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file.')
    parser.add_argument('-k', '--keep_mols', metavar='molnames.txt', required=False, default=None,
                        help='file with mol names to keep. Molecules which are not in the list will be ignored.')
    parser.add_argument('-r', '--radius', metavar='NUMBER', required=False, default=1,
                        help='radius of molecular context (in bonds) which will be taken into account. Default: 1.')
    parser.add_argument('-a', '--max_heavy_atoms', metavar='NUMBER', required=False, default=20,
                        help='maximum number of heavy atoms in cores. If the number of atoms exceeds the limit '
                             'fragment will be discarded. Default: 20.')
    parser.add_argument('-s', '--keep_stereo', action='store_true', default=False,
                        help='set this flag if you want to keep stereo in context and core parts.')
    parser.add_argument('-c', '--ncpu', metavar='NUMBER', required=False, default=1,
                        help='number of cpus used for computation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "out": output_fname = v
        if o == "verbose": verbose = v
        if o == "radius": radius = int(v)
        if o == "keep_stereo": keep_stereo = v
        if o == "ncpu": ncpu = int(v)
        if o == "max_heavy_atoms": max_heavy_atoms = int(v)
        if o == "keep_mols": keep_mols = v

    main(input_fname=input_fname,
         output_fname=output_fname,
         keep_mols=keep_mols,
         radius=radius,
         keep_stereo=keep_stereo,
         max_heavy_atoms=max_heavy_atoms,
         ncpu=ncpu,
         verbose=verbose)
