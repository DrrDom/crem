__author__ = 'pavel'

import argparse
import sys
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from rdkit.Chem import rdMMPA


def fragment_mol(smi, smi_id=''):

    mol = Chem.MolFromSmiles(smi)

    outlines = set()

    if mol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % smi)
    else:
        # heavy atoms
        frags = rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4, resultsAsMols=False, maxCutBonds=30)
        frags += rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=3, resultsAsMols=False, maxCutBonds=30)
        frags = set(frags)
        for core, chains in frags:
            output = '%s,%s,%s,%s\n' % (smi, smi_id, core, chains)
            outlines.add(output)
        # hydrogen splitting
        mol = Chem.AddHs(mol)
        n = mol.GetNumAtoms() - mol.GetNumHeavyAtoms()
        if n < 60:
            frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=False, maxCutBonds=100)
            for core, chains in frags:
                output = '%s,%s,%s,%s\n' % (smi, smi_id, core, chains)
                outlines.add(output)
    return outlines


def process_line(line):
    tmp = line.strip().split(',')
    if len(tmp) == 1:
        return fragment_mol(tmp[0])
    else:
        return fragment_mol(tmp[0], tmp[1])


def main(input_fname, output_fname, ncpu, verbose):

    ncpu = min(cpu_count(), max(ncpu, 1))
    p = Pool(ncpu)

    with open(output_fname, 'wt') as out:

        with open(input_fname) as f:

            for i, res in enumerate(p.imap_unordered(process_line, f, chunksize=100), 1):

                out.write(''.join(res))

                if verbose and i % 1000 == 0:
                    sys.stderr.write('\r%i molecules fragmented' % i)
                    sys.stderr.flush()

    p.close()


def entry_point():
    parser = argparse.ArgumentParser(description='Fragment input compounds by cutting bonds matching bond SMARTS.')
    parser.add_argument('-i', '--input', metavar='input.smi', required=True,
                        help='input SMILES with optional comma-separated ID).')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='fragmented molecules.')
    parser.add_argument('-c', '--ncpu', metavar='NUMBER', required=False, default=1,
                        help='number of cpus used for computation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "out": output_fname = v
        if o == "verbose": verbose = v
        if o == "ncpu": ncpu = int(v)

    main(input_fname=input_fname,
         output_fname=output_fname,
         ncpu=ncpu,
         verbose=verbose)


if __name__ == '__main__':
    entry_point()
