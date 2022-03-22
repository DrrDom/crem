__author__ = 'pavel'

import argparse
import sys
from functools import partial
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from rdkit.Chem import rdMMPA


def fragment_mol(smi, smi_id='', mode=0):

    mol = Chem.MolFromSmiles(smi)

    outlines = set()

    if mol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % smi)
    else:
        # heavy atoms
        if mode == 0 or mode == 1:
            frags = rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4, resultsAsMols=False, maxCutBonds=30)
            frags += rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=3, resultsAsMols=False, maxCutBonds=30)
            frags = set(frags)
            for core, chains in frags:
                output = '%s,%s,%s,%s\n' % (smi, smi_id, core, chains)
                outlines.add(output)
        # hydrogen splitting
        if mode == 1 or mode == 2:
            mol = Chem.AddHs(mol)
            n = mol.GetNumAtoms() - mol.GetNumHeavyAtoms()
            if n < 60:
                frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=False, maxCutBonds=100)
                for core, chains in frags:
                    output = '%s,%s,%s,%s\n' % (smi, smi_id, core, chains)
                    outlines.add(output)
    return outlines


def process_line(line, sep, mode):
    tmp = line.strip().split(sep)
    if tmp:
        if len(tmp) == 1:
            return fragment_mol(tmp[0], mode=mode)
        else:
            return fragment_mol(tmp[0], tmp[1], mode=mode)
    else:
        return None


def main(input_fname, output_fname, mode, sep, ncpu, verbose):

    ncpu = min(cpu_count(), max(ncpu, 1))
    p = Pool(ncpu)

    with open(output_fname, 'wt') as out:

        with open(input_fname) as f:

            for i, res in enumerate(p.imap_unordered(partial(process_line, sep=sep, mode=mode), f, chunksize=100), 1):

                if res:
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
    parser.add_argument('-s', '--sep', metavar='STRING', required=False, default=None,
                        help='separator in input file. Default: Tab.')
    parser.add_argument('-m', '--mode', metavar='INTEGER', required=False, default=0, choices=[0, 1, 2], type=int,
                        help='fragmentation mode: 0 - all atoms constitute a fragment, 1- heavy atoms only, '
                             '2 - hydrogen atoms only. Default: 0.')
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
        if o == "sep": sep = v
        if o == "mode": mode = v

    main(input_fname=input_fname,
         output_fname=output_fname,
         sep=sep,
         mode=mode,
         ncpu=ncpu,
         verbose=verbose)


if __name__ == '__main__':
    entry_point()
