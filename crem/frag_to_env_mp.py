__author__ = 'pavel'

import argparse
import sys

from itertools import permutations
from multiprocessing import Pool, cpu_count
from rdkit import Chem

from .mol_context import get_std_context_core_permutations


def process_line(line):
    # returns env_smi, core_smi, heavy_atoms_num, core_smarts

    output = []
    smi, id, core, context = line.strip().split(',')

    if (not core and not context) or (_keep_mols and id not in _keep_mols):
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
                    if num_heavy_atoms <= _max_heavy_atoms:
                        env, cores = get_std_context_core_permutations(context, core, _radius, _keep_stereo)
                        if env and cores:
                            # for 1 cut cores will always contain 1 item
                            if not _store_comp_id:
                                output.append((env, cores[0], num_heavy_atoms))
                            else:
                                output.append((env, cores[0], num_heavy_atoms, id))
            else:
                sys.stderr.write('more than two fragments in context (%s) where core is empty' % context)
                sys.stderr.flush()
        # two or more splits
        else:
            mm = Chem.MolFromSmiles(core, sanitize=False)
            num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float('inf')
            if num_heavy_atoms <= _max_heavy_atoms:
                env, cores = get_std_context_core_permutations(context, core, _radius, _keep_stereo)
                if env and cores:
                    for c in cores:
                        if not _store_comp_id:
                            output.append((env, c, num_heavy_atoms))
                        else:
                            output.append((env, c, num_heavy_atoms, id))
        return output


def init(keep_mols, radius, keep_stereo, max_heavy_atoms, store_comp_id):
    global _keep_mols
    global _radius
    global _keep_stereo
    global _max_heavy_atoms
    global _store_comp_id
    _keep_mols = set([line.strip() for line in open(keep_mols).readlines()]) if keep_mols else set()
    _radius = radius
    _keep_stereo = keep_stereo
    _max_heavy_atoms = max_heavy_atoms
    _store_comp_id = store_comp_id


def main(input_fname, output_fname, keep_mols, radius, keep_stereo, max_heavy_atoms, ncpu, store_comp_id, verbose):

    # radius and remove_stereo are supplied to process_context_core via global environment (ugly but working solution)

    ncpu = min(cpu_count(), max(ncpu, 1))
    p = Pool(ncpu, initializer=init, initargs=(keep_mols, radius, keep_stereo, max_heavy_atoms, store_comp_id))

    try:
        with open(output_fname, 'wt') as out:

            with open(input_fname) as f:

                for i, res in enumerate(p.imap_unordered(process_line, f, chunksize=1000), 1):

                    for item in res:
                        if item:
                            out.write(','.join(map(str, item)) + '\n')

                    if verbose and i % 1000 == 0:
                        sys.stderr.write('\r%i lines passed' % i)
                        sys.stderr.flush()

    finally:
        p.close()


def entry_point():
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
    parser.add_argument('--store_comp_id', action='store_true', default=False,
                        help='store compound id in output (only for debug).')
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
        if o == "store_comp_id": store_comp_id = v

    main(input_fname=input_fname,
         output_fname=output_fname,
         keep_mols=keep_mols,
         radius=radius,
         keep_stereo=keep_stereo,
         max_heavy_atoms=max_heavy_atoms,
         ncpu=ncpu,
         store_comp_id=store_comp_id,
         verbose=verbose)


if __name__ == '__main__':
    entry_point()
