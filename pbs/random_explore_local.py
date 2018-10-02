# Pavel Polishchuk, 2017

import random
import argparse
import os
import datetime

from crem import mutate_mol

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

from multiprocessing import Pool


def divide(n, iterable):
    # https://github.com/erikrose/more-itertools/blob/dfde6c75eb520605ac6a2d374bcdb72e7ce2d894/more_itertools/more.py
    """Divide the elements from *iterable* into *n* parts, maintaining
    order.
        >>> group_1, group_2 = divide(2, [1, 2, 3, 4, 5, 6])
        >>> list(group_1)
        [1, 2, 3]
        >>> list(group_2)
        [4, 5, 6]
    If the length of *iterable* is not evenly divisible by *n*, then the
    length of the returned iterables will not be identical:
        >>> children = divide(3, [1, 2, 3, 4, 5, 6, 7])
        >>> [list(c) for c in children]
        [[1, 2, 3], [4, 5], [6, 7]]
    If the length of the iterable is smaller than n, then the last returned
    iterables will be empty:
        >>> children = divide(5, [1, 2, 3])
        >>> [list(c) for c in children]
        [[1], [2], [3], [], []]
    This function will exhaust the iterable before returning and may require
    significant storage. If order is not important, see :func:`distribute`,
    which does not first pull the iterable into memory.
    """
    if n < 1:
        raise ValueError('n must be at least 1')

    seq = tuple(iterable)
    q, r = divmod(len(seq), n)

    ret = []
    for i in range(n):
        start = (i * q) + (i if i < r else r)
        stop = ((i + 1) * q) + (i + 1 if i + 1 < r else r)
        ret.append(iter(seq[start:stop]))

    return ret


def process(smi, mol_id, db_fname, radius, min_freq, output_dir, iteration, batch_id):

    res = []

    with open(os.path.join(output_dir, 'out_i%s_j%s.txt' % (str(iteration).zfill(3), str(batch_id).zfill(2))), 'wt') as w:

        w.write('\t'.join(['SMILES', 'ID', 'Parent', 'Iteration', "MW", 'transformation']) + '\n')

        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        iterator = mutate_mol(m, db_fname, radius=radius, min_freq=min_freq,
                              min_size=0, max_size=8,
                              min_rel_size=0, max_rel_size=0.3,
                              min_inc=-1, max_inc=1, replace_cycles=False, ncores=1,
                              max_replacements=25000)
        for k, (new_smi, rxn) in enumerate(iterator):
            new_id = '%i-%i-%i' % (iteration, batch_id, k)
            mw = MolWt(Chem.MolFromSmiles(new_smi))
            if mw <= 500:
                w.write('\t'.join((new_smi,
                                   new_id,
                                   mol_id,
                                   str(iteration),
                                   str(round(mw, 1)),
                                   rxn)) + '\n')
                res.append([new_smi, new_id, mw])

    return res


def process_mp(args):
    return process(*args)


def get_data(input_mols, db_fname, radius, min_freq, output_dir, iteration):
    for j, (smi, id) in enumerate(input_mols):
        yield smi, id, db_fname, radius, min_freq, output_dir, iteration, j


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Stochastically explore chemical space.')
    parser.add_argument('-o', '--out', metavar='OUTPUT_DIR', required=True,
                        help='output dir to store results.')
    parser.add_argument('-d', '--db', metavar='REPLACEMENT _DB', required=True,
                        help='database with replacement.')
    parser.add_argument('-r', '--radius', metavar='INTEGER', required=False, default=1,
                        type=int,
                        help='radius of environment. Default: 1.')
    parser.add_argument('-f', '--min_freq', metavar='INTEGER', required=False, default=0,
                        type=int,
                        help='minimum occurrence of fragments in DB. Default: 0.')
    parser.add_argument('-m', '--nmols', metavar='INTEGER', required=False, default=20,
                        type=int,
                        help='number of compounds selected on each step. Default: 20.')
    parser.add_argument('-s', '--nsteps', metavar='INTEGER', required=False, default=150,
                        type=int,
                        help='number of steps. Default: 150.')

    args = vars(parser.parse_args())
    radius = args['radius']
    min_freq = args['min_freq']
    db_fname = os.path.abspath(args['db'])
    output_dir = os.path.abspath(args['out'])
    nmols = args['nmols']
    niter = args['nsteps']

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    input_mols = [('c1ccccc1', '0')]

    p = Pool(32)

    for i in range(niter):

        uniq_smi = dict()  # smi: (id/name, mw)

        for res in p.imap_unordered(process_mp, get_data(input_mols, db_fname, radius, min_freq, output_dir, i)):
            for smi, new_id, mw in res:
                if smi not in uniq_smi:
                    uniq_smi[smi] = (new_id, mw)


        # for j, (smi, mol_id) in enumerate(reversed(input_mols)):
        #
        #     with open(os.path.join(output_dir, 'out_i%s_j%s.txt' % (str(i).zfill(3), str(j).zfill(2))), 'wt') as w:
        #
        #         w.write('\t'.join(['SMILES', 'ID', 'Parent', 'Iteration', "MW", 'transformation']) + '\n')
        #
        #         m = Chem.AddHs(Chem.MolFromSmiles(smi))
        #         iterator = mutate_mol(m, db_fname, radius=radius, min_freq=min_freq,
        #                               min_size=0, max_size=8,
        #                               min_rel_size=0, max_rel_size=0.3,
        #                               min_inc=-1, max_inc=1, replace_cycles=False, ncores=32,
        #                               max_replacements=25000)
        #         for k, (new_smi, rxn) in enumerate(iterator):
        #             new_id = '%i-%i-%i' % (i, j, k)
        #             mw = MolWt(Chem.MolFromSmiles(new_smi))
        #             if mw <= 500:
        #                 w.write('\t'.join((new_smi,
        #                                    new_id,
        #                                    mol_id,
        #                                    str(i),
        #                                    str(round(mw, 1)),
        #                                    rxn)) + '\n')
        #             if new_smi not in uniq_smi and mw <= 500:
        #                 uniq_smi[new_smi] = (new_id, mw)
        #
        #         del iterator

        if len(uniq_smi) > nmols:
            mw = sorted((v[1], (k, v[0])) for k, v in uniq_smi.items())  # v[0] - name, v[1] - mw
            input_mols = []
            for c in divide(nmols, mw):
                input_mols.append(random.choice(list(c))[1])
        else:
            input_mols = [(k, v[0]) for k, v in uniq_smi.items()]

        print("[%s] Iteration: %i, mols: %i" % (str(datetime.datetime.now()), i + 1, len(uniq_smi)))
