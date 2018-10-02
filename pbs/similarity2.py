#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 17-09-2018
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2018
# license         : 
#==============================================================================

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from multiprocessing import Pool
from itertools import combinations


def process_line(line):
    smi = line.strip()
    fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2)
    # s = DataStructs.TanimotoSimilarity(ref_fp, fp)
    return fp


def process_pair(fps):
    return DataStructs.TanimotoSimilarity(fps[0], fps[1])


with open(sys.argv[1]) as f:
    ref = f.readline().strip()
    ref_fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(ref), 2)
    # sys.stdout.write('%s\t1.000\n' % ref)

    p = Pool(32)
    res = [r for r in p.imap_unordered(process_line, f)]
    s = sum(p.imap_unordered(process_pair, combinations(res, 2)))
    p.close()
    print((s * 2 + len(res)) / len(res) ** 2)
