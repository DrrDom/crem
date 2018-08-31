#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 31-08-2018
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2018
# license         : 
#==============================================================================

import re
import sys
from rdkit import Chem


def mol_to_smarts(mol):

    mol.UpdatePropertyCache()

    # change the isotope to 42
    for atom in mol.GetAtoms():
        atom.SetIsotope(42)

    # print out the smiles - all the atom attributes will be fully specified
    smarts = Chem.MolToSmiles(mol, isomericSmiles=True)
    # remove the 42 isotope labels
    smarts = re.sub(r'\[42', "[", smarts)

    return smarts


def smiles_to_smarts(smi):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % smi)
        return None
    return mol_to_smarts(mol)
