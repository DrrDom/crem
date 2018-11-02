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


def mol_to_smarts(mol, keep_h=True):
    # keep_h - will increase the count of H atoms for atoms with attached hydrogens to create a valid smarts
    # e.g. [H]-[CH2]-[*] -> [H]-[CH3]-[*]

    mol = Chem.Mol(mol)
    mol.UpdatePropertyCache()

    # change the isotope to 42
    for atom in mol.GetAtoms():
        if keep_h:
            s = sum(na.GetAtomicNum() == 1 for na in atom.GetNeighbors())
            if s:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + s)
        atom.SetIsotope(42)

    # print out the smiles - all the atom attributes will be fully specified
    smarts = Chem.MolToSmiles(mol, isomericSmiles=True, allBondsExplicit=True)
    # remove the 42 isotope labels
    smarts = re.sub(r'\[42', "[", smarts)

    return smarts


def smiles_to_smarts(smi, keep_h=True):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % smi)
        return None
    return mol_to_smarts(mol, keep_h)
