import sys
from collections import OrderedDict
from multiprocessing import cpu_count

import joblib

from .crem import grow_mol2, mutate_mol2


def __get_child_added_atom_ids(child_mol):
    '''
    Returns ids of atoms in a current mol which were added upon grow/mutation procedure
    # After RDKit reaction procedure there is a field <react_atom_idx> with initial parent atom idx in a child mol
    '''
    added_mol_ids = []
    for a in child_mol.GetAtoms():
        if not a.HasProp('react_atom_idx'):
            added_mol_ids.append(a.GetIdx())
    return sorted(added_mol_ids)


def __get_child_protected_atom_ids(mol, protected_parent_ids):
    '''

    :param mol:
    :param protected_parent_ids: ids of a parent molecule which were protected and should be transferred to the
                                 current molecule
    :type  protected_parent_ids: list[int]
    :return: sorted list of integers
    '''
    # After RDKit reaction procedure there is a field <react_atom_idx> with initial parent atom idx in product mol
    protected_mol_ids = []
    for a in mol.GetAtoms():
        if a.HasProp('react_atom_idx') and int(a.GetProp('react_atom_idx')) in protected_parent_ids:
            protected_mol_ids.append(a.GetIdx())
    return sorted(protected_mol_ids)


def __mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol


def enumerate_compounds(mol, db_fname, mode='scaffold', n_iterations=1, radius=3, max_replacements=None,
                        protected_ids=None, replace_ids=None, min_freq=0, protect_added_frag=False, return_smi=False,
                        ncpu=None, **kwargs):
    '''
    Convenience function to perform scaffold decoration or enumeration of analog series. This performs in multiple
    iterations by modification of compounds enumerated on the previous iteration. May result in combinatorial explosion.
    The function returns the list of distinct molecules generated over all iterations.

    Scaffold decoration uses grow procedure. Hydrogens will be added to the supplied molecule and replaced with
    fragments from the database. A user can protect particular positions from expansion by setting protect_ids argument.
    Note: The value of protect_added_frag parameter is ignored. New fragments cannot be attached to previously
    added fragments.

    Enumeration of analog series uses mutate procedure. The molecule should be supplied with explicit hydrogens if one
    wants to replace them as well, otherwise only heavy atoms will be considered for replacements. It is recommended to
    set an altered part of a molecule not too large to obtain reasonable suggestions.

    :param mol:
    :param db_fname: path to DB file with fragment replacements.
    :param mode: 'scaffold' decoration or 'analogs' enumeration. In 'scaffold' mode the supplied molecule will be
                 substituted with fragments from a database. In 'analogs' mode the supplied molecule wil be mutated.
                 Default: scaffold.
    :param n_iterations: the number of rounds of generation. Molecules generated on the previous round are supplied
                         to the next iteration. Be careful setting this parameter as it may cause combinatorial
                         explosion. Default: 1.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param max_replacements: maximum number of replacements to make for each molecule on each iteration. None will
                             result in all possible replacements or it is possible to set a randomly chosen number of
                             replacements. Default: None.
    :param protected_ids: iterable with atom ids which will not be altered. In 'scaffold' mode this can be ids of
                          hydrogens or heavy atoms whose hydrogens should be protected from expansion. In 'analogs' mode
                          these are ids of hydrogens and heavy atoms. Please note, a molecule should be supplied with
                          explicit hydrogens in 'analogs' mode if one wants to replace them. Default: None.
    :param replace_ids: iterable with atom ids to replace, it has lower priority then `protected_ids`. Default: None.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param protect_added_frag: True or False. If set True new fragments cannot be attached/replace fragments added on
                               previous iterations. Applicable only in 'analogs' mode. In 'scaffold' mode user input is
                               ignored and the argument internally set to True. Default: False
    :param return_smi: if True will return the list of SMILES instead of Mol objects. Default: False.
    :param ncpu: number of cores. None means all cpus.

    :param kwargs: these arguments will be passed to grow_mol (for 'scaffold' mode) and mutate_mol ('analogs' mode)
    functions.

    arguments relevant to 'scaffold' mode
    min_atoms: minimum number of atoms in the fragment which will replace H
    max_atoms: maximum number of atoms in the fragment which will replace H

    arguments relevant to 'analogs' mode
    :min_size: minimum number of heavy atoms in a fragment to replace. If 0 - hydrogens will be replaced
                     (if they are explicit). Default: 0.
    :max_size: maximum number of heavy atoms in a fragment to replace. Default: 10.
    :min_inc: minimum change of a number of heavy atoms in replacing fragments to a number of heavy atoms in
                    replaced one. Negative value means that the replacing fragments would be smaller than the replaced
                    one on a specified number of heavy atoms. Default: -2.
    :max_inc: maximum change of a number of heavy atoms in replacing fragments to a number of heavy atoms in
                    replaced one. Default: 2.
    :replace_cycles: looking for replacement of a fragment containing cycles irrespectively of the fragment size.
                    Default: False.

    '''

    if mode not in ['scaffold', 'analogs']:
        raise ValueError('Wrong mode. Please choose one from the list - "analogs","scaffold"')

    if ncpu is None:
        ncpu = cpu_count()
    else:
        ncpu = max(1, min(int(ncpu), cpu_count()))
    pool = joblib.Parallel(n_jobs=ncpu)

    # to check if the statical arguments are in the kwargs dict
    for kw in ['return_mol', 'return_rxn', 'return_rxn_freq', 'ncores']:
        if kw in kwargs:
            kwargs.pop(kw)

    if mode == 'scaffold':
        protect_added_frag = True

    if protected_ids is None and replace_ids is not None:
        protected_ids = list(set(a.GetIdx() for a in mol.GetAtoms()).difference(replace_ids))
    if protected_ids is None:
        protected_ids = []

    start_mols = {mol: protected_ids}
    # to get results ordered by iterations
    generated_mols = OrderedDict()
    n = 0

    for n in range(n_iterations):
        new_mols = ()
        if mode == 'scaffold':
            new_mols = pool(joblib.delayed(grow_mol2)(m, db_name=db_fname, protected_ids=prot_ids,
                                                      min_freq=min_freq, radius=radius,
                                                      max_replacements=max_replacements,
                                                      return_mol=True, return_rxn=False, return_rxn_freq=False,
                                                      ncores=1 if len(start_mols) > 1 else ncpu, **kwargs)
                            for m, prot_ids in start_mols.items())
        if mode == 'analogs':
            new_mols = pool(joblib.delayed(mutate_mol2)(m, db_name=db_fname, protected_ids=prot_ids,
                                                        min_freq=0, radius=radius, max_replacements=max_replacements,
                                                        return_mol=True, return_rxn=False, return_rxn_freq=False,
                                                        ncores=1 if len(start_mols) > 1 else ncpu, **kwargs)
                            for m, prot_ids in start_mols.items())

        parent_protected_ids_list = start_mols.values()
        start_mols = OrderedDict()
        for childs, parent_protected_ids in zip(new_mols, parent_protected_ids_list):
            for items in childs:
                if items[0] not in generated_mols:
                    protected_ids = __get_child_protected_atom_ids(items[1], parent_protected_ids)
                    if protect_added_frag:
                        protect_added_ids = __get_child_added_atom_ids(items[1])
                        protected_ids = set(protected_ids + protect_added_ids)

                    generated_mols[items[0]] = items[1]
                    start_mols[items[1]] = protected_ids

        if not start_mols:
            break

    if n + 1 < n_iterations:
        sys.stderr.write(f'INFO. Procedure is finished after {n + 1} iterations instead of {n_iterations}\n')

    if not return_smi:
        return list(generated_mols.values())
    else:
        return list(generated_mols.keys())
