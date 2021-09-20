import sys
from collections import OrderedDict

import joblib

from .crem import grow_mol2, mutate_mol2


def __get_child_added_atom_ids(child_mol):
    '''
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
    :param protected_parent_ids: list[int]
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


def enumerate_compounds(mol, db_fname, num_iteration, mode='scaffold', radius=3, max_replacements=None,
                        min_freq=0, protected_ids=None, replase_ids=None, ncpu=1, protect_added_frag=False, **kwargs):
    '''
    :param mol:
    :param db_fname: path to DB file with fragment replacements.
    :param num_iteration: int. The number of rounds of substitution
    :param protected_ids: iterable with atom ids which will not be growed/mutated. If the molecule was supplied with explicit
                          hydrogen the ids of protected hydrogens should be supplied as well, otherwise they will be
                          replaced. Default: None.
    :param replased_ids: iterable with atom ids to replace, it has lower priority over `protected_ids` (replace_ids
                        which are present in protected_ids would be protected).
                        Ids of hydrogen atoms (if any) connected to the specified heavy atoms will be automatically
                        labeled as replaceable. Default: None.
    :param ncpu: number of cores
    :param mode: scaffold or analogs
    *scaffold
    To decorate user-chosen scaffold.
    Scaffold decoration uses grow procedure. Set user-defined atoms to grow/protect.
    Hydrogens will be added and replaced with fragments from the database.
    Note: The value of protect_added_frag parameter is ignored. New fragments cannot be attached to added fragments.

    *analogs
    To enumerate analog series.
    Enumeration of analog series can use grow/mutate procedures.
    There are possible next variants:
    1) To replace only hydrogens and grow molecule set min_size and max_size parameters as 0.
    2) To mutate non-hydrogens not-protected atoms of molecule set min_size and max_size more than 0
    3) To mutate all not-protected atoms of molecule (include hydrogens) set min_size as 0 (Default: 0) and max_size more than 0 (Default: 10).
    Note: Hydrogens should be added by user before running this option. Hydrogens ids should be numerate in protected_ids or replase_ids as well.

    :param protect_added_frag: True or False. If set True new fragments cannot be attached to previously added fragments.
    For scaffold mode protect_added_frag arguments set True and user set will be ignored. Default: False

    *Common grow/mutate arguments*
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied. Default: None.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param protected_ids: iterable with atom ids which will not be mutated. If the molecule was supplied with explicit
                          hydrogen the ids of protected hydrogens should be supplied as well, otherwise they will be
                          replaced.
                          Ids of all equivalent atoms should be supplied (e.g. to protect meta-position in toluene
                          ids of both carbons in meta-positions should be supplied)
                          This argument has a higher priority over `replace_ids`. Default: None.
    :param replace_ids: iterable with atom ids to replace, it has lower priority over `protected_ids` (replace_ids
                        which are present in protected_ids would be protected).
                        Ids of hydrogen atoms (if any) connected to the specified heavy atoms will be automatically
                        labeled as replaceable. Default: None.

    :param kwargs:
    *grow arguments*
    :min_atoms: minimum number of atoms in the fragment which will replace H
    :max_atoms: maximum number of atoms in the fragment which will replace H

    *mutate arguments*
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

    pool = joblib.Parallel(n_jobs=ncpu)

    if mode not in ['scaffold', 'analogs']:
        print('Wrong mode. Please choose one from the list - "analogs","scaffold"')
        return None

    # to check if the statical arguments are in the kwargs dict
    for kw in ['return_mol', 'return_rxn', 'return_rxn_freq', 'ncores']:
        if kw in kwargs:
            kwargs.pop(kw)

    if mode == 'scaffold':
        protect_added_frag = True

    if protected_ids is None and replase_ids is not None:
        protected_ids = list(set(a.GetIdx() for a in mol.GetAtoms()).difference(replase_ids))

    start_mols = {mol: protected_ids}
    # to get results ordered by iterations
    generated_mols = OrderedDict()
    n = 0

    for n in range(num_iteration):
        new_mols = ()
        if mode == 'scaffold':
            new_mols = pool(joblib.delayed(grow_mol2)(m, db_name=db_fname, protected_ids=prot_ids,
                                                      min_freq=min_freq, radius=radius,
                                                      max_replacements=max_replacements,
                                                      return_mol=True, return_rxn=False, return_rxn_freq=False,
                                                      ncores=1 if len(start_mols) != 1 else ncpu, **kwargs)
                            for m, prot_ids in start_mols.items())
        if mode == 'analogs':
            new_mols = pool(joblib.delayed(mutate_mol2)(m, db_name=db_fname, protected_ids=prot_ids,
                                                        min_freq=0, radius=radius, max_replacements=max_replacements,
                                                        return_mol=True, return_rxn=False, return_rxn_freq=False,
                                                        ncores=1 if len(start_mols) != 1 else ncpu, **kwargs)
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

    sys.stderr.write(f'Procedure is finished.\n{n + 1} iterations were completed successfully.\n'
                     f'Totally {len(generated_mols)} new compounds were generated.\n')
    return generated_mols