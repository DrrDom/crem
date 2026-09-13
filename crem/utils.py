import joblib
import numpy as np
import random
import sys

from collections import OrderedDict, defaultdict
from multiprocessing import cpu_count

from crem.crem import CREM_MARKER_PROP, _get_replacements
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .crem import grow_mol2, mutate_mol2


# Marks an atom which must not be altered. Set on the starting molecule and inherited by
# every product which keeps that atom, because atom properties travel through fragmentation
# and product assembly. That is what lets an iteration protect the same positions round after
# round without mapping parent atom ids onto child atom ids - the child arrives protected.
PROTECTED_ATOM_PROP = '__crem_protected'


def __get_child_added_atom_ids(child_mol):
    '''
    Returns ids of atoms which the generation step that produced this molecule inserted.

    crem marks them with `CREM_MARKER_PROP` and clears the markers of earlier calls, so this
    is the fragment added in the last round rather than every fragment ever added.
    '''
    return sorted(a.GetIdx() for a in child_mol.GetAtoms() if a.HasProp(CREM_MARKER_PROP))


def __protect_atoms(mol, atom_ids):
    '''
    Returns a copy of `mol` whose atoms in `atom_ids` are marked as protected.

    :param mol: RDKit Mol
    :param atom_ids: iterable of atom ids to protect
    :return: a new Mol; the input is left untouched
    '''
    mol = Chem.Mol(mol)
    for i in atom_ids:
        mol.GetAtomWithIdx(int(i)).SetBoolProp(PROTECTED_ATOM_PROP, True)
    return mol


def __get_protected_atom_ids(mol):
    '''Returns ids of atoms which carry the protection mark, inherited or newly set.'''
    return sorted(a.GetIdx() for a in mol.GetAtoms() if a.HasProp(PROTECTED_ATOM_PROP))


def __clear_protection(mol):
    '''Removes the protection mark, so molecules handed back look like ordinary output.'''
    for a in mol.GetAtoms():
        if a.HasProp(PROTECTED_ATOM_PROP):
            a.ClearProp(PROTECTED_ATOM_PROP)


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

    :param kwargs: additional keyword arguments forwarded to the underlying generator -
                   ``grow_mol`` in 'scaffold' mode (e.g. ``min_atoms``, ``max_atoms``) and
                   ``mutate_mol`` in 'analogs' mode (e.g. ``min_size``, ``max_size``, ``min_inc``,
                   ``max_inc``, ``replace_cycles``). See those functions for the meaning of each argument.

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

    # The protection is stored on the atoms themselves, so every product inherits the marks of
    # the atoms it kept and no parent-to-child id mapping is needed between iterations.
    start_mols = [__protect_atoms(mol, protected_ids)]
    # to get results ordered by iterations
    generated_mols = OrderedDict()
    n = 0

    for n in range(n_iterations):
        new_mols = ()
        if mode == 'scaffold':
            new_mols = pool(joblib.delayed(grow_mol2)(m, db_name=db_fname,
                                                      protected_ids=__get_protected_atom_ids(m),
                                                      min_freq=min_freq, radius=radius,
                                                      max_replacements=max_replacements,
                                                      return_mol=True, return_rxn=False, return_rxn_freq=False,
                                                      ncores=1 if len(start_mols) > 1 else ncpu, **kwargs)
                            for m in start_mols)
        if mode == 'analogs':
            new_mols = pool(joblib.delayed(mutate_mol2)(m, db_name=db_fname,
                                                        protected_ids=__get_protected_atom_ids(m),
                                                        min_freq=0, radius=radius, max_replacements=max_replacements,
                                                        return_mol=True, return_rxn=False, return_rxn_freq=False,
                                                        ncores=1 if len(start_mols) > 1 else ncpu, **kwargs)
                            for m in start_mols)

        start_mols = []
        for childs in new_mols:
            for items in childs:
                if items[0] not in generated_mols:
                    child = items[1]
                    if protect_added_frag:
                        for i in __get_child_added_atom_ids(child):
                            child.GetAtomWithIdx(i).SetBoolProp(PROTECTED_ATOM_PROP, True)
                    generated_mols[items[0]] = child
                    start_mols.append(child)

        if not start_mols:
            break

    if n + 1 < n_iterations:
        sys.stderr.write(f'INFO. Procedure is finished after {n + 1} iterations instead of {n_iterations}\n')

    for child in generated_mols.values():
        __clear_protection(child)

    if not return_smi:
        return list(generated_mols.values())
    else:
        return list(generated_mols.keys())


def sample_csp3(row_ids, cur, radius, n):
    """
    Performs random selection of fragments proportionally to a squared fraction of sp3 carbon atoms.
    :param row_ids: the list of row ids of fragments to consider
    :param cur: cursor to the fragment database
    :param radius: context radius
    :param n: the number of fragments to select
    :return: the list of row ids of selected fragments
    """
    if n >= len(row_ids):
        return row_ids
    d = defaultdict(list)
    for rowid, core_smi, _, _ in _get_replacements(cur, radius, row_ids):
        d[core_smi].append(rowid)
    smis = list(d.keys())
    values = [rdMolDescriptors.CalcFractionCSP3(Chem.MolFromSmiles(smi)) ** 2 for smi in smis]
    values = [v + 0.01 for v in values]
    values = np.array(values) / sum(values)
    selected_smiles = np.random.choice(smis, n, replace=False, p=values).tolist()
    ids = []
    for smi in selected_smiles:
        ids.extend(d[smi])
    if len(ids) < n:
        ids = random.sample(ids, n)
    return ids


def filter_max_ring_size(row_ids, cur, radius, max_size=6):
    """
    Remove fragments having a ring size greater than a maximum threshold value
    :param row_ids: the list of row ids of fragments to consider
    :param cur: cursor to the fragment database
    :param radius: context radius
    :param max_size: maximum allowed ring size
    :return: the list of row ids of selected fragments
    """
    d = defaultdict(list)
    for rowid, core_smi, _, _ in _get_replacements(cur, radius, row_ids):
        d[core_smi].append(rowid)
    smis = list(d.keys())
    for smi in smis:
        mol = Chem.MolFromSmiles(smi)
        w = mol.GetRingInfo()
        rings = w.AtomRings()
        if rings and max(len(atom_ids) for atom_ids in rings) > max_size:
            del d[smi]
    ids = []
    for v in d.values():
        ids.extend(v)
    return ids


def filter_acyclic_attachment_points(row_ids, cur, radius):
    """
    Keep only fragments where each attachment point [*:n] is attached to an acyclic atom.

    :param row_ids: the list of row ids of fragments to consider
    :param cur: cursor to the fragment database
    :param radius: context radius
    :return: the list of row ids of selected fragments
    """
    d = defaultdict(list)
    for rowid, core_smi, _, _ in _get_replacements(cur, radius, row_ids):
        d[core_smi].append(rowid)

    for smi in list(d.keys()):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            del d[smi]
            continue

        keep = True
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 0:
                continue
            # Attachment dummy atoms in CReM cores should have a single neighbor.
            if not atom.GetNeighbors():
                keep = False
                break
            neighbor = atom.GetNeighbors()[0]
            if neighbor.IsInRing():
                keep = False
                break

        if not keep:
            del d[smi]

    ids = []
    for v in d.values():
        ids.extend(v)
    return ids
