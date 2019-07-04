# Pavel Polishchuk, 2017

import sys
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMMPA
from mol_context import get_canon_context_core, combine_core_env_to_rxn_smarts
from multiprocessing import Pool, cpu_count
import sqlite3
import random
from itertools import product

cycle_pattern = re.compile("[a-zA-Z\]][1-9]+")


def __fragment_mol(mol, radius=3, return_ids=True, keep_stereo=False, protected_ids=None):
    """
    INPUT:
        mol - Mol
        radius - integer, number of bonds to cut context
        keep_stereo - bool, keep or discard information about stereoconfiguration
        protected_ids - set/list/tuple os atom ids which cannot be present in core fragments

    OUTPUT:
        list of tuples (env_smi, core_smi, tuple of core atom ids)
        ('C[*:1].C[*:2]', 'CC(C(=O)O)c1ccc(CC([*:1])[*:2])c(Br)c1', (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
        ('Cc(c)c(cc)[*:1]', '[H][*:1]', (25,))

    If input mol has explicit hydrogens the output will contain also fragments where core = [H][*:1].
    Smiles of fragments with heavy atoms will contain only heavy atoms
    """

    def get_atom_prop(molecule, prop="Index"):
        res = []
        for a in molecule.GetAtoms():
            if a.GetAtomicNum():
                res.append(a.GetIntProp(prop))
        return tuple(sorted(res))

    if protected_ids:
        return_ids = True

    output = []

    # set original atom idx to keep them in fragmented mol
    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    # heavy atoms
    frags = rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4, resultsAsMols=True, maxCutBonds=30)
    # hydrogen atoms
    frags += rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)

    for i, (core, chains) in enumerate(frags):
        if core is None:  # single cut
            components = list(Chem.GetMolFrags(chains, asMols=True))
            ids_0 = get_atom_prop(components[0]) if return_ids else tuple()
            ids_1 = get_atom_prop(components[1]) if return_ids else tuple()
            if Chem.MolToSmiles(components[0]) != '[H][*:1]':  # context cannot be H
                env, frag = get_canon_context_core(components[0], components[1], radius, keep_stereo)
                output.append((env, frag, ids_1))
            if Chem.MolToSmiles(components[1]) != '[H][*:1]':  # context cannot be H
                env, frag = get_canon_context_core(components[1], components[0], radius, keep_stereo)
                output.append((env, frag, ids_0))
        else:   # multiple cuts
            # there are no checks for H needed because H can be present only in single cuts
            env, frag = get_canon_context_core(chains, core, radius, keep_stereo)
            output.append((env, frag, get_atom_prop(core) if return_ids else tuple()))

    if protected_ids:
        protected_ids = set(protected_ids)
        output = [item for item in output if protected_ids.isdisjoint(item[2])]

    return output  # list of tuples (env smiles, core smiles, list of atom ids)


def __fragment_mol_link(mol1, mol2, radius=3, keep_stereo=False, protected_ids_1=None, protected_ids_2=None,
                        return_ids=True):

    def filter_frags(frags, protected_ids):
        output = []
        protected_ids = set(protected_ids)
        for _, chains in frags:
            for atom in chains.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    for d in atom.GetNeighbors():
                        if d.GetAtomicNum() != 1 and d.GetIdx() not in protected_ids:
                            output.append((None, chains))
        return output

    def prep_frags(frags, keep_stereo=False):
        # frags is a list of tuples [(None, frag_mol_1), (None, frag_mol_2), ...]
        ls = []
        for _, chains in frags:
            for atom in chains.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    for d in atom.GetNeighbors():
                        if d.GetAtomicNum() == 1:
                            ids = [d.GetIntProp('Index')]
            a, b = Chem.MolToSmiles(chains, isomericSmiles=keep_stereo).split('.')
            if a == '[H][*:1]':
                ls.append([b, ids])
            else:
                ls.append([a, ids])
        return ls

    if protected_ids_1 or protected_ids_2:
        return_ids = True

    if return_ids:
        for atom in mol1.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())
        for atom in mol2.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    frags_1 = rdMMPA.FragmentMol(mol1, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    frags_2 = rdMMPA.FragmentMol(mol2, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)

    # frags_1 = rdMMPA.FragmentMol(Chem.AddHs(mol1), pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    # frags_2 = rdMMPA.FragmentMol(Chem.AddHs(mol2), pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    # frags_1 = tuple((_, Chem.RemoveHs(m)) for _, m in frags_1)
    # frags_2 = tuple((_, Chem.RemoveHs(m)) for _, m in frags_2)

    if protected_ids_1:
        frags_1 = filter_frags(frags_1, protected_ids_1)

    if protected_ids_2:
        frags_2 = filter_frags(frags_2, protected_ids_2)

    frags_1 = prep_frags(frags_1, keep_stereo)
    frags_2 = prep_frags(frags_2, keep_stereo)

    for i in range(len(frags_1)):
        frags_1[i][0] = frags_1[i][0].replace('*:1', '*:2')

    q = []
    for (fr1, ids1), (fr2, ids2) in product(frags_1, frags_2):
        q.append(['%s.%s' % (fr1, fr2), ids1, ids2])

    fake_core = '[*:1]C[*:2]'
    output = []

    for (chains, ids_1, ids_2) in q:
        env, frag = get_canon_context_core(chains, fake_core, radius=radius, keep_stereo=keep_stereo)
        output.append((env, '[H][*:1].[H][*:2]', ids_1, ids_2))

    return output  # list of tuples (env smiles, core smiles, list of atom ids)


def __frag_replace(mol1, mol2, frag_sma, replace_sma, radius, frag_ids_1=None, frag_ids_2=None):
    """
    INPUT
        mol1:        mol for mutate, grow or link,
        mol2:        for link,
        frag_sma:    SMARTS of a fragment,
        replace_sma: SMARTS of a replacement (from DB),
        radius:      context radius considered,
        frag_ids_1, frag_ids_2: atom ids of a fragment if you need to make exact replacement, if None all possible
                                matches will be replaced
    OUTPUT
        generator returns canonical isomeric SMILES, RDKit Mol and rxn rule
    """

    def set_protected_atoms(mol, ids, radius):

        def extend_ids(mol, atom_id, r, ids):
            if r:
                for a in mol.GetAtomWithIdx(atom_id).GetNeighbors():
                    a_id = a.GetIdx()
                    if a_id not in ids:
                        ids.add(a_id)
                    extend_ids(mol, a_id, r-1, ids)

        if ids:
            ids_ext = set(ids)
            # extend atom ids on neighbour atoms
            for i in ids:
                extend_ids(mol, i, radius + 1, ids_ext)
            # protect untouched atoms
            for a in mol.GetAtoms():
                if a.GetAtomicNum() > 1 and a.GetIdx() not in ids_ext:
                    a.SetProp('_protected', '1')
                else:
                    a.ClearProp('_protected')

    link = False
    if not isinstance(mol1, Chem.Mol):
        raise StopIteration("The first molecule in __gen_replacement always must be specified")
    if isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol):
        link = True

    frag_sma = frag_sma.replace('*', '!#1')  # to avoid map H in mol with explicit H (lead to wrong replacement)
    rxn_sma = "%s>>%s" % (frag_sma, replace_sma)
    rxn = AllChem.ReactionFromSmarts(rxn_sma)

    set_protected_atoms(mol1, frag_ids_1, radius)
    if link:
        set_protected_atoms(mol2, frag_ids_2, radius)

    if link:
        reactants = [mol1, mol2]
    else:
        reactants = [mol1]
    ps = rxn.RunReactants(reactants)

    products = set()
    for y in ps:
        for p in y:
            e = Chem.SanitizeMol(p, catchErrors=True)
            if e:
                sys.stderr.write("Molecule %s caused sanitization error %i" % (Chem.MolToSmiles(p, isomericSmiles=True), e))
                sys.stderr.flush()
            else:
                smi = Chem.MolToSmiles(Chem.RemoveHs(p), isomericSmiles=True)
                if smi not in products:
                    products.add(smi)
                    yield smi, p, rxn_sma


def __get_replacements(db_cur, env, dist, min_atoms, max_atoms, radius, min_freq=0):
    sql = """SELECT core_smi, core_sma, freq
             FROM radius%i
             WHERE env = ? AND 
                   freq >= ? AND
                   core_num_atoms BETWEEN ? AND ?""" % radius
    if isinstance(dist, int):
        sql += " AND dist2 = %i" % dist
    elif isinstance(dist, tuple) and len(dist) == 2:
        sql += " AND dist2 BETWEEN %i AND %i" % dist
    db_cur.execute(sql, (env, min_freq, min_atoms, max_atoms))
    return db_cur.fetchall()


def __gen_replacements(mol1, mol2, db_name, radius, dist=None, min_size=0, max_size=8, min_rel_size=0, max_rel_size=1,
                       min_inc=-2, max_inc=2, max_replacements=None, replace_cycles=False,
                       protected_ids_1=None, protected_ids_2=None, min_freq=10):

    def func():
        for env, core, *ids in f:  # if link = True ids is two tuples, if link = False ids is a single tuple

            num_heavy_atoms = Chem.MolFromSmiles(core).GetNumHeavyAtoms()
            hac_ratio = num_heavy_atoms / mol_hac

            if (min_size <= num_heavy_atoms <= max_size and min_rel_size <= hac_ratio <= max_rel_size) \
                    or (replace_cycles and cycle_pattern.search(core)):

                frag_sma = combine_core_env_to_rxn_smarts(core, env)

                min_atoms = num_heavy_atoms + min_inc
                max_atoms = num_heavy_atoms + max_inc

                rep = __get_replacements(cur, env, dist, min_atoms, max_atoms, radius, min_freq)
                for core_smi, core_sma, freq in rep:
                    if core_smi != core:
                        if link:
                            yield frag_sma, core_sma, freq, ids[0], ids[1]
                        else:
                            yield frag_sma, core_sma, freq, ids[0]

    link = False
    if not isinstance(mol1, Chem.Mol):
        raise StopIteration("The first molecule in __gen_replacement always must be specified")
    if isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol):
        link = True

    if link:
        f = __fragment_mol_link(mol1=mol1, mol2=mol2, radius=radius, protected_ids_1=protected_ids_1,
                                protected_ids_2=protected_ids_2)
        mol = Chem.CombineMols(mol1, mol2)
    else:
        mol = mol1
        f = __fragment_mol(mol, radius, protected_ids=protected_ids_1)

    mol_hac = mol.GetNumHeavyAtoms()

    con = sqlite3.connect(db_name)
    cur = con.cursor()

    if max_replacements is not None:
        res = list(func())
        random.shuffle(res)
        for items in res[:max_replacements]:
            yield items
    else:
        for items in func():
            yield items


def __frag_replace_mp(items):
    # return smi, rxn_smarts, rxn_smarts_freq, mol
    return [(*item, items[-1]) for item in __frag_replace(*items[:-1])]


def __get_data(mol, db_name, radius, min_size, max_size, min_rel_size, max_rel_size, min_inc, max_inc,
               replace_cycles, protected_ids, min_freq, max_replacements):
    for frag_sma, core_sma, freq, ids in __gen_replacements(mol1=mol, mol2=None, db_name=db_name, radius=radius,
                                                            min_size=min_size, max_size=max_size,
                                                            min_rel_size=min_rel_size, max_rel_size=max_rel_size,
                                                            min_inc=min_inc, max_inc=max_inc,
                                                            max_replacements=max_replacements,
                                                            replace_cycles=replace_cycles,
                                                            protected_ids_1=protected_ids, protected_ids_2=None,
                                                            min_freq=min_freq):
        yield mol, None, frag_sma, core_sma, radius, ids, None, freq


def __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms, protected_ids_1, protected_ids_2, min_freq,
                    max_replacements):
    for frag_sma, core_sma, freq, ids_1, ids_2 in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name,
                                                                     radius=radius, dist=dist,
                                                                     min_size=0, max_size=0,
                                                                     min_rel_size=0, max_rel_size=1,
                                                                     min_inc=min_atoms, max_inc=max_atoms,
                                                                     max_replacements=max_replacements,
                                                                     replace_cycles=False,
                                                                     protected_ids_1=protected_ids_1,
                                                                     protected_ids_2=protected_ids_2,
                                                                     min_freq=min_freq):
        yield mol1, mol2, frag_sma, core_sma, radius, ids_1, ids_2, freq


def mutate_mol(mol, db_name, radius=3, min_size=0, max_size=10, min_rel_size=0, max_rel_size=1, min_inc=-2, max_inc=2,
               max_replacements=None, replace_cycles=False, replace_ids=None, protected_ids=None, min_freq=10,
               return_rxn=True, return_rxn_freq=False, return_mol=False, ncores=1):
    """
    Generator of new molecules by replacement of fragments of the supplied molecule with fragments from DB having
    the same chemical context

    :param mol: RDKit Mol object
    :param db_name: path to DB with replacements
    :param radius: radius of context which will be considered for replacement
    :param min_size, max_size: number of heavy atoms in a fragment to replace
    :param min_rel_size, max_rel_size: Relative size of a fragment to the whole mol
                                       (in terms of a number of heavy atoms)
    :param min_inc, max_inc: Relative minimum and maximum size of new fragments which will replace
                             the existed one. -2 and 2 mean that the existed fragment with N atoms will be replaced
                             with fragments from a DB having from N-2 to N+2 atoms.
    :param max_replacements: maximum number of replacements to make. If the number of available replacements is more
                             than the specified threshold only the specified number of randomly chosen replacements
                             will be applied.
    :param replace_cycles: looking for replacement of a fragment containing cycles irrespectively of the fragment size
    :param replace_ids: iterable with atom ids to replace, it has lower priority over protected_ids (replace_ids
                        available in protected_ids would be ignored).
                        Ids of hydrogen atoms (if any) connected to the specified heavy atoms will be automatically
                        labeled as replaceable.
    :param protected_ids: iterable with atom ids which cannot be mutated. If the molecule was supplied with explicit
                          hydrogen the ids of protected hydrogens should be supplied as well, otherwise they will be
                          replaced.
    :param min_freq: minimum occurrence of fragments in DB for replacement
    :param return_rxn: control whether to additionally return rxn of a transformation or return only generated SMILES
    :param return_rxn_freq: return the frequency of a transformation in the DB
    :param return_mol return RDKit Mol object of a generated molecule
    :param ncores: number of cores
    :return: generator over new molecules. Each entry is a list. The first item is SMILES. If return_rxn = True
             the next items is SMARTS of a transformation. If return_rxn_freq = True the next item is frequency of
             this transformation in the DB (will only added if return_rxn = True). The next item is RDKit Mol
             if return_mol == True.
             Only entries with distinct SMILES will be returned.

    Note: supply mol with explicit Hs if H replacement is desired
    """

    products = set()

    protected_ids = set(protected_ids) if protected_ids else set()

    if replace_ids:
        ids = set()
        for i in replace_ids:
            ids.update(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum() == 1)
        ids = set(a.GetIdx() for a in mol.GetAtoms()).difference(ids).difference(replace_ids)  # ids which should be protected
        protected_ids.update(ids)  # since protected_ids has a higher priority add them anyway

    protected_ids = sorted(protected_ids)

    if ncores == 1:

        for frag_sma, core_sma, freq, ids in __gen_replacements(mol1=mol, mol2=None, db_name=db_name, radius=radius,
                                                                min_size=min_size, max_size=max_size,
                                                                min_rel_size=min_rel_size, max_rel_size=max_rel_size,
                                                                min_inc=min_inc, max_inc=max_inc,
                                                                max_replacements=max_replacements,
                                                                replace_cycles=replace_cycles,
                                                                protected_ids_1=protected_ids, protected_ids_2=None,
                                                                min_freq=min_freq):
            for smi, m, rxn in __frag_replace(mol, None, frag_sma, core_sma, radius, ids, None):
                if max_replacements is None or (max_replacements is not None and len(products) < max_replacements):
                    if smi not in products:
                        products.add(smi)
                        res = [smi]
                        if return_rxn:
                            res.append(rxn)
                            if return_rxn_freq:
                                res.append(freq)
                        if return_mol:
                            res.append(m)
                        yield res
    else:

        p = Pool(min(ncores, cpu_count()))
        for items in p.imap(__frag_replace_mp, __get_data(mol, db_name, radius, min_size, max_size, min_rel_size,
                                                          max_rel_size, min_inc, max_inc, replace_cycles,
                                                          protected_ids, min_freq, max_replacements),
                            chunksize=100):
            for smi, m, rxn, freq in items:
                if max_replacements is None or (max_replacements is not None and len(products) < max_replacements):
                    if smi not in products:
                        products.add(smi)
                        res = [smi]
                        if return_rxn:
                            res.append(rxn)
                            if return_rxn_freq:
                                res.append(freq)
                        if return_mol:
                            res.append(m)
                        yield res


def grow_mol(mol, db_name, radius=3, min_atoms=1, max_atoms=2, max_replacements=None, replace_ids=None,
             protected_ids=None, min_freq=10, return_rxn=True, return_rxn_freq=False, ncores=1):
    """
    Replace hydrogens with fragments from the database having the same chemical context

    :param mol: RDKit Mol object
    :param db_name: path to DB with replacements
    :param radius: radius of context which will be considered for replacement
    :param min_atoms: minimum number of atoms in the fragment which will replace H
    :param max_atoms: maximum number of atoms in the fragment which will replace H
    :param max_replacements: maximum number of replacements to make. If the number of available replacements is more
                             than the specified threshold only the specified number of randomly chosen replacements
                             will be applied.
    :param replace_ids: iterable with ids of heavy atom with replaceable Hs or/and ids of H atoms to replace,
                        it has lower priority over protected_ids (replace_ids available in protected_ids would be
                        ignored). If ids of heavy atoms were supplied all attached hydrogens will be replaced.
    :param protected_ids: iterable with ids of heavy atoms at which no H replacement should be made or ids of
                          protected Hs.
                          Ids of all equivalent atoms should be supplied (e.g. to protect meta-position in toluene
                          ids of both carbons in meta-positions should be supplied)
    :param min_freq: minimum occurrence of fragments in DB for replacement
    :param return_rxn: control whether to additionally return rxn of a transformation or return only generated SMILES
    :param ncores: number of cores
    :return: generator over new molecules. Each entry is a tuple of SMILES of a new molecule and
             SMARTS of an applied transformation. Only entries with unique SMILES will be returned.
    """

    m = Chem.AddHs(mol)

    # create the list of ids of protected Hs only would be enough, however in the first case (replace_ids) the full list
    # of protected atom ids is created
    if protected_ids:

        ids = []
        for i in protected_ids:
            if m.GetAtomWithIdx(i).GetAtomicNum() == 1:
                ids.append(i)
            else:
                for a in m.GetAtomWithIdx(i).GetNeighbors():
                    if a.GetAtomicNum() == 1:
                        ids.append(a.GetIdx())
        protected_ids = set(ids)  # ids of protected Hs

    else:
        protected_ids = set()

    if replace_ids:

        ids = set()  # ids if replaceable Hs
        for i in replace_ids:
            if m.GetAtomWithIdx(i).GetAtomicNum() == 1:
                ids.add(i)
            else:
                ids.update(a.GetIdx() for a in m.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum() == 1)
        ids = set(a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 1).difference(ids)  # ids of Hs to protect
        protected_ids.update(ids)  # since protected_ids has a higher priority add them anyway

    return mutate_mol(m, db_name, radius, min_size=0, max_size=0, min_inc=min_atoms, max_inc=max_atoms,
                      max_replacements=max_replacements, replace_ids=None, protected_ids=protected_ids,
                      min_freq=min_freq, return_rxn=return_rxn, return_rxn_freq=return_rxn_freq, ncores=ncores)


def link_mols(mol1, mol2, db_name, radius=3, dist=None, min_atoms=1, max_atoms=2, max_replacements=None,
              replace_ids_1=None, replace_ids_2=None, protected_ids_1=None, protected_ids_2=None,
              min_freq=10, return_rxn=True, return_rxn_freq=False, return_mol=False, ncores=1):
    """
    Link two molecules by a linker from the database

    :param mol1: RDKit Mol object
    :param mol2: RDKit Mol object
    :param db_name: path to DB with replacements
    :param radius: radius of context which will be considered for replacement
    :param dist: topological distance between two attachment points in the fragment which will link molecules.
                 Can be a single integer or a tuple of lower and upper values.
    :param min_atoms: minimum number of heavy atoms in the fragment which will link molecules
    :param max_atoms: maximum number of heavy atoms in the fragment which will link molecules
    :param max_replacements: maximum number of replacements to make. If the number of available replacements is more
                             than the specified threshold only the specified number of randomly chosen replacements
                             will be applied.
    :param replace_ids_1, replace_ids_2: iterable with atom ids to replace, it has lower priority over protected_ids
                                         (replace_ids available in protected_ids would be ignored).
                                         If ids of heavy atoms were supplied the attached hydrogens will be replaced.
    :param protected_ids_1, protected_ids_2: iterable with ids of heavy atoms at which no H replacement should be made
                                             and/or ids of protected H.
                                             Ids of all equivalent atoms should be supplied (e.g. to protect
                                             meta-position in toluene ids of both carbons in meta-positions should
                                             be supplied)
    :param min_freq: minimum occurrence of fragments in DB for replacement
    :param return_rxn: control whether to additionally return rxn of a transformation or return only generated SMILES
    :param return_rxn_freq: return the frequency of a transformation in the DB
    :param return_mol return RDKit Mol object of a generated molecule
    :param ncores: number of cores
    :return: generator over new molecules. Each entry is a list. The first item is SMILES. If return_rxn = True
             the next items is SMARTS of a transformation. If return_rxn_freq = True the next item is frequency of
             this transformation in the DB (will only added if return_rxn = True). The next item is RDKit Mol
             if return_mol == True.
             Only entries with distinct SMILES will be returned.
    """

    def __get_protected_ids(m, replace_ids, protected_ids):
        # the list of ids of heavy atom with protected hydrogens should be returned

        if protected_ids:

            ids = set()
            for i in protected_ids:
                if m.GetAtomWithIdx(i).GetAtomicNum() == 1:
                    ids.update(a.GetIdx() for a in m.GetAtomWithIdx(i).GetNeighbors())
                else:
                    ids.add(i)
            protected_ids = ids

        else:
            protected_ids = set()

        if replace_ids:

            ids = set()
            for i in replace_ids:
                if m.GetAtomWithIdx(i).GetAtomicNum() == 1:
                    ids.update(a.GetIdx() for a in m.GetAtomWithIdx(i).GetNeighbors())
                else:
                    ids.add(i)
            heavy_atom_ids = set(a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1)
            ids = heavy_atom_ids.difference(ids)  # ids of heavy atoms which should be protected
            protected_ids.update(ids)  # since protected_ids has a higher priority add them anyway

        return protected_ids

    products = set()

    mol1 = Chem.AddHs(mol1)
    mol2 = Chem.AddHs(mol2)

    protected_ids_1 = __get_protected_ids(mol1, replace_ids_1, protected_ids_1)
    protected_ids_2 = __get_protected_ids(mol2, replace_ids_2, protected_ids_2)

    if ncores == 1:

        for frag_sma, core_sma, freq, ids_1, ids_2 in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name,
                                                                         radius=radius, dist=dist,
                                                                         min_size=0, max_size=0, min_rel_size=0,
                                                                         max_rel_size=1, min_inc=min_atoms,
                                                                         max_inc=max_atoms, replace_cycles=False,
                                                                         max_replacements=max_replacements,
                                                                         protected_ids_1=protected_ids_1,
                                                                         protected_ids_2=protected_ids_2,
                                                                         min_freq=min_freq):
            for smi, m, rxn in __frag_replace(mol1, mol2, frag_sma, core_sma, radius, ids_1, ids_2):
                if max_replacements is None or (max_replacements is not None and len(products) < max_replacements):
                    if smi not in products:
                        products.add(smi)
                        res = [smi]
                        if return_rxn:
                            res.append(rxn)
                            if return_rxn_freq:
                                res.append(freq)
                        if return_mol:
                            res.append(m)
                        yield res

    else:

        p = Pool(min(ncores, cpu_count()))
        for items in p.imap(__frag_replace_mp, __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms,
                                                               protected_ids_1, protected_ids_2, min_freq,
                                                               max_replacements),
                            chunksize=100):
            for smi, m, rxn, freq in items:
                if max_replacements is None or (max_replacements is not None and len(products) < max_replacements):
                    if smi not in products:
                        products.add(smi)
                        res = [smi]
                        if return_rxn:
                            res.append(rxn)
                            if return_rxn_freq:
                                res.append(freq)
                        if return_mol:
                            res.append(m)
                        yield res
