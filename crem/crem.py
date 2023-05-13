# Pavel Polishchuk, 2017

import sys
import re
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMMPA
from crem.mol_context import get_canon_context_core, combine_core_env_to_rxn_smarts
from multiprocessing import Pool, cpu_count
import sqlite3
import random
from itertools import product
from crem.mol_context import patt_remove_map

cycle_pattern = re.compile("[a-zA-Z\]][1-9]+")
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
patt_remove_brackets = re.compile('\(\)')


def __extend_output_by_equivalent_atoms(mol, output):
    """
    Generate additional fragments which cover equivalent atoms to extend the output and make replacements for
    equivalent atoms as well

    :param mol:
    :param output:
    :return:
    """

    atom_ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False, includeIsotopes=False))
    tmp = defaultdict(list)
    for i, rank in enumerate(atom_ranks):
        tmp[rank].append(i)
    atom_eq = dict()  # dict of equivalent atoms
    for ids in tmp.values():
        if len(ids) > 1:
            for i in ids:
                atom_eq[i] = [j for j in ids if j != i]

    extended_output = []
    for item in output:
        if all(i in atom_eq.keys() for i in item[2]):  # if all atoms of a fragment have equivalent atoms
            smi = patt_remove_map.sub('', item[1])
            smi = patt_remove_brackets.sub('', smi)
            ids_list = [set(i) for i in mol.GetSubstructMatches(Chem.MolFromSmarts(smi))]
            for ids_matched in ids_list:
                for ids_eq in product(*(atom_eq[i] for i in item[2])):  # enumerate all combinations of equivalent atoms
                    if ids_matched == set(ids_eq):
                        extended_output.append((item[0], item[1], tuple(sorted(ids_eq))))
    return extended_output


def __fragment_mol(mol, radius=3, return_ids=True, keep_stereo=False, protected_ids=None, symmetry_fixes=False):
    """
    INPUT:
        mol - Mol
        radius - integer, number of bonds to cut context
        keep_stereo - bool, keep or discard information about stereoconfiguration
        protected_ids - set/list/tuple os atom ids which cannot be present in core fragments
        symmetry_fixes - if set, then duplicated fragments having different ids will be returned (useful when one
                         wants to alter only small part of a molecule and it is important atoms with which ids will be
                         replaced)

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

    # due to the bug https://github.com/rdkit/rdkit/issues/3040
    # outputs of rdMMPA.FragmentMol calls will contain duplicated fragments
    # their are removed by using this set
    output = set()

    # set original atom idx to keep them in fragmented mol
    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    # heavy atoms
    frags = rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4, resultsAsMols=True, maxCutBonds=30)
    frags += rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=3, resultsAsMols=True, maxCutBonds=30)
    # hydrogen atoms
    frags += rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)

    for i, (core, chains) in enumerate(frags):
        if core is None:  # single cut
            components = list(Chem.GetMolFrags(chains, asMols=True))
            ids_0 = get_atom_prop(components[0]) if return_ids else tuple()
            ids_1 = get_atom_prop(components[1]) if return_ids else tuple()
            if Chem.MolToSmiles(components[0]) != '[H][*:1]':  # context cannot be H
                env, frag = get_canon_context_core(components[0], components[1], radius, keep_stereo)
                output.add((env, frag, ids_1))
            if Chem.MolToSmiles(components[1]) != '[H][*:1]':  # context cannot be H
                env, frag = get_canon_context_core(components[1], components[0], radius, keep_stereo)
                output.add((env, frag, ids_0))
        else:  # multiple cuts
            # there are no checks for H needed because H can be present only in single cuts
            env, frag = get_canon_context_core(chains, core, radius, keep_stereo)
            output.add((env, frag, get_atom_prop(core) if return_ids else tuple()))

    if symmetry_fixes:
        extended_output = __extend_output_by_equivalent_atoms(mol, output)
        if extended_output:
            output.update(extended_output)

    if protected_ids:
        protected_ids = set(protected_ids)
        output = [item for item in output if protected_ids.isdisjoint(item[2])]

    return list(output)  # list of tuples (env smiles, core smiles, list of atom ids)


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
            ids = []
            for atom in chains.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    for d in atom.GetNeighbors():
                        if d.GetAtomicNum() == 1:
                            ids = [d.GetIntProp('Index')]
                if ids:
                    break   # only one such occurrence can be
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
        reactants = [[mol1, mol2], [mol2, mol1]]   # separate parts in SMARTS match only the corresponding molecule
                                                   # in C.O C will match only the first reactnt and O - only the second
    else:
        reactants = [[mol1]]

    products = set()
    for r in reactants:
        ps = rxn.RunReactants(r)
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


def __get_replacements_rowids(db_cur, env, dist, min_atoms, max_atoms, radius, min_freq=0, **kwargs):
    sql = f"""SELECT rowid
              FROM radius{radius}
              WHERE env = '{env}' AND 
                    freq >= {min_freq} AND
                    core_num_atoms BETWEEN {min_atoms} AND {max_atoms}"""
    if isinstance(dist, int):
        sql += f" AND dist2 = {dist}"
    elif isinstance(dist, tuple) and len(dist) == 2:
        sql += f" AND dist2 BETWEEN {dist[0]} AND {dist[1]}"
    for k, v in kwargs.items():
        if isinstance(v, tuple) and len(v) == 2:
            sql += f" AND {k} BETWEEN {v[0]} AND {v[1]}"
        else:
            sql += f" AND {k} = {v}"
    db_cur.execute(sql)
    return set(i[0] for i in db_cur.fetchall())


def __get_replacements(db_cur, radius, row_ids):
    sql = f"""SELECT rowid, core_smi, core_sma, freq
             FROM radius{radius}
             WHERE rowid IN ({','.join(map(str, row_ids))})"""
    db_cur.execute(sql)
    return db_cur.fetchall()


def __gen_replacements(mol1, mol2, db_name, radius, dist=None, min_size=0, max_size=8, min_rel_size=0, max_rel_size=1,
                       min_inc=-2, max_inc=2, max_replacements=None, replace_cycles=False,
                       protected_ids_1=None, protected_ids_2=None, min_freq=10, symmetry_fixes=False, filter_func=None,
                       return_frag_smi_only=False, **kwargs):

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
        f = __fragment_mol(mol, radius, protected_ids=protected_ids_1, symmetry_fixes=symmetry_fixes)

    if f:
        mol_hac = mol.GetNumHeavyAtoms()

        con = sqlite3.connect(db_name)
        cur = con.cursor()

        replacements = dict()  # to store unused   row_id: (frag_sma, core, ids)
        returned_values = 0
        preliminary_return = 0
        if max_replacements is not None:
            random.shuffle(f)
            preliminary_return = max_replacements // len(f)
            if preliminary_return == 0:
                preliminary_return = 1

        for env, core, *ids in f:  # if link = True ids is two tuples, if link = False ids is a single tuple

            num_heavy_atoms = Chem.MolFromSmiles(core).GetNumHeavyAtoms()
            hac_ratio = num_heavy_atoms / mol_hac

            if (min_size <= num_heavy_atoms <= max_size and min_rel_size <= hac_ratio <= max_rel_size) \
                    or (replace_cycles and cycle_pattern.search(core)):

                frag_sma = combine_core_env_to_rxn_smarts(core, env)

                min_atoms = num_heavy_atoms + min_inc
                max_atoms = num_heavy_atoms + max_inc

                row_ids = __get_replacements_rowids(cur, env, dist, min_atoms, max_atoms, radius, min_freq, **kwargs)

                if filter_func:
                    row_ids = filter_func(row_ids, cur, radius)

                if max_replacements is None:
                    res = __get_replacements(cur, radius, row_ids)
                else:
                    selected_row_ids = random.sample(list(row_ids), min(len(row_ids), preliminary_return))
                    row_ids.difference_update(selected_row_ids)
                    replacements.update({i: (frag_sma, core, ids) for i in row_ids})
                    res = __get_replacements(cur, radius, selected_row_ids)

                for row_id, core_smi, core_sma, freq in res:
                    if core_smi != core:
                        if return_frag_smi_only:
                            yield core_smi
                        else:
                            if link:
                                yield frag_sma, core_sma, freq, ids[0], ids[1]
                            else:
                                yield frag_sma, core_sma, freq, ids[0]
                        if max_replacements is not None:
                            returned_values += 1
                            if returned_values >= max_replacements:
                                return

        if max_replacements is not None:
            selected_row_ids = random.sample(list(replacements.keys()), min(len(replacements), max_replacements - returned_values))
            res = __get_replacements(cur, radius, selected_row_ids)
            for row_id, core_smi, core_sma, freq in res:
                if core_smi != replacements[row_id][1]:
                    if return_frag_smi_only:
                        yield core_smi
                    else:
                        if link:
                            yield replacements[row_id][0], core_sma, freq, replacements[row_id][2][0], replacements[row_id][2][1]
                        else:
                            yield replacements[row_id][0], core_sma, freq, replacements[row_id][2][0]


def __frag_replace_mp(items):
    # return smi, rxn_smarts, rxn_smarts_freq, mol
    return [(*item, items[-1]) for item in __frag_replace(*items[:-1])]


def __get_data(mol, db_name, radius, min_size, max_size, min_rel_size, max_rel_size, min_inc, max_inc,
               replace_cycles, protected_ids, min_freq, max_replacements, symmetry_fixes, filter_func=None, **kwargs):
    for frag_sma, core_sma, freq, ids in __gen_replacements(mol1=mol, mol2=None, db_name=db_name, radius=radius,
                                                            min_size=min_size, max_size=max_size,
                                                            min_rel_size=min_rel_size, max_rel_size=max_rel_size,
                                                            min_inc=min_inc, max_inc=max_inc,
                                                            max_replacements=max_replacements,
                                                            replace_cycles=replace_cycles,
                                                            protected_ids_1=protected_ids, protected_ids_2=None,
                                                            min_freq=min_freq, symmetry_fixes=symmetry_fixes,
                                                            filter_func=filter_func, return_frag_smi_only=False,
                                                            **kwargs):
        yield mol, None, frag_sma, core_sma, radius, ids, None, freq


def __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms, protected_ids_1, protected_ids_2, min_freq,
                    max_replacements, filter_func=None, **kwargs):
    for frag_sma, core_sma, freq, ids_1, ids_2 in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name,
                                                                     radius=radius, dist=dist,
                                                                     min_size=0, max_size=0,
                                                                     min_rel_size=0, max_rel_size=1,
                                                                     min_inc=min_atoms, max_inc=max_atoms,
                                                                     max_replacements=max_replacements,
                                                                     replace_cycles=False,
                                                                     protected_ids_1=protected_ids_1,
                                                                     protected_ids_2=protected_ids_2,
                                                                     min_freq=min_freq, filter_func=filter_func,
                                                                     return_frag_smi_only=False, **kwargs):
        yield mol1, mol2, frag_sma, core_sma, radius, ids_1, ids_2, freq


def mutate_mol(mol, db_name, radius=3, min_size=0, max_size=10, min_rel_size=0, max_rel_size=1, min_inc=-2, max_inc=2,
               max_replacements=None, replace_cycles=False, replace_ids=None, protected_ids=None, symmetry_fixes=False,
               min_freq=0, return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
               **kwargs):
    """
    Generator of new molecules by replacement of fragments in the supplied molecule with fragments from DB.

    :param mol: RDKit Mol object. If hydrogens are explicit they will be replaced as well, otherwise not.
    :param db_name: path to DB file with fragment replacements.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param min_size: minimum number of heavy atoms in a fragment to replace. If 0 - hydrogens will be replaced
                     (if they are explicit). Default: 0.
    :param max_size: maximum number of heavy atoms in a fragment to replace. Default: 10.
    :param min_rel_size: minimum relative size of a replaced fragment to the whole molecule
                         (in terms of a number of heavy atoms)
    :param max_rel_size: maximum relative size of a replaced fragment to the whole molecule
                         (in terms of a number of heavy atoms)
    :param min_inc: minimum change of a number of heavy atoms in replacing fragments to a number of heavy atoms in
                    replaced one. Negative value means that the replacing fragments would be smaller than the replaced
                    one on a specified number of heavy atoms. Default: -2.
    :param max_inc: maximum change of a number of heavy atoms in replacing fragments to a number of heavy atoms in
                    replaced one. Default: 2.
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied. Default: None.
    :param replace_cycles: looking for replacement of a fragment containing cycles irrespectively of the fragment size.
                           Default: False.
    :param replace_ids: iterable with atom ids to replace, it has lower priority over `protected_ids` (replace_ids
                        which are present in protected_ids would be protected).
                        Ids of hydrogen atoms (if any) connected to the specified heavy atoms will be automatically
                        labeled as replaceable. Default: None.
    :param protected_ids: iterable with atom ids which will not be mutated. If the molecule was supplied with explicit
                          hydrogen the ids of protected hydrogens should be supplied as well, otherwise they will be
                          replaced.
                          Ids of all equivalent atoms should be supplied (e.g. to protect meta-position in toluene
                          ids of both carbons in meta-positions should be supplied)
                          This argument has a higher priority over `replace_ids`. Default: None.
    :param symmetry_fixes: if set True duplicated fragments with equivalent atoms having different ids will be
                           enumerated. This makes sense if one wants to replace particular atom(s) which have
                           equivalent ones. By default, among equivalent atoms only those with the lowest ids
                           are replaced. This will result in generation of duplicated molecules if several equivalent
                           atoms are selected which will be filtered later nevertheless. So, it is not very reasonable
                           to use this argument and select several equivalent atoms to replace.
                           This solves the issue of rdkit MMPA fragmenter which removes duplicates internally.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param return_rxn: whether to additionally return rxn of a transformation. Default: False.
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB.  Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule.  Default: False.
    :param ncores: number of cores. Default: 1.
    :param filter_func: a function which will filter selected fragments by additional rules
                        (in this way one may add arbitrary selection constrains). The function takes necessary first
                        three arguments: row_ids (list or set of row_ids from the fragment database supplied to
                        the mutate_mol function), cursor of that fragment database and radius (int). This is required
                        access the selected fragments. Other arguments are custom and user-defined.
                        It is the most convenient to define a filtering function, implement specific logic inside and
                        pass it to mutate_mol using functools.partial. The filtering function should return a list/set
                        of row ids which are a subset of the input row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. Name is a name of a column in
                     the radiusX table of the fragment database. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional), RDKit Mol object (optional).
             Only entries with distinct SMILES will be returned.

    Note: supply RDKit Mol object with explicit hydrogens if H replacement is required

    """

    products = {Chem.MolToSmiles(mol)}

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
                                                                min_freq=min_freq, symmetry_fixes=symmetry_fixes,
                                                                filter_func=filter_func, return_frag_smi_only=False,
                                                                **kwargs):
            for smi, m, rxn in __frag_replace(mol, None, frag_sma, core_sma, radius, ids, None):
                if max_replacements is None or len(products) < (max_replacements + 1):  # +1 because we added source mol to output smiles
                    if smi not in products:
                        products.add(smi)
                        res = [smi]
                        if return_rxn:
                            res.append(rxn)
                            if return_rxn_freq:
                                res.append(freq)
                        if return_mol:
                            res.append(m)
                        if len(res) == 1:
                            yield res[0]
                        else:
                            yield res
    else:

        p = Pool(min(ncores, cpu_count()))
        for items in p.imap(__frag_replace_mp, __get_data(mol, db_name, radius, min_size, max_size, min_rel_size,
                                                          max_rel_size, min_inc, max_inc, replace_cycles,
                                                          protected_ids, min_freq, max_replacements, symmetry_fixes,
                                                          filter_func=filter_func, **kwargs),
                            chunksize=100):
            for smi, m, rxn, freq in items:
                if max_replacements is None or len(products) < (max_replacements + 1):  # +1 because we added source mol to output smiles
                    if smi not in products:
                        products.add(smi)
                        res = [smi]
                        if return_rxn:
                            res.append(rxn)
                            if return_rxn_freq:
                                res.append(freq)
                        if return_mol:
                            res.append(m)
                        if len(res) == 1:
                            yield res[0]
                        else:
                            yield res
        p.close()


def grow_mol(mol, db_name, radius=3, min_atoms=1, max_atoms=2, max_replacements=None, replace_ids=None,
             protected_ids=None, symmetry_fixes=False, min_freq=0, return_rxn=False, return_rxn_freq=False,
             return_mol=False, ncores=1, filter_func=None, **kwargs):
    """
    Replace hydrogens with fragments from the database.

    :param mol: RDKit Mol object.
    :param db_name: path to DB file with fragment replacements.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param min_atoms: minimum number of atoms in the fragment which will replace H
    :param max_atoms: maximum number of atoms in the fragment which will replace H
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied. Default: None.
    :param replace_ids: iterable with ids of heavy atom with replaceable Hs or/and ids of H atoms to replace,
                        it has lower priority over `protected_ids` (replace_ids
                        which are present in protected_ids would be protected). Default: None.
    :param protected_ids: iterable with hydrogen atom ids or ids of heavy atoms at which hydrogens will not be replaced.
                          Ids of all equivalent atoms should be supplied (e.g. to protect meta-position in toluene
                          ids of both carbons in meta-positions should be supplied).
                          This argument has a higher priority over `replace_ids`. Default: None.
    :param symmetry_fixes: if set True duplicated fragments with equivalent atoms having different ids will be
                           enumerated. This makes sense if one wants to replace particular atom(s) which have
                           equivalent ones. By default, among equivalent atoms only those with the lowest ids
                           are replaced. This will result in generation of duplicated molecules if several equivalent
                           atoms are selected which will be filtered later nevertheless. So, it is not very reasonable
                           to use this argument and select several equivalent atoms to replace.
                           This solves the issue of rdkit MMPA fragmenter which removes duplicates internally.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param return_rxn: whether to additionally return rxn of a transformation. Default: False.
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB.  Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule.  Default: False.
    :param ncores: number of cores. Default: 1.
    :param filter_func: a function which will filter selected fragments by additional rules
                        (in this way one may add arbitrary selection constrains). The function takes necessary first
                        three arguments: row_ids (list or set of row_ids from the fragment database supplied to
                        the grow_mol function), cursor of that fragment database and radius (int). This is required
                        access the selected fragments. Other arguments are custom and user-defined.
                        It is the most convenient to define a filtering function, implement specific logic inside and
                        pass it to grow_mol using functools.partial. The filtering function should return a list/set
                        of row ids which are a subset of the input row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. Name is a name of a column in
                     the radiusX table of the fragment database. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional), RDKit Mol object (optional).
             Only entries with distinct SMILES will be returned.

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
                      min_freq=min_freq, return_rxn=return_rxn, return_rxn_freq=return_rxn_freq, return_mol=return_mol,
                      ncores=ncores, symmetry_fixes=symmetry_fixes, filter_func=filter_func, **kwargs)


def link_mols(mol1, mol2, db_name, radius=3, dist=None, min_atoms=1, max_atoms=2, max_replacements=None,
              replace_ids_1=None, replace_ids_2=None, protected_ids_1=None, protected_ids_2=None,
              min_freq=0, return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
              **kwargs):
    """
    Link two molecules by a linker from the database.

    :param mol1: the first RDKit Mol object
    :param mol2: the second RDKit Mol object
    :param db_name: path to DB file with fragment replacements.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param dist: topological distance between two attachment points in the fragment which will link molecules.
                 Can be a single integer or a tuple of lower and upper bound values.
    :param min_atoms: minimum number of heavy atoms in the fragment which will link molecules
    :param max_atoms: maximum number of heavy atoms in the fragment which will link molecules
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied. Default: None.
    :param replace_ids_1: iterable with ids of heavy atom of the first molecule with replaceable Hs or/and ids of H
                          atoms to replace,
                          it has lower priority over `protected_ids_1` (replace_ids
                          which are present in protected_ids would be protected). Default: None.
    :param replace_ids_2: iterable with ids of heavy atom of the second molecule with replaceable Hs or/and ids of H
                          atoms to replace,
                          it has lower priority over `protected_ids_2` (replace_ids
                          which are present in protected_ids would be protected). Default: None.
    :param protected_ids_1: iterable with ids of heavy atoms of the first molecule at which no H replacement should
                            be made and/or ids of protected hydrogens.
                            This argument has a higher priority over `replace_ids_1`. Default: None.
    :param protected_ids_2: iterable with ids of heavy atoms of the second molecule at which no H replacement should
                            be made and/or ids of protected hydrogens.
                            This argument has a higher priority over `replace_ids_2`. Default: None.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param return_rxn: whether to additionally return rxn of a transformation. Default: False.
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB.  Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule.  Default: False.
    :param ncores: number of cores. Default: 1.
    :param filter_func: a function which will filter selected fragments by additional rules
                        (in this way one may add arbitrary selection constrains). The function takes necessary first
                        three arguments: row_ids (list or set of row_ids from the fragment database supplied to
                        the link_mols function), cursor of that fragment database and radius (int). This is required
                        access the selected fragments. Other arguments are custom and user-defined.
                        It is the most convenient to define a filtering function, implement specific logic inside and
                        pass it to link_mols using functools.partial. The filtering function should return a list/set
                        of row ids which are a subset of the input row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. Name is a name of a column in
                     the radiusX table of the fragment database. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional), RDKit Mol object (optional).
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
                                                                         min_freq=min_freq, filter_func=filter_func,
                                                                         return_frag_smi_only=False, **kwargs):
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
                        if len(res) == 1:
                            yield res[0]
                        else:
                            yield res

    else:

        p = Pool(min(ncores, cpu_count()))
        for items in p.imap(__frag_replace_mp, __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms,
                                                               protected_ids_1, protected_ids_2, min_freq,
                                                               max_replacements, filter_func=filter_func, **kwargs),
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
                        if len(res) == 1:
                            yield res[0]
                        else:
                            yield res
        p.close()


def mutate_mol2(*args, **kwargs):
    """
    Convenience function which can be used to process molecules in parallel using multiprocessing module.
    It calls mutate_mol which cannot be used directly in multiprocessing because it is a generator

    :param args: positional arguments, the same as in mutate_mol function
    :param kwargs: keyword arguments, the same as in mutate_mol function
    :return: list with output molecules

    """
    return list(mutate_mol(*args, **kwargs))


def grow_mol2(*args, **kwargs):
    """
    Convenience function which can be used to process molecules in parallel using multiprocessing module.
    It calls grow_mol which cannot be used directly in multiprocessing because it is a generator

    :param args: positional arguments, the same as in grow_mol function
    :param kwargs: keyword arguments, the same as in grow_mol function
    :return: list with output molecules

    """
    return list(grow_mol(*args, **kwargs))


def link_mols2(*args, **kwargs):
    """
    Convenience function which can be used to process molecules in parallel using multiprocessing module.
    It calls link_mols which cannot be used directly in multiprocessing because it is a generator

    :param args: positional arguments, the same as in link_mols function
    :param kwargs: keyword arguments, the same as in link_mols function
    :return: list with output molecules

    """
    return list(link_mols(*args, **kwargs))


def get_replacements(mol1, db_name, radius, mol2=None, dist=None, min_size=0, max_size=8, min_rel_size=0,
                     max_rel_size=1, min_inc=-2, max_inc=2, max_replacements=None, replace_cycles=False,
                     protected_ids_1=None, protected_ids_2=None, replace_ids_1=None, replace_ids_2=None, min_freq=0,
                     symmetry_fixes=False, filter_func=None, **kwargs):
    """
    An auxiliary function, which returns smiles of fragments in a DB which satisfy given criteria
    :param mol1: RDKit Mol object
    :param db_name: path to DB file with fragment replacements.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param mol2: a second RDKit Mol object if searching for linking fragments
    :param dist: topological distance between two attachment points in the fragment which will link molecules.
                 Can be a single integer or a tuple of lower and upper bound values.
    :param min_size: minimum number of heavy atoms in a fragment to replace. If 0 - hydrogens will be replaced
                     (if they are explicit).
    :param max_size: maximum number of heavy atoms in a fragment to replace.
    :param min_rel_size: minimum relative size of a replaced fragment to the whole molecule
                         (in terms of a number of heavy atoms)
    :param max_rel_size: maximum relative size of a replaced fragment to the whole molecule
                         (in terms of a number of heavy atoms)
    :param min_inc: minimum change of a number of heavy atoms in replacing fragments to a number of heavy atoms in
                    replaced one. Negative value means that the replacing fragments would be smaller than the replaced
                    one on a specified number of heavy atoms.
    :param max_inc: maximum change of a number of heavy atoms in replacing fragments to a number of heavy atoms in
                    replaced one.
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied.
    :param replace_cycles: looking for replacement of a fragment containing cycles irrespectively of the fragment size.
                           Default: False.

    :param protected_ids_1: iterable with atom ids which will not be mutated in mol1. If the molecule was supplied with
                            explicit hydrogen the ids of protected hydrogens should be supplied as well, otherwise they
                            will be replaced.
                            Ids of all equivalent atoms should be supplied (e.g. to protect meta-position in toluene
                            ids of both carbons in meta-positions should be supplied)
                            This argument has a higher priority over `replace_ids_1`.
    :param protected_ids_2: iterable with atom ids which will not be mutated in mol2. If the molecule was supplied with
                            explicit hydrogen the ids of protected hydrogens should be supplied as well, otherwise they
                            will be replaced.
                            Ids of all equivalent atoms should be supplied (e.g. to protect meta-position in toluene
                            ids of both carbons in meta-positions should be supplied)
                            This argument has a higher priority over `replace_ids_2`.
    :param replace_ids_1: iterable with atom ids to replace in mol1, it has lower priority over `protected_ids`
                          (replace_ids which are present in protected_ids would be protected).
    :param replace_ids_2: iterable with atom ids to replace in mol2, it has lower priority over `protected_ids`
                          (replace_ids which are present in protected_ids would be protected).
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param symmetry_fixes: if set True duplicated fragments with equivalent atoms having different ids will be
                           enumerated. This makes sense if one wants to replace particular atom(s) which have
                           equivalent ones. By default, among equivalent atoms only those with the lowest ids
                           are replaced. This will result in generation of duplicated molecules if several equivalent
                           atoms are selected which will be filtered later nevertheless. So, it is not very reasonable
                           to use this argument and select several equivalent atoms to replace.
                           This solves the issue of rdkit MMPA fragmenter which removes duplicates internally.
    :param filter_func: a function which will filter selected fragments by additional rules
                        (in this way one may add arbitrary selection constrains). The function takes necessary first
                        three arguments: row_ids (list or set of row_ids from the fragment database supplied to
                        the mutate_mol function), cursor of that fragment database and radius (int). This is required
                        access the selected fragments. Other arguments are custom and user-defined.
                        It is the most convenient to define a filtering function, implement specific logic inside and
                        pass it to mutate_mol using functools.partial. The filtering function should return a list/set
                        of row ids which are a subset of the input row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. Name is a name of a column in
                     the radiusX table of the fragment database. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over smiles of fragments in a DB which satisfy given criteria
    """

    protected_ids_1 = set(protected_ids_1) if protected_ids_1 else set()
    if replace_ids_1:
        replace_ids_1 = set(replace_ids_1) if replace_ids_1 else set()
        protected_ids_1 = set(protected_ids_1) | set(range(mol1.GetNumAtoms())).difference(replace_ids_1)
    if isinstance(mol2, Chem.Mol):
        protected_ids_2 = set(protected_ids_2) if protected_ids_2 else set()
        if replace_ids_2:
            replace_ids_2 = set(replace_ids_2) if replace_ids_2 else set()
            protected_ids_2 = set(protected_ids_2) | set(range(mol2.GetNumAtoms())).difference(replace_ids_2)
    else:
        protected_ids_2 = None

    for frag_smi in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name, radius=radius, dist=dist,
                                       min_size=min_size, max_size=max_size, min_rel_size=min_rel_size,
                                       max_rel_size=max_rel_size, min_inc=min_inc, max_inc=max_inc,
                                       max_replacements=max_replacements, replace_cycles=replace_cycles,
                                       protected_ids_1=protected_ids_1, protected_ids_2=protected_ids_2,
                                       min_freq=min_freq, symmetry_fixes=symmetry_fixes, filter_func=filter_func,
                                       return_frag_smi_only=True, **kwargs):
        yield frag_smi
