# Pavel Polishchuk, 2017

import os
import sys
import re
from collections import defaultdict
from rdkit import Chem, RDLogger
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMMPA
from crem.mol_context import get_canon_context_core
from multiprocessing import Pool, cpu_count
import sqlite3
import random
from itertools import product, combinations
from crem.mol_context import patt_remove_map

__cycle_pattern = re.compile("[a-zA-Z\]][1-9]+")
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
__patt_remove_brackets = re.compile('\(\)')

__molzip_params = rdmolops.MolzipParams()
__molzip_params.label = rdmolops.MolzipLabel.AtomMapNumber


def __check_db_existence(fname):
    if not os.path.exists(fname):
        raise FileNotFoundError(f'Database {fname} does not exist. ')


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
        core_ids = item[2]
        tail = item[3:] if len(item) > 3 else tuple()
        if all(i in atom_eq.keys() for i in core_ids):  # if all atoms of a fragment have equivalent atoms
            smi = patt_remove_map.sub('', item[1])
            smi = __patt_remove_brackets.sub('', smi)
            ids_list = [set(i) for i in mol.GetSubstructMatches(Chem.MolFromSmarts(smi))]
            for ids_matched in ids_list:
                for ids_eq in product(*(atom_eq[i] for i in core_ids)):  # enumerate all combinations of equivalent atoms
                    if ids_matched == set(ids_eq):
                        extended_output.append((item[0], item[1], tuple(sorted(ids_eq)), *tail))
    return extended_output


def __renumber_attachment_points(mol, old_to_new_map):
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() in old_to_new_map:
            a.SetAtomMapNum(old_to_new_map[a.GetAtomMapNum()])
    return m


def __context_to_std_smi(context_mol, old_to_new_map):
    context_std = __renumber_attachment_points(context_mol, old_to_new_map)
    return Chem.MolToSmiles(context_std, isomericSmiles=True)


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
        list of tuples (env_smi, core_smi, context_mol)
        `context_mol` has dummy atom map numbers aligned to `env_smi`.

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

    protected_ids = set(protected_ids) if protected_ids else set()

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

    def add_output(context_mol, core_mol):
        if Chem.MolToSmiles(context_mol) == '[H][*:1]':  # context cannot be H
            return
        core_ids = get_atom_prop(core_mol) if return_ids else tuple()
        if protected_ids and not protected_ids.isdisjoint(core_ids):
            return
        env, frag, old_to_new_map = get_canon_context_core(context_mol, core_mol, radius, keep_stereo,
                                                           return_att_map=True)
        context_std = __context_to_std_smi(context_mol, old_to_new_map)
        output.add((env, frag, core_ids, context_std))  # context_std should be SMILES, otherwise substantial slow down

    for core, chains in frags:
        if core is None:  # single cut
            components = list(Chem.GetMolFrags(chains, asMols=True))
            add_output(components[0], components[1])
            add_output(components[1], components[0])
        else:  # multiple cuts
            add_output(chains, core)

    if symmetry_fixes:
        extended_output = __extend_output_by_equivalent_atoms(mol, output)
        if extended_output:
            output.update(extended_output)

    res = []
    for env, frag, _, context_mol in output:
        context_mol = Chem.MolFromSmiles(context_mol)
        res.append((env, frag, context_mol))
    return res


def __fragment_mol_link(mol1, mol2, radius=3, keep_stereo=False, protected_ids_1=None, protected_ids_2=None,
                        return_ids=True):

    def _prepare_single_cut_contexts(frags, protected_ids):
        output = []
        for core, chains in frags:
            if core is not None:
                continue
            components = list(Chem.GetMolFrags(chains, asMols=True))
            if len(components) != 2:
                continue
            if Chem.MolToSmiles(components[0]) == '[H][*:1]':
                context = components[1]
            else:
                context = components[0]

            valid = True
            for a in context.GetAtoms():
                if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                    for n in a.GetNeighbors():
                        if n.GetAtomicNum() != 0 and n.HasProp("Index"):
                            anchor_id = n.GetIntProp("Index")
                            if protected_ids and anchor_id in protected_ids:   # TODO: why we need this?
                                valid = False
                            break
            if valid:
                output.append(context)
        return output

    if protected_ids_1 or protected_ids_2:
        return_ids = True

    if return_ids:
        for atom in mol1.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())
        for atom in mol2.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    frags_1 = rdMMPA.FragmentMol(mol1, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    frags_2 = rdMMPA.FragmentMol(mol2, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)

    frags_1 = _prepare_single_cut_contexts(frags_1, protected_ids_1)
    frags_2 = _prepare_single_cut_contexts(frags_2, protected_ids_2)
    fake_core = '[*:1]C[*:2]'
    output = set()

    for ctx_1, ctx_2 in product(frags_1, frags_2):
        # keep historical convention: attachment from mol1 starts with map 2
        ctx_1 = __renumber_attachment_points(ctx_1, {1: 2})

        chains = Chem.CombineMols(ctx_1, ctx_2)
        env, frag, old_to_new_map = get_canon_context_core(chains, fake_core, radius=radius, keep_stereo=keep_stereo,
                                                            return_att_map=True)
        context_smi = __context_to_std_smi(chains, old_to_new_map)
        output.add((env, '[H][*:1].[H][*:2]', context_smi))

    res = []
    for env, core, context_smi in output:
        context_mol = Chem.MolFromSmiles(context_smi)
        res.append((env, core, context_mol))
    return res


def __fragment_mol_macrocycle(mol, radius=3, keep_stereo=False, protected_ids=None, return_ids=True):

    def _prepare_single_cut_contexts(frags, protected_ids):
        output = {}  # anchor_idx -> context_mol
        for core, chains in frags:
            if core is not None:
                continue
            components = list(Chem.GetMolFrags(chains, asMols=True))
            if len(components) != 2:
                continue
            if Chem.MolToSmiles(components[0]) == '[H][*:1]':
                context = components[1]
            else:
                context = components[0]

            anchor_id = None
            valid = True
            for a in context.GetAtoms():
                if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                    for n in a.GetNeighbors():
                        if n.GetAtomicNum() != 0 and n.HasProp("Index"):
                            anchor_id = n.GetIntProp("Index")
                            if anchor_id in protected_ids:
                                valid = False
                            break
                    break
            if valid and anchor_id is not None and anchor_id not in output:
                output[anchor_id] = context
        return sorted(output.items())

    if protected_ids:
        return_ids = True

    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    contexts = _prepare_single_cut_contexts(frags, protected_ids)
    fake_core = '[*:1]C[*:2]'
    output = set()

    for (anchor_1, ctx_1), (anchor_2, ctx_2) in combinations(contexts, 2):
        if anchor_1 == anchor_2:
            continue
        # keep convention consistent with linker generation
        ctx_1 = __renumber_attachment_points(ctx_1, {1: 2})
        chains = Chem.CombineMols(ctx_1, ctx_2)
        env, frag, old_to_new_map = get_canon_context_core(chains, fake_core, radius=radius, keep_stereo=keep_stereo,
                                                            return_att_map=True)

        # transfer attachemnt point from the second chain to the first one to get a single Mol object
        # with two attachment points suitable for molzip
        # replace a hydrogen in the first chain with a dummy atom having the attachement map number from the second part
        anchor_id = None
        att_map = None
        chains = __renumber_attachment_points(chains, old_to_new_map)
        chains = Chem.GetMolFrags(chains, asMols=True)
        for a in chains[1].GetAtoms():
            if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                att_map = a.GetAtomMapNum()
                for n in a.GetNeighbors():
                    if n.GetAtomicNum() != 0 and n.HasProp("Index"):
                        anchor_id = n.GetIntProp("Index")
                        break
                break
        a = chains[0].GetAtomWithIdx(anchor_id)
        for n in a.GetNeighbors():
            if n.GetAtomicNum() == 1:
                n.SetAtomicNum(0)
                n.SetAtomMapNum(att_map)
                break
        context_smi = Chem.MolToSmiles(chains[0], isomericSmiles=True)
        output.add((env, '[H][*:1].[H][*:2]', context_smi))

    res = []
    for env, core, context_smi in output:
        context_mol = Chem.MolFromSmiles(context_smi)
        res.append((env, core, context_mol))
    return res


def __frag_replace(mol1, mol2, old_frag_smi, new_frag_smi, radius, context_mol=None, frag_ids_1=None, frag_ids_2=None,
                   intramol=False):
    """
    INPUT
        mol1:         mol for mutate, grow or link
        mol2:         for link (unused in fast path)
        old_frag_smi: SMILES of a source fragment core
        new_frag_smi: SMILES of a replacement core (from DB)
        radius:       context radius considered
        context_mol:  RDKit Mol containing pre-built whole context part with mapped dummies
    OUTPUT
        generator returns canonical isomeric SMILES, RDKit Mol and rxn rule
    """

    if not isinstance(mol1, Chem.Mol):
        raise StopIteration("The first molecule in __gen_replacement always must be specified")

    if isinstance(context_mol, str):
        context_mol = Chem.MolFromSmiles(context_mol)
    if not isinstance(context_mol, Chem.Mol):
        return

    repl_core_mol = Chem.MolFromSmiles(new_frag_smi)

    transformation_smi = f"{old_frag_smi}>>{new_frag_smi}"
    p = rdmolops.molzip(Chem.CombineMols(context_mol, repl_core_mol), __molzip_params)
    e = Chem.SanitizeMol(p, catchErrors=True)
    if e:
        sys.stderr.write(f"Molecule {Chem.MolToSmiles(p, isomericSmiles=True)} caused "
                         f"sanitization error {e}. "
                         f"Transformation {transformation_smi} was applied to "
                         f"{Chem.MolToSmiles(mol1, isomericSmiles=True)} and "
                         f"{Chem.MolToSmiles(mol2, isomericSmiles=True) if isinstance(mol2, Chem.Mol) else 'None'}\n")
        sys.stderr.flush()
        return

    smi = Chem.MolToSmiles(Chem.RemoveHs(p), isomericSmiles=True)
    yield smi, p, transformation_smi


def __get_replacements_rowids(db_cur, env, dist, min_atoms, max_atoms, radius, min_freq=0, set_names=None, **kwargs):

    user_version = db_cur.execute("PRAGMA user_version").fetchone()[0]

    if user_version == 0:
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
            elif isinstance(v, (int, float, complex)) and not isinstance(v, bool):
                sql += f" AND {k} = {v}"

    elif user_version == 1:
        def _sql_value(value):
            if isinstance(value, str):
                return "'" + value.replace("'", "''") + "'"
            if isinstance(value, bool):
                return "1" if value else "0"
            if value is None:
                return "NULL"
            return str(value)

        frags_columns = {row[1] for row in db_cur.execute("PRAGMA table_info(frags)").fetchall()}
        frags_h_columns = {row[1] for row in db_cur.execute("PRAGMA table_info(frags_h)").fetchall()}

        # Discover available set columns in the radius table.
        # core_num_atoms and dist2 are denormalized into radius{N} alongside
        # env_id / core_smi_id and must be excluded when listing set columns.
        radius_columns = {row[1] for row in db_cur.execute(f"PRAGMA table_info(radius{radius})").fetchall()}
        reserved = {'env_id', 'core_smi_id', 'core_num_atoms', 'dist2'}
        available = sorted(radius_columns - reserved)

        # Normalize set_names: None → all available set columns
        if set_names is None:
            set_names_list = available
        elif isinstance(set_names, str):
            set_names_list = [set_names]
        else:
            set_names_list = list(set_names)

        # Validate requested names
        missing = [sn for sn in set_names_list if sn not in radius_columns]
        if missing:
            raise ValueError(f"Column(s) {missing} not found in radius{radius}. Available set names: {available}")

        # Build OR condition for frequency filter
        freq_clause = " OR ".join(f"r.{sn} >= {_sql_value(min_freq)}" for sn in set_names_list)
        if len(set_names_list) > 1:
            freq_clause = f"({freq_clause})"

        # Decide which side tables we actually need to join. The hot path
        # filters (env, core_num_atoms, dist2, set columns) are all on
        # radius{N} or envs. frags / frags_h are joined only when a kwarg
        # targets one of their columns.
        kwarg_target = {}  # column name -> 'f' or 'h'
        for k in kwargs:
            if k in frags_columns:
                kwarg_target[k] = 'f'
            elif k in frags_h_columns:
                kwarg_target[k] = 'h'
            else:
                raise ValueError(f"Column {k} not found in frags or frags_h")

        need_frags = any(t == 'f' for t in kwarg_target.values())
        need_frags_h = any(t == 'h' for t in kwarg_target.values())

        joins = ["JOIN envs e ON r.env_id = e.env_id"]
        if need_frags or need_frags_h:
            joins.append("JOIN frags f ON r.core_smi_id = f.core_smi_id")
        if need_frags_h:
            joins.append("JOIN frags_h h ON f.core_smi_h_id = h.core_smi_h_id")

        sql = f"""SELECT r.rowid
                  FROM radius{radius} r
                  {' '.join(joins)}
                  WHERE e.env = {_sql_value(env)} AND
                        {freq_clause} AND
                        r.core_num_atoms BETWEEN {_sql_value(min_atoms)} AND {_sql_value(max_atoms)}"""

        if dist is not None:
            if isinstance(dist, tuple):
                if len(dist) != 2:
                    raise ValueError("dist must be a single value or a tuple of two values")
                sql += f" AND r.dist2 BETWEEN {_sql_value(dist[0])} AND {_sql_value(dist[1])}"
            else:
                sql += f" AND r.dist2 = {_sql_value(dist)}"

        for k, v in kwargs.items():
            column = f"{kwarg_target[k]}.{k}"

            if isinstance(v, tuple):
                if len(v) != 2:
                    raise ValueError(f"Value for {k} must be a single value or a tuple of two values")
                sql += f" AND {column} BETWEEN {_sql_value(v[0])} AND {_sql_value(v[1])}"
            else:
                sql += f" AND {column} = {_sql_value(v)}"
    else:
        raise NotImplementedError('Not implemented for database version other than 0 and 1')
    db_cur.execute(sql)
    return set(i[0] for i in db_cur.fetchall())


def _get_replacements(db_cur, radius, row_ids):
    user_version = db_cur.execute("PRAGMA user_version").fetchone()[0]
    if user_version == 0:
        sql = f"""SELECT rowid, core_smi, core_sma, freq
                      FROM radius{radius}
                      WHERE rowid IN ({','.join(map(str, row_ids))})"""
    elif user_version == 1:
        # Note: freq was removed from DB, therefore 0 is returned (maybe None is better)
        sql = f"""SELECT r.rowid, f.core_smi, e.env
                  FROM radius{radius} r
                  JOIN frags f ON r.core_smi_id = f.core_smi_id
                  JOIN envs e ON r.env_id = e.env_id
                  WHERE r.rowid IN ({','.join(map(str, row_ids))})"""
    else:
        raise NotImplementedError('Not implemented for database version other than 0 and 1')
    db_cur.execute(sql)
    if user_version == 1:
        # Keep tuple shape identical to user_version 0 for compatibility.
        return [(row_id, core_smi, core_smi, 0) for row_id, core_smi, env in db_cur.fetchall()]
    return db_cur.fetchall()


def __gen_replacements(mol1, mol2, db_name, radius, dist=None, min_size=0, max_size=8, min_rel_size=0, max_rel_size=1,
                       min_inc=-2, max_inc=2, max_replacements=None, replace_cycles=False,
                       protected_ids_1=None, protected_ids_2=None, min_freq=10, set_names=None,
                       symmetry_fixes=False, filter_func=None, sample_func=None, return_frag_smi_only=False,
                       macrocycle=False, **kwargs):

    link = False
    if not isinstance(mol1, Chem.Mol):
        raise StopIteration("The first molecule in __gen_replacement always must be specified")
    if isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol):
        link = True

    if macrocycle:
        if link:
            raise ValueError("macrocycle mode expects a single molecule")
        mol = mol1
        f = __fragment_mol_macrocycle(mol=mol, radius=radius, protected_ids=protected_ids_1)
    elif link:
        f = __fragment_mol_link(mol1=mol1, mol2=mol2, radius=radius, protected_ids_1=protected_ids_1,
                                protected_ids_2=protected_ids_2)
        mol = Chem.CombineMols(mol1, mol2)
    else:
        mol = mol1
        f = __fragment_mol(mol, radius, protected_ids=protected_ids_1, symmetry_fixes=symmetry_fixes)

    if not f:
        return

    mol_hac = mol.GetNumHeavyAtoms()

    with sqlite3.connect(db_name) as con:
        # Read-side tuning. The fragment DB is queried, never written, in this
        # path. mmap_size lets SQLite read pages directly from the kernel page
        # cache without copying into its own cache, which is the dominant win
        # on Linux for large DBs. The other PRAGMAs are cheap and reduce
        # per-call overhead for the many short queries this function issues.
        con.execute("PRAGMA query_only = ON")
        con.execute("PRAGMA mmap_size = 2147483648")   # 2 GiB
        con.execute("PRAGMA cache_size = -262144")     # 256 MiB page cache
        con.execute("PRAGMA temp_store = MEMORY")
        cur = con.cursor()

        replacements = dict()  # to store unused row_id: (source_core_smi, context_mol)
        returned_values = 0
        preliminary_return = 0
        if max_replacements is not None:
            random.shuffle(f)
            preliminary_return = max_replacements // len(f)
            if preliminary_return == 0:
                preliminary_return = 1

        for env, core, context_mol in f:

            RDLogger.DisableLog('rdApp.warning')
            num_heavy_atoms = Chem.MolFromSmiles(core).GetNumHeavyAtoms()
            hac_ratio = num_heavy_atoms / mol_hac
            RDLogger.EnableLog('rdApp.warning')

            if (min_size <= num_heavy_atoms <= max_size and min_rel_size <= hac_ratio <= max_rel_size) \
                    or (replace_cycles and __cycle_pattern.search(core)):

                min_atoms = num_heavy_atoms + min_inc
                max_atoms = num_heavy_atoms + max_inc

                row_ids = __get_replacements_rowids(cur, env, dist, min_atoms, max_atoms, radius, min_freq,
                                                    set_names=set_names, **kwargs)

                if filter_func:
                    row_ids = set(filter_func(row_ids, cur, radius))

                if max_replacements is None:
                    res = _get_replacements(cur, radius, row_ids)
                else:
                    n = min(len(row_ids), preliminary_return)
                    if sample_func is not None:
                        selected_row_ids = sample_func(list(row_ids), cur, radius, n)
                    else:
                        selected_row_ids = random.sample(list(row_ids), n)
                    row_ids.difference_update(selected_row_ids)
                    replacements.update({i: (core, context_mol) for i in row_ids})
                    res = _get_replacements(cur, radius, selected_row_ids)

                for row_id, core_smi, _, freq in res:
                    if core_smi != core:
                        if return_frag_smi_only:
                            yield core_smi
                        else:
                            yield core, core_smi, freq, context_mol
                        if max_replacements is not None:
                            returned_values += 1
                            if returned_values >= max_replacements:
                                return

        if max_replacements is not None:
            n = min(len(replacements), max_replacements - returned_values)
            if sample_func is not None:
                selected_row_ids = sample_func(list(replacements.keys()), cur, radius, n)
            else:
                selected_row_ids = random.sample(list(replacements.keys()), n)
            res = _get_replacements(cur, radius, selected_row_ids)
            for row_id, core_smi, _, freq in res:
                src_core, src_context = replacements[row_id]
                if core_smi != src_core:
                    if return_frag_smi_only:
                        yield core_smi
                    else:
                        yield src_core, core_smi, freq, src_context


def __frag_replace_mp(items):
    # return smi, transformation, transformation_freq, mol
    if len(items) >= 2 and isinstance(items[-1], bool):
        freq = items[-2]
        intramol = items[-1]
        args = items[:-2]
    else:
        freq = items[-1]
        intramol = False
        args = items[:-1]
    return [(*item, freq) for item in __frag_replace(*args, intramol=intramol)]


def __get_data(mol, db_name, radius, min_size, max_size, min_rel_size, max_rel_size, min_inc, max_inc,
               replace_cycles, protected_ids, min_freq, set_names, max_replacements, symmetry_fixes, filter_func=None,
               sample_func=None, **kwargs):
    for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol, mol2=None, db_name=db_name,
                                                                     radius=radius, min_size=min_size,
                                                                     max_size=max_size, min_rel_size=min_rel_size,
                                                                     max_rel_size=max_rel_size, min_inc=min_inc,
                                                                     max_inc=max_inc,
                                                                     max_replacements=max_replacements,
                                                                     replace_cycles=replace_cycles,
                                                                     protected_ids_1=protected_ids,
                                                                     protected_ids_2=None, min_freq=min_freq,
                                                                     set_names=set_names,
                                                                     symmetry_fixes=symmetry_fixes,
                                                                     filter_func=filter_func,
                                                                     sample_func=sample_func,
                                                                     return_frag_smi_only=False, **kwargs):
        yield mol, None, frag_sma, core_sma, radius, context_mol, freq


def __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms, protected_ids_1, protected_ids_2, min_freq,
                    set_names, max_replacements, filter_func=None, sample_func=None, **kwargs):
    for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name,
                                                                     radius=radius, dist=dist,
                                                                     min_size=0, max_size=0,
                                                                     min_rel_size=0, max_rel_size=1,
                                                                     min_inc=min_atoms, max_inc=max_atoms,
                                                                     max_replacements=max_replacements,
                                                                     replace_cycles=False,
                                                                     protected_ids_1=protected_ids_1,
                                                                     protected_ids_2=protected_ids_2,
                                                                     min_freq=min_freq, set_names=set_names,
                                                                     filter_func=filter_func,
                                                                     sample_func=sample_func,
                                                                     return_frag_smi_only=False, **kwargs):
        yield mol1, mol2, frag_sma, core_sma, radius, context_mol, freq


def __get_data_macrocycle(mol, db_name, radius, dist, min_size, max_size, protected_ids, min_freq, set_names,
                          max_replacements, filter_func=None, sample_func=None, **kwargs):
    for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol, mol2=None, db_name=db_name,
                                                                     radius=radius, dist=dist,
                                                                     min_size=0, max_size=0,
                                                                     min_rel_size=0, max_rel_size=1,
                                                                     min_inc=min_size, max_inc=max_size,
                                                                     max_replacements=max_replacements,
                                                                     replace_cycles=False,
                                                                     protected_ids_1=protected_ids,
                                                                     protected_ids_2=None,
                                                                     min_freq=min_freq, set_names=set_names,
                                                                     filter_func=filter_func,
                                                                     sample_func=sample_func,
                                                                     return_frag_smi_only=False,
                                                                     macrocycle=True, **kwargs):
        yield mol, None, frag_sma, core_sma, radius, context_mol, freq


def mutate_mol(mol, db_name, radius=3, min_size=0, max_size=10, min_rel_size=0, max_rel_size=1, min_inc=-2, max_inc=2,
               max_replacements=None, replace_cycles=False, replace_ids=None, protected_ids=None, symmetry_fixes=False,
               min_freq=0, return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
               sample_func=None, set_names=None, **kwargs):
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
    :param set_names: column name or list of column names in radius tables defining the set(s) of fragments to use
                      with min_freq in v1 databases. A fragment is included if at least one of the named sets
                      satisfies the min_freq threshold (OR logic). If None (default), all available set columns
                      are used. If a column name is not found, a ValueError is raised listing available set names.
                      Ignored for v0 databases. Default: None.
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
    :param sample_func: a function which will sample selected fragments if max_replacements is supplied. If omitted
                        uniform sampling will be used. The function takes necessary first four arguments: row_ids
                        (list or set of row_ids from the fragment database), cursor of that fragment database,
                        radius (int) and the number of returned items (int). This is required to access the selected
                        fragments. Other arguments can be custom and user-defined. The function should return
                        a list/set of selected row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
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

    __check_db_existence(db_name)
    products = {Chem.MolToSmiles(Chem.RemoveHs(mol))}

    protected_ids = set(protected_ids) if protected_ids else set()

    if replace_ids:
        ids = set()
        for i in replace_ids:
            ids.update(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum() == 1)
        ids = set(a.GetIdx() for a in mol.GetAtoms()).difference(ids).difference(replace_ids)  # ids which should be protected
        protected_ids.update(ids)  # since protected_ids has a higher priority add them anyway

    # protected_ids = sorted(protected_ids)  # why we made sorted?

    if ncores == 1:

        for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol, mol2=None, db_name=db_name,
                                                                         radius=radius, min_size=min_size,
                                                                         max_size=max_size,
                                                                         min_rel_size=min_rel_size,
                                                                         max_rel_size=max_rel_size,
                                                                         min_inc=min_inc, max_inc=max_inc,
                                                                         max_replacements=max_replacements,
                                                                         replace_cycles=replace_cycles,
                                                                         protected_ids_1=protected_ids,
                                                                         protected_ids_2=None, min_freq=min_freq,
                                                                         set_names=set_names,
                                                                         symmetry_fixes=symmetry_fixes,
                                                                         filter_func=filter_func,
                                                                         sample_func=sample_func,
                                                                         return_frag_smi_only=False, **kwargs):
            for smi, m, rxn in __frag_replace(mol, None, frag_sma, core_sma, radius, context_mol):
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
        try:
            for items in p.imap(__frag_replace_mp, __get_data(mol, db_name, radius, min_size, max_size, min_rel_size,
                                                              max_rel_size, min_inc, max_inc, replace_cycles,
                                                              protected_ids, min_freq, set_names, max_replacements,
                                                              symmetry_fixes, filter_func=filter_func,
                                                              sample_func=sample_func, **kwargs),
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
        finally:
            p.close()
            p.join()


def grow_mol(mol, db_name, radius=3, min_atoms=1, max_atoms=2, max_replacements=None, replace_ids=None,
             protected_ids=None, symmetry_fixes=False, min_freq=0, return_rxn=False, return_rxn_freq=False,
             return_mol=False, ncores=1, filter_func=None, sample_func=None, set_names=None, **kwargs):
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
    :param symmetry_fixes: if Sset True duplicated fragments with equivalent atoms having different ids will be
                           enumerated. This makes sense if one wants to replace particular atom(s) which have
                           equivalent ones. By default, among equivalent atoms only those with the lowest ids
                           are replaced. This will result in generation of duplicated molecules if several equivalent
                           atoms are selected which will be filtered later nevertheless. So, it is not very reasonable
                           to use this argument and select several equivalent atoms to replace.
                           This solves the issue of rdkit MMPA fragmenter which removes duplicates internally.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param set_names: column name or list of column names in radius tables defining the set(s) of fragments to use
                      with min_freq in v1 databases. A fragment is included if at least one of the named sets
                      satisfies the min_freq threshold (OR logic). If None (default), all available set columns
                      are used. If a column name is not found, a ValueError is raised listing available set names.
                      Ignored for v0 databases. Default: None.
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
    :param sample_func: a function which will sample selected fragments if max_replacements is supplied. If omitted
                        uniform sampling will be used. The function takes necessary first four arguments: row_ids
                        (list or set of row_ids from the fragment database), cursor of that fragment database,
                        radius (int) and the number of returned items (int). This is required to access the selected
                        fragments. Other arguments can be custom and user-defined. The function should return
                        a list/set of selected row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional), RDKit Mol object (optional).
             Only entries with distinct SMILES will be returned.

    """

    __check_db_existence(db_name)
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
                      min_freq=min_freq, set_names=set_names, return_rxn=return_rxn, return_rxn_freq=return_rxn_freq,
                      return_mol=return_mol, ncores=ncores, symmetry_fixes=symmetry_fixes, filter_func=filter_func,
                      sample_func=sample_func, **kwargs)


def link_mols(mol1, mol2, db_name, radius=3, dist=None, min_atoms=1, max_atoms=2, max_replacements=None,
              replace_ids_1=None, replace_ids_2=None, protected_ids_1=None, protected_ids_2=None,
              min_freq=0, return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
              sample_func=None, set_names=None, **kwargs):
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
    :param set_names: column name or list of column names in radius tables defining the set(s) of fragments to use
                      with min_freq in v1 databases. A fragment is included if at least one of the named sets
                      satisfies the min_freq threshold (OR logic). If None (default), all available set columns
                      are used. If a column name is not found, a ValueError is raised listing available set names.
                      Ignored for v0 databases. Default: None.
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
    :param sample_func: a function which will sample selected fragments if max_replacements is supplied. If omitted
                        uniform sampling will be used. The function takes necessary first four arguments: row_ids
                        (list or set of row_ids from the fragment database), cursor of that fragment database,
                        radius (int) and the number of returned items (int). This is required to access the selected
                        fragments. Other arguments can be custom and user-defined. The function should return
                        a list/set of selected row ids.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
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

    __check_db_existence(db_name)
    products = set()

    mol1 = Chem.AddHs(mol1)
    mol2 = Chem.AddHs(mol2)

    protected_ids_1 = __get_protected_ids(mol1, replace_ids_1, protected_ids_1)
    protected_ids_2 = __get_protected_ids(mol2, replace_ids_2, protected_ids_2)

    if ncores == 1:

        for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol1, mol2=mol2,
                                                                         db_name=db_name, radius=radius,
                                                                         dist=dist, min_size=0,
                                                                         max_size=0, min_rel_size=0,
                                                                         max_rel_size=1,
                                                                         min_inc=min_atoms,
                                                                         max_inc=max_atoms,
                                                                         replace_cycles=False,
                                                                         max_replacements=max_replacements,
                                                                         protected_ids_1=protected_ids_1,
                                                                         protected_ids_2=protected_ids_2,
                                                                         min_freq=min_freq,
                                                                         set_names=set_names,
                                                                         filter_func=filter_func,
                                                                         sample_func=sample_func,
                                                                         return_frag_smi_only=False, **kwargs):
            for smi, m, rxn in __frag_replace(mol1, mol2, frag_sma, core_sma, radius, context_mol):
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
        try:
            for items in p.imap(__frag_replace_mp, __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms,
                                                                   protected_ids_1, protected_ids_2, min_freq,
                                                                   set_names, max_replacements, filter_func=filter_func,
                                                                   sample_func=sample_func, **kwargs),
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
        finally:
            p.close()
            p.join()


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
                     symmetry_fixes=False, filter_func=None, sample_func=None, return_frag_smi_only=True,
                     set_names=None, **kwargs):
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
    :param set_names: column name or list of column names in radius tables defining the set(s) of fragments to use
                      with min_freq in v1 databases. A fragment is included if at least one of the named sets
                      satisfies the min_freq threshold (OR logic). If None (default), all available set columns
                      are used. If a column name is not found, a ValueError is raised listing available set names.
                      Ignored for v0 databases. Default: None.
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
    :param sample_func: a function which will sample selected fragments if max_replacements is supplied. If omitted
                        uniform sampling will be used. The function takes necessary first four arguments: row_ids
                        (list or set of row_ids from the fragment database), cursor of that fragment database,
                        radius (int) and the number of returned items (int). This is required to access the selected
                        fragments. Other arguments can be custom and user-defined. The function should return
                        a list/set of selected row ids.
    :param return_frag_smi_only: control whether to return only SMILES of fragments selected from a database or return
                                 a tuple `(source_core_smi, replacement_core_smi, freq, context_mol)` which can be
                                 further passed to `get_mols_from_replacements`.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
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

    for res in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name, radius=radius, dist=dist,
                                       min_size=min_size, max_size=max_size, min_rel_size=min_rel_size,
                                       max_rel_size=max_rel_size, min_inc=min_inc, max_inc=max_inc,
                                       max_replacements=max_replacements, replace_cycles=replace_cycles,
                                       protected_ids_1=protected_ids_1, protected_ids_2=protected_ids_2,
                                       min_freq=min_freq, set_names=set_names, symmetry_fixes=symmetry_fixes,
                                       filter_func=filter_func, sample_func=sample_func,
                                       return_frag_smi_only=return_frag_smi_only, **kwargs):
        yield res


def get_mols_from_replacements(mol1, radius, replacements, mol2=None, return_rxn=False, return_rxn_freq=False,
                               return_mol=False):

    if isinstance(mol2, Chem.Mol):
        products = set()
    else:
        products = {Chem.MolToSmiles(Chem.RemoveHs(mol1))}

    for items in replacements:

        if len(items) == 4:
            frag_sma, core_sma, freq, context_mol = items
        else:
            raise ValueError('Each replacement tuple should have 4 items: '
                             '(source_core_smi, replacement_core_smi, freq, context_mol)\n')

        for smi, m, rxn in __frag_replace(mol1, mol2, frag_sma, core_sma, radius, context_mol):
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


def make_macrocycle(mol, db_name, radius=3, dist=None, min_atoms=1, max_atoms=10, max_replacements=None,
                    replace_cycles=False, replace_ids=None, protected_ids=None, symmetry_fixes=False, min_freq=0,
                    return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
                    sample_func=None, set_names=None, **kwargs):
    """
    Generate macrocycles by linking two atoms in the same molecule with a linker from DB.

    :param mol: RDKit Mol object.
    :param db_name: path to DB file with fragment replacements.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param dist: topological distance between two attachment points in the linking fragment.
                 Can be a single integer or a tuple of lower and upper bound values.
    :param min_atoms: minimum number of heavy atoms in the linker. Default: 1.
    :param max_atoms: maximum number of heavy atoms in the linker. Default: 10.
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied. Default: None.
    :param replace_cycles: accepted for API compatibility with mutate/grow functions but not used for macrocyclization.
    :param replace_ids: iterable with ids of heavy atom with replaceable Hs or/and ids of H atoms to replace,
                        it has lower priority over `protected_ids` (replace_ids
                        which are present in protected_ids would be protected). Default: None.
    :param protected_ids: iterable with ids of heavy atoms at which no H replacement should be made and/or ids of
                          protected hydrogens. This argument has a higher priority over `replace_ids`. Default: None.
    :param symmetry_fixes: accepted for API compatibility with mutate/grow functions but not used for macrocyclization.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param return_rxn: whether to additionally return rxn of a transformation. Default: False.
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB. Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule. Default: False.
    :param ncores: number of cores. Default: 1.
    :param filter_func: a function which will filter selected fragments by additional rules
                        (in this way one may add arbitrary selection constrains). The function takes necessary first
                        three arguments: row_ids (list or set of row_ids from the fragment database supplied to
                        make_macrocycle), cursor of that fragment database and radius (int). This is required to
                        access the selected fragments. Other arguments are custom and user-defined.
                        It is the most convenient to define a filtering function, implement specific logic inside and
                        pass it using functools.partial. The filtering function should return a list/set
                        of row ids which are a subset of the input row ids.
    :param sample_func: a function which will sample selected fragments if max_replacements is supplied. If omitted
                        uniform sampling will be used. The function takes necessary first four arguments: row_ids
                        (list or set of row_ids from the fragment database), cursor of that fragment database,
                        radius (int) and the number of returned items (int). This is required to access the selected
                        fragments. Other arguments can be custom and user-defined. The function should return
                        a list/set of selected row ids.
    :param set_names: column name or list of column names in radius tables defining the set(s) of fragments to use
                      with min_freq in v1 databases. A fragment is included if at least one of the named sets
                      satisfies the min_freq threshold (OR logic). If None (default), all available set columns
                      are used. If a column name is not found, a ValueError is raised listing available set names.
                      Ignored for v0 databases. Default: None.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment.
    :return: generator over new molecules. If no additional return arguments were requested this is a generator over
             SMILES of new molecules. If additional return values were requested, the function yields a list where
             the first item is SMILES, then rxn string (optional), frequency (optional), RDKit Mol object (optional).
             Only entries with distinct SMILES will be returned.
    """

    def __get_protected_ids(m, replace_ids, protected_ids):
        # the list of ids of heavy atoms with protected hydrogens should be returned

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

    __check_db_existence(db_name)
    products = set()

    mol = Chem.AddHs(mol)
    source_smi = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    protected_ids = __get_protected_ids(mol, replace_ids, protected_ids)

    if ncores == 1:

        for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol, mol2=None, db_name=db_name,
                                                                         radius=radius, dist=dist,
                                                                         min_size=0, max_size=0,
                                                                         min_rel_size=0, max_rel_size=1,
                                                                         min_inc=min_atoms, max_inc=max_atoms,
                                                                         max_replacements=max_replacements,
                                                                         replace_cycles=False,
                                                                         protected_ids_1=protected_ids,
                                                                         protected_ids_2=None,
                                                                         min_freq=min_freq, set_names=set_names,
                                                                         filter_func=filter_func,
                                                                         sample_func=sample_func,
                                                                         return_frag_smi_only=False,
                                                                         macrocycle=True, **kwargs):
            for smi, m, rxn in __frag_replace(mol, None, frag_sma, core_sma, radius, context_mol):
                if max_replacements is None or (max_replacements is not None and len(products) < max_replacements):
                    if smi != source_smi and smi not in products:
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
        try:
            for items in p.imap(__frag_replace_mp, __get_data_macrocycle(mol, db_name, radius, dist,
                                                                         min_atoms, max_atoms, protected_ids,
                                                                         min_freq, set_names, max_replacements,
                                                                         filter_func=filter_func,
                                                                         sample_func=sample_func, **kwargs),
                                chunksize=100):
                for smi, m, rxn, freq in items:
                    if max_replacements is None or (max_replacements is not None and len(products) < max_replacements):
                        if smi != source_smi and smi not in products:
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
        finally:
            p.close()
            p.join()
