# Pavel Polishchuk, 2017

import os
import sys
import re
import math
import warnings
from collections import defaultdict
from rdkit import Chem, rdBase
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMMPA
from crem.mol_context import RADIUS0_ENV_CLASSES, get_canon_context_core
from multiprocessing import Pool, cpu_count
import sqlite3
import random
from itertools import product, combinations
from crem.mol_context import patt_remove_map
from crem.ring_fragments import (ATOM_INDEX_PROP, RING_CUT_DUMMY_ISOTOPE, _ensure_atom_indices,
                                 iter_partial_ring_fragments)

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
__patt_remove_brackets = re.compile(r'\(\)')

# Ring-cut attachment point are labeled with isotope labels to separate the from acyclic cuts.
# During DB creation isotope lebels of soutce structures are srtipped. During structure
# generation isotope labels in parent structures are transferred through intermediate storate
# to output structures. In all cases it is expeected that only ring-arc attachment points will
# have isotope labels.
__RING_CUT_DUMMY_PATT = f'[{RING_CUT_DUMMY_ISOTOPE}*'

#: Columns of radius{N} that are schema metadata rather than per-set occurrence counts.
#: `is_ring_closure` is absent from v2 tables; naming it here is harmless because set
#: discovery is a set difference.
_RESERVED_RADIUS_COLUMNS = frozenset(
    {'env_id', 'core_smi_id', 'core_num_atoms', 'dist2', 'is_ring_closure'}
)

__molzip_params = rdmolops.MolzipParams()
__molzip_params.label = rdmolops.MolzipLabel.AtomMapNumber
__explicit_h_query = Chem.MolFromSmarts("[#1]")
__remove_hs_params = Chem.RemoveHsParameters()
__remove_hs_params.removeDefiningBondStereo = True
# (canonical SMILES, canonical ranking) of an isolated aromatic ring system -> its fixed
# pairwise distances, used by the ring-geometry filter (see __rigid_pair_distances)
__rigid_geometry_cache = {}
__atom_properties_to_backup = ("isotope",)
__atom_property_backup_handlers = {
    "isotope": (
        lambda atom: atom.GetIsotope(),
        lambda atom, value: atom.SetIsotope(value),
        0,
    ),
}


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


def __standardize_context_mol(context_mol, old_to_new_map):
    # SMILES is only the dedup/query key; the Mol value carries atom props.
    context_std = __renumber_attachment_points(context_mol, old_to_new_map)
    context_smi = Chem.MolToSmiles(context_std, isomericSmiles=True)
    return context_smi, context_std


def __clear_atom_prop(mol, prop):
    for atom in mol.GetAtoms():
        if atom.HasProp(prop):
            atom.ClearProp(prop)


def __atom_backup_prop_name(name):
    return f"__crem_{name}"


def __backup_atom_properties(mol, names):
    mol = Chem.Mol(mol)
    names = tuple(names)
    unsupported = set(names) - set(__atom_property_backup_handlers)
    if unsupported:
        raise ValueError(f"Unsupported atom property backup(s): {sorted(unsupported)}")
    for atom in mol.GetAtoms():
        for name in names:
            backup_name = __atom_backup_prop_name(name)
            if atom.HasProp(backup_name):
                atom.ClearProp(backup_name)
            getter, setter, default_value = __atom_property_backup_handlers[name]
            value = getter(atom)
            if value != default_value:
                atom.SetIntProp(backup_name, value)
                setter(atom, default_value)
    return mol


def __restore_atom_properties(mol):
    for atom in mol.GetAtoms():
        for name, (_, setter, _) in __atom_property_backup_handlers.items():
            backup_name = __atom_backup_prop_name(name)
            if atom.HasProp(backup_name):
                setter(atom, atom.GetIntProp(backup_name))
                atom.ClearProp(backup_name)


def __prepare_context_mol_for_output(context_mol):
    __restore_atom_properties(context_mol)
    __clear_atom_prop(context_mol, ATOM_INDEX_PROP)
    return context_mol


def __has_ring(mol):
    """Return True if `mol` contains at least one ring, without requiring
    sanitization or ring perception.

    The cyclomatic number (bonds - atoms + connected_components) is a purely 
    topological ring count that needs neither sanitization nor ring perception; 
    dummy attachment atoms are pendant and so do not change it.
    """
    num_components = len(Chem.GetMolFrags(mol))
    return mol.GetNumBonds() - mol.GetNumAtoms() + num_components > 0


def __fragment_mol(mol, radius=3, return_ids=True, keep_stereo=False, protected_ids=None, symmetry_fixes=False,
                   min_core_atoms=None, max_core_atoms=None, include_cyclic_cores=False):
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

    def get_atom_prop(molecule, prop=ATOM_INDEX_PROP):
        res = []
        for a in molecule.GetAtoms():
            if a.GetAtomicNum():
                res.append(a.GetIntProp(prop))
        return tuple(sorted(res))

    if protected_ids:
        return_ids = True

    protected_ids = set(protected_ids) if protected_ids else set()

    # due to the bug https://github.com/rdkit/rdkit/issues/3040
    # outputs of rdMMPA.FragmentMol calls will contain duplicated fragments;
    # SMILES keys remove duplicates while Mol values preserve atom properties.
    output = {}

    # set original atom idx to keep them in fragmented mol
    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp(ATOM_INDEX_PROP, atom.GetIdx())

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
        num_heavy_atoms = core_mol.GetNumHeavyAtoms()
        if (min_core_atoms is not None and num_heavy_atoms < min_core_atoms) or \
                (max_core_atoms is not None and num_heavy_atoms > max_core_atoms):
            if not include_cyclic_cores or not __has_ring(core_mol):
                return
        env, frag, old_to_new_map = get_canon_context_core(context_mol, core_mol, radius, keep_stereo,
                                                           return_att_map=True)
        context_smi, context_std = __standardize_context_mol(context_mol, old_to_new_map)
        output[(env, frag, core_ids, context_smi, num_heavy_atoms)] = context_std

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
            context_by_smi = {key[3]: context_mol for key, context_mol in output.items()}
            for item in extended_output:
                output[item] = context_by_smi[item[3]]

    res = []
    for (env, frag, _, context_smi, num_heavy_atoms), context_mol in output.items():
        res.append((env, frag, context_mol, num_heavy_atoms))
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
                        if n.GetAtomicNum() != 0 and n.HasProp(ATOM_INDEX_PROP):
                            anchor_id = n.GetIntProp(ATOM_INDEX_PROP)
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
            atom.SetIntProp(ATOM_INDEX_PROP, atom.GetIdx())
        for atom in mol2.GetAtoms():
            atom.SetIntProp(ATOM_INDEX_PROP, atom.GetIdx())

    frags_1 = rdMMPA.FragmentMol(mol1, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    frags_2 = rdMMPA.FragmentMol(mol2, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)

    frags_1 = _prepare_single_cut_contexts(frags_1, protected_ids_1)
    frags_2 = _prepare_single_cut_contexts(frags_2, protected_ids_2)
    fake_core = '[*:1]C[*:2]'
    output = {}

    for ctx_1, ctx_2 in product(frags_1, frags_2):
        # keep historical convention: attachment from mol1 starts with map 2
        ctx_1 = __renumber_attachment_points(ctx_1, {1: 2})

        chains = Chem.CombineMols(ctx_1, ctx_2)
        env, frag, old_to_new_map = get_canon_context_core(chains, fake_core, radius=radius, keep_stereo=keep_stereo,
                                                            return_att_map=True)
        context_smi, context_std = __standardize_context_mol(chains, old_to_new_map)
        output[(env, '[H][*:1].[H][*:2]', context_smi, 0)] = context_std

    res = []
    for (env, core, context_smi, num_heavy_atoms), context_mol in output.items():
        res.append((env, core, context_mol, num_heavy_atoms))
    return res


def __fragment_mol_macrocycle(mol, radius=3, ring_size=None, keep_stereo=False, protected_ids=None,
                              return_ids=True, label_variants=(False,)):
    """Anchor-pair fragmentation via two independent single H cuts (broad make_cycle).

    Unlike __fragment_mol_ring_closure the two contexts stay separate, so the env is
    disconnected. `label_variants` has the same meaning as there; the default is the plain
    (v1) form because this fragmenter exists to reach acyclic-cut rows.
    """

    def _is_h_cut_component(component):
        if component.GetNumAtoms() != 2 or component.GetNumBonds() != 1:
            return False
        atoms = list(component.GetAtoms())
        return any(a.GetAtomicNum() == 1 for a in atoms) and any(
            a.GetAtomicNum() == 0 and a.GetAtomMapNum() for a in atoms
        )

    def _prepare_single_cut_contexts(frags, protected_ids):
        output = {}  # anchor_idx -> context_mol
        for core, chains in frags:
            if core is not None:
                continue
            components = list(Chem.GetMolFrags(chains, asMols=True))
            if len(components) != 2:
                continue
            if _is_h_cut_component(components[0]):
                context = components[1]
            elif _is_h_cut_component(components[1]):
                context = components[0]
            else:
                continue

            anchor_id = None
            valid = True
            for a in context.GetAtoms():
                if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                    for n in a.GetNeighbors():
                        if n.GetAtomicNum() != 0 and n.HasProp(ATOM_INDEX_PROP):
                            anchor_id = n.GetIntProp(ATOM_INDEX_PROP)
                            if anchor_id in protected_ids:
                                valid = False
                            break
                    break
            if valid and anchor_id is not None and anchor_id not in output:
                output[anchor_id] = context
        return sorted(output.items())

    if protected_ids:
        return_ids = True
    protected_ids = set(protected_ids) if protected_ids else set()

    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp(ATOM_INDEX_PROP, atom.GetIdx())

    if ring_size is None:
        rs_min = rs_max = None
    elif isinstance(ring_size, int):
        rs_min = rs_max = ring_size
    else:
        if len(ring_size) != 2:
            raise ValueError("ring_size must be an int or a (min, max) tuple")
        rs_min, rs_max = ring_size

    distance_matrix = Chem.GetDistanceMatrix(mol)

    frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    contexts = _prepare_single_cut_contexts(frags, protected_ids)
    # see __fragment_mol_ring_closure: the stand-in core must match the context's classes
    fake_cores = {True: '[1*:1]C[1*:2]', False: '[*:1]C[*:2]'}
    output = {}

    for (anchor_1, ctx_1), (anchor_2, ctx_2) in combinations(contexts, 2):
        if anchor_1 == anchor_2:
            continue

        # ring_size = d_in + dist2_frag → derive per-pair dist2 window.
        d_in = int(distance_matrix[anchor_1, anchor_2])
        if rs_min is not None:
            lo = max(1, rs_min - d_in)
            hi = rs_max - d_in
            if hi < lo:
                continue  # target ring smaller than the in-mol shortest path
            frag_dist = (lo, hi) if lo != hi else lo
        else:
            frag_dist = None

        # keep convention consistent with linker generation
        ctx_1_renumbered = __renumber_attachment_points(ctx_1, {1: 2})

        for labelled in label_variants:
            # See __fragment_mol_ring_closure: the anchors are the prospective
            # ring-closure ends, so on a v2 database they must carry the ring-cut label.
            # Rebuilding chains per variant is required, not an accident: it is what keeps
            # the label out of the plain variant. MMPA-created dummies always arrive with
            # isotope 0, so only the labelled variant needs to touch them.
            chains = Chem.CombineMols(ctx_1_renumbered, ctx_2)
            if labelled:
                for a in chains.GetAtoms():
                    if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                        a.SetIsotope(RING_CUT_DUMMY_ISOTOPE)
            env, frag, old_to_new_map = get_canon_context_core(chains, fake_cores[labelled],
                                                              radius=radius,
                                                              keep_stereo=keep_stereo,
                                                              return_att_map=True,
                                                              preserve_dummy_isotopes=True)

            # Move the attachment point from the second chain into the first one.
            # "__crem_index" stores the original molecule atom id; fragment-local atom
            # indices can differ after GetMolFrags().
            anchor_id = None
            att_map = None
            chains = __renumber_attachment_points(chains, old_to_new_map)
            chains = Chem.GetMolFrags(chains, asMols=True)
            if len(chains) != 2:
                raise RuntimeError("Expected two macrocycle context fragments after attachment renumbering")
            for a in chains[1].GetAtoms():
                if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                    att_map = a.GetAtomMapNum()
                    for n in a.GetNeighbors():
                        if n.GetAtomicNum() != 0 and n.HasProp(ATOM_INDEX_PROP):
                            anchor_id = n.GetIntProp(ATOM_INDEX_PROP)
                            break
                    break
            if anchor_id is None or att_map is None:
                raise RuntimeError("Could not find mapped attachment point and original anchor id")
            a = None
            for atom in chains[0].GetAtoms():
                if atom.HasProp(ATOM_INDEX_PROP) and atom.GetIntProp(ATOM_INDEX_PROP) == anchor_id:
                    a = atom
                    break
            if a is None:
                raise RuntimeError(f"Could not find original anchor atom {anchor_id} in macrocycle context")
            replaced_hydrogen = False
            for n in a.GetNeighbors():
                if n.GetAtomicNum() == 1:
                    n.SetAtomicNum(0)
                    n.SetAtomMapNum(att_map)
                    # This dummy is a hydrogen converted in place, so unlike the
                    # MMPA-created ones above it can arrive with an isotope of its own -
                    # an explicit [2H] becomes [2*:n] unless it is overwritten here.
                    n.SetIsotope(RING_CUT_DUMMY_ISOTOPE if labelled else 0)
                    replaced_hydrogen = True
                    break
            if not replaced_hydrogen:
                raise RuntimeError(f"Could not find replaceable hydrogen on original anchor atom {anchor_id}")
            context_mol = Chem.Mol(chains[0])
            context_smi = Chem.MolToSmiles(context_mol, isomericSmiles=True)
            output[(env, '[H][*:1].[H][*:2]', context_smi, 0, frag_dist)] = context_mol

    res = []
    for (env, core, context_smi, num_heavy_atoms, frag_dist), context_mol in output.items():
        res.append((env, core, context_mol, num_heavy_atoms, frag_dist))
    return res


def __fragment_mol_ring_closure(mol, radius=3, ring_size=None, keep_stereo=False, protected_ids=None,
                                 return_ids=True, label_variants=(True,)):
    """Anchor-pair fragmentation for `make_cycle(ring_closures=True)`.

    Cuts pairs of H bonds in one MMPA call (maxCuts=2) so the resulting
    context is a single connected Mol with two attachment points — matching
    the connected-env shape that ring-closure rows in the DB carry.

    `ring_size` (int or (min,max)) is translated per-pair into a `dist2`
    window using d_in = topological distance between the two anchor heavy
    atoms in the input molecule. Pairs whose feasible dist2 window is empty
    (target ring smaller than the in-mol path) are dropped.

    The two anchors are the prospective ring-closure points, so on a v2 database they
    must carry the ring-cut isotope label that v2 stores on every ring cut. That is what
    `label_variants` selects, per env variant to emit from the single MMPA pass:

    * ``(True,)``  - labelled envs only; matches v2 ring-arc rows (strict make_cycle).
    * ``(False,)`` - plain envs only; the v1 convention, where provenance instead comes
      from the ``is_ring_closure`` predicate.
    * ``(True, False)`` - both, for broad ``make_cycle(ring_closures=False)`` on v2, where
      "either provenance" has to be expressed as a union over the two env forms.

    Returns 5-tuples: (env, '[H][*:1].[H][*:2]', context_mol, 0, dist2_filter).
    """
    protected_ids = set(protected_ids) if protected_ids else set()

    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp(ATOM_INDEX_PROP, atom.GetIdx())

    if ring_size is None:
        rs_min = rs_max = None
    elif isinstance(ring_size, int):
        rs_min = rs_max = ring_size
    else:
        if len(ring_size) != 2:
            raise ValueError("ring_size must be an int or a (min, max) tuple")
        rs_min, rs_max = ring_size

    distance_matrix = Chem.GetDistanceMatrix(mol)
    # `fake_core` stands in for the real core ('[H][*:1].[H][*:2]') during standardisation.
    # Its attachment points must carry the same class labels as the context: at radius 0 the
    # env is derived from the core, so an unlabelled fake core would report a ring closure
    # as an acyclic class.
    fake_cores = {True: '[1*:1]C[1*:2]', False: '[*:1]C[*:2]'}
    seen_pairs = set()
    output = {}

    frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=2,
                               resultsAsMols=True, maxCutBonds=100)
    for mmpa_core, mmpa_chains in frags:
        # In MMPA's convention for maxCuts=2 the *core* is the connected
        # rest-of-mol carrying both * attachment points, and *chains* is the
        # disconnected pair of cut substituents (here `[H][*:1].[H][*:2]`).
        # Single-cut results (core is None) come along for the ride — skip
        # them, ring-closure semantics need both anchors at once.
        if mmpa_core is None:
            continue
        chain_components = list(Chem.GetMolFrags(mmpa_chains, asMols=True))
        if len(chain_components) != 2:
            continue
        chain_smis = {Chem.MolToSmiles(c) for c in chain_components}
        if chain_smis != {'[H][*:1]', '[H][*:2]'}:
            continue  # not a pure 2-H cut (one of the cuts was on a heavy bond)
        context = mmpa_core
        # Identify the two anchor heavy atoms and their map numbers.
        anchors = {}  # map_num -> __crem_index property of anchor heavy atom
        for a in context.GetAtoms():
            if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                map_num = a.GetAtomMapNum()
                for n in a.GetNeighbors():
                    if n.GetAtomicNum() != 0 and n.HasProp(ATOM_INDEX_PROP):
                        anchors[map_num] = n.GetIntProp(ATOM_INDEX_PROP)
                        break
        if len(anchors) != 2:
            continue
        anchor_ids = sorted(anchors.values())
        if anchor_ids[0] == anchor_ids[1]:
            continue
        if protected_ids and not protected_ids.isdisjoint(anchor_ids):
            continue
        pair = tuple(anchor_ids)
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)

        # ring_size = d_in + dist2_frag → derive per-pair dist2 window.
        d_in = int(distance_matrix[anchor_ids[0], anchor_ids[1]])
        if rs_min is not None:
            lo = max(1, rs_min - d_in)
            hi = rs_max - d_in
            if hi < lo:
                continue  # target ring smaller than the in-mol path
            frag_dist = (lo, hi) if lo != hi else lo
        else:
            frag_dist = None

        for labelled in label_variants:
            # The anchors are the prospective ring-closure ends. Labelling them makes the
            # env match v2's ring-arc rows; leaving them plain makes it match v1's (where
            # provenance comes from the is_ring_closure predicate instead).
            # The per-variant copy is required, not an accident: it is what keeps the label
            # out of the plain variant. Every dummy here was created by MMPA and so arrives
            # with isotope 0, which is why only the labelled variant touches them.
            variant = Chem.Mol(context)
            if labelled:
                for a in variant.GetAtoms():
                    if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
                        a.SetIsotope(RING_CUT_DUMMY_ISOTOPE)

            env, frag, old_to_new_map = get_canon_context_core(
                variant, fake_cores[labelled], radius=radius, keep_stereo=keep_stereo,
                return_att_map=True, preserve_dummy_isotopes=True,
            )
            if env is None or not frag:
                continue

            # Renumber the context Mol so its * map numbers match the
            # standardised env / DB-side core numbering.
            context_smi, std_context = __standardize_context_mol(variant, old_to_new_map)
            output[(env, '[H][*:1].[H][*:2]', context_smi, 0, frag_dist)] = std_context

    res = []
    for (env, core, context_smi, num_heavy_atoms, frag_dist), context_mol in output.items():
        res.append((env, core, context_mol, num_heavy_atoms, frag_dist))
    return res


def __fragment_mol_partial_cycles(mol, radius=3, keep_stereo=False, protected_ids=None, symmetry_fixes=False,
                                  return_ids=True, min_core_atoms=None, max_core_atoms=None,
                                  include_cyclic_cores=False, side_cut_mode="all",
                                  label_all_ring_cuts=True):
    """Fragment partial cycles for mutate_mol(..., replace_cycles='partial_*').

    The source core is one connected ring arc made by two non-aromatic single
    ring-bond cuts plus 0-2 acyclic side cuts. Explicit hydrogens are stripped
    before enumeration so they do not create distinct ring-replacement contexts.

    `label_all_ring_cuts` must follow the schema version of the DB being queried
    (True for v2, False for v1); see crem.ring_fragments.iter_partial_ring_fragments.
    """
    protected_ids = set(protected_ids) if protected_ids else set()

    work_mol = _ensure_atom_indices(mol) if return_ids else Chem.Mol(mol)
    work_mol = Chem.RemoveHs(work_mol)

    atom_specific_output = {}
    iter_min_core_atoms = None if include_cyclic_cores else min_core_atoms
    iter_max_core_atoms = None if include_cyclic_cores else max_core_atoms
    for core_mol, context_mol, core_atom_ids in iter_partial_ring_fragments(
        work_mol,
        label_all_ring_cuts=label_all_ring_cuts,
        max_acyclic_cuts=2,
        min_core_atoms=iter_min_core_atoms,
        max_core_atoms=iter_max_core_atoms,
        side_cut_mode=side_cut_mode,
    ):
        core_atom_ids = tuple(sorted(core_atom_ids))
        if protected_ids and not protected_ids.isdisjoint(core_atom_ids):
            continue

        num_heavy_atoms = core_mol.GetNumHeavyAtoms()
        if (min_core_atoms is not None and num_heavy_atoms < min_core_atoms) or \
                (max_core_atoms is not None and num_heavy_atoms > max_core_atoms):
            if not include_cyclic_cores or not __has_ring(core_mol):
                continue
        env, frag, old_to_new_map = get_canon_context_core(
            context_mol,
            core_mol,
            radius,
            keep_stereo,
            return_att_map=True,
            preserve_dummy_isotopes=True,
        )
        if env is None or not frag:
            continue
        context_smi, context_std = __standardize_context_mol(context_mol, old_to_new_map)
        atom_specific_output[(env, frag, core_atom_ids, context_smi, num_heavy_atoms)] = context_std

    if symmetry_fixes:
        output = atom_specific_output
    else:
        output = {}
        for (env, frag, _, context_smi, num_heavy_atoms), context_mol in atom_specific_output.items():
            output[(env, frag, tuple(), context_smi, num_heavy_atoms)] = context_mol

    res = []
    for (env, frag, _, context_smi, num_heavy_atoms), context_mol in output.items():
        res.append((env, frag, context_mol, num_heavy_atoms))
    return res


# def __molzip_would_duplicate_bond(mol):
#     neighbors_by_map = defaultdict(list)
#     for atom in mol.GetAtoms():
#         if atom.GetAtomicNum() != 0 or not atom.GetAtomMapNum():
#             continue
#         heavy_neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum()]
#         if len(heavy_neighbors) == 1:
#             neighbors_by_map[atom.GetAtomMapNum()].append(heavy_neighbors[0])
#
#     for neighbor_ids in neighbors_by_map.values():
#         if len(neighbor_ids) != 2:
#             continue
#         if neighbor_ids[0] == neighbor_ids[1]:
#             return True
#         if mol.GetBondBetweenAtoms(neighbor_ids[0], neighbor_ids[1]) is not None:
#             return True
#     return False


def __bfs_path(mol, start, end, allowed):
    """Shortest path (list of atom ids) from start to end traversing only atoms in `allowed`
    (a set that must include start and end). Returns None if no such path exists."""
    if start == end:
        return [start]
    from collections import deque
    prev = {start: None}
    q = deque([start])
    while q:
        cur = q.popleft()
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            idx = nb.GetIdx()
            if idx in allowed and idx not in prev:
                prev[idx] = cur
                if idx == end:
                    path = [end]
                    while prev[path[-1]] is not None:
                        path.append(prev[path[-1]])
                    return path[::-1]
                q.append(idx)
    return None


def __identify_new_ring(p):
    """Detect a 2-attachment-point intramolecular ring closure in a product `p` of
    __frag_replace, using the ``__crem`` labels carried by inserted (replacement) atoms.

    A ring closure is a transformation whose inserted bridge connects two atoms that are
    already part of the same connected parent fragment (make_cycle and partial-ring mutate).
    grow / single-cut mutate (one attachment) and link (two attachments across two separate
    molecules) are not ring closures and return None.

    Returns a dict {bridge, a1, a2, parent_arc, ring_atoms, ring_size} describing the smallest
    newly formed ring, or None when `p` is not a simple two-point ring closure.
    """
    crem = set(a.GetIdx() for a in p.GetAtoms() if a.HasProp('__crem'))
    if not crem:
        return None
    noncrem = set(a.GetIdx() for a in p.GetAtoms() if a.GetIdx() not in crem and a.GetAtomicNum() > 1)
    # junction bonds connect an inserted (crem) heavy atom to a parent (non-crem) heavy atom;
    # record both the parent-side anchor and the inserted-side neighbour
    junctions = []
    for b in p.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        ic, jc = i in crem, j in crem
        if ic != jc:
            anchor, cside = (j, i) if ic else (i, j)
            if p.GetAtomWithIdx(anchor).GetAtomicNum() <= 1:
                continue
            junctions.append((anchor, cside))
    if len(junctions) != 2 or junctions[0][0] == junctions[1][0]:
        return None
    (a1, b1), (a2, b2) = junctions
    parent_arc = __bfs_path(p, a1, a2, noncrem)
    if parent_arc is None:
        return None  # anchors not connected through the parent (e.g. link) -> not a closure
    # bridge core must traverse the inserted atoms, not shortcut through a parent bond
    bridge_core = [b1] if b1 == b2 else __bfs_path(p, b1, b2, crem)
    if bridge_core is None:
        return None
    ring_atoms = set(parent_arc) | set(bridge_core) | {a1, a2}
    # the cycle in traversal order: along the parent from a1 to a2, then back through the bridge
    ring_path = list(parent_arc) + list(reversed(bridge_core))
    return {'bridge': set(bridge_core), 'crem': crem, 'a1': a1, 'a2': a2, 'parent_arc': parent_arc,
            'ring_atoms': ring_atoms, 'ring_size': len(ring_atoms), 'ring_path': ring_path}


def __aromatic_arc(p, ring, ring_atoms):
    """The contiguous run of `ring` (an aromatic ring given as an ordered tuple of atom ids)
    that lies inside the new ring `ring_atoms`, as an ordered list of atom ids. Returns None
    when the overlap is empty, covers the whole aromatic ring, or is broken into several runs
    (the new ring then enters and leaves the aromatic system more than once and the simple
    span argument below does not apply)."""
    n = len(ring)
    inside = set(i for i, idx in enumerate(ring) if idx in ring_atoms)
    if not inside or len(inside) == n:
        return None
    start = next(i for i in inside if (i - 1) % n not in inside)
    if any((start + k) % n not in inside for k in range(len(inside))):
        return None  # overlap is not a single contiguous arc
    return [ring[(start + k) % n] for k in range(len(inside))]


def __ring_topology_ok(p, ringinfo):
    """Stage 0 geometry check (calibrated). Returns False only for proven-impossible small
    ring closures across rigid aromatic systems; True (accept) otherwise.

    The test is keyed on the rigid arc *inside* the new ring rather than on the atoms the
    inserted fragment happens to be bonded to. A ring is impossible whenever it contains a
    contiguous arc of an aromatic ring whose two ends are held further apart than the rest of
    the ring can span - which is a property of the arc, no matter whether the attachment
    points are the arc ends themselves (make_cycle between two ring atoms) or lie one or more
    bonds outside it (a closure starting on a side chain).

    Reject bands, each ending at the smallest bridge that closes in practice:
      * arc of 3 atoms of a six-membered ring (meta positions) -> new ring 4-8; ring 9 is
        [6]metacyclophane, which exists, and bare templates embed there;
      * arc of 4 atoms of a six-membered ring (para positions) -> new ring 5-9; ring 10 is
        [6]paracyclophane, isolable though strained;
      * arc of 3 atoms of a five-membered aromatic ring (the 1,3-span across its middle atom)
        -> new ring 4-8. This span is the exact analogue of benzene-meta: two divergent
        exocyclic vectors. Bare furan, pyrrole and thiophene templates fail at every bridge
        length up to ring 7 whatever the middle atom is; at ring 8 only the span across an O
        closes, and it closes only with the ring bent by 33 deg, so the band covers it too.
    Everything else passes: ortho/adjacent arcs (indane, acenaphthene, fluorene), sp3 spans,
    and rings large enough to be genuine cyclophanes.
    """
    ring_size = ringinfo['ring_size']
    if not (4 <= ring_size <= 9):
        return True
    ring_atoms = ringinfo['ring_atoms']
    for ring in p.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6):
            continue
        if not all(p.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        arc = __aromatic_arc(p, ring, ring_atoms)
        if arc is None:
            continue
        if len(ring) == 6:
            if len(arc) == 3 and ring_size <= 8:       # meta positions bridged
                return False
            if len(arc) == 4 and 5 <= ring_size <= 9:  # para positions bridged
                return False
        elif len(arc) == 3 and ring_size <= 8:  # 1,3-span across a five-membered ring
            return False
    return True


def __aromatic_ring_systems(p, atoms=None):
    """Maximal groups of aromatic ring atoms fused through shared atoms - one rigid planar unit
    each. Rings joined only by a rotatable bond (biaryls) stay separate units, since their
    mutual geometry is not fixed. Restricted to units touching `atoms` when given."""
    units = []
    for ring in p.GetRingInfo().AtomRings():
        if not all(p.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        rset = set(ring)
        overlapping = [u for u in units if u & rset]
        units = [u for u in units if not (u & rset)]
        units.append(rset.union(*overlapping) if overlapping else rset)
    if atoms is not None:
        units = [u for u in units if u & atoms]
    return units


def __rigid_pair_distances(p, unit):
    """Interatomic distances (Angstrom) inside one rigid aromatic unit, taken from an embedded
    copy of the *isolated* ring system, keyed by atom index in `p`. A fused aromatic system is
    planar and rigid, so a single conformer fixes every distance in it - which is exactly the
    information `GetMoleculeBoundsMatrix` lacks (it bounds distant pairs by van der Waals
    contact only, which is why the plain triangle-smoothing test never fires).

    Results are cached so a scaffold recurring across thousands of products is embedded once.
    The key is the isolated system's canonical SMILES *together with its canonical ranking*:
    symmetry-equivalent atoms of the bare system (the two ortho positions of a biphenyl ring,
    say) are interchangeable in the stripped fragment but not in the product, and their rank
    assignment can differ between two extractions of the same scaffold. Keying on the SMILES
    alone would then read a cached distance under a mismatched labelling, making the verdict
    depend on which products happened to be processed first. Returns None when the system
    cannot be extracted or embedded."""
    idx = sorted(unit)
    em = Chem.RWMol(p)
    for i in sorted(set(range(p.GetNumAtoms())) - unit, reverse=True):
        em.RemoveAtom(i)
    sub = em.GetMol()
    try:
        Chem.SanitizeMol(sub)
        smi = Chem.MolToSmiles(sub)
        ranks = list(Chem.CanonicalRankAtoms(sub))
    except Exception:
        return None
    key = (smi, tuple(ranks))
    if key not in __rigid_geometry_cache:
        __rigid_geometry_cache[key] = __embed_rigid_system(sub, ranks)
    by_rank = __rigid_geometry_cache[key]
    if by_rank is None:
        return None
    out = {}
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            span = by_rank.get((min(ranks[a], ranks[b]), max(ranks[a], ranks[b])))
            if span is not None:
                out[(idx[a], idx[b])] = span
    return out


def __embed_rigid_system(sub, ranks):
    """Embed an isolated aromatic ring system and return its pairwise heavy-atom distances
    keyed by canonical rank pairs (so the result is reusable for any copy of the system).

    A fused system is rigid, so one conformer gives every distance. A biaryl-style assembly has
    one rotatable bond joining its two halves; there the torsion is scanned and the *minimum*
    distance over the scan is kept, which is the sound lower bound. Assemblies with more than
    one such bond are not analysed (returns None) rather than guessed at."""
    try:
        from rdkit.Chem import rdDistGeom, rdMolTransforms
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        params.maxIterations = 4  # give up quickly; a failure only skips the injection
        m = Chem.Mol(sub)         # heavy atoms are enough and embed several times faster
        if rdDistGeom.EmbedMolecule(m, params) != 0:
            return None
        conf = m.GetConformer()
        hinges = [b for b in sub.GetBonds()
                  if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE]
        if len(hinges) > 1:
            return None
        torsion = None
        if hinges:
            b, c = hinges[0].GetBeginAtom(), hinges[0].GetEndAtom()
            a = next((x.GetIdx() for x in b.GetNeighbors() if x.GetIdx() != c.GetIdx()), None)
            d = next((x.GetIdx() for x in c.GetNeighbors() if x.GetIdx() != b.GetIdx()), None)
            if a is None or d is None:
                return None
            torsion = (a, b.GetIdx(), c.GetIdx(), d)
        n = sub.GetNumAtoms()
        out = {}
        for angle in (range(0, 360, 15) if torsion else (None,)):
            if angle is not None:
                rdMolTransforms.SetDihedralDeg(conf, *torsion, float(angle))
            pos = [conf.GetAtomPosition(i) for i in range(n)]
            for x in range(n):
                for y in range(x + 1, n):
                    key = (min(ranks[x], ranks[y]), max(ranks[x], ranks[y]))
                    dist = pos[x].Distance(pos[y])
                    lo, hi = out.get(key, (dist, dist))
                    out[key] = (min(lo, dist), max(hi, dist))
        if torsion:
            # a hinged assembly is only a model of the real geometry: closing a strained ring
            # over it splays the two halves further apart than the relaxed scan ever gets, so
            # only the lower bound (the reach limit) may be trusted
            out = {k: (lo, None) for k, (lo, _) in out.items()}
        return out
    except Exception:
        return None


def __rigid_groups(p, atoms):
    """Rigid (or nearly rigid) aromatic groups touching `atoms`: each fused ring system, plus
    every pair of fused systems joined by a single bond (a biaryl, rigid except for one
    torsion)."""
    units = __aromatic_ring_systems(p, atoms)
    groups = list(units)
    for a in range(len(units)):
        for b in range(a + 1, len(units)):
            joins = [bond for bond in p.GetBonds()
                     if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} & units[a]
                     and {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} & units[b]]
            if len(joins) == 1 and joins[0].GetBondType() == Chem.BondType.SINGLE:
                groups.append(units[a] | units[b])
    return groups


def __ring_reach_ok(p, ringinfo, slack=0.10, max_reach=10.0):
    """Stage 0b geometry check: can the rest of the new ring even reach across the rigid
    aromatic stretch it contains? Returns False only when it provably cannot.

    Stage 0's arc rules cover spans inside a single six- or five-membered ring, where the
    obstruction is angular. The complementary failure is one of pure distance: a ring closing
    over a *fused* system (naphthalene 2,6 is ~5 A apart) or over a biaryl (4,4' is ~8.6 A)
    needs a bridge long enough to get there. The span is read off the real geometry of the
    isolated ring system (see __rigid_pair_distances) and compared against the bridge stretched
    perfectly straight - the sum of its maximum bond lengths, which no conformation can exceed.
    Nothing here depends on bond angles, so strained-but-real rings cannot be rejected.
    """
    path = ringinfo.get('ring_path')
    if not path or len(path) < 4:
        return True
    n = len(path)
    pt = Chem.GetPeriodicTable()

    def max_bond(x, y):
        return (pt.GetRcovalent(p.GetAtomWithIdx(x).GetAtomicNum())
                + pt.GetRcovalent(p.GetAtomWithIdx(y).GetAtomicNum()) + 0.10)

    for group in __rigid_groups(p, set(path)):
        inside = [k for k in range(n) if path[k] in group]
        if len(inside) < 2 or len(inside) == n:
            continue
        starts = [k for k in inside if (k - 1) % n not in inside]
        if len(starts) != 1:
            continue  # the ring dips in and out of the system - no single span to measure
        start = starts[0]
        run = [(start + k) % n for k in range(len(inside))]
        if any(k not in inside for k in run):
            continue
        a, b = path[run[0]], path[run[-1]]
        # the rest of the cycle, stretched straight, is all the bridge can ever span
        rest = [path[(run[-1] + k) % n] for k in range(n - len(run) + 1)] + [a]
        reach = sum(max_bond(rest[k], rest[k + 1]) for k in range(len(rest) - 1))
        if reach >= max_reach:
            continue
        dists = __rigid_pair_distances(p, group)
        if not dists:
            continue
        span = dists.get((min(a, b), max(a, b)))
        if span is not None and span[0] * (1.0 - slack) > reach:
            return False
    return True


def __frag_replace(mol1, mol2, old_frag_smi, new_frag_smi, radius, context_mol=None, frag_ids_1=None, frag_ids_2=None,
                   discard_ring_geometry=True):
    """
    INPUT
        mol1:         mol for mutate, grow or link
        mol2:         for link (unused in fast path)
        old_frag_smi: SMILES of a source fragment core
        new_frag_smi: SMILES of a replacement core (from DB)
        radius:       context radius considered
        context_mol:  RDKit Mol containing pre-built whole context part with mapped dummies
        discard_ring_geometry: if True (default), silently discard ring-forming products whose
                      newly created ring is geometrically impossible to embed in 3D
                      (aromatic arc rules plus the rigid-span reach test).
    OUTPUT
        generator returns canonical isomeric SMILES, RDKit Mol and rxn rule
    """

    if not isinstance(mol1, Chem.Mol):
        raise StopIteration("The first molecule in __gen_replacement always must be specified")

    if isinstance(context_mol, str):
        context_mol = Chem.MolFromSmiles(context_mol)
    if not isinstance(context_mol, Chem.Mol):
        return

    with rdBase.BlockLogs():
        repl_core_mol = Chem.MolFromSmiles(new_frag_smi)
    # label new atoms in generated structure with bool prop
    for atom in repl_core_mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            atom.SetBoolProp("__crem", True)

    transformation_smi = f"{old_frag_smi}>>{new_frag_smi}"
    try:
        p = rdmolops.molzip(Chem.CombineMols(context_mol, repl_core_mol), __molzip_params)
    except RuntimeError as exc:
        try:
            context_smi = Chem.MolToSmiles(context_mol, isomericSmiles=True)
        except Exception:
            context_smi = '<unavailable>'
        sys.stderr.write(
            f"molzip RuntimeError: {exc}. "
            f"old fragment: {old_frag_smi}; "
            f"new fragment: {new_frag_smi}; "
            f"context mol: {context_smi}; "
            f"mol1: {Chem.MolToSmiles(mol1, isomericSmiles=True)}; "
            f"mol2: {Chem.MolToSmiles(mol2, isomericSmiles=True) if isinstance(mol2, Chem.Mol) else 'None'}\n"
        )
        sys.stderr.flush()
        return
    e = Chem.SanitizeMol(p, catchErrors=True)
    if e:
        sys.stderr.write(f"Molecule {Chem.MolToSmiles(p, isomericSmiles=True)} caused "
                         f"sanitization error {e}. "
                         f"Transformation {transformation_smi} was applied to "
                         f"{Chem.MolToSmiles(mol1, isomericSmiles=True)} and "
                         f"{Chem.MolToSmiles(mol2, isomericSmiles=True) if isinstance(mol2, Chem.Mol) else 'None'}\n")
        sys.stderr.flush()
        return

    __restore_atom_properties(p)
    __clear_atom_prop(p, ATOM_INDEX_PROP)
    if p.HasSubstructMatch(__explicit_h_query):
        p = Chem.RemoveHs(p, __remove_hs_params)

    if discard_ring_geometry:
        # discard ring-forming products with an impossible-to-embed new ring; any error in the
        # geometry analysis must never discard a structure (would be a false negative)
        try:
            ringinfo = __identify_new_ring(p)
            if ringinfo is not None and (
                    not __ring_topology_ok(p, ringinfo) or not __ring_reach_ok(p, ringinfo)):
                return
        except Exception:
            pass

    smi = Chem.MolToSmiles(p, isomericSmiles=True)
    yield smi, p, transformation_smi


def _load_schema_meta(db_cur, radius):
    """Read the static schema info needed by query helpers in one pass.

    Each `PRAGMA user_version` / `PRAGMA table_info(...)` is a separate
    round-trip into SQLite; the values don't change for the life of the
    connection. Calling this once per `__gen_replacements` invocation and
    passing the result down to the per-fragment helpers eliminates 4–5
    PRAGMA round-trips that would otherwise fire for every fragment.
    """
    user_version = db_cur.execute("PRAGMA user_version").fetchone()[0]
    meta = {'user_version': user_version}
    if user_version == 1:
        # FutureWarning rather than DeprecationWarning: the latter is suppressed by Python's
        # default filters when raised from a library module rather than __main__, and the
        # stack depth from here to user code varies by entry point, so stacklevel cannot be
        # relied on to surface it.
        warnings.warn(
            "This database uses the v1 schema, which will be deprecated: cremdb_create now "
            "only writes v2 and support for reading v1 may be dropped in a future release, "
            "leaving v0 and v2. Note there is no v1 -> v2 converter - the current upgrade "
            "is a rebuild from the source SMILES with cremdb_create. By request we may "
            "develop v1 -> v2 converter, however, this is non-trivial.",
            FutureWarning,
        )
    # Fragment convention implied by the schema version. `label_all_ring_cuts` must be
    # handed to every ring fragmenter so the env/core strings it emits match the ones
    # stored in this database; `provenance_in_env` says whether ring provenance is already
    # encoded in the env string (v2) or still needs the is_ring_closure predicate (v0/v1).
    # Both are derived here, from one place, because setting them independently produces
    # silently empty results rather than an error.
    meta['label_all_ring_cuts'] = user_version >= 2
    meta['provenance_in_env'] = user_version >= 2
    if user_version in (1, 2):
        meta['frags_columns'] = {
            row[1] for row in db_cur.execute("PRAGMA table_info(frags)").fetchall()
        }
        meta['frags_h_columns'] = {
            row[1] for row in db_cur.execute("PRAGMA table_info(frags_h)").fetchall()
        }
        meta['radius_columns'] = {
            row[1] for row in db_cur.execute(f"PRAGMA table_info(radius{radius})").fetchall()
        }
    return meta


def _resolve_set_columns(schema_meta, set_names, radius):
    """Resolve `set_names` to the radius{N} count columns it selects.

    Shared by the row-id filter and the fragment fetch so that the frequency a query
    filters on and the frequency it reports are derived from the same column list; keeping
    two copies of this normalisation is what let them drift apart in the first place.

    Returns `(set_names_list, is_full_set)`, where `is_full_set` says whether the selection
    covers every set column the table has. The filter needs that distinction - it drops the
    frequency predicate entirely for a full selection at min_freq <= 0, and requires
    membership (count >= 1) for a strict subset - and it is answered here rather than by the
    caller so that it stays correct after the deduplication below.
    """
    radius_columns = schema_meta['radius_columns']
    available = sorted(radius_columns - _RESERVED_RADIUS_COLUMNS)
    if not available:
        # Every supported build path creates at least one set column, so this means the
        # database is malformed rather than merely unusual. Raising here also replaces the
        # dangling "AND" that an empty predicate list used to produce downstream.
        raise ValueError(
            f"radius{radius} has no per-set count columns, so fragment frequencies cannot "
            f"be evaluated. The database appears malformed - rebuild it with cremdb_create."
        )

    # Normalize set_names: None → all available set columns. Duplicates are dropped while
    # preserving order; a repeated name would otherwise emit the column twice in the SQL and
    # make the full-set test below fail on an otherwise complete selection.
    if set_names is None:
        set_names_list = list(available)
    elif isinstance(set_names, str):
        set_names_list = [set_names]
    else:
        set_names_list = list(dict.fromkeys(set_names))

    # Validate against the set columns rather than every column: schema columns such as
    # dist2 or core_num_atoms are in radius_columns but are not per-set counts, and
    # accepting one would build a silently meaningless filter and frequency.
    missing = [sn for sn in set_names_list if sn not in available]
    if missing:
        raise ValueError(f"Column(s) {missing} not found in radius{radius}. "
                         f"Available set names: {available}")
    return set_names_list, set(set_names_list) == set(available)


def __get_count_sql_expr(set_names_list):
    """SQL expression for a row's frequency: the largest count among the selected sets.

    Maximum rather than sum, because the filter accepts a row when *any* selected set
    clears the threshold (`col_a >= t OR col_b >= t`), and `MAX(...) >= t` is exactly that
    disjunction — so the reported number is the same quantity the row was selected on.

    Note the single-column case must not use MAX(): with one argument SQLite parses it as
    the aggregate, which would collapse the whole result set to one row.
    """
    if len(set_names_list) == 1:
        return f"r.{set_names_list[0]}"
    return f"MAX({', '.join(f'r.{sn}' for sn in set_names_list)})"


def __load_convention(db_name, radius):
    """Read only the schema-implied fragment convention, without holding a connection.

    `__gen_replacements` must know the convention *before* it fragments the input, and
    fragmentation happens outside the main query connection.
    """
    __check_db_existence(db_name)
    with sqlite3.connect(f'file:{db_name}?mode=ro', uri=True) as con:
        return _load_schema_meta(con.cursor(), radius)


def __check_env_provenance(env, is_ring_closure, radius):
    """Guard against a fragmenter/schema convention desync on a v2 database.

    On v2 the env string carries the provenance, so a query whose env disagrees with the
    requested provenance would silently return either nothing or rows of the wrong kind.
    Raising here converts that into an explicit failure. `is_ring_closure=None` means
    "either provenance" and is always acceptable - broad make_cycle queries both variants.
    """
    if is_ring_closure is None:
        return
    if env in RADIUS0_ENV_CLASSES:
        # A radius-0 env is the attachment-point class rather than SMILES: the ring cuts
        # are stated by the leading R<x> count, not by an isotope inside the string.
        labelled = not env.startswith('R0')
    else:
        labelled = __RING_CUT_DUMMY_PATT in env
    if bool(is_ring_closure) != labelled:
        raise RuntimeError(
            f"fragment convention mismatch on a v2 database (radius{radius}): "
            f"is_ring_closure={is_ring_closure} but env {'carries' if labelled else 'lacks'} "
            f"ring-cut attachment points: {env!r}. The fragmenter must be called with "
            f"label_all_ring_cuts taken from _load_schema_meta()."
        )


def __get_replacements_rowids(db_cur, env, dist, min_atoms, max_atoms, radius, min_freq=0, set_names=None,
                              schema_meta=None, is_ring_closure=None, **kwargs):

    if schema_meta is None:
        schema_meta = _load_schema_meta(db_cur, radius)
    user_version = schema_meta['user_version']

    if user_version == 0:
        if is_ring_closure:
            raise ValueError(
                "ring_closures=True requires a v1 database with the "
                "is_ring_closure column. Rebuild the DB with cremdb_create."
            )
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

    elif user_version in (1, 2):
        def _sql_value(value):
            if isinstance(value, str):
                return "'" + value.replace("'", "''") + "'"
            if isinstance(value, bool):
                return "1" if value else "0"
            if value is None:
                return "NULL"
            return str(value)

        frags_columns = schema_meta['frags_columns']
        frags_h_columns = schema_meta['frags_h_columns']

        # Discover available set columns in the radius table.
        # core_num_atoms and dist2 are denormalized into radius{N} alongside
        # env_id / core_smi_id and must be excluded when listing set columns.
        set_names_list, is_full_set = _resolve_set_columns(schema_meta, set_names, radius)

        # Build the frequency / membership filter.
        #
        # Set columns are INTEGER NOT NULL DEFAULT 0; a row carries the
        # per-set count for each set it belongs to (0 if absent from that
        # set). The user-facing semantics:
        #
        #   * If `set_names_list` covers every available set column the user
        #     is asking for "any membership" — and any row in radius{N}
        #     belongs to at least one set, so the predicate is trivially true
        #     when min_freq <= 0 and can be dropped entirely. Dropping it is
        #     important: the denormalized lookup index becomes a covering
        #     index for the query, eliminating a per-matched-row probe into
        #     the radius{N} heap.
        #
        #   * If `set_names_list` is a strict subset, the user wants rows
        #     that are members of *those specific* sets. Membership means
        #     count > 0, so the per-set threshold is `max(1, min_freq)` even
        #     when min_freq is 0. This costs the heap probe but is the
        #     correct semantics; recovering the covering optimization for
        #     this case would require per-set partial indices.
        mf = min_freq if min_freq is not None else 0
        if is_full_set and mf <= 0:
            freq_clause = None
        else:
            # threshold = mf for the all-sets case (mf > 0 here),
            #             max(1, mf) for explicit subsets.
            threshold = mf if is_full_set else max(1, mf)
            freq_clause = " OR ".join(
                f"r.{sn} >= {_sql_value(threshold)}" for sn in set_names_list
            )
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
                        r.core_num_atoms BETWEEN {_sql_value(min_atoms)} AND {_sql_value(max_atoms)}"""
        if freq_clause is not None:
            sql += f" AND {freq_clause}"

        if dist is not None:
            if isinstance(dist, tuple):
                if len(dist) != 2:
                    raise ValueError("dist must be a single value or a tuple of two values")
                sql += f" AND r.dist2 BETWEEN {_sql_value(dist[0])} AND {_sql_value(dist[1])}"
            else:
                sql += f" AND r.dist2 = {_sql_value(dist)}"

        # Filter by fragment provenance. is_ring_closure=1 restricts to
        # ring-cut rows; 0 restricts to acyclic-cut rows; None disables the
        # filter, which is used by broad cycle generation.
        #
        # v2 needs no predicate at all: every ring cut carries isotope 1, so the env
        # string itself is provenance-specific (a labelled env can only have come from
        # a ring cut, a plain one only from an acyclic cut). The caller still passes the
        # intent, and it is checked against the env below rather than pushed into SQL.
        if schema_meta['provenance_in_env']:
            __check_env_provenance(env, is_ring_closure, radius)
        elif is_ring_closure is not None:
            radius_columns = schema_meta['radius_columns']
            if 'is_ring_closure' not in radius_columns:
                if is_ring_closure:
                    raise ValueError(
                        f"radius{radius} has no is_ring_closure column - "
                        f"this DB predates ring-closure support. "
                        f"Rebuild with cremdb_create to use ring_closures=True."
                    )
                # Legacy v1 DB without the column: skip the predicate entirely
                # (every existing row was acyclic-cut, matching is_ring_closure=0).
            else:
                sql += f" AND r.is_ring_closure = {_sql_value(is_ring_closure)}"

        for k, v in kwargs.items():
            column = f"{kwarg_target[k]}.{k}"

            if isinstance(v, tuple):
                if len(v) != 2:
                    raise ValueError(f"Value for {k} must be a single value or a tuple of two values")
                sql += f" AND {column} BETWEEN {_sql_value(v[0])} AND {_sql_value(v[1])}"
            else:
                sql += f" AND {column} = {_sql_value(v)}"
    else:
        raise NotImplementedError('Not implemented for database version other than 0, 1 and 2')
    db_cur.execute(sql)
    return set(i[0] for i in db_cur.fetchall())


def _get_replacements(db_cur, radius, row_ids, schema_meta=None, set_names=None):
    """Fetch (rowid, core_smi, core_sma, freq) for the given radius{N} row ids.

    `freq` is the occurrence count reported to users via `return_rxn_freq`. On v0 it is the
    stored `freq` column. On v1/v2 the single column was replaced by one count column per
    set, and the reported value is the largest count among the sets named by `set_names`
    (all of them when it is None) — see `__get_count_sql_expr` for why maximum and not sum.

    `set_names` must be the same value passed to `__get_replacements_rowids`, so that the
    reported frequency is the quantity the row was selected on.
    """
    if schema_meta is None:
        schema_meta = _load_schema_meta(db_cur, radius)
    user_version = schema_meta['user_version']
    if user_version == 0:
        sql = f"""SELECT rowid, core_smi, core_sma, freq
                      FROM radius{radius}
                      WHERE rowid IN ({','.join(map(str, row_ids))})"""
    elif user_version in (1, 2):
        # v2 differs from v1 only by the absence of is_ring_closure, which is not selected
        # here, so the same statement serves both. The set columns live on radius{N}, which
        # is already the driving table, so reporting the count costs no extra join.
        set_names_list, _is_full_set = _resolve_set_columns(schema_meta, set_names, radius)
        sql = f"""SELECT r.rowid, f.core_smi, {__get_count_sql_expr(set_names_list)}
                  FROM radius{radius} r
                  JOIN frags f ON r.core_smi_id = f.core_smi_id
                  WHERE r.rowid IN ({','.join(map(str, row_ids))})"""
    else:
        raise NotImplementedError('Not implemented for database version other than 0, 1 and 2')
    db_cur.execute(sql)
    if user_version in (1, 2):
        # Keep tuple shape identical to user_version 0 for compatibility: v1/v2 store no
        # SMARTS for the core, so core_smi stands in for core_sma.
        return [(row_id, core_smi, core_smi, freq)
                for row_id, core_smi, freq in db_cur.fetchall()]
    return db_cur.fetchall()


def __fragment_tuple_sort_key(item):
    context_mol = item[2]
    if isinstance(context_mol, Chem.Mol):
        context_smi = Chem.MolToSmiles(context_mol, isomericSmiles=True)
    else:
        raise TypeError(f'context_mol can only be RDKit.Mol object: {type(context_mol)}')
    return item[0], item[1], context_smi, item[3], repr(item[4]), -1 if item[5] is None else item[5]


def _normalize_replace_cycles(replace_cycles):
    if replace_cycles is False or replace_cycles == "no":
        return "no"
    if replace_cycles is True or replace_cycles == "forced":
        return "forced"
    if replace_cycles in ("partial_all", "partial_exo"):
        return replace_cycles
    raise ValueError(
        "replace_cycles must be one of False, True, 'no', 'forced', "
        "'partial_all', or 'partial_exo'"
    )


def __gen_replacements(mol1, mol2, db_name, radius, dist=None, min_size=0, max_size=8, min_rel_size=0, max_rel_size=1,
                       min_inc=-2, max_inc=2, max_replacements=None, replace_cycles=False,
                       protected_ids_1=None, protected_ids_2=None, min_freq=10, set_names=None,
                       symmetry_fixes=False, filter_func=None, sample_func=None, return_frag_smi_only=False,
                       operation="mutate", ring_closures=False, ring_size=None,
                       seed=None, **kwargs):

    rng = random.Random(seed)

    link = False
    if not isinstance(mol1, Chem.Mol):
        raise StopIteration("The first molecule in __gen_replacement always must be specified")
    if isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol):
        link = True

    if operation not in {"mutate", "link", "cycle"}:
        raise ValueError("operation must be one of: 'mutate', 'link', 'cycle'")

    replace_cycle_mode = _normalize_replace_cycles(replace_cycles)

    # The fragment convention is a property of the database schema version, so it has to be
    # known *before* the input is fragmented - the env/core strings emitted here must match
    # the ones the DB stores, and a mismatch yields silently empty results. The metadata is
    # read once here and reused for the queries below instead of being re-read inside the
    # connection.
    schema_meta = __load_convention(db_name, radius)
    label_all_ring_cuts = schema_meta['label_all_ring_cuts']
    # v0/v1 keep the historical link behaviour of not constraining provenance at all; on v2
    # the env already encodes it, so state the intent explicitly and let the guard check it.
    link_is_ring_closure = 0 if schema_meta['provenance_in_env'] else None

    # fragmentation output f should be a tuple of
    # (env: str,
    #  core: str,
    #  context_mol: Chem.Mol,
    #  num_heavy_atoms: int,
    #  query_dist2: None | int | tuple[int, int],
    #  is_ring_closure: None | int)
    # individual fragmentation functions return different data shape, that are reshaped afterwards as needed
    if operation == "mutate":
        if link:
            raise ValueError("mutate operation expects a single molecule")
        mol = mol1
        mol_hac = mol.GetNumHeavyAtoms()
        lower_core_atoms = max(min_size, math.ceil(min_rel_size * mol_hac))
        upper_core_atoms = min(max_size, math.floor(max_rel_size * mol_hac))
        f = [
            (*frag, None, 0)
            for frag in __fragment_mol(
                mol,
                radius,
                protected_ids=protected_ids_1,
                symmetry_fixes=symmetry_fixes,
                min_core_atoms=lower_core_atoms,
                max_core_atoms=upper_core_atoms,
                include_cyclic_cores=(replace_cycle_mode == "forced"),
            )
        ]
        if replace_cycle_mode in {"partial_all", "partial_exo"}:
            side_cut_mode = "all" if replace_cycle_mode == "partial_all" else "exo"
            f.extend(
                (*frag, None, 1)
                for frag in __fragment_mol_partial_cycles(
                    mol,
                    radius=radius,
                    protected_ids=protected_ids_1,
                    symmetry_fixes=symmetry_fixes,
                    min_core_atoms=lower_core_atoms,
                    max_core_atoms=upper_core_atoms,
                    include_cyclic_cores=False,
                    side_cut_mode=side_cut_mode,
                    label_all_ring_cuts=label_all_ring_cuts,
                )
            )
    elif operation == "cycle":
        if link:
            raise ValueError("cycle operation expects a single molecule")
        mol = mol1
        if ring_closures:
            # Strict cycle mode: arc-cut fragmenter only, in the DB's own convention.
            f = [
                (*frag, 1)
                for frag in __fragment_mol_ring_closure(
                    mol=mol,
                    radius=radius,
                    ring_size=ring_size,
                    protected_ids=protected_ids_1,
                    label_variants=(label_all_ring_cuts,),
                )
            ]
        else:
            # Broad cycle mode: union of both fragmenters so the DB can return
            # both arc-cut (connected-env) and acyclic-cut (disconnected-env)
            # rows. On v0/v1 the SQL-side is_ring_closure filter is disabled; on v2
            # provenance lives in the env, so "either provenance" becomes a union over
            # the labelled and plain env variants of each fragmenter.
            broad_variants = (True, False) if label_all_ring_cuts else (False,)
            f = [
                (*frag, None)
                for frag in __fragment_mol_macrocycle(
                    mol=mol,
                    radius=radius,
                    ring_size=ring_size,
                    protected_ids=protected_ids_1,
                    label_variants=broad_variants,
                )
            ]
            f.extend(
                (*frag, None)
                for frag in __fragment_mol_ring_closure(
                    mol=mol,
                    radius=radius,
                    ring_size=ring_size,
                    protected_ids=protected_ids_1,
                    label_variants=broad_variants,
                )
            )
    elif operation == "link":
        if not link:
            raise ValueError("link operation expects two molecules")
        f = __fragment_mol_link(mol1=mol1, mol2=mol2, radius=radius, protected_ids_1=protected_ids_1,
                                protected_ids_2=protected_ids_2)
        f = [(*frag, None, link_is_ring_closure) for frag in f]
        mol = Chem.CombineMols(mol1, mol2)

    if not f:
        return

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
            f.sort(key=__fragment_tuple_sort_key)  # canonical order before seeded shuffle
            rng.shuffle(f)
            preliminary_return = max_replacements // len(f)
            if preliminary_return == 0:
                preliminary_return = 1

        for frag_tuple in f:
            env, core, context_mol, num_heavy_atoms, query_dist2, is_ring_closure = frag_tuple
            effective_dist = query_dist2 if query_dist2 is not None else dist

            min_atoms = num_heavy_atoms + min_inc
            max_atoms = num_heavy_atoms + max_inc

            row_ids = __get_replacements_rowids(cur, env, effective_dist, min_atoms, max_atoms, radius, min_freq,
                                                set_names=set_names, schema_meta=schema_meta,
                                                is_ring_closure=is_ring_closure, **kwargs)

            if filter_func:
                row_ids = set(filter_func(row_ids, cur, radius))

            if max_replacements is None:
                res = _get_replacements(cur, radius, row_ids, schema_meta=schema_meta,
                                        set_names=set_names)
            else:
                n = min(len(row_ids), preliminary_return)
                if sample_func is not None:
                    selected_row_ids = sample_func(list(row_ids), cur, radius, n)
                else:
                    selected_row_ids = rng.sample(sorted(row_ids), n)
                row_ids.difference_update(selected_row_ids)
                replacements.update({i: (core, context_mol) for i in row_ids})
                res = _get_replacements(cur, radius, selected_row_ids, schema_meta=schema_meta,
                                        set_names=set_names)

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
                selected_row_ids = rng.sample(sorted(replacements.keys()), n)
            res = _get_replacements(cur, radius, selected_row_ids, schema_meta=schema_meta,
                                    set_names=set_names)
            for row_id, core_smi, _, freq in res:
                src_core, src_context = replacements[row_id]
                if core_smi != src_core:
                    if return_frag_smi_only:
                        yield core_smi
                    else:
                        yield src_core, core_smi, freq, src_context


def __frag_replace_mp(items):
    # return smi, transformation, transformation_freq, mol
    # data generators yield (mol1, mol2, old_frag, new_frag, radius, context_mol, freq,
    # discard_ring_geometry); the last two items are unpacked here
    *args, freq, discard_ring_geometry = items
    return [(*item, freq) for item in __frag_replace(*args, discard_ring_geometry=discard_ring_geometry)]


def __get_data(mol, db_name, radius, min_size, max_size, min_rel_size, max_rel_size, min_inc, max_inc,
               replace_cycles, protected_ids, min_freq, set_names, max_replacements, symmetry_fixes, filter_func=None,
               sample_func=None, seed=None, discard_ring_geometry=True, **kwargs):
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
                                                                    return_frag_smi_only=False,
                                                                    operation="mutate",
                                                                    seed=seed, **kwargs):
        yield mol, None, frag_sma, core_sma, radius, context_mol, freq, discard_ring_geometry


def __get_data_link(mol1, mol2, db_name, radius, dist, min_atoms, max_atoms, protected_ids_1, protected_ids_2, min_freq,
                    set_names, max_replacements, filter_func=None, sample_func=None, seed=None,
                    discard_ring_geometry=True, **kwargs):
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
                                                                     operation="link",
                                                                     return_frag_smi_only=False, seed=seed, **kwargs):
        yield mol1, mol2, frag_sma, core_sma, radius, context_mol, freq, discard_ring_geometry


def __get_data_cycle(mol, db_name, radius, ring_size, ring_closures, min_size, max_size, protected_ids,
                     min_freq, set_names, max_replacements, filter_func=None, sample_func=None,
                     seed=None, discard_ring_geometry=True, **kwargs):
    for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol, mol2=None, db_name=db_name,
                                                                    radius=radius,
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
                                                                    operation="cycle",
                                                                    ring_closures=ring_closures,
                                                                    ring_size=ring_size,
                                                                    seed=seed, **kwargs):
        yield mol, None, frag_sma, core_sma, radius, context_mol, freq, discard_ring_geometry


def mutate_mol(mol, db_name, radius=3, min_size=0, max_size=10, min_rel_size=0, max_rel_size=1, min_inc=-2, max_inc=2,
               max_replacements=None, replace_cycles="no", replace_ids=None, protected_ids=None,
               symmetry_fixes=False, min_freq=0, return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1,
               filter_func=None, sample_func=None, set_names=None, seed=None, discard_ring_geometry=True, **kwargs):
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
    :param replace_cycles: controls replacement of cyclic source fragments.
                           ``"no"``/False uses ordinary acyclic-cut mutation
                           only. ``"forced"``/True allows cyclic cores from
                           ordinary fragmentation to be replaced ignoring the size filters.
                           ``"partial_all"`` additionally replaces partial
                           ring arcs using exhaustive side cuts.
                           ``"partial_exo"`` additionally replaces partial
                           ring arcs using only exo side cuts adjacent to the
                           selected ring arc. Default: ``"no"``.
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
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB,
                            that is the number of occurrences of the replacement fragment. On a database
                            with several sets this is the maximum count across the sets selected by
                            `set_names` (all of them by default), which matches the `min_freq` filter -
                            a fragment is selected when any one set reaches the threshold. Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule.  Default: False.
                       In the returned Mol, atoms inserted from the replacement fragment carry a boolean
                       property ``__crem`` set to True (other atoms do not).
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
    :param seed: random seed for reproducible fragment selection when max_replacements is set. Default: None.
    :param discard_ring_geometry: if True (default), silently discard ring-forming products (here, partial
                     ring-arc replacements selected via replace_cycles) whose newly created ring cannot be
                     embedded in 3D, e.g. a too-short bridge across the meta/para positions of a six-membered
                     aromatic ring or across the sulfur of a thiophene. Set False to return all generated
                     structures without this geometric filtering. Non ring-forming transformations are
                     unaffected. Default: True.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional; maximum count across the selected sets), RDKit Mol
             object (optional).
             Only entries with distinct SMILES will be returned.

    Note: supply RDKit Mol object with explicit hydrogens if H replacement is required

    """

    replace_cycles = _normalize_replace_cycles(replace_cycles)

    __check_db_existence(db_name)
    products = {Chem.MolToSmiles(Chem.RemoveHs(mol))}
    mol = __backup_atom_properties(mol, __atom_properties_to_backup)

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
                                                                        return_frag_smi_only=False,
                                                                        operation="mutate",
                                                                        seed=seed, **kwargs):
            for smi, m, rxn in __frag_replace(mol, None, frag_sma, core_sma, radius, context_mol,
                                              discard_ring_geometry=discard_ring_geometry):
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
                                                              sample_func=sample_func,
                                                              seed=seed,
                                                              discard_ring_geometry=discard_ring_geometry, **kwargs),
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
             return_mol=False, ncores=1, filter_func=None, sample_func=None, set_names=None, seed=None, **kwargs):
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
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB,
                            that is the number of occurrences of the replacement fragment. On a database
                            with several sets this is the maximum count across the sets selected by
                            `set_names` (all of them by default), which matches the `min_freq` filter -
                            a fragment is selected when any one set reaches the threshold. Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule.  Default: False.
                       In the returned Mol, atoms inserted from the replacement fragment carry a boolean
                       property ``__crem`` set to True (other atoms do not).
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
    :param seed: random seed for reproducible fragment selection when max_replacements is set. Default: None.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional; maximum count across the selected sets), RDKit Mol
             object (optional).
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
                      sample_func=sample_func, seed=seed, **kwargs)


def link_mols(mol1, mol2, db_name, radius=3, dist=None, min_atoms=1, max_atoms=2, max_replacements=None,
              replace_ids_1=None, replace_ids_2=None, protected_ids_1=None, protected_ids_2=None,
              min_freq=0, return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
              sample_func=None, set_names=None, seed=None, **kwargs):
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
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB,
                            that is the number of occurrences of the replacement fragment. On a database
                            with several sets this is the maximum count across the sets selected by
                            `set_names` (all of them by default), which matches the `min_freq` filter -
                            a fragment is selected when any one set reaches the threshold. Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule.  Default: False.
                       In the returned Mol, atoms inserted from the replacement fragment carry a boolean
                       property ``__crem`` set to True (other atoms do not).
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
    :param seed: random seed for reproducible fragment selection when max_replacements is set. Default: None.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over new molecules. If no additional return arguments were called this would be a generator over
             SMILES of new molecules. If any of additional return values were asked the function will return a list
             of list where the first item is SMILES, then rxn string of a transformation (optional), frequency of
             fragment occurrence in the DB (optional; maximum count across the selected sets), RDKit Mol
             object (optional).
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

    mol1 = __backup_atom_properties(Chem.AddHs(mol1), __atom_properties_to_backup)
    mol2 = __backup_atom_properties(Chem.AddHs(mol2), __atom_properties_to_backup)

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
                                                                         return_frag_smi_only=False,
                                                                         operation="link",
                                                                         seed=seed, **kwargs):
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
                                                                   sample_func=sample_func, seed=seed, **kwargs),
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


def make_cycle2(*args, **kwargs):
    """
    Convenience function which can be used to process molecules in parallel using multiprocessing module.
    It calls make_cycle which cannot be used directly in multiprocessing because it is a generator

    :param args: positional arguments, the same as in make_cycle function
    :param kwargs: keyword arguments, the same as in make_cycle function
    :return: list with output molecules

    """
    return list(make_cycle(*args, **kwargs))


def get_replacements(mol1, db_name, radius, mol2=None, dist=None, min_size=0, max_size=8, min_rel_size=0,
                     max_rel_size=1, min_inc=-2, max_inc=2, max_replacements=None, replace_cycles="no",
                     protected_ids_1=None, protected_ids_2=None, replace_ids_1=None,
                     replace_ids_2=None, min_freq=0, symmetry_fixes=False, filter_func=None, sample_func=None,
                     return_frag_smi_only=True,
                     set_names=None, seed=None, **kwargs):
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
    :param replace_cycles: controls replacement of cyclic source fragments for
                           single-molecule searches. ``"no"``/False uses
                           ordinary acyclic-cut mutation only.
                           ``"forced"``/True allows cyclic cores from
                           ordinary fragmentation to be replaced ignoring the size filters.
                           ``"partial_all"`` additionally searches partial
                           ring arcs with exhaustive side cuts.
                           ``"partial_exo"`` additionally searches partial
                           ring arcs with only exo side cuts. Ignored for
                           link searches. Default: ``"no"``.

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
                                 further passed to `get_mols_from_replacements`. `freq` is the number of occurrences
                                 of the replacement fragment in the DB - the maximum count across the sets selected
                                 by `set_names`.
    :param seed: random seed for reproducible fragment selection when max_replacements is set. Default: None.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment. This can be useful to annotate
                     fragments with additional custom properties (e.g. number of particular pharmacophore features,
                     lipophilicity, etc) and use these parameters to additionally restrict selected fragments.
    :return: generator over smiles of fragments in a DB which satisfy given criteria
    """

    replace_cycles = _normalize_replace_cycles(replace_cycles)

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

    mol1 = __backup_atom_properties(mol1, __atom_properties_to_backup)
    if isinstance(mol2, Chem.Mol):
        mol2 = __backup_atom_properties(mol2, __atom_properties_to_backup)

    for res in __gen_replacements(mol1=mol1, mol2=mol2, db_name=db_name, radius=radius, dist=dist,
                                  min_size=min_size, max_size=max_size, min_rel_size=min_rel_size,
                                  max_rel_size=max_rel_size, min_inc=min_inc, max_inc=max_inc,
                                  max_replacements=max_replacements, replace_cycles=replace_cycles,
                                  protected_ids_1=protected_ids_1, protected_ids_2=protected_ids_2,
                                  min_freq=min_freq, set_names=set_names, symmetry_fixes=symmetry_fixes,
                                  filter_func=filter_func, sample_func=sample_func,
                                  return_frag_smi_only=return_frag_smi_only,
                                  operation=("link" if isinstance(mol2, Chem.Mol) else "mutate"),
                                  seed=seed, **kwargs):
        if return_frag_smi_only:
            yield res
        else:
            src_core, repl_core, freq, context_mol = res
            yield src_core, repl_core, freq, __prepare_context_mol_for_output(context_mol)


def get_mols_from_replacements(mol1, radius, replacements, mol2=None, return_rxn=False, return_rxn_freq=False,
                               return_mol=False, discard_ring_geometry=True):
    """
    Build product molecules from replacements previously obtained with
    ``get_replacements(..., return_frag_smi_only=False)``.

    :param mol1: the first RDKit Mol object (the same molecule passed to get_replacements).
    :param radius: context radius; must match the value used in get_replacements.
    :param replacements: iterable of (source_core_smi, replacement_core_smi, freq, context_mol) tuples.
    :param mol2: the second RDKit Mol object for link operations, otherwise None.
    :param return_rxn: whether to additionally return the rxn string of a transformation. Default: False.
    :param return_rxn_freq: whether to additionally return the frequency of a transformation, taken from the
                            supplied replacements tuples (maximum count across sets, see `get_replacements`).
                            Default: False.
    :param return_mol: whether to additionally return the RDKit Mol object of a generated molecule. Default: False.
                       In the returned Mol, atoms inserted from the replacement fragment carry a boolean
                       property ``__crem`` set to True (other atoms do not).
    :param discard_ring_geometry: if True (default), silently discard ring-closure products whose newly created
                     ring cannot be embedded in 3D (a too-short bridge across the meta/para positions of a
                     six-membered aromatic ring or across the sulfur of a thiophene). Set False to keep all
                     products without geometric filtering. Default: True.
    :return: generator over new molecules - SMILES, optionally followed by the rxn string, frequency,
             and/or the RDKit Mol object. Only entries with distinct SMILES are returned.
    """

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

        for smi, m, rxn in __frag_replace(mol1, mol2, frag_sma, core_sma, radius, context_mol,
                                          discard_ring_geometry=discard_ring_geometry):
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


def make_cycle(mol, db_name, radius=3, ring_size=None, ring_closures=True,
               min_atoms=1, max_atoms=10, max_replacements=None,
               replace_ids=None, protected_ids=None, symmetry_fixes=False, min_freq=0,
               return_rxn=False, return_rxn_freq=False, return_mol=False, ncores=1, filter_func=None,
               sample_func=None, set_names=None, seed=None, discard_ring_geometry=True, **kwargs):
    """
    Generate new rings (macrocycles or smaller native cycles) by linking two
    atoms in the same molecule with a 2-attachment-point fragment from the DB.

    Two complementary modes:

    * ``ring_closures=False`` (broad): query **any** linker fragment.
      Internally both fragmenters are run on the input molecule (the
      connected-env arc-cut fragmenter and the disconnected-env macrocycle
      fragmenter) and the ``is_ring_closure`` provenance column is **not**
      filtered, so DB rows of either provenance can match.
    * ``ring_closures=True`` (strict): only the connected-env arc-cut
      fragmenter runs and the query is restricted to ``is_ring_closure=1``
      rows (populated by ``--frag-mode ring`` / ``both`` or the corresponding
      ``*_optimal`` modes at DB build time).
      Useful for closing native (typically aliphatic) rings.

    :param mol: RDKit Mol object.
    :param db_name: path to DB file with fragment replacements.
    :param radius: radius of context which will be considered for replacement. Default: 3.
    :param ring_size: size of the *new* ring being formed (in atoms = bonds).
                      ``int`` for a single size, ``(min, max)`` tuple for a
                      window. ``None`` imposes no ring-size constraint. The
                      per-anchor-pair ``dist2`` filter is derived as
                      ``ring_size − d_in`` where ``d_in`` is the topological
                      distance between the two anchor heavy atoms in the
                      input molecule.
    :param ring_closures: if True, query ring-closure (arc) fragments in DB
                          (rows with ``is_ring_closure = 1``). If False
                          (default) query acyclic-cut linker fragments.
    :param min_atoms: minimum number of heavy atoms in the linker fragment. Default: 1.
    :param max_atoms: maximum number of heavy atoms in the linker fragment. Default: 10.
    :param max_replacements: maximum number of replacements to make. If the number of replacements available in DB is
                             greater than the specified value the specified number of randomly chosen replacements
                             will be applied. Default: None.
    :param replace_ids: iterable with ids of heavy atom with replaceable Hs or/and ids of H atoms to replace,
                        it has lower priority over `protected_ids` (replace_ids
                        which are present in protected_ids would be protected). Default: None.
    :param protected_ids: iterable with ids of heavy atoms at which no H replacement should be made and/or ids of
                          protected hydrogens. This argument has a higher priority over `replace_ids`. Default: None.
    :param symmetry_fixes: accepted for API compatibility with mutate/grow functions but not used here.
    :param min_freq: minimum occurrence of fragments in DB for replacement. Default: 0.
    :param return_rxn: whether to additionally return rxn of a transformation. Default: False.
    :param return_rxn_freq: whether to additionally return the frequency of a transformation in the DB,
                            that is the number of occurrences of the replacement fragment. On a database
                            with several sets this is the maximum count across the sets selected by
                            `set_names` (all of them by default), which matches the `min_freq` filter -
                            a fragment is selected when any one set reaches the threshold. Default: False.
    :param return_mol: whether to additionally return RDKit Mol object of a generated molecule. Default: False.
                       In the returned Mol, atoms inserted from the replacement fragment carry a boolean
                       property ``__crem`` set to True (other atoms do not).
    :param ncores: number of cores. Default: 1.
    :param filter_func: a function which will filter selected fragments by additional rules
                        (in this way one may add arbitrary selection constrains). The function takes necessary first
                        three arguments: row_ids (list or set of row_ids from the fragment database supplied to
                        make_cycle), cursor of that fragment database and radius (int). This is required to
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
    :param seed: random seed for reproducible fragment selection when max_replacements is set. Default: None.
    :param discard_ring_geometry: if True (default), silently discard ring-forming products (here, the new ring formed
                     by make_cycle) whose newly created ring cannot be
                     embedded in 3D, e.g. a too-short bridge across the meta/para positions of a six-membered
                     aromatic ring or across the sulfur of a thiophene. Set False to return all generated
                     structures without this geometric filtering. Non ring-forming transformations are
                     unaffected. Default: True.
    :param **kwargs: named arguments to additionally filter replacing fragments. For v0 DB use columns from radiusX,
                     for v1 DB use columns from frags or frags_h. Values are a single value or 2-item tuple with lower
                     and upper bound of the corresponding parameter of a fragment.
    :return: generator over new molecules. If no additional return arguments were requested this is a generator over
             SMILES of new molecules. If additional return values were requested, the function yields a list where
             the first item is SMILES, then rxn string (optional), frequency (optional; maximum count
             across the selected sets), RDKit Mol object (optional).
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
    mol = __backup_atom_properties(mol, __atom_properties_to_backup)

    if ncores == 1:

        for frag_sma, core_sma, freq, context_mol in __gen_replacements(mol1=mol, mol2=None, db_name=db_name,
                                                                        radius=radius,
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
                                                                        operation="cycle",
                                                                        ring_closures=ring_closures,
                                                                        ring_size=ring_size,
                                                                        seed=seed, **kwargs):
            for smi, m, rxn in __frag_replace(mol, None, frag_sma, core_sma, radius, context_mol,
                                              discard_ring_geometry=discard_ring_geometry):
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
            for items in p.imap(__frag_replace_mp, __get_data_cycle(mol, db_name, radius, ring_size,
                                                                    ring_closures, min_atoms, max_atoms,
                                                                    protected_ids, min_freq, set_names,
                                                                    max_replacements,
                                                                    filter_func=filter_func,
                                                                    sample_func=sample_func, seed=seed,
                                                                    discard_ring_geometry=discard_ring_geometry,
                                                                    **kwargs),
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
