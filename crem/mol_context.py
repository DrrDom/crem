import re
from itertools import product, permutations, combinations
from collections import defaultdict
from rdkit import Chem
from .functions import mol_to_smarts
from .ring_fragments import RING_CUT_DUMMY_ISOTOPE

__author__ = 'pavel'

patt_remove_map = re.compile(r"\[(?:\d+)?\*:[0-9]+\]")    # to change CC([*:1])O to CC([*])O
patt_remove_h = re.compile(r"(?<!\[)H[1-9]*(?=:[0-9])")   # to remove H after atoms with maps: [CH2:1] to [C:1], but not touching [H] or [nH]


# --------------------------------------------------------------------------- #
# radius-0 attachment-point classes
#
# At radius 0 no context is recorded, so all that remains of an environment is how many
# attachment points a fragment has and which of them close a ring. That is encoded as
# `R<x>A<y>`: x ring-cut points, y acyclic ones. Ring cuts always come in a pair (two ring
# bonds are always cut together), so x is 0 or 2 - verified against all 32 769 176 rows of
# the chembl36 `frags` table, where n_ring is never anything else. With n_att capped at 4
# that leaves exactly seven classes.
#
# The class string is the radius-0 `env`. It is deliberately not valid SMILES: it cannot be
# mistaken for a real env, and unlike a synthetic `[1*:1].[*:2]` form it states the class
# counts without implying a map numbering.
# --------------------------------------------------------------------------- #
RADIUS0_ENV_CLASSES = ('R0A1', 'R0A2', 'R0A3', 'R0A4', 'R2A0', 'R2A1', 'R2A2')


def __radius0_env(n_ring, n_acyclic):
    """`env` string for a radius-0 attachment-point class."""
    return f'R{n_ring}A{n_acyclic}'


def __attachment_point_classes(mol):
    """Split a fragment's attachment points into (ring_cut_ids, acyclic_ids).

    Ring cuts are the dummies carrying RING_CUT_DUMMY_ISOTOPE - the v2 convention applied
    in crem.ring_fragments._materialize_fragment.
    """
    ring, acyclic = [], []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
            (ring if a.GetIsotope() == RING_CUT_DUMMY_ISOTOPE else acyclic).append(a.GetIdx())
    return ring, acyclic


def __radius0_fake_context(core):
    """A stand-in context for radius-0 standardisation: one aliphatic carbon per attachment
    point, carrying that point's class label.

    Radius 0 keeps no context, but the standardisation pipeline is driven by one - it
    derives the attachment-point numbering, the symmetry orbits and hence the set of
    equivalent labellings from the env. Feeding it a context in which every point of a
    class sits in an identical component makes all points of that class one orbit, so
    get_std_context_core_permutations returns exactly the within-class relabellings,
    canonically numbered, with no radius-0-specific standardisation code. Ring and acyclic
    points stay in separate orbits because the isotope is part of each component's
    canonical SMILES.

    Returns `(context_mol, env)`, where `env` is the class string that replaces the fake
    context's own env in the result.
    """
    ring, acyclic = __attachment_point_classes(core)
    parts = []
    for idx in ring + acyclic:
        a = core.GetAtomWithIdx(idx)
        iso = RING_CUT_DUMMY_ISOTOPE if a.GetIsotope() == RING_CUT_DUMMY_ISOTOPE else ''
        parts.append(f'C[{iso}*:{a.GetAtomMapNum()}]')
    return Chem.MolFromSmiles('.'.join(parts)), __radius0_env(len(ring), len(acyclic))


def __get_submol(mol, atom_ids):
    bond_ids = []
    for pair in combinations(atom_ids, 2):
        b = mol.GetBondBetweenAtoms(*pair)
        if b:
            bond_ids.append(b.GetIdx())
    m = Chem.PathToSubmol(mol, bond_ids)
    m.UpdatePropertyCache()
    return m


def __bonds_to_atoms(mol, bond_ids):
    output = []
    for i in bond_ids:
        b = mol.GetBondWithIdx(i)
        output.append(b.GetBeginAtom().GetIdx())
        output.append(b.GetEndAtom().GetIdx())
    return tuple(set(output))


def __get_context_env(mol, radius):
    """
    INPUT:
        mol - Mol object containing chain(s) of molecular context
        radius - integer, number of bonds to cut context
    OUTPUT:
        Mol containing only atoms within the specified radius from the attachment point(s).
        All explicit Hs will be stripped.
    """
    # mol is context consisting of one or more groups with single attachment point

    m = Chem.RemoveHs(mol)
    m = Chem.RWMol(m)

    bond_ids = set()
    for a in m.GetAtoms():
        if a.GetSymbol() == "*":
            i = radius
            b = Chem.FindAtomEnvironmentOfRadiusN(m, i, a.GetIdx())
            while not b and i > 0:
                i -= 1
                b = Chem.FindAtomEnvironmentOfRadiusN(m, i, a.GetIdx())
            bond_ids.update(b)

    atom_ids = set(__bonds_to_atoms(m, bond_ids))

    dummy_atoms = []

    for a in m.GetAtoms():
        if a.GetIdx() not in atom_ids:
            nei_ids = set(na.GetIdx() for na in a.GetNeighbors())
            intersect = nei_ids & atom_ids
            if intersect:
                dummy_atom_bonds = []
                for ai in intersect:
                    dummy_atom_bonds.append((ai, m.GetBondBetweenAtoms(a.GetIdx(), ai).GetBondType()))
                dummy_atoms.append(dummy_atom_bonds)

    for data in dummy_atoms:
        dummy_id = m.AddAtom(Chem.Atom(0))
        for atom_id, bond_type in data:
            m.AddBond(dummy_id, atom_id, bond_type)
        atom_ids.add(dummy_id)

    m = __get_submol(m, atom_ids)

    return m


def __replace_att(mol, repl_dict):
    for a in mol.GetAtoms():
        map_num = a.GetAtomMapNum()
        if map_num in repl_dict:
            a.SetAtomMapNum(repl_dict[map_num])


def __clear_ignored_isotopes(mol, keep_stereo=False, preserve_dummy_isotopes=False):
    if keep_stereo:
        return
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0 or not preserve_dummy_isotopes:
            atom.SetIsotope(0)


def __mol_to_smiles(mol, keep_stereo=False, preserve_dummy_isotopes=False, **kwargs):
    if keep_stereo or not preserve_dummy_isotopes:
        return Chem.MolToSmiles(mol, isomericSmiles=keep_stereo, **kwargs)

    # RDKit omits isotopes when isomericSmiles=False. This branch removes
    # stereochemistry explicitly, then uses isomeric output only to keep
    # isotope labels on dummy atoms such as [1*:2].
    tmp = Chem.Mol(mol)
    Chem.RemoveStereochemistry(tmp)
    __clear_ignored_isotopes(
        tmp,
        keep_stereo=False,
        preserve_dummy_isotopes=preserve_dummy_isotopes,
    )
    return Chem.MolToSmiles(tmp, isomericSmiles=True, **kwargs)


def __compute_att_keys(env, keep_stereo=False, preserve_dummy_isotopes=False):
    """
    Compute a sort key for every attachment point (* atom with non-zero map)
    in env. Returns a list of (key, original_map) tuples.

    Each key is `(component_canonical_smiles, intra_component_atom_rank)`.

    For the historical disconnected-env case (every component has exactly
    one * atom), the intra-component rank is irrelevant — the * is the only
    candidate inside its component — so ordering reduces to sorting by
    component canonical SMILES. This matches the behaviour of the original
    `__get_maps_and_ranks` implementation and preserves compatibility with
    pre-existing fragment databases that were built against that ordering.

    For the connected case (one component with multiple * atoms — produced
    by ring-closure DB rows or by close-anchor input fragmentation), all *s
    share the component SMILES, so the intra-component canonical atom rank
    breaks the order. Symmetry-equivalent *s within a component receive
    identical keys and form a single orbit.
    """
    tmp = Chem.Mol(env)
    out = []
    for comp in Chem.GetMolFrags(tmp, asMols=True, sanitizeFrags=False):
        __clear_ignored_isotopes(
            comp,
            keep_stereo,
            preserve_dummy_isotopes=preserve_dummy_isotopes,
        )
        star_atoms = [a for a in comp.GetAtoms()
                      if a.GetAtomicNum() == 0 and a.GetAtomMapNum()]
        if not star_atoms:
            continue
        # Zero out maps before computing canonical SMILES / ranks so the
        # map numbers themselves don't bias the ordering.
        recorded = [(a.GetIdx(), a.GetAtomMapNum()) for a in star_atoms]
        for a in star_atoms:
            a.SetAtomMapNum(0)
        comp_smi = __mol_to_smiles(
            comp,
            keep_stereo,
            preserve_dummy_isotopes=preserve_dummy_isotopes,
        )
        if len(star_atoms) > 1:
            if comp.NeedsUpdatePropertyCache():
                comp.UpdatePropertyCache()
            atom_ranks = list(Chem.CanonicalRankAtoms(comp, breakTies=False))
            for idx, original_map in recorded:
                out.append(((comp_smi, atom_ranks[idx]), original_map))
        else:
            # Single-* component — atom rank is degenerate; pin to 0 so the
            # key reduces to (comp_smi, 0). Identical components produce
            # identical keys → orbit, matching the historical behaviour.
            for _, original_map in recorded:
                out.append(((comp_smi, 0), original_map))
    return out


def __standardize_att_by_env(env, core, keep_stereo=False, preserve_dummy_isotopes=False):
    """
    Set attachment point map numbers in core and env based on canonical sort
    keys of the * atoms in env. Map numbers are assigned 1..N in sorted-key
    order; within a key (symmetry orbit) ties are broken deterministically
    by sorting on the original map number. Modifies env and core in place.
    Returns {original_map: new_map}.
    """
    pairs = __compute_att_keys(
        env,
        keep_stereo,
        preserve_dummy_isotopes=preserve_dummy_isotopes,
    )
    sorted_pairs = sorted(pairs, key=lambda kp: (kp[0], kp[1]))
    new_att = {original_map: i + 1 for i, (_, original_map) in enumerate(sorted_pairs)}
    __replace_att(core, new_att)
    __replace_att(env, new_att)
    return new_att


def __get_att_permutations(env, keep_stereo=False, preserve_dummy_isotopes=False):
    """
    Return possible permutations of attachment-point map numbers consistent
    with env's symmetry. Each returned dict maps current-map to permuted-map.
    Suitable for use after __standardize_att_by_env.
    """
    pairs = __compute_att_keys(
        env,
        keep_stereo,
        preserve_dummy_isotopes=preserve_dummy_isotopes,
    )
    if not pairs:
        return ({},)
    by_key = defaultdict(list)
    for key, m in pairs:
        by_key[key].append(m)
    orbits = [sorted(group) for _, group in sorted(by_key.items())]
    per_orbit = [[dict(zip(orbit, perm)) for perm in permutations(orbit)]
                 for orbit in orbits]
    return tuple(__merge_dicts(*item) for item in product(*per_orbit))


def __permute_att(mol, d):
    new_mol = Chem.Mol(mol)
    for a in new_mol.GetAtoms():
        i = a.GetAtomMapNum()
        if i in d:
            a.SetAtomMapNum(d[i])
    return new_mol


def __merge_dicts(*dicts):
    res = dicts[0].copy()
    for item in dicts[1:]:
        res.update(item)
    return res


def __att_point_smiles(isotope, atom_map):
    if isotope:
        return f"[{isotope}*:{atom_map}]"
    return f"[*:{atom_map}]"


def __standardize_smiles_with_att_points(mol, keep_stereo=False, preserve_dummy_isotopes=False):
    """
    to avoid different order of atoms in SMILES with different map number of attachment points

    smi = ["ClC1=C([*:1])C(=S)C([*:2])=C([*:3])N1",
           "ClC1=C([*:1])C(=S)C([*:3])=C([*:2])N1",
           "ClC1=C([*:2])C(=S)C([*:1])=C([*:3])N1",
           "ClC1=C([*:2])C(=S)C([*:3])=C([*:1])N1",
           "ClC1=C([*:3])C(=S)C([*:1])=C([*:2])N1",
           "ClC1=C([*:3])C(=S)C([*:2])=C([*:1])N1"]

    these will produce different output with RDKit MolToSmiles():
        S=c1c([*:1])c(Cl)[nH]c([*:3])c1[*:2]
        S=c1c([*:1])c(Cl)[nH]c([*:2])c1[*:3]
        S=c1c([*:1])c([*:3])[nH]c(Cl)c1[*:2]
        S=c1c([*:2])c(Cl)[nH]c([*:1])c1[*:3]
        S=c1c([*:1])c([*:2])[nH]c(Cl)c1[*:3]
        S=c1c([*:2])c([*:1])[nH]c(Cl)c1[*:3]

    output of this function
        S=c1c([*:2])c([*:3])[nH]c(Br)c1[*:1]
        S=c1c([*:3])c([*:2])[nH]c(Br)c1[*:1]
        S=c1c([*:1])c([*:3])[nH]c(Br)c1[*:2]
        S=c1c([*:3])c([*:1])[nH]c(Br)c1[*:2]
        S=c1c([*:1])c([*:2])[nH]c(Br)c1[*:3]
        S=c1c([*:2])c([*:1])[nH]c(Br)c1[*:3]

    https://sourceforge.net/p/rdkit/mailman/message/35862258/
    """

    # update property cache if needed
    if mol.NeedsUpdatePropertyCache():
        mol.UpdatePropertyCache()

    __clear_ignored_isotopes(
        mol,
        keep_stereo,
        preserve_dummy_isotopes=preserve_dummy_isotopes,
    )

    # store original maps and remove map numbers from mol
    backup_atom_map = "backupAtomMap"
    for a in mol.GetAtoms():
        atom_map = a.GetAtomMapNum()
        if atom_map:
            a.SetIntProp(backup_atom_map, atom_map)
            a.SetAtomMapNum(0)

    # get canonical ranks for atoms for a mol without maps
    atoms = list(zip(list(Chem.CanonicalRankAtoms(mol)), [a.GetIdx() for a in mol.GetAtoms()]))
    atoms.sort()

    # set new atom maps based on canonical order
    rep = {}
    atom_map = 1
    for pos, atom_idx in atoms:
        a = mol.GetAtomWithIdx(atom_idx)
        if a.HasProp(backup_atom_map):
            a.SetAtomMapNum(atom_map)
            isotope = a.GetIsotope() if keep_stereo or preserve_dummy_isotopes else 0
            original_map = a.GetIntProp(backup_atom_map)
            rep[__att_point_smiles(isotope, atom_map)] = __att_point_smiles(
                isotope, original_map
            )
            atom_map += 1

    # get SMILES and relabel with original map numbers
    s = __mol_to_smiles(
        mol,
        keep_stereo,
        preserve_dummy_isotopes=preserve_dummy_isotopes,
    )
    if rep:
        rep = dict((re.escape(k), v) for k, v in rep.items())
        patt = re.compile("|".join(rep.keys()))
        s = patt.sub(lambda m: rep[re.escape(m.group(0))], s)

    return s


def get_std_context_core_permutations(context, core, radius, keep_stereo, return_att_map=False,
                                      preserve_dummy_isotopes=True):
    """
    INPUT:
        context - Mol or SMILES containing full chain(s) of a context with labeled attachment point(s),
                  if context is absent (e.g.for radius 0) specify empty string or empty Mol
        core    - Mol or SMILES of a core fragment with labeled attachment point(s)
        keep_stereo - boolean to keep stereo information in output
        radius  - integer (0, 1, 2, etc), number of bonds to cut context
        preserve_dummy_isotopes - keep isotope labels on dummy atoms when keep_stereo is False.
                  Ring-cut attachment points carry isotope 1 (see crem.ring_fragments), so
                  stripping them merges ring-arc fragments into the acyclic string space.
                  Defaults to True because that is the non-lossy behaviour; pass False only
                  to deliberately obtain the unlabelled variant of an env/core (used by broad
                  make_cycle, which must match both provenances). Has no effect when the input
                  carries no dummy isotopes, and is ignored entirely when keep_stereo is True
                  (isotopes are preserved in that case regardless).
    OUTPUT:
        SMILES of a context environment of a specified radius,
        list of SMILES of a core fragment with possible permutations of attachment point numbers
        env_smi, (core_smi_1, core_smi_2, ...)

        env_smi will not contain any Hs

        for radius 0 attachment point numbers will be stripped, but the string will correspond to core SMILES with
        radius > 0 if remove all map numbers from SMILES

    Output SMILES are standardized
    """

    if isinstance(context, str):
        context = Chem.MolFromSmiles(context)
    if isinstance(core, str):
        core = Chem.MolFromSmiles(core)

    # remove Hs from context and core
    if context:  # context cannot be H (no check needed), if so the user will obtain meaningless output
        context = Chem.RemoveHs(context)
    if core and Chem.MolToSmiles(core) != '[H][*:1]':
        core = Chem.RemoveHs(core)

    # Radius 0 keeps no context, so the environment reduces to the attachment-point
    # classes, `R<x>A<y>` for x ring cuts and y acyclic points. Rather than a separate
    # standardisation path, substitute a fake context in which every point of a class sits
    # in an identical component and let the ordinary pipeline below do the work: it then
    # numbers the points canonically and returns exactly the within-class relabellings as
    # its permutation tuple. Only the env string is swapped for the class on the way out.
    radius0_env_str = None
    if radius == 0 and core:
        context, radius0_env_str = __radius0_fake_context(core)
        radius = 1

    if core and context:

        # Count attachment points directly — counting context components
        # undercounts whenever two attachment points sit inside the same
        # connected component (close-anchor / ring-closure case).
        att_num = sum(a.GetAtomicNum() == 0 and a.GetAtomMapNum() > 0 for a in core.GetAtoms())

        if not keep_stereo:
            Chem.RemoveStereochemistry(context)
            Chem.RemoveStereochemistry(core)

        env = __get_context_env(context, radius)   # cut context to radius
        old_to_std = __standardize_att_by_env(
            env,
            core,
            keep_stereo,
            preserve_dummy_isotopes=preserve_dummy_isotopes,
        )
        env_smi = __mol_to_smiles(
            env,
            keep_stereo,
            preserve_dummy_isotopes=preserve_dummy_isotopes,
            allBondsExplicit=True,
        )
        if radius0_env_str is not None:
            # the fake context has served its purpose (numbering + orbits); report the
            # attachment-point class instead of its SMILES
            env_smi = radius0_env_str

        if att_num == 1:
            core_smi = __standardize_smiles_with_att_points(
                core,
                keep_stereo,
                preserve_dummy_isotopes=preserve_dummy_isotopes,
            )
            if return_att_map:
                return env_smi, (core_smi, ), {core_smi: old_to_std}
            return env_smi, (core_smi, )

        else:
            p = __get_att_permutations(
                env,
                keep_stereo,
                preserve_dummy_isotopes=preserve_dummy_isotopes,
            )
            smi_to_map = {}

            # permute attachment point numbering only in core,
            # since permutations in env will give the same canonical smiles
            if len(p) > 1:
                for d in p:
                    c = __permute_att(core, d)
                    smi = __standardize_smiles_with_att_points(
                        c,
                        keep_stereo,
                        preserve_dummy_isotopes=preserve_dummy_isotopes,
                    )
                    if smi not in smi_to_map:
                        """
                        Because d and the required mapping are in different coordinate systems, 
                        this assignment will not work - smi_to_map[smi] = d
                        - old_to_std is mapping from original input map numbers to standardized map numbers.
                        - d is mapping from standardized map numbers to permuted standardized map numbers.
                        The function needs to return mapping from original input maps to the final map numbers used in that smi.
                        So it must compose them:
                        old -> std -> permuted_std
                        """
                        smi_to_map[smi] = {old: d.get(std, std) for old, std in old_to_std.items()}
            else:
                smi = __standardize_smiles_with_att_points(
                    core,
                    keep_stereo,
                    preserve_dummy_isotopes=preserve_dummy_isotopes,
                )
                smi_to_map[smi] = old_to_std

            d = tuple(smi_to_map.keys())

            if return_att_map:
                return env_smi, d, smi_to_map
            return env_smi, d

    if return_att_map:
        return None
    return None, None


def get_radius0_rows(core, keep_stereo=False):
    """Radius-0 rows for one fragment: `(env, (core_smi, ...))`.

    `env` is the `RxAy` attachment-point class. The cores are every within-class
    relabelling of the fragment, with labellings that differ only by a symmetry of the
    fragment itself collapsed to one entry - those would splice to the identical product,
    so storing them all would weight symmetric fragments up under uniform row sampling
    (`C([1*:1])[1*:2]` and `C([1*:2])[1*:1]` are one row, not two).

    Both radius-0 builders go through this, so the from-scratch and the frags-derived
    tables cannot disagree about which rows a fragment produces.
    """
    if isinstance(core, str):
        core = Chem.MolFromSmiles(core)
    if core is None:
        return None, ()

    env, cores = get_std_context_core_permutations(
        '', core, 0, keep_stereo, preserve_dummy_isotopes=True,
    )
    if not cores:
        return env, ()

    # Symmetries of the fragment itself are exactly the permutations its own attachment
    # points admit when the fragment is used as its own env, so the existing orbit
    # machinery identifies them - no separate automorphism code is needed.
    #
    # The permutations must be derived from the labelling they are applied to. A dict from
    # __get_att_permutations maps map number to map number, so it is only meaningful in the
    # coordinate system of the molecule it was computed from - the same trap documented for
    # `d` versus `old_to_std` above. Computing them once from the *input* core and applying
    # them to the canonically renumbered strings silently permutes the wrong atoms, which
    # makes the returned set depend on the labelling the fragment arrived with. Both
    # builders rely on that set being labelling-independent, so each orientation gets its
    # own automorphisms here; two orientations related by a symmetry then share an orbit
    # and collapse onto the same representative.
    representatives = set()
    for smi in cores:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            representatives.add(smi)
            continue
        automorphisms = __get_att_permutations(mol, keep_stereo, preserve_dummy_isotopes=True)
        equivalent = {
            __standardize_smiles_with_att_points(
                __permute_att(mol, d), keep_stereo, preserve_dummy_isotopes=True)
            for d in automorphisms
        }
        representatives.add(min(equivalent))
    return env, tuple(sorted(representatives))


def get_canon_context_core(context, core, radius, keep_stereo=False, return_att_map=False,
                           preserve_dummy_isotopes=True):
    # context and core are Mols or SMILES
    # returns SMILES by default
    # preserve_dummy_isotopes defaults to True - see get_std_context_core_permutations
    res = get_std_context_core_permutations(
        context,
        core,
        radius,
        keep_stereo,
        return_att_map=return_att_map,
        preserve_dummy_isotopes=preserve_dummy_isotopes,
    )
    if res is not None:
        if return_att_map:
            env, cores, smi_to_map = res
            core_smi = sorted(cores)[0]
            return env, core_smi, smi_to_map[core_smi]
        env, cores = res
        return env, sorted(cores)[0]
    else:
        if return_att_map:
            return None, None, {}
        return None, None


def combine_core_env_to_rxn_smarts(core, env, keep_h=True):

    if isinstance(env, str):
        m_env = Chem.MolFromSmiles(env, sanitize=False)
    if isinstance(core, str):
        m_frag = Chem.MolFromSmiles(core, sanitize=False)

    backup_atom_map = "backupAtomMap"

    # put all atom maps to atom property and remove them
    for a in m_env.GetAtoms():
        atom_map = a.GetAtomMapNum()
        if atom_map:
            a.SetIntProp(backup_atom_map, atom_map)
            a.SetAtomMapNum(0)

    for a in m_frag.GetAtoms():
        atom_map = a.GetAtomMapNum()
        if atom_map:
            a.SetIntProp(backup_atom_map, atom_map)
            a.SetAtomMapNum(0)

    # set canonical ranks for atoms in env without maps
    m_env.UpdatePropertyCache()
    for atom_id, rank in zip([a.GetIdx() for a in m_env.GetAtoms()], list(Chem.CanonicalRankAtoms(m_env))):
        a = m_env.GetAtomWithIdx(atom_id)
        if not a.HasProp(backup_atom_map):
            a.SetAtomMapNum(rank + 1)  # because ranks start from 0

    m = Chem.RWMol(Chem.CombineMols(m_frag, m_env))

    links = defaultdict(list)  # pairs of atom ids to create bonds
    att_to_remove = []  # ids of att points to remove
    for a in m.GetAtoms():
        if a.HasProp(backup_atom_map):
            i = a.GetIntProp(backup_atom_map)
            links[i].append(a.GetNeighbors()[0].GetIdx())
            att_to_remove.append(a.GetIdx())

    for i, j in links.values():
        m.AddBond(i, j, Chem.BondType.SINGLE)

    for i in sorted(att_to_remove, reverse=True):
        m.RemoveAtom(i)

    comb_sma = mol_to_smarts(m, keep_h)
    if not keep_h:  # remove H only in mapped env part
        comb_sma = patt_remove_h.sub('', comb_sma)
    return comb_sma
