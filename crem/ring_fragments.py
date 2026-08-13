from itertools import combinations

from rdkit import Chem


ATOM_INDEX_PROP = "__crem_index"
RING_CUT_DUMMY_ISOTOPE = 1


def _ensure_atom_indices(mol):
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() and not atom.HasProp(ATOM_INDEX_PROP):
            atom.SetIntProp(ATOM_INDEX_PROP, atom.GetIdx())
    return mol


def _bond_atom_ids(mol, bond_id):
    bond = mol.GetBondWithIdx(bond_id)
    return bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()


def _components_after_cuts(mol, cut_bond_ids, allowed_atoms=None):
    if allowed_atoms is None:
        allowed_atoms = set(range(mol.GetNumAtoms()))
    else:
        allowed_atoms = set(allowed_atoms)

    adjacency = {idx: [] for idx in allowed_atoms}
    for bond in mol.GetBonds():
        bond_id = bond.GetIdx()
        if bond_id in cut_bond_ids:
            continue
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in allowed_atoms and end in allowed_atoms:
            adjacency[begin].append(end)
            adjacency[end].append(begin)

    components = []
    seen = set()
    for atom_id in sorted(allowed_atoms):
        if atom_id not in seen:
            stack = [atom_id]
            seen.add(atom_id)
            component = set()
            while stack:
                current = stack.pop()
                component.add(current)
                for neighbor in adjacency[current]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        stack.append(neighbor)
            components.append(component)
    return components  # tuple of sets of atom ids


def _ring_atoms_for_bonds(mol, bond_ids):
    atom_ids = set()
    for bond_id in bond_ids:
        begin, end = _bond_atom_ids(mol, bond_id)
        atom_ids.add(begin)
        atom_ids.add(end)
    return atom_ids


def _acyclic_side_cut_candidates(mol, base_atoms, ring_arc_atoms, ring_cut_bond_ids, side_cut_mode="all"):
    candidates = []
    ring_cut_bond_ids = set(ring_cut_bond_ids)
    for bond in mol.GetBonds():
        bond_id = bond.GetIdx()
        if bond_id in ring_cut_bond_ids:
            continue
        if bond.IsInRing() or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin not in base_atoms or end not in base_atoms:
            continue
        if not bond.GetBeginAtom().GetAtomicNum() or not bond.GetEndAtom().GetAtomicNum():
            continue
        if side_cut_mode == "exo":
            begin_in_arc = begin in ring_arc_atoms
            end_in_arc = end in ring_arc_atoms
            if begin_in_arc == end_in_arc:
                continue

        components = _components_after_cuts(
            mol,
            ring_cut_bond_ids | {bond_id},
            allowed_atoms=base_atoms,
        )
        core_side = [component for component in components if ring_arc_atoms <= component]
        if len(core_side) != 1:
            continue
        detached = [component for component in components if component is not core_side[0]]
        if len(detached) != 1:
            continue
        candidates.append((bond_id, detached[0]))
    return candidates


def _valid_side_cut_combo(combo):
    detached_sets = [detached for _, detached in combo]
    for i, detached in enumerate(detached_sets):
        for other in detached_sets[i + 1:]:
            if detached & other:
                return False
    return True


def _atom_index_set(mol):
    atom_ids = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() and atom.HasProp(ATOM_INDEX_PROP):
            atom_ids.add(atom.GetIntProp(ATOM_INDEX_PROP))
    return atom_ids


def _convert_dummy_isotopes_to_maps(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            isotope = atom.GetIsotope()
            if isotope:
                atom.SetIsotope(0)
                atom.SetAtomMapNum(isotope)


def _label_ring_cut_dummies(mol, ring_cut_maps):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() in ring_cut_maps:
            atom.SetIsotope(RING_CUT_DUMMY_ISOTOPE)


def _combine_mols(mols):
    if not mols:
        return None
    combined = Chem.Mol(mols[0])
    for mol in mols[1:]:
        combined = Chem.CombineMols(combined, mol)
    return combined


def _materialize_fragment(mol, cut_bond_ids, core_atom_ids, label_all_ring_cuts,
                          ring_cut_maps=None):
    cut_bond_ids = list(cut_bond_ids)
    ring_cut_maps = set(ring_cut_maps or ())
    dummy_labels = [(i + 1, i + 1) for i in range(len(cut_bond_ids))]
    try:
        cut = Chem.FragmentOnBonds(mol, cut_bond_ids, dummyLabels=dummy_labels)
        pieces = list(Chem.GetMolFrags(cut, asMols=True))
    except Exception:
        return None

    core_mol = None
    context_mols = []
    for piece in pieces:
        _convert_dummy_isotopes_to_maps(piece)
        atom_ids = _atom_index_set(piece)
        if atom_ids == core_atom_ids:
            if core_mol is not None:
                return None
            core_mol = piece
        else:
            context_mols.append(piece)

    if core_mol is None or not context_mols:
        return None

    n_core_attachments = sum(
        atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0
        for atom in core_mol.GetAtoms()
    )
    if not 2 <= n_core_attachments <= 4:
        return None

    context_mol = _combine_mols(context_mols)
    n_context_attachments = sum(
        atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0
        for atom in context_mol.GetAtoms()
    )
    if n_context_attachments != n_core_attachments:
        return None

    # Mark the ring-cut ends with isotope 1 so that later canonicalization can distinguish
    # them from acyclic (exo) attachment points. Under the v2 fragment convention every
    # ring cut is labelled; under v1 only fragments with a third/fourth attachment point
    # were labelled, because with exactly two points both are ring cuts by construction
    # and there was nothing to disambiguate.
    if (label_all_ring_cuts or n_core_attachments > 2) and ring_cut_maps:
        _label_ring_cut_dummies(core_mol, ring_cut_maps)
        _label_ring_cut_dummies(context_mol, ring_cut_maps)

    return core_mol, context_mol, core_atom_ids


def iter_partial_ring_fragments(mol, *, label_all_ring_cuts, max_acyclic_cuts=2,
                                min_core_atoms=None, max_core_atoms=None,
                                side_cut_mode="all"):
    """Yield connected partial-ring fragments with 2-4 context attachments.

    Two non-aromatic single bonds from the same ring are always cut. Up to
    `max_acyclic_cuts` additional acyclic single heavy-atom bonds may be cut
    on distinct detached side components of the selected ring arc. With
    ``side_cut_mode="exo"``, those side cuts are limited to acyclic bonds
    adjacent to an atom of the selected ring arc.

    ``label_all_ring_cuts`` selects the fragment convention and is deliberately
    keyword-only with no default: it must follow the schema version of the database
    the fragments will be matched against, and silently inheriting the wrong value
    yields zero matches rather than an error. Pass True for v2 databases (every
    ring cut carries isotope 1) and False for v1 (only 3-/4-attachment fragments
    are labelled).
    """
    if side_cut_mode not in {"all", "exo"}:
        raise ValueError("side_cut_mode must be one of: 'all', 'exo'")

    mol = _ensure_atom_indices(mol)
    try:
        bond_rings = mol.GetRingInfo().BondRings()
    except (RuntimeError, IndexError):
        return

    if not bond_rings:
        return

    seen_ring_pairs = set()
    for ring_bond_ids in bond_rings:
        ring_bond_ids = tuple(sorted(ring_bond_ids))
        ring_atom_ids = _ring_atoms_for_bonds(mol, ring_bond_ids)
        for bond_1, bond_2 in combinations(ring_bond_ids, 2):
            ring_pair = (bond_1, bond_2)
            if ring_pair in seen_ring_pairs:
                continue
            seen_ring_pairs.add(ring_pair)

            b1 = mol.GetBondWithIdx(bond_1)
            b2 = mol.GetBondWithIdx(bond_2)
            if b1.GetBondType() != Chem.BondType.SINGLE or b1.GetIsAromatic():
                continue
            if b2.GetBondType() != Chem.BondType.SINGLE or b2.GetIsAromatic():
                continue

            ring_cut_bond_ids = {bond_1, bond_2}
            base_components = _components_after_cuts(mol, ring_cut_bond_ids)
            if len(base_components) != 2:
                continue

            for base_atoms in base_components:
                ring_arc_atoms = base_atoms & ring_atom_ids
                if not ring_arc_atoms:
                    continue
                side_candidates = _acyclic_side_cut_candidates(
                    mol,
                    base_atoms,
                    ring_arc_atoms,
                    ring_cut_bond_ids,
                    side_cut_mode=side_cut_mode,
                )
                max_cuts = min(max_acyclic_cuts, len(side_candidates))
                for n_cuts in range(max_cuts + 1):
                    for side_combo in combinations(side_candidates, n_cuts):
                        if not _valid_side_cut_combo(side_combo):
                            continue
                        detached = set()
                        side_bond_ids = []
                        for bond_id, detached_atoms in side_combo:
                            side_bond_ids.append(bond_id)
                            detached.update(detached_atoms)
                        core_atoms = set(base_atoms) - detached
                        core_atom_ids = {
                            mol.GetAtomWithIdx(atom_id).GetIntProp(ATOM_INDEX_PROP)
                            for atom_id in core_atoms
                        }
                        if min_core_atoms is not None and len(core_atom_ids) < min_core_atoms:
                            continue
                        if max_core_atoms is not None and len(core_atom_ids) > max_core_atoms:
                            continue
                        fragment = _materialize_fragment(
                            mol,
                            sorted(ring_cut_bond_ids) + sorted(side_bond_ids),
                            core_atom_ids,
                            label_all_ring_cuts,
                            ring_cut_maps=range(1, len(ring_cut_bond_ids) + 1),
                        )
                        if fragment is not None:
                            yield fragment
