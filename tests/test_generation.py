import pytest
from rdkit import Chem

import crem.crem as crem_module
from crem.crem import (
    get_mols_from_replacements,
    get_replacements,
    grow_mol,
    link_mols,
    make_cycle,
    mutate_mol,
)


def _valid(smi):
    return Chem.MolFromSmiles(smi) is not None


def _mark_parent_atoms(mol, prefix="parent"):
    for atom in mol.GetAtoms():
        atom.SetProp("parent_marker", f"{prefix}:{atom.GetIdx()}")


def _parent_marker_values(mol):
    return [
        atom.GetProp("parent_marker")
        for atom in mol.GetAtoms()
        if atom.HasProp("parent_marker") and not atom.HasProp("__crem")
    ]


def _assert_internal_index_removed(mol):
    assert all(not atom.HasProp("__crem_index") for atom in mol.GetAtoms())
    assert all(not atom.HasProp("__crem_isotope") for atom in mol.GetAtoms())


# ---------------------------------------------------------------------------
# mutate_mol
# ---------------------------------------------------------------------------

def test_mutate_returns_results(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert res


def test_mutate_valid_smiles(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert all(_valid(s) for s in res)


def test_mutate_no_duplicates(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert len(res) == len(set(res))


def test_mutate_source_not_in_output(db, mol_aniline):
    src = Chem.MolToSmiles(mol_aniline)
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert src not in res


@pytest.mark.parametrize("n", [1, 3])
def test_mutate_max_replacements_cap(db, mol_aniline, n):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                          max_size=8, max_replacements=n))
    assert len(res) <= n


def test_mutate_max_replacements_zero(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                          max_size=8, max_replacements=0))
    assert res == []


def test_mutate_min_freq_reduces_output(db, mol_aniline):
    r0 = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    r5 = list(mutate_mol(mol_aniline, db, radius=3, min_freq=5, max_size=8))
    assert len(r5) <= len(r0)


def test_mutate_high_min_freq_empty(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=10_000_000))
    assert res == []


def test_mutate_set_names_subset(db, mol_aniline):
    r_all = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    r_set = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8,
                            set_names=["test"]))
    assert len(r_set) <= len(r_all)


def test_mutate_set_names_invalid_raises(db, mol_aniline):
    with pytest.raises(ValueError, match="not found"):
        list(mutate_mol(mol_aniline, db, radius=3, set_names=["nonexistent"]))


def test_mutate_return_rxn(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                          max_size=8, return_rxn=True))
    assert res
    assert isinstance(res[0], list) and len(res[0]) == 2


def test_mutate_h_replacement_clears_invalid_alkene_stereo_and_removes_h():
    mol = Chem.MolFromSmiles(r"CC\C=C\[12CH3]")
    replace_ids = []
    for atom in mol.GetAtoms():
        if atom.GetIsotope():
            replace_ids.append(atom.GetIdx())
            atom.SetIsotope(0)
    protected_ids = set(range(mol.GetNumAtoms())).difference(replace_ids)

    fragments = crem_module.__fragment_mol(
        mol,
        radius=3,
        protected_ids=protected_ids,
        min_core_atoms=1,
        max_core_atoms=1,
    )
    assert len(fragments) == 1

    _, old_core, context_mol, num_heavy_atoms = fragments[0]
    assert old_core == "C[*:1]"
    assert num_heavy_atoms == 1

    products = list(
        crem_module.__frag_replace(
            mol,
            None,
            old_core,
            "[H][*:1]",
            3,
            context_mol,
        )
    )
    assert len(products) == 1
    product_smi, product_mol, transformation = products[0]
    assert product_smi == "C=CCC"
    assert transformation == "C[*:1]>>[H][*:1]"
    assert not product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#1]"))
    assert all(
        bond.GetStereo() == Chem.BondStereo.STEREONONE
        for bond in product_mol.GetBonds()
        if bond.GetBondType() == Chem.BondType.DOUBLE
    )


@pytest.mark.parametrize("ncores", [1, 2])
def test_mutate_return_mol(db, mol_aniline, ncores):
    _mark_parent_atoms(mol_aniline)
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                          max_size=8, return_mol=True, ncores=ncores))
    assert res
    assert isinstance(res[0], list)
    assert isinstance(res[0][1], Chem.Mol)
    assert any(
        any(atom.HasProp("__crem") for atom in item[1].GetAtoms())
        for item in res
    )
    assert any(_parent_marker_values(item[1]) for item in res)
    for item in res:
        _assert_internal_index_removed(item[1])


@pytest.mark.parametrize("ncores", [1, 2])
def test_mutate_return_mol_restores_parent_isotopes(db, mol_aniline, ncores):
    mol = Chem.Mol(mol_aniline)
    mol.GetAtomWithIdx(0).SetIsotope(13)
    _mark_parent_atoms(mol)
    res = list(mutate_mol(mol, db, radius=3, min_freq=0,
                          max_size=8, replace_ids=[6],
                          return_mol=True, ncores=ncores))
    assert res
    assert mol.GetAtomWithIdx(0).GetIsotope() == 13
    assert all(not atom.HasProp("__crem_isotope") for atom in mol.GetAtoms())
    assert any(
        any(atom.GetIsotope() == 13 and atom.HasProp("parent_marker")
            for atom in item[1].GetAtoms())
        for item in res
    )
    for item in res:
        _assert_internal_index_removed(item[1])


def test_mutate_return_rxn_and_freq(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8,
                          return_rxn=True, return_rxn_freq=True))
    assert res
    assert isinstance(res[0], list) and len(res[0]) == 3


def test_mutate_all_protected_empty(db, mol_aniline):
    all_ids = list(range(mol_aniline.GetNumAtoms()))
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                          protected_ids=all_ids))
    assert res == []


def test_mutate_replace_ids_restricts_output(db, mol_aniline):
    # Replace only atom 6 (the N); result set must be a subset of unrestricted
    full = set(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    restricted = set(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                                max_size=8, replace_ids=[6]))
    assert restricted.issubset(full)


def test_mutate_nonexistent_db():
    mol = Chem.MolFromSmiles("CCO")
    with pytest.raises(FileNotFoundError):
        list(mutate_mol(mol, "/nonexistent/path.db"))


def test_mutate_replace_cycles_no_queries_acyclic_rows_only(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")

    def assert_acyclic(row_ids, cur, radius):
        # On v2 provenance is the isotope-1 label on the core, not a column.
        if row_ids:
            bad = cur.execute(
                f"SELECT count(*) FROM radius{radius} r "
                f"JOIN frags f ON r.core_smi_id = f.core_smi_id "
                f"WHERE r.rowid IN ({','.join(map(str, row_ids))}) "
                f"AND f.core_smi LIKE '%[1*%'"
            ).fetchone()[0]
            assert bad == 0
        return row_ids

    list(mutate_mol(mol, db_rc, radius=1, min_freq=0, min_size=1,
                    max_size=8, replace_cycles="no",
                    filter_func=assert_acyclic))


def test_mutate_replace_cycles_partial_all_adds_partial_cycle_products(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    kw = dict(radius=1, min_freq=0, min_size=1, max_size=8,
              min_inc=-2, max_inc=4)
    base = set(mutate_mol(mol, db_rc, replace_cycles="no", **kw))
    ring = set(mutate_mol(mol, db_rc, replace_cycles="partial_all", **kw))
    assert base.issubset(ring)
    assert ring - base
    assert all(_valid(s) for s in ring)


def test_mutate_partial_cycles_return_mol_preserves_parent_atom_props(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    _mark_parent_atoms(mol)
    res = list(mutate_mol(mol, db_rc, radius=1, min_freq=0,
                          min_size=1, max_size=8, min_inc=-2, max_inc=4,
                          replace_cycles="partial_all", return_mol=True))
    assert res
    assert any(_parent_marker_values(item[1]) for item in res)
    for item in res:
        _assert_internal_index_removed(item[1])


def test_mutate_replace_cycles_partial_exo_is_subset_of_partial_all(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    kw = dict(radius=1, min_freq=0, min_size=1, max_size=8,
              min_inc=-2, max_inc=4)
    partial_all = set(mutate_mol(mol, db_rc, replace_cycles="partial_all", **kw))
    partial_exo = set(mutate_mol(mol, db_rc, replace_cycles="partial_exo", **kw))
    assert partial_exo
    assert partial_exo <= partial_all


def test_mutate_replace_cycles_boolean_aliases(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    kw = dict(radius=1, min_freq=0, min_size=1, max_size=8)
    assert set(mutate_mol(mol, db_rc, replace_cycles=False, **kw)) == \
        set(mutate_mol(mol, db_rc, replace_cycles="no", **kw))
    assert set(mutate_mol(mol, db_rc, replace_cycles=True, **kw)) == \
        set(mutate_mol(mol, db_rc, replace_cycles="forced", **kw))


def test_mutate_replace_cycles_forced_keeps_oversized_ring_core(db_rc):
    # Regression: forced mode must replace a cyclic source core even when it is
    # larger than max_size. This is the only path that hits the
    # include_cyclic_cores branch in __fragment_mol, which previously called
    # GetRingInfo() on an unsanitized rdMMPA.FragmentMol core and raised
    # "RingInfo not initialized".
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")  # cyclohexane ring (6 atoms) > max_size
    kw = dict(radius=1, min_freq=0, min_size=1, max_size=2, min_inc=-4, max_inc=4)
    no = set(mutate_mol(mol, db_rc, replace_cycles="no", **kw))
    forced = set(mutate_mol(mol, db_rc, replace_cycles="forced", **kw))
    assert no.issubset(forced)              # forced only ever adds products
    assert forced - no                      # oversized ring core actually got replaced
    assert all(_valid(s) for s in forced)


def test_mutate_rejects_unknown_replace_cycles_mode(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    with pytest.raises(ValueError, match="replace_cycles must be one of"):
        list(mutate_mol(mol, db_rc, radius=1, replace_cycles="partial"))


# ---------------------------------------------------------------------------
# grow_mol
# ---------------------------------------------------------------------------

def test_grow_returns_results(db, mol_aniline):
    res = list(grow_mol(mol_aniline, db, radius=3, min_freq=0,
                        min_atoms=1, max_atoms=3))
    assert res


def test_grow_valid_smiles(db, mol_aniline):
    res = list(grow_mol(mol_aniline, db, radius=3, min_freq=0,
                        min_atoms=1, max_atoms=3))
    assert all(_valid(s) for s in res)


def test_grow_outputs_larger_molecules(db, mol_aniline):
    hac = mol_aniline.GetNumHeavyAtoms()
    res = list(grow_mol(mol_aniline, db, radius=3, min_freq=0,
                        min_atoms=1, max_atoms=3))
    for smi in res:
        assert Chem.MolFromSmiles(smi).GetNumHeavyAtoms() > hac


def test_grow_max_replacements_cap(db, mol_aniline):
    res = list(grow_mol(mol_aniline, db, radius=3, min_freq=0,
                        min_atoms=1, max_atoms=3, max_replacements=1))
    assert len(res) <= 1


def test_grow_no_duplicates(db, mol_aniline):
    res = list(grow_mol(mol_aniline, db, radius=3, min_freq=0,
                        min_atoms=1, max_atoms=3))
    assert len(res) == len(set(res))


# ---------------------------------------------------------------------------
# link_mols
# ---------------------------------------------------------------------------

def test_link_returns_results(db, mol_link1, mol_link2):
    res = list(link_mols(mol_link1, mol_link2, db, radius=1, min_freq=0,
                         min_atoms=1, max_atoms=5))
    assert res


def test_link_valid_smiles(db, mol_link1, mol_link2):
    res = list(link_mols(mol_link1, mol_link2, db, radius=1, min_freq=0,
                         min_atoms=1, max_atoms=5))
    assert all(_valid(s) for s in res)


def test_link_no_duplicates(db, mol_link1, mol_link2):
    res = list(link_mols(mol_link1, mol_link2, db, radius=1, min_freq=0,
                         min_atoms=1, max_atoms=5))
    assert len(res) == len(set(res))


def test_link_return_mol_preserves_parent_atom_props(db, mol_link1, mol_link2):
    _mark_parent_atoms(mol_link1, "left")
    _mark_parent_atoms(mol_link2, "right")
    res = list(link_mols(mol_link1, mol_link2, db, radius=1, min_freq=0,
                         min_atoms=1, max_atoms=5, return_mol=True))
    assert res
    assert any(
        any(value.startswith("left:") for value in _parent_marker_values(item[1])) and
        any(value.startswith("right:") for value in _parent_marker_values(item[1]))
        for item in res
    )
    for item in res:
        _assert_internal_index_removed(item[1])


def test_link_max_replacements_cap(db, mol_link1, mol_link2):
    res = list(link_mols(mol_link1, mol_link2, db, radius=1, min_freq=0,
                         min_atoms=1, max_atoms=5, max_replacements=2))
    assert len(res) <= 2


# ---------------------------------------------------------------------------
# make_cycle
# ---------------------------------------------------------------------------

def test_macrocycle_returns_results(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8,
                          ring_closures=False))
    assert res


def test_macrocycle_valid_smiles(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8,
                          ring_closures=False))
    assert all(_valid(s) for s in res)


def test_macrocycle_no_duplicates(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8,
                          ring_closures=False))
    assert len(res) == len(set(res))


def test_macrocycle_return_mol_preserves_parent_atom_props(db, mol_macrocycle):
    _mark_parent_atoms(mol_macrocycle)
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8, return_mol=True,
                          ring_closures=False))
    assert res
    assert any(_parent_marker_values(item[1]) for item in res)
    for item in res:
        _assert_internal_index_removed(item[1])


def test_macrocycle_max_replacements_cap(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8,
                          max_replacements=1, ring_closures=False))
    assert len(res) <= 1


# ---------------------------------------------------------------------------
# make_cycle ring_closures=True (ring-bond arc fragments)
# ---------------------------------------------------------------------------

def test_ring_closure_returns_cycloalkane_from_butane(db_rc):
    # Butane(1,4) anchors are terminal CH3s with d_in=3. ring_size=5
    # combined with cyclopentane-arc fragments (1-carbon arc, dist2=2)
    # in db_rc must produce cyclopentane.
    m = Chem.MolFromSmiles("CCCC")
    res = list(make_cycle(m, db_rc, radius=1, ring_size=5,
                          ring_closures=True, min_freq=0,
                          min_atoms=1, max_atoms=8))
    assert "C1CCCC1" in res


def test_ring_closure_ring_size_window(db_rc):
    # ring_size=(5,7) on butane should generate cyclopentane through cycloheptane.
    m = Chem.MolFromSmiles("CCCC")
    res = list(make_cycle(m, db_rc, radius=1, ring_size=(5, 7),
                          ring_closures=True, min_freq=0,
                          min_atoms=1, max_atoms=8))
    assert "C1CCCC1" in res
    assert "C1CCCCC1" in res


def test_ring_closure_excludes_ring_size_below_d_in(db_rc):
    # ring_size=2 is impossible: butane(1,4) has d_in=3, so any new ring
    # must be at least 3+1=4 atoms. ring_size=2 yields no products.
    m = Chem.MolFromSmiles("CCCC")
    res = list(make_cycle(m, db_rc, radius=1, ring_size=2,
                          ring_closures=True, min_freq=0,
                          min_atoms=1, max_atoms=8))
    assert res == []


def test_ring_closure_provenance_separation(db_acyclic, db_rc):
    # ring_closures=True must query is_ring_closure=1 rows only. Against a
    # DB built with --frag-mode acyclic (no ring rows), a ring-closure call
    # must yield zero products. The same call against db_rc (which has both
    # provenances) must yield products for the same input molecule.
    m = Chem.MolFromSmiles("CCCC")
    res_acyclic_only = list(make_cycle(m, db_acyclic, radius=1, ring_size=5,
                                        ring_closures=True, min_freq=0,
                                        min_atoms=1, max_atoms=8))
    res_with_rings = list(make_cycle(m, db_rc, radius=1, ring_size=5,
                                      ring_closures=True, min_freq=0,
                                      min_atoms=1, max_atoms=8))
    assert res_acyclic_only == []
    assert res_with_rings  # non-empty


def test_ring_closure_close_anchors_connected_env(db_rc):
    # Close-anchor case: anchors 2 bonds apart on a small substrate exercise
    # the connected-env code path in mol_context. Use propane: anchors at
    # C0/C2 → d_in=2 → ring_size=4 needs dist2=2.
    m = Chem.MolFromSmiles("CCC")
    res = list(make_cycle(m, db_rc, radius=1, ring_size=4,
                          ring_closures=True, min_freq=0,
                          min_atoms=1, max_atoms=8))
    # Don't assert a specific product — corpus may not have a 2-bond arc
    # whose env matches propane. Just assert the call doesn't raise and the
    # result is a list of valid SMILES (or empty).
    for s in res:
        assert _valid(s)


def test_ring_closure_legacy_db_raises(db_acyclic):
    # The legacy-DB error path belongs to the v1 reader: a v1 database predating
    # ring-closure support has no is_ring_closure column, so ring_closures=True cannot be
    # satisfied and must raise rather than silently return nothing. db_acyclic is built by
    # the current code and is therefore v2 (where the column is absent by design and
    # provenance lives in the env), so stamp the copy back to v1 to exercise that path.
    import shutil, sqlite3, tempfile
    legacy_db = tempfile.mkstemp(suffix=".db")[1]
    shutil.copyfile(db_acyclic, legacy_db)
    with sqlite3.connect(legacy_db) as c:
        c.execute("PRAGMA user_version = 1")

    m = Chem.MolFromSmiles("CCCC")
    with pytest.raises(ValueError, match="is_ring_closure"):
        list(make_cycle(m, legacy_db, radius=1, ring_size=5,
                        ring_closures=True, min_freq=0,
                        min_atoms=1, max_atoms=8))

    # Same legacy DB with ring_closures=False must still work.
    res = list(make_cycle(m, legacy_db, radius=1, ring_size=5,
                          ring_closures=False, min_freq=0,
                          min_atoms=1, max_atoms=8))
    # Either products or none — the call must complete without raising.
    assert isinstance(res, list)


# ---------------------------------------------------------------------------
# get_replacements + get_mols_from_replacements
# ---------------------------------------------------------------------------

def test_get_replacements_frag_smi_only(db, mol_aniline):
    res = list(get_replacements(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert res
    assert all(isinstance(s, str) for s in res)


def test_get_replacements_tuple_mode(db, mol_aniline):
    res = list(get_replacements(mol_aniline, db, radius=3, min_freq=0,
                                max_size=8, return_frag_smi_only=False))
    assert res
    assert all(len(t) == 4 for t in res)


def test_get_replacements_context_restores_parent_isotopes(db, mol_aniline):
    mol = Chem.Mol(mol_aniline)
    mol.GetAtomWithIdx(0).SetIsotope(13)
    res = list(get_replacements(mol, db, radius=3, min_freq=0,
                                max_size=8, replace_ids_1=[6],
                                return_frag_smi_only=False))
    assert res
    context_mols = [item[3] for item in res]
    assert any(
        any(atom.GetIsotope() == 13 for atom in context_mol.GetAtoms())
        for context_mol in context_mols
    )
    for context_mol in context_mols:
        _assert_internal_index_removed(context_mol)


def test_get_mols_from_replacements_consistent_with_mutate(db, mol_aniline):
    replacements = list(get_replacements(mol_aniline, db, radius=3, min_freq=0,
                                         max_size=8, return_frag_smi_only=False))
    via_repl = set(get_mols_from_replacements(mol_aniline, 3, replacements))
    direct = set(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert via_repl == direct


def test_get_mols_from_replacements_return_mol_removes_internal_index(db, mol_aniline):
    _mark_parent_atoms(mol_aniline)
    replacements = list(get_replacements(mol_aniline, db, radius=3, min_freq=0,
                                         max_size=8, return_frag_smi_only=False))
    res = list(get_mols_from_replacements(mol_aniline, 3, replacements, return_mol=True))
    assert res
    assert any(_parent_marker_values(item[1]) for item in res)
    for item in res:
        _assert_internal_index_removed(item[1])


def test_get_replacements_max_replacements_cap(db, mol_aniline):
    res = list(get_replacements(mol_aniline, db, radius=3, min_freq=0,
                                max_size=8, max_replacements=2))
    assert len(res) <= 2


def test_get_replacements_replace_cycles_partial_all_matches_mutate(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    kw = dict(radius=1, min_freq=0, min_size=1, max_size=8,
              min_inc=-2, max_inc=4, replace_cycles="partial_all")
    replacements = list(get_replacements(mol, db_rc, return_frag_smi_only=False, **kw))
    via_repl = set(get_mols_from_replacements(mol, 1, replacements))
    direct = set(mutate_mol(mol, db_rc, **kw))
    assert via_repl == direct
