import pytest
from rdkit import Chem

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


def test_mutate_return_mol(db, mol_aniline):
    res = list(mutate_mol(mol_aniline, db, radius=3, min_freq=0,
                          max_size=8, return_mol=True))
    assert res
    assert isinstance(res[0], list)
    assert isinstance(res[0][1], Chem.Mol)


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


def test_link_max_replacements_cap(db, mol_link1, mol_link2):
    res = list(link_mols(mol_link1, mol_link2, db, radius=1, min_freq=0,
                         min_atoms=1, max_atoms=5, max_replacements=2))
    assert len(res) <= 2


# ---------------------------------------------------------------------------
# make_cycle
# ---------------------------------------------------------------------------

def test_macrocycle_returns_results(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8))
    assert res


def test_macrocycle_valid_smiles(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8))
    assert all(_valid(s) for s in res)


def test_macrocycle_no_duplicates(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8))
    assert len(res) == len(set(res))


def test_macrocycle_max_replacements_cap(db, mol_macrocycle):
    res = list(make_cycle(mol_macrocycle, db, radius=1, min_freq=0,
                          min_atoms=1, max_atoms=8, max_replacements=1))
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
    # db_acyclic was built with --frag-mode acyclic, so it has no ring rows.
    # ring_closures=True should still WORK against it (just yield zero
    # results) because the column exists and just has no rows. To exercise
    # the legacy-DB error path, we simulate a pre-feature DB by dropping
    # the column from a fresh copy.
    import shutil, sqlite3, tempfile
    legacy_db = tempfile.mkstemp(suffix=".db")[1]
    shutil.copyfile(db_acyclic, legacy_db)
    with sqlite3.connect(legacy_db) as c:
        # Drop is_ring_closure column from radius1 (test only that radius).
        cols = [r[1] for r in c.execute("PRAGMA table_info(radius1)")]
        keep = [c2 for c2 in cols if c2 != "is_ring_closure"]
        c.execute(f"CREATE TABLE radius1_new ({','.join(keep)})")
        c.execute(f"INSERT INTO radius1_new SELECT {','.join(keep)} FROM radius1")
        c.execute("DROP TABLE radius1")
        c.execute("ALTER TABLE radius1_new RENAME TO radius1")

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


def test_get_mols_from_replacements_consistent_with_mutate(db, mol_aniline):
    replacements = list(get_replacements(mol_aniline, db, radius=3, min_freq=0,
                                         max_size=8, return_frag_smi_only=False))
    via_repl = set(get_mols_from_replacements(mol_aniline, 3, replacements))
    direct = set(mutate_mol(mol_aniline, db, radius=3, min_freq=0, max_size=8))
    assert via_repl == direct


def test_get_replacements_max_replacements_cap(db, mol_aniline):
    res = list(get_replacements(mol_aniline, db, radius=3, min_freq=0,
                                max_size=8, max_replacements=2))
    assert len(res) <= 2
