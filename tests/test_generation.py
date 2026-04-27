import pytest
from rdkit import Chem

from crem.crem import (
    get_mols_from_replacements,
    get_replacements,
    grow_mol,
    link_mols,
    make_macrocycle,
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
# make_macrocycle
# ---------------------------------------------------------------------------

def test_macrocycle_returns_results(db, mol_macrocycle):
    res = list(make_macrocycle(mol_macrocycle, db, radius=1, min_freq=0,
                               min_atoms=1, max_atoms=8))
    assert res


def test_macrocycle_valid_smiles(db, mol_macrocycle):
    res = list(make_macrocycle(mol_macrocycle, db, radius=1, min_freq=0,
                               min_atoms=1, max_atoms=8))
    assert all(_valid(s) for s in res)


def test_macrocycle_no_duplicates(db, mol_macrocycle):
    res = list(make_macrocycle(mol_macrocycle, db, radius=1, min_freq=0,
                               min_atoms=1, max_atoms=8))
    assert len(res) == len(set(res))


def test_macrocycle_max_replacements_cap(db, mol_macrocycle):
    res = list(make_macrocycle(mol_macrocycle, db, radius=1, min_freq=0,
                               min_atoms=1, max_atoms=8, max_replacements=1))
    assert len(res) <= 1


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
