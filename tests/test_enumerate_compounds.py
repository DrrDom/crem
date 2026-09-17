"""Iterative enumeration with crem.utils.enumerate_compounds.

Every round works on the molecules produced by the previous one, so it has to know two
things about a product: which atoms were just added, and which atoms the caller wanted
left alone. Both are read off the molecule itself.

  * added atoms carry `CREM_MARKER_PROP`, which crem sets on the fragment it inserts and
    clears for the markers of any earlier call;
  * protected atoms carry `PROTECTED_ATOM_PROP`, set by enumerate_compounds on the starting
    molecule. Atom properties survive fragmentation and product assembly, so a child arrives
    with the protection of every atom it kept - no parent-to-child id mapping is involved.

RDKit's own <react_atom_idx> is not available here: crem clears the index bookkeeping of any
molecule it hands back.
"""
from rdkit import Chem

from crem import utils
from crem.crem import CREM_MARKER_PROP, grow_mol2
from crem.utils import PROTECTED_ATOM_PROP, enumerate_compounds

added_atom_ids = getattr(utils, "__get_child_added_atom_ids")
protect_atoms = getattr(utils, "__protect_atoms")
protected_atom_ids = getattr(utils, "__get_protected_atom_ids")

ACETANILIDE = "CC(=O)Nc1ccccc1"
ACETYL = Chem.MolFromSmarts("[CH3]C(=O)N")


# ---------------------------------------------------------------------------
# the two marks the iteration depends on
# ---------------------------------------------------------------------------

def test_added_atom_ids_are_the_inserted_ones(db, mol_aniline):
    _, product = grow_mol2(mol_aniline, db_name=db, radius=1, max_replacements=1,
                           return_mol=True)[0]
    added = added_atom_ids(product)
    assert added == sorted(a.GetIdx() for a in product.GetAtoms() if a.HasProp(CREM_MARKER_PROP))
    assert 0 < len(added) < product.GetNumAtoms()


def test_protection_is_inherited_by_products(db, mol_aniline):
    """The mark travels with the atom, which is what removes the need for an id mapping."""
    protected = protect_atoms(mol_aniline, [0, 1])
    assert protected_atom_ids(protected) == [0, 1]
    assert not any(a.HasProp(PROTECTED_ATOM_PROP) for a in mol_aniline.GetAtoms()), \
        "the input molecule must not be modified"

    products = grow_mol2(protected, db_name=db, radius=1, max_replacements=3,
                         protected_ids=[0, 1], return_mol=True)
    assert products
    for _, product in products:
        inherited = protected_atom_ids(product)
        assert len(inherited) == 2
        # an inherited mark never lands on an atom this call inserted
        assert not set(inherited) & set(added_atom_ids(product))


def test_protect_atoms_accepts_an_empty_selection(db, mol_aniline):
    untouched = protect_atoms(mol_aniline, [])
    assert protected_atom_ids(untouched) == []


# ---------------------------------------------------------------------------
# iteration
# ---------------------------------------------------------------------------

def test_scaffold_decoration_iterates(db):
    mol = Chem.MolFromSmiles(ACETANILIDE)
    one = enumerate_compounds(mol, db, mode="scaffold", n_iterations=1, radius=1,
                              return_smi=True, ncpu=1)
    two = enumerate_compounds(mol, db, mode="scaffold", n_iterations=2, radius=1,
                              return_smi=True, ncpu=1)
    assert len(one) > 0
    assert len(two) > len(one), "a second round must decorate the remaining positions"
    assert set(one) <= set(two), "results of earlier rounds are kept"


def test_scaffold_decoration_keeps_the_scaffold(db):
    mol = Chem.MolFromSmiles(ACETANILIDE)
    products = enumerate_compounds(mol, db, mode="scaffold", n_iterations=2, radius=1,
                                   return_smi=True, ncpu=1)
    assert all(Chem.MolFromSmiles(smi).HasSubstructMatch(mol) for smi in products)


def test_protected_positions_survive_every_round(db):
    mol = Chem.MolFromSmiles(ACETANILIDE)
    protected = enumerate_compounds(mol, db, mode="scaffold", n_iterations=2, radius=1,
                                    protected_ids=[0], return_smi=True, ncpu=1)
    unprotected = enumerate_compounds(mol, db, mode="scaffold", n_iterations=2, radius=1,
                                      return_smi=True, ncpu=1)
    assert protected, "protecting one position must not empty the result"
    assert len(protected) < len(unprotected)
    # the protected methyl is untouched in every product of both rounds
    assert all(Chem.MolFromSmiles(smi).HasSubstructMatch(ACETYL) for smi in protected)
    assert any(not Chem.MolFromSmiles(smi).HasSubstructMatch(ACETYL) for smi in unprotected)


def test_replace_ids_restricts_every_round(db):
    """replace_ids is the complement of protected_ids and must be just as persistent."""
    mol = Chem.MolFromSmiles(ACETANILIDE)
    ring = list(mol.GetSubstructMatch(Chem.MolFromSmarts("c1ccccc1")))
    products = enumerate_compounds(mol, db, mode="scaffold", n_iterations=2, radius=1,
                                   replace_ids=ring, return_smi=True, ncpu=1)
    assert products
    assert all(Chem.MolFromSmiles(smi).HasSubstructMatch(ACETYL) for smi in products), \
        "only the ring was selected, so the acetyl must survive both rounds"


def test_analog_enumeration_iterates_with_frozen_fragments(db, mol_macrocycle):
    one = enumerate_compounds(mol_macrocycle, db, mode="analogs", n_iterations=1, radius=1,
                              max_replacements=5, protect_added_frag=True,
                              return_smi=True, ncpu=1)
    two = enumerate_compounds(mol_macrocycle, db, mode="analogs", n_iterations=2, radius=1,
                              max_replacements=5, protect_added_frag=True,
                              return_smi=True, ncpu=1)
    assert one
    assert len(two) > len(one)


def test_returned_molecules_carry_no_protection_marks(db):
    mol = Chem.MolFromSmiles(ACETANILIDE)
    products = enumerate_compounds(mol, db, mode="scaffold", n_iterations=2, radius=1,
                                   max_replacements=3, protected_ids=[0], ncpu=1)
    assert products
    assert all(not a.HasProp(PROTECTED_ATOM_PROP)
               for product in products for a in product.GetAtoms())


def test_parallel_iteration_agrees_with_serial(db):
    """The marks are atom properties, so they have to survive being pickled to workers."""
    mol = Chem.MolFromSmiles(ACETANILIDE)
    options = dict(mode="scaffold", n_iterations=2, radius=1, protected_ids=[0],
                   return_smi=True)
    assert (sorted(enumerate_compounds(mol, db, ncpu=1, **options))
            == sorted(enumerate_compounds(mol, db, ncpu=2, **options)))
