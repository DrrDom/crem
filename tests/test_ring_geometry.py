"""Tests for the discard_ring_geometry filter (geometrically impossible ring closures).

Two layers:
  * unit/integration through the private __frag_replace assembly function, which exercises the
    real molzip + sanitize + geometry-check path on hand-built context/bridge pairs;
  * public-API invariants on make_cycle / mutate_mol (filtering only removes structures; the
    flag forwards through the picklable wrappers; non ring-forming transforms are unaffected).
"""
import pytest
from rdkit import Chem

from crem import crem as C
from crem.crem import make_cycle, make_cycle2, mutate_mol

# module-level "__"-prefixed functions are plain attributes (no class name-mangling)
_frag_replace = C.__dict__["__frag_replace"]


def _assemble(ctx_smi, bridge_smi, discard):
    """Run __frag_replace on a benzene/azole context with two mapped dummies and a mapped
    bridge fragment; return the list of product SMILES."""
    ctx = Chem.MolFromSmiles(ctx_smi)
    parent = Chem.MolFromSmiles(ctx_smi.replace("[*:1]", "[H]").replace("[*:2]", "[H]"))
    return [t[0] for t in _frag_replace(parent, None, "[*:1][H]", bridge_smi, 1,
                                        context_mol=ctx, discard_ring_geometry=discard)]


# (context, bridge, description) — newly formed ring is geometrically impossible
IMPOSSIBLE = [
    ("[*:1]c1cccc([*:2])c1", "[*:1]CCC[*:2]", "benzene meta -> ring 6"),
    ("[*:1]c1cccc([*:2])c1", "[*:1]CC[*:2]", "benzene meta -> ring 5"),
    ("[*:1]c1ccc([*:2])s1", "[*:1]CCC[*:2]", "thiophene across S -> ring 6"),
    # the arc the new ring contains decides, not the atoms the bridge attaches to: these close
    # from a side-chain atom onto a ring atom, so neither anchor pair is meta/para itself
    ("CC(=O)N([*:1])c1cccc([*:2])c1", "[*:1]C[*:2]", "anilide N + meta ring C -> ring 5"),
    ("CC(=O)N([*:1])c1cccc([*:2])c1", "[*:1]CC[*:2]", "anilide N + meta ring C -> ring 6"),
    ("CC(=O)N([*:1])c1ccc([*:2])cc1", "[*:1]C[*:2]", "anilide N + para ring C -> ring 6"),
    ("OC([*:1])c1cccc([*:2])c1", "[*:1]C[*:2]", "benzylic C + meta ring C -> ring 5"),
]

# (context, bridge, description) — embeddable; must never be discarded (no false negatives)
POSSIBLE = [
    ("[*:1]c1cccc([*:2])c1", "[*:1]CCCCCC[*:2]", "benzene meta -> ring 9"),
    ("[*:1]c1ccccc1[*:2]", "[*:1]CCC[*:2]", "benzene ortho -> indane ring 5"),
    ("[*:1]c1ccccc1[*:2]", "[*:1]CC[*:2]", "benzene ortho -> benzocyclobutene ring 4"),
    ("[*:1]c1ccc([*:2])o1", "[*:1]CCC[*:2]", "furan across O -> ring 6"),
    ("[*:1]c1ccc([*:2])s1", "[*:1]CC[*:2]", "thiophene across S -> ring 5"),
    ("CC(=O)N([*:1])c1ccccc1[*:2]", "[*:1]C[*:2]", "anilide N + ortho ring C -> ring 4"),
    ("CC(=O)N([*:1])c1ccccc1[*:2]", "[*:1]CC[*:2]", "anilide N + ortho ring C -> ring 5"),
    ("OC([*:1])c1ccccc1[*:2]", "[*:1]CC[*:2]", "benzylic C + ortho ring C -> ring 5"),
    ("[*:1]c1cccc2cccc([*:2])c12", "[*:1]CC[*:2]", "naphthalene peri -> acenaphthene"),
]


@pytest.mark.parametrize("ctx,bridge,desc", IMPOSSIBLE)
def test_impossible_ring_discarded_when_on(ctx, bridge, desc):
    assert _assemble(ctx, bridge, discard=True) == [], desc


@pytest.mark.parametrize("ctx,bridge,desc", IMPOSSIBLE)
def test_impossible_ring_kept_when_off(ctx, bridge, desc):
    # with filtering disabled the structure is returned unchanged
    assert len(_assemble(ctx, bridge, discard=False)) == 1, desc


@pytest.mark.parametrize("ctx,bridge,desc", POSSIBLE)
def test_embeddable_ring_never_discarded(ctx, bridge, desc):
    # the critical no-false-negative property: on and off give the same (non-empty) result
    on = _assemble(ctx, bridge, discard=True)
    off = _assemble(ctx, bridge, discard=False)
    assert on == off and len(on) == 1, desc


def test_make_cycle_false_is_superset(db_rc):
    mol = Chem.MolFromSmiles("Cc1ccccc1CN")
    on = set(make_cycle(mol, db_rc, radius=1, ring_closures=False, discard_ring_geometry=True))
    off = set(make_cycle(mol, db_rc, radius=1, ring_closures=False, discard_ring_geometry=False))
    assert on <= off  # filtering can only ever remove structures


def test_make_cycle_single_and_multicore_agree(db_rc):
    mol = Chem.MolFromSmiles("CC(N)CCCO")
    for discard in (True, False):
        s1 = set(make_cycle(mol, db_rc, radius=1, ring_closures=False,
                            discard_ring_geometry=discard, ncores=1))
        s2 = set(make_cycle(mol, db_rc, radius=1, ring_closures=False,
                            discard_ring_geometry=discard, ncores=2))
        assert s1 == s2


def test_flag_forwards_through_make_cycle2(db_rc):
    mol = Chem.MolFromSmiles("CC(N)CCCO")
    on = set(make_cycle2(mol, db_rc, radius=1, ring_closures=False, discard_ring_geometry=True))
    off = set(make_cycle2(mol, db_rc, radius=1, ring_closures=False, discard_ring_geometry=False))
    assert on <= off


def test_non_ring_transform_unaffected(db):
    # mutate without ring formation: the flag must not change the output
    mol = Chem.MolFromSmiles("c1ccccc1N")
    on = set(mutate_mol(mol, db, radius=1, discard_ring_geometry=True))
    off = set(mutate_mol(mol, db, radius=1, discard_ring_geometry=False))
    assert on == off
