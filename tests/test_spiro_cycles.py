"""make_cycle(spiro=...) - closing a ring through a single atom.

Both ends of the linker attach to one ring atom, so the new ring meets the
existing one there and nowhere else: a spiro centre. `spiro` is False (off,
the default), True (added to the ordinary two-atom closures) or 'only' (those
closures alone). An acyclic atom is not eligible - a closure through one of
those is an ordinary core replacement that mutate_mol and grow_mol already
produce.

Nothing on the database side is special: `iter_partial_ring_fragments` cuts
adjacent ring-bond pairs like any other pair, so a spiro compound in the corpus
already contributes the env with two attachment points on one atom. See the
`db_spiro` fixture for why that corpus cannot be the ring_closures one.
"""
import pytest
from rdkit import Chem

from crem import crem as crem_mod
from crem.crem import make_cycle

dist2_sql = getattr(crem_mod, "__dist2_sql")

# every call in this file uses the same tiny-corpus settings
KW = dict(radius=1, ring_closures=True, min_freq=0, min_atoms=1, max_atoms=8)


def cycles(smi, db, **kwargs):
    """make_cycle over `smi` as a set of product SMILES."""
    return set(make_cycle(Chem.MolFromSmiles(smi), db, **{**KW, **kwargs}))


def spiro_atoms(smi):
    """Indices of atoms shared by two rings that meet at that atom alone."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return []
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    out = []
    for idx in range(mol.GetNumAtoms()):
        member = [r for r in rings if idx in r]
        if len(member) >= 2 and all(len(a & b) == 1
                                    for i, a in enumerate(member)
                                    for b in member[i + 1:]):
            out.append(idx)
    return out


def smallest_ring(smi):
    mol = Chem.MolFromSmiles(smi)
    rings = mol.GetRingInfo().AtomRings()
    return min((len(r) for r in rings), default=None)


# ---------------------------------------------------------------------------
# the three settings of the flag
# ---------------------------------------------------------------------------

def test_spiro_is_off_by_default(db_spiro):
    # Cyclopentane has no two-atom closure to a three-membered ring, so with the
    # feature off there is nothing at all; the spiro closure is opt-in.
    assert cycles("C1CCCC1", db_spiro, ring_size=3) == set()
    assert cycles("C1CCCC1", db_spiro, ring_size=3, spiro=False) == set()
    assert cycles("C1CCCC1", db_spiro, ring_size=3, spiro=True)


def test_spiro_true_adds_to_the_two_atom_closures(db_spiro):
    # Methylcyclohexane at ring_size=3 has both kinds: the methyl and its ring
    # neighbour close a cyclopropane between them, and three ring CH2 carry a
    # spiro one. True is exactly the union - nothing is substituted.
    off = cycles("CC1CCCCC1", db_spiro, ring_size=3)
    only = cycles("CC1CCCCC1", db_spiro, ring_size=3, spiro="only")
    both = cycles("CC1CCCCC1", db_spiro, ring_size=3, spiro=True)
    assert off and only
    assert not (off & only)
    assert both == off | only


def test_spiro_only_returns_the_same_atom_closures_alone(db_spiro):
    only = cycles("CC1CCCCC1", db_spiro, ring_size=3, spiro="only")
    assert Chem.CanonSmiles("CC1CCCCC12CC2") in only
    assert all(spiro_atoms(s) for s in only)
    # the two-atom product of the same call is gone
    assert Chem.CanonSmiles("C1CCC2(CC1)CC2") not in only


def test_an_unknown_spiro_value_is_rejected(db_spiro):
    with pytest.raises(ValueError, match="spiro must be"):
        list(make_cycle(Chem.MolFromSmiles("C1CCCC1"), db_spiro, spiro="yes"))


# ---------------------------------------------------------------------------
# which atoms are eligible
# ---------------------------------------------------------------------------

def test_same_atom_closure_builds_a_spiro_ring(db_spiro):
    # Cyclopentane, both linker ends on one CH2, ring_size=3: the two-carbon arc
    # in db_spiro must give spiro[2.4]heptane, which is not itself in the corpus.
    res = cycles("C1CCCC1", db_spiro, ring_size=3, spiro=True, replace_ids=[0])
    assert Chem.CanonSmiles("C1CC2(CC2)CC1") in res
    assert all(spiro_atoms(s) for s in res)


def test_an_acyclic_atom_is_not_eligible(db_spiro):
    # Butane's C2 has two hydrogens but no ring to be spiro to, so the closure
    # through it is an ordinary core replacement and make_cycle leaves it to
    # mutate_mol. The two-atom closures of the same call are unaffected.
    res = cycles("CCCC", db_spiro, ring_size=5, spiro=True)
    assert Chem.CanonSmiles("C1CCCC1") in res, "two-atom closure missing"
    assert Chem.CanonSmiles("CCC1(C)CCCC1") not in res
    with pytest.warns(RuntimeWarning):
        assert cycles("CCCC", db_spiro, ring_size=5, spiro="only") == set()


def test_a_ring_atom_with_one_hydrogen_is_not_eligible(db_spiro):
    # C1 of methylcyclohexane is in a ring but carries the methyl, so it has a
    # single hydrogen and cannot host both linker ends.
    assert cycles("CC1CCCCC1", db_spiro, ring_size=3, spiro=True,
                  replace_ids=[1]) == set()
    assert cycles("CC1CCCCC1", db_spiro, ring_size=3, spiro=True,
                  replace_ids=[2]) == {Chem.CanonSmiles("CC1CCCCC12CC2")}


def test_spiro_only_warns_when_nothing_is_eligible(db_spiro):
    # 'only' drops every other source of products, so an ineligible request is
    # empty by construction - the user hears about it rather than seeing silence.
    with pytest.warns(RuntimeWarning, match="spiro='only'"):
        res = cycles("CCCC", db_spiro, ring_size=3, spiro="only")
    assert res == set()


# ---------------------------------------------------------------------------
# a ring needs three atoms
# ---------------------------------------------------------------------------

def test_same_atom_closure_refuses_a_two_membered_ring(db_spiro):
    # ring_size=2 would ask for a one-atom linker on one anchor, which is a
    # second bond between two already-bonded atoms rather than a ring.
    assert cycles("C1CCCC1", db_spiro, ring_size=2, spiro=True,
                  replace_ids=[0]) == set()


def test_unconstrained_same_atom_closure_still_has_no_tiny_rings(db_spiro):
    # ring_size=None leaves the size open, so the three-atom floor is carried by
    # the one-sided dist2 window rather than by the caller.
    res = cycles("C1CCCC1", db_spiro, ring_size=None, spiro="only",
                 replace_ids=[0])
    assert res
    assert all(smallest_ring(s) >= 3 for s in res)


# ---------------------------------------------------------------------------
# symmetry
# ---------------------------------------------------------------------------

def test_symmetry_equivalent_anchors_give_the_same_products(db_spiro):
    # MMPA offers each fragmentation at one representative of a symmetry class, so
    # a same-atom cut read off MMPA would exist for only one of two equivalent
    # atoms and naming the other would return nothing. Building the pairs directly
    # removes that asymmetry: piperidine's C3 and C5 must agree, and neither may
    # be empty.
    m = Chem.MolFromSmiles("C1CCNCC1")
    ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    classes = {}
    for idx, rank in enumerate(ranks):
        classes.setdefault(rank, []).append(idx)
    pairs = [ids for ids in classes.values() if len(ids) == 2]
    assert pairs, "piperidine must expose at least one symmetric pair of atoms"
    for first, second in pairs:
        a = cycles("C1CCNCC1", db_spiro, ring_size=3, spiro=True, replace_ids=[first])
        b = cycles("C1CCNCC1", db_spiro, ring_size=3, spiro=True, replace_ids=[second])
        assert a == b
        assert a


# ---------------------------------------------------------------------------
# arc closures are untouched
# ---------------------------------------------------------------------------

def test_arc_closures_are_unaffected(db_rc):
    # Two-atom anchor pairs must behave exactly as before, on the corpus the
    # existing ring-closure tests use: butane(1,4) still closes to cyclopentane.
    res = cycles("CCCC", db_rc, ring_size=5)
    assert Chem.CanonSmiles("C1CCCC1") in res


# ---------------------------------------------------------------------------
# what 'only' skips
# ---------------------------------------------------------------------------

def test_spiro_only_does_not_run_the_two_cut_pass(db_spiro, monkeypatch):
    # Same-atom contexts are built from the molecule directly, so 'only' has no
    # use for the maxCuts=2 MMPA pass - not running it is most of what makes the
    # mode cheaper than filtering a full run.
    class Boom:
        def FragmentMol(self, *args, **kwargs):
            raise AssertionError("MMPA must not run for spiro='only'")

    monkeypatch.setattr(crem_mod, "rdMMPA", Boom())
    assert cycles("C1CCCC1", db_spiro, ring_size=3, spiro="only")
    with pytest.raises(AssertionError, match="MMPA must not run"):
        cycles("C1CCCC1", db_spiro, ring_size=3, spiro=True)


def test_spiro_only_skips_the_disconnected_env_fragmenter(db_spiro, monkeypatch):
    # The broad fragmenter pairs two independent single cuts, one per atom, so a
    # same-atom pair is not expressible in it and it has nothing to contribute.
    def boom(*args, **kwargs):
        raise AssertionError("the macrocycle fragmenter must not run")

    monkeypatch.setattr(crem_mod, "__fragment_mol_macrocycle", boom)
    assert cycles("C1CCCC1", db_spiro, ring_size=3, spiro="only", ring_closures=False)


def test_spiro_only_is_indifferent_to_ring_closures(db_spiro):
    # A same-atom context can only be matched by a fragment that was itself cut
    # out of a ring - the two cut bonds and the core close a cycle through the
    # atom - so the strict mode's is_ring_closure filter removes nothing that the
    # broad mode would have added.
    strict = cycles("CC1CCCCC1", db_spiro, ring_size=None, spiro="only", ring_closures=True)
    broad = cycles("CC1CCCCC1", db_spiro, ring_size=None, spiro="only", ring_closures=False)
    assert strict
    assert strict == broad


def test_spiro_survives_multiprocessing(db_spiro):
    # the flag has to reach the fragmenter through __get_data_cycle as well
    one = cycles("CC1CCCCC1", db_spiro, ring_size=None, spiro=True)
    two = cycles("CC1CCCCC1", db_spiro, ring_size=None, spiro=True, ncores=2)
    assert one
    assert one == two


# ---------------------------------------------------------------------------
# the dist2 predicate builder
# ---------------------------------------------------------------------------

def test_dist2_sql_windows():
    assert dist2_sql(None) == ""
    assert dist2_sql(3) == " AND dist2 = 3"
    assert dist2_sql((3, 6)) == " AND dist2 BETWEEN 3 AND 6"
    assert dist2_sql((3, None)) == " AND dist2 >= 3"
    assert dist2_sql((None, 6)) == " AND dist2 <= 6"
    assert dist2_sql((None, None)) == ""


def test_dist2_sql_honours_column_and_formatter():
    quoted = lambda v: f"<{v}>"
    assert dist2_sql((3, None), column="r.dist2", fmt=quoted) == " AND r.dist2 >= <3>"
    assert dist2_sql(4, column="r.dist2", fmt=quoted) == " AND r.dist2 = <4>"


def test_dist2_sql_rejects_a_bad_tuple():
    with pytest.raises(ValueError):
        dist2_sql((1, 2, 3))


# ---------------------------------------------------------------------------
# replace_ids is a restriction, not a hint
# ---------------------------------------------------------------------------

def test_replace_ids_restricts_to_the_named_atom(db_spiro):
    # A two-atom closure needs both of its anchors named, so naming one ring atom
    # leaves its spiro closure as the only possibility. On methylcyclohexane the
    # unrestricted call also closes the methyl onto C1; that product must go.
    kw = dict(ring_size=3, spiro=True)
    everything = cycles("CC1CCCCC1", db_spiro, **kw)
    only_c3 = cycles("CC1CCCCC1", db_spiro, replace_ids=[3], **kw)
    fused = Chem.CanonSmiles("C1CCC2(CC1)CC2")
    assert fused in everything
    assert fused not in only_c3
    assert only_c3
    assert only_c3 < everything


def test_replace_ids_restricts_on_a_fully_symmetric_molecule(db_spiro):
    # Every carbon of cyclohexane is one symmetry class, so widening the request
    # to that class would leave nothing protected and lift the restriction
    # altogether. The named atom has to be taken exactly: only closures through
    # it survive, and they are strictly fewer than the unrestricted result.
    kw = dict(ring_size=6, spiro=True)
    everything = cycles("C1CCCCC1", db_spiro, **kw)
    only_zero = cycles("C1CCCCC1", db_spiro, replace_ids=[0], **kw)
    assert only_zero
    assert only_zero < everything, "replace_ids did not restrict anything"
    assert all(spiro_atoms(s) for s in only_zero)


def test_naming_two_atoms_allows_more_than_naming_one(db_spiro):
    # [i, j] permits the closure through i, the one through j, and the arc between
    # them, so it must be a proper superset of [i].
    kw = dict(ring_size=None, spiro=True)
    one = cycles("C1CCCCC1", db_spiro, replace_ids=[0], **kw)
    two = cycles("C1CCCCC1", db_spiro, replace_ids=[0, 1], **kw)
    assert one
    assert one < two


def test_duplicate_ids_are_discarded(db_spiro):
    # [i, i] carries no extra meaning: the duplicate is dropped and the call is
    # the [i] call.
    kw = dict(ring_size=None, spiro=True)
    assert (cycles("C1CCCCC1", db_spiro, replace_ids=[0, 0], **kw) ==
            cycles("C1CCCCC1", db_spiro, replace_ids=[0], **kw) != set())
