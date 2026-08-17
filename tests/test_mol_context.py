import pytest

from crem.mol_context import (RADIUS0_ENV_CLASSES, get_radius0_rows,
                              get_std_context_core_permutations)


CONTEXT = "C([1*:1])[1*:2].C[*:3]"
CORE = "CC(CCC([1*:1])[*:3])C[1*:2]"
UNLABELED_CONTEXT = "C([*:1])[*:2].C[*:3]"
UNLABELED_CORE = "CC(CCC([*:1])[*:3])C[*:2]"


def test_dummy_isotopes_are_stripped_when_explicitly_disabled():
    env, cores = get_std_context_core_permutations(
        CONTEXT, CORE, 2, False, preserve_dummy_isotopes=False,
    )
    expected_env, expected_cores = get_std_context_core_permutations(
        UNLABELED_CONTEXT, UNLABELED_CORE, 2, False, preserve_dummy_isotopes=False,
    )
    assert "[1*" not in env
    assert all("[1*" not in c for c in cores)
    assert (env, cores) == (expected_env, expected_cores)


def test_dummy_isotopes_are_preserved_by_default():
    # The default is preserve, because stripping merges ring-arc fragments into the
    # acyclic string space and a caller that forgets the flag would silently query the
    # wrong env class on a v2 database.
    for kwargs in ({}, {"preserve_dummy_isotopes": True}):
        env, cores = get_std_context_core_permutations(CONTEXT, CORE, 2, False, **kwargs)
        assert "[1*" in env
        assert any("[1*" in c for c in cores)


def test_stripping_is_a_no_op_for_unlabelled_input():
    # Verified property that lets the build path preserve unconditionally: for input
    # carrying no dummy isotopes the flag makes no difference at all.
    a = get_std_context_core_permutations(
        UNLABELED_CONTEXT, UNLABELED_CORE, 2, False, preserve_dummy_isotopes=False,
    )
    b = get_std_context_core_permutations(
        UNLABELED_CONTEXT, UNLABELED_CORE, 2, False, preserve_dummy_isotopes=True,
    )
    assert a == b


# ---------------------------------------------------------------------------
# radius 0
# ---------------------------------------------------------------------------

RADIUS0_CORES = [
    "C[*:1]",                             # R0A1
    "C([*:1])[*:2]",                      # R0A2, symmetric
    "O(C[*:1])[*:2]",                     # R0A2, asymmetric
    "C([1*:1])[1*:2]",                    # R2A0, symmetric arc
    "C(C([*:1])[1*:2])[1*:3]",            # R2A1
    "C([*:1])([1*:2])[1*:3]",             # R2A1, equivalent ring points
    "C(C(C[1*:2])C[1*:3])[*:1]",          # R2A1, equivalent ring branches
    "C(C([*:1])[1*:2])([*:4])[1*:3]",     # R2A2
    "CC([*:1])([*:2])[*:3]",              # R0A3, all points equivalent
    "C([*:1])([*:2])([*:3])[*:4]",        # R0A4, all points equivalent
    "C(C[*:2])C([*:3])([*:4])C[*:1]",     # R0A4, asymmetric
]


def _relabellings(core_smi):
    """Every class-preserving relabelling of a fragment's attachment points."""
    from itertools import permutations
    from rdkit import Chem
    m = Chem.MolFromSmiles(core_smi)
    ring = [a.GetIdx() for a in m.GetAtoms()
            if a.GetAtomicNum() == 0 and a.GetAtomMapNum() and a.GetIsotope() == 1]
    acyclic = [a.GetIdx() for a in m.GetAtoms()
               if a.GetAtomicNum() == 0 and a.GetAtomMapNum() and not a.GetIsotope()]
    slots = list(range(1, len(ring) + len(acyclic) + 1))
    out = []
    for ring_slots in permutations(slots, len(ring)):
        rest = [s for s in slots if s not in ring_slots]
        for acyclic_slots in permutations(rest):
            mm = Chem.Mol(m)
            for idx, slot in zip(ring, ring_slots):
                mm.GetAtomWithIdx(idx).SetAtomMapNum(slot)
            for idx, slot in zip(acyclic, acyclic_slots):
                mm.GetAtomWithIdx(idx).SetAtomMapNum(slot)
            out.append(Chem.MolToSmiles(mm))
    return out


@pytest.mark.parametrize("core", RADIUS0_CORES)
def test_radius0_rows_are_idempotent(core):
    """Feeding a returned row back in must reproduce the same row set.

    Both radius-0 builders rely on a fragment producing one fixed set of rows. This caught
    a coordinate-system bug where the automorphism permutations were computed on the input
    labelling and applied to the canonically renumbered strings, which made the result
    depend on the labelling the fragment arrived with.
    """
    env, rows = get_radius0_rows(core)
    assert env is not None and rows
    for row in rows:
        assert get_radius0_rows(row) == (env, rows)


@pytest.mark.parametrize("core", RADIUS0_CORES)
def test_radius0_rows_are_labelling_invariant(core):
    """Every relabelling of one fragment must give the identical (env, rows)."""
    expected = get_radius0_rows(core)
    for variant in _relabellings(core):
        assert get_radius0_rows(variant) == expected, variant


@pytest.mark.parametrize("core", RADIUS0_CORES)
def test_radius0_env_matches_attachment_point_classes(core):
    """The env must report the actual ring / acyclic attachment-point counts."""
    from rdkit import Chem
    m = Chem.MolFromSmiles(core)
    n_ring = sum(1 for a in m.GetAtoms()
                 if a.GetAtomicNum() == 0 and a.GetAtomMapNum() and a.GetIsotope() == 1)
    n_att = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 0 and a.GetAtomMapNum())
    env, _rows = get_radius0_rows(core)
    assert env == f"R{n_ring}A{n_att - n_ring}"
    assert env in RADIUS0_ENV_CLASSES


def test_radius0_symmetric_relabellings_collapse_to_one_row():
    """A fragment whose attachment points are symmetry-equivalent yields a single row.

    Storing both C([1*:1])[1*:2] and C([1*:2])[1*:1] would weight symmetric fragments up
    under uniform row sampling, since the two splice to the identical product.
    """
    assert get_radius0_rows("C([1*:1])[1*:2]")[1] == ("C([1*:1])[1*:2]",)
    assert len(get_radius0_rows("CC([*:1])([*:2])[*:3]")[1]) == 1
    assert len(get_radius0_rows("C([*:1])([*:2])([*:3])[*:4]")[1]) == 1
    # asymmetric points must NOT be collapsed
    assert len(get_radius0_rows("O(C[*:1])[*:2]")[1]) == 2
    assert len(get_radius0_rows("C(C([*:1])[1*:2])[1*:3]")[1]) == 2
