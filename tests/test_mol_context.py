from crem.mol_context import get_std_context_core_permutations


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
