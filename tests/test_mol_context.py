from crem.mol_context import get_std_context_core_permutations


def test_dummy_isotopes_are_preserved_only_when_requested():
    context = "C([1*:1])[1*:2].C[*:3]"
    core = "CC(CCC([1*:1])[*:3])C[1*:2]"
    unlabeled_context = "C([*:1])[*:2].C[*:3]"
    unlabeled_core = "CC(CCC([*:1])[*:3])C[*:2]"

    env, cores = get_std_context_core_permutations(context, core, 2, False)
    expected_env, expected_cores = get_std_context_core_permutations(
        unlabeled_context,
        unlabeled_core,
        2,
        False,
    )
    assert "[1*" not in env
    assert all("[1*" not in c for c in cores)
    assert (env, cores) == (expected_env, expected_cores)

    env, cores = get_std_context_core_permutations(
        context,
        core,
        2,
        False,
        preserve_dummy_isotopes=True,
    )
    assert "[1*" in env
    assert any("[1*" in c for c in cores)
