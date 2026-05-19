import os
import subprocess
import sys

import pytest
from rdkit import Chem

from crem.crem import grow_mol, link_mols, mutate_mol


def test_mutate_same_seed_five_runs(db):
    mol = Chem.MolFromSmiles("c1ccccc1N")
    kw = dict(radius=3, min_freq=0, max_size=8, max_replacements=4, seed=42)
    results = [list(mutate_mol(mol, db, **kw)) for _ in range(5)]
    assert all(r == results[0] for r in results[1:])


def test_mutate_ring_closures_same_seed_five_runs(db_rc):
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    kw = dict(radius=1, min_freq=0, min_size=1, max_size=8,
              min_inc=-2, max_inc=4, ring_closures=True,
              max_replacements=6, seed=42)
    results = [list(mutate_mol(mol, db_rc, **kw)) for _ in range(5)]
    assert all(r == results[0] for r in results[1:])


def test_grow_same_seed_five_runs(db):
    mol = Chem.MolFromSmiles("c1ccccc1N")
    kw = dict(radius=3, min_freq=0, min_atoms=1, max_atoms=3, max_replacements=2, seed=7)
    results = [list(grow_mol(mol, db, **kw)) for _ in range(5)]
    assert all(r == results[0] for r in results[1:])


def test_link_same_seed_five_runs(db):
    m1 = Chem.MolFromSmiles("CCO")
    m2 = Chem.MolFromSmiles("CN")
    kw = dict(radius=1, min_freq=0, min_atoms=1, max_atoms=5, max_replacements=3, seed=13)
    results = [list(link_mols(m1, m2, db, **kw)) for _ in range(5)]
    assert all(r == results[0] for r in results[1:])


def test_different_seeds_differ(db):
    mol = Chem.MolFromSmiles("c1ccccc1N")
    kw = dict(radius=3, min_freq=0, max_size=8, max_replacements=4)
    r42 = list(mutate_mol(mol, db, seed=42, **kw))
    r99 = list(mutate_mol(mol, db, seed=99, **kw))
    assert r42 != r99


def test_seed_none_vs_seed_set(db):
    mol = Chem.MolFromSmiles("c1ccccc1N")
    kw = dict(radius=3, min_freq=0, max_size=8, max_replacements=4)
    # Seeded run must be reproducible; unseeded has no such guarantee
    r_seeded = list(mutate_mol(mol, db, seed=42, **kw))
    r_seeded2 = list(mutate_mol(mol, db, seed=42, **kw))
    assert r_seeded == r_seeded2


def test_seed_independent_of_pythonhashseed(db):
    """Results with seed=42 must be identical regardless of PYTHONHASHSEED."""
    cmd = [
        sys.executable, "-c",
        f"import sys; sys.path.insert(0, '.')\n"
        f"from rdkit import Chem\n"
        f"from crem.crem import mutate_mol\n"
        f"mol = Chem.MolFromSmiles('c1ccccc1N')\n"
        f"res = list(mutate_mol(mol, '{db}', radius=3, min_freq=0, max_size=8,\n"
        f"                      max_replacements=4, seed=42))\n"
        f"print(sorted(res))\n",
    ]
    outputs = set()
    for hs in ["1", "2", "42", "999", "random"]:
        env = {**os.environ, "PYTHONHASHSEED": hs}
        out = subprocess.run(cmd, capture_output=True, text=True,
                             env=env, check=True)
        outputs.add(out.stdout.strip())
    assert len(outputs) == 1, (
        f"seed=42 output differs across PYTHONHASHSEED values: {outputs}"
    )


def test_unseeded_nondeterministic(db):
    """Without a seed, repeated calls should not always return the same subset.

    We draw max_replacements=4 from 6 available results, so the probability that
    all 20 independent runs pick the exact same 4 molecules is (1/C(6,4))^19 ≈ 0.
    """
    mol = Chem.MolFromSmiles("c1ccccc1N")
    kw = dict(radius=3, min_freq=0, max_size=8, max_replacements=4)
    seen = set()
    for _ in range(20):
        r = tuple(sorted(mutate_mol(mol, db, **kw)))
        seen.add(r)
    assert len(seen) > 1, "unseeded calls should not always return the same subset"
