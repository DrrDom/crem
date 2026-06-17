import subprocess
import sys

import pytest
from rdkit import Chem

# 20 molecules chosen to yield non-trivial fragment coverage at radii 1–3
# and to produce results for every public API function.
CORPUS = [
    "CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccccc1N",
    "c1ccccc1O", "c1ccccc1F", "CC(=O)O", "CCNCC", "c1ccc(N)cc1",
    "c1ccc(O)cc1", "c1ccncc1", "CC(C)O", "CCOCC", "CC(=O)Nc1ccccc1",
    "CCCC(=O)O", "c1ccc2ccccc2c1", "CCC(N)=O", "NCCCCN", "CC(N)CCC(=O)O",
]


@pytest.fixture(scope="session")
def db(tmp_path_factory):
    """Build a tiny v1 DB from CORPUS once per session; shared by all tests."""
    d = tmp_path_factory.mktemp("db")
    smi_file = d / "input.smi"
    db_path = str(d / "test.db")
    smi_file.write_text(
        "\n".join(f"{s} mol{i}" for i, s in enumerate(CORPUS))
    )
    subprocess.run(
        [sys.executable, "-m", "crem.scripts.cremdb_create",
         "-i", str(smi_file), "-o", db_path,
         "-s", "test", "--radii", "1", "2", "3", "--ncpu", "1",
         ],
        check=True, capture_output=True,
    )
    return db_path


@pytest.fixture
def mol_aniline():
    return Chem.MolFromSmiles("c1ccccc1N")


@pytest.fixture
def mol_link1():
    return Chem.MolFromSmiles("CCO")


@pytest.fixture
def mol_link2():
    return Chem.MolFromSmiles("CN")


@pytest.fixture
def mol_macrocycle():
    # CC(N)CCCO yields macrocycle results at radius=1 against CORPUS
    return Chem.MolFromSmiles("CC(N)CCCO")


@pytest.fixture(scope="session")
def db_rc(tmp_path_factory):
    """DB built with --frag-mode both from a corpus rich in saturated rings
    and 1,2-disubstituted benzenes; used by ring-closure tests."""
    import os
    here = os.path.dirname(__file__)
    smi_file = os.path.join(here, "data", "ring_closures.smi")
    d = tmp_path_factory.mktemp("db_rc")
    db_path = str(d / "ring_closures.db")
    subprocess.run(
        [sys.executable, "-m", "crem.scripts.cremdb_create",
         "-i", smi_file, "-o", db_path,
         "-s", "test", "--radii", "1", "2", "3", "--ncpu", "1",
         "--frag-mode", "both"],
        check=True, capture_output=True,
    )
    return db_path


@pytest.fixture(scope="session")
def db_acyclic(tmp_path_factory):
    """Same corpus as db_rc but built with --frag-mode acyclic — used to
    exercise the legacy-DB error path of make_cycle(ring_closures=True)
    against a DB that has the column but no ring rows (and to cross-check
    that the error path triggers cleanly when the column is absent in
    older DBs that predate this feature)."""
    import os
    here = os.path.dirname(__file__)
    smi_file = os.path.join(here, "data", "ring_closures.smi")
    d = tmp_path_factory.mktemp("db_acyclic")
    db_path = str(d / "acyclic_only.db")
    subprocess.run(
        [sys.executable, "-m", "crem.scripts.cremdb_create",
         "-i", smi_file, "-o", db_path,
         "-s", "test", "--radii", "1", "2", "3", "--ncpu", "1",
         "--frag-mode", "acyclic"],
        check=True, capture_output=True,
    )
    return db_path
