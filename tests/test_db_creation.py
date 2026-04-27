import sqlite3

from rdkit import Chem


def test_user_version(db):
    with sqlite3.connect(db) as c:
        assert c.execute("PRAGMA user_version").fetchone()[0] == 1


def test_required_tables(db):
    with sqlite3.connect(db) as c:
        tables = {r[0] for r in c.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
    assert {"envs", "frags_h", "frags", "radius1", "radius2", "radius3"}.issubset(tables)


def test_envs_not_empty(db):
    with sqlite3.connect(db) as c:
        assert c.execute("SELECT count(*) FROM envs").fetchone()[0] > 0


def test_frags_h_valid_smiles(db):
    with sqlite3.connect(db) as c:
        smiles = [r[0] for r in c.execute("SELECT smi FROM frags_h")]
    assert smiles
    assert all(Chem.MolFromSmiles(s) is not None for s in smiles)


def test_frags_core_smi_valid_and_atom_count_matches(db):
    with sqlite3.connect(db) as c:
        rows = c.execute("SELECT core_smi, core_num_atoms FROM frags").fetchall()
    assert rows
    for smi, nha in rows:
        m = Chem.MolFromSmiles(smi)
        assert m is not None, f"invalid SMILES in frags: {smi}"
        assert m.GetNumHeavyAtoms() == nha, (
            f"core_num_atoms mismatch for {smi}: stored {nha}, "
            f"actual {m.GetNumHeavyAtoms()}"
        )


def test_radius3_columns(db):
    with sqlite3.connect(db) as c:
        cols = {r[1] for r in c.execute("PRAGMA table_info(radius3)")}
    assert {"env_id", "core_smi_id", "core_num_atoms", "dist2", "test"} == cols


def test_radius3_not_empty(db):
    with sqlite3.connect(db) as c:
        assert c.execute("SELECT count(*) FROM radius3").fetchone()[0] > 0


def test_covering_index_exists(db):
    with sqlite3.connect(db) as c:
        indices = {r[0] for r in c.execute(
            "SELECT name FROM sqlite_master WHERE type='index'"
        )}
    assert "idx_radius3_lookup" in indices


def test_frags_h_foreign_key_integrity(db):
    with sqlite3.connect(db) as c:
        orphans = c.execute("""
            SELECT count(*) FROM frags f
            LEFT JOIN frags_h h ON f.core_smi_h_id = h.core_smi_h_id
            WHERE h.core_smi_h_id IS NULL
        """).fetchone()[0]
    assert orphans == 0


def test_radius_test_counts_positive(db):
    # Single set DB: every row must have been counted at least once
    with sqlite3.connect(db) as c:
        zeros = c.execute(
            "SELECT count(*) FROM radius3 WHERE test <= 0"
        ).fetchone()[0]
    assert zeros == 0


def test_env_ids_consistent(db):
    # radius3.env_id must reference a real envs row
    with sqlite3.connect(db) as c:
        orphans = c.execute("""
            SELECT count(*) FROM radius3 r
            LEFT JOIN envs e ON r.env_id = e.env_id
            WHERE e.env_id IS NULL
        """).fetchone()[0]
    assert orphans == 0


def test_core_smi_ids_consistent(db):
    with sqlite3.connect(db) as c:
        orphans = c.execute("""
            SELECT count(*) FROM radius3 r
            LEFT JOIN frags f ON r.core_smi_id = f.core_smi_id
            WHERE f.core_smi_id IS NULL
        """).fetchone()[0]
    assert orphans == 0


def test_denormalized_core_num_atoms_matches_frags(db):
    # core_num_atoms in radius3 must equal the value stored in frags
    with sqlite3.connect(db) as c:
        mismatches = c.execute("""
            SELECT count(*) FROM radius3 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
            WHERE r.core_num_atoms != f.core_num_atoms
        """).fetchone()[0]
    assert mismatches == 0
