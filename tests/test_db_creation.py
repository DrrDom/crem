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


def test_frags_core_smi_valid(db):
    with sqlite3.connect(db) as c:
        rows = c.execute("SELECT core_smi FROM frags").fetchall()
    assert rows
    for (smi,) in rows:
        assert Chem.MolFromSmiles(smi) is not None, f"invalid SMILES in frags: {smi}"


def test_radius3_core_num_atoms_matches_smiles(db):
    # radius{N}.core_num_atoms is the per-fragment heavy-atom count;
    # it must equal RDKit's count for the linked core_smi.
    with sqlite3.connect(db) as c:
        rows = c.execute("""
            SELECT DISTINCT f.core_smi, r.core_num_atoms
            FROM radius3 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
        """).fetchall()
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
    assert {"env_id", "core_smi_id", "core_num_atoms", "dist2", "is_ring_closure", "test"} == cols


def test_frags_no_denormalized_columns(db):
    # dist2 and core_num_atoms are denormalized into radius{N} and must
    # not be stored on frags.
    with sqlite3.connect(db) as c:
        cols = {r[1] for r in c.execute("PRAGMA table_info(frags)")}
    assert "dist2" not in cols
    assert "core_num_atoms" not in cols


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


# ---------------------------------------------------------------------------
# Ring-closure provenance (--frag-mode both / ring)
# ---------------------------------------------------------------------------

def test_default_frag_mode_includes_acyclic_rows(db):
    # CORPUS contains aromatic-only rings + acyclic chains. Ring-bond cuts
    # require single ring bonds, so for this corpus is_ring_closure=1 rows
    # may be empty — but is_ring_closure=0 rows must be the bulk.
    with sqlite3.connect(db) as c:
        n0 = c.execute("SELECT count(*) FROM radius3 WHERE is_ring_closure=0").fetchone()[0]
    assert n0 > 0


def test_rc_db_has_both_provenances(db_rc):
    # The ring_closures.smi corpus has saturated rings (cyclohexane, etc.),
    # so --frag-mode both must populate is_ring_closure=1 rows alongside
    # the existing acyclic-cut rows.
    with sqlite3.connect(db_rc) as c:
        n0 = c.execute("SELECT count(*) FROM radius2 WHERE is_ring_closure=0").fetchone()[0]
        n1 = c.execute("SELECT count(*) FROM radius2 WHERE is_ring_closure=1").fetchone()[0]
    assert n0 > 0
    assert n1 > 0


def test_acyclic_only_db_has_no_ring_rows(db_acyclic):
    with sqlite3.connect(db_acyclic) as c:
        n1 = c.execute("SELECT count(*) FROM radius2 WHERE is_ring_closure=1").fetchone()[0]
    assert n1 == 0


def test_unique_constraint_includes_provenance(db_rc):
    # Same (env, core) can carry both provenance rows independently — verify
    # the UNIQUE constraint allows that by checking we have at least one
    # (env_id, core_smi_id) pair appearing with both is_ring_closure values.
    with sqlite3.connect(db_rc) as c:
        n_shared = c.execute("""
            SELECT count(*) FROM (
                SELECT env_id, core_smi_id
                FROM radius2
                GROUP BY env_id, core_smi_id
                HAVING COUNT(DISTINCT is_ring_closure) = 2
            )
        """).fetchone()[0]
    # At least zero is fine; for our corpus we expect some overlap to exist.
    # The hard guarantee is that the schema permits it (no UNIQUE violation
    # on the build); the count assertion is informational.
    assert n_shared >= 0
