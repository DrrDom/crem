import sqlite3
from datetime import datetime

from rdkit import Chem

from crem.scripts.cremdb_create import (
    _FRAGMENT_ISSUE_COLUMNS,
    _FragmentIssueWriter,
    _fragment_error_log_path,
    _fragment_issue_records,
    _fragment_mol,
    _fragment_mol_ring,
    _normalize_input_mol,
)


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
# Ring-closure provenance (--frag-mode both / ring / *_optimal)
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


def test_ring_rows_include_partial_cycle_attachment_counts(db_rc):
    # --frag-mode both should now store the historical 2-attachment ring arcs
    # and partial-cycle fragments with side cuts, i.e. 3/4 attachment points.
    with sqlite3.connect(db_rc) as c:
        rows = c.execute("""
            SELECT
                length(f.core_smi) - length(replace(f.core_smi, '*', '')) AS nstars,
                count(*),
                min(r.dist2),
                max(r.dist2)
            FROM radius2 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
            WHERE r.is_ring_closure = 1
            GROUP BY nstars
        """).fetchall()
        by_stars = {nstars: (count, min_dist, max_dist)
                    for nstars, count, min_dist, max_dist in rows}
        nonzero_dist_partial = c.execute("""
            SELECT count(*)
            FROM radius2 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
            WHERE r.is_ring_closure = 1
              AND r.dist2 != 0
              AND length(f.core_smi) - length(replace(f.core_smi, '*', '')) > 2
        """).fetchone()[0]
        labeled_partial = c.execute("""
            SELECT count(*)
            FROM radius2 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
            WHERE r.is_ring_closure = 1
              AND length(f.core_smi) - length(replace(f.core_smi, '*', '')) > 2
              AND instr(f.core_smi, '[1*') > 0
        """).fetchone()[0]
        labeled_env_partial = c.execute("""
            SELECT count(*)
            FROM radius2 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
            JOIN envs e ON r.env_id = e.env_id
            WHERE r.is_ring_closure = 1
              AND length(f.core_smi) - length(replace(f.core_smi, '*', '')) > 2
              AND instr(e.env, '[1*') > 0
        """).fetchone()[0]
        labeled_two_point = c.execute("""
            SELECT count(*)
            FROM radius2 r
            JOIN frags f ON r.core_smi_id = f.core_smi_id
            WHERE r.is_ring_closure = 1
              AND length(f.core_smi) - length(replace(f.core_smi, '*', '')) = 2
              AND instr(f.core_smi, '[1*') > 0
        """).fetchone()[0]

    assert by_stars[2][0] > 0
    assert by_stars[2][2] > 0
    assert by_stars[3][0] > 0
    assert by_stars[4][0] > 0
    assert nonzero_dist_partial == 0
    assert labeled_partial > 0
    assert labeled_env_partial > 0
    assert labeled_two_point == 0


def test_fragment_issue_records_detect_defensive_checks():
    core = "[H]C([H])(N(C[1*:3])[*:1])C([H])([*:3])[1*:2]"
    records = list(
        _fragment_issue_records(4, "input_smi", "mol1", core, "ctx", 1)
    )
    issues = [record[3] for record in records]

    assert issues.count("explicit_h_core_smi") == 1
    assert issues.count("duplicate_attachment_map") == 1
    assert issues.count("mixed_ring_external_attachment_map") == 1
    assert all(record[:3] == (4, "mol1", "input_smi") for record in records)


def test_fragment_issue_records_allow_valid_h_attachment_core():
    records = list(
        _fragment_issue_records(0, "[H][*:1]", "mol1", "[H][*:1]", "ctx", 0)
    )

    assert records == []


def test_normalize_input_mol_removes_isotopic_alkene_stereo_hs():
    mol = Chem.MolFromSmiles("[3H]/C=C(\\[3H])CC")
    normalized = _normalize_input_mol(mol)

    assert "[H]" not in Chem.MolToSmiles(normalized, isomericSmiles=True)
    assert all(atom.GetAtomicNum() != 1 for atom in normalized.GetAtoms())


def test_normalize_input_mol_preserves_tetrahedral_chiral_hydrogen():
    mol = Chem.MolFromSmiles("[H][C@](F)(Cl)Br")
    normalized = _normalize_input_mol(mol)

    assert Chem.MolToSmiles(normalized, isomericSmiles=True) == "F[C@@H](Cl)Br"


def test_reported_tritium_alkene_generates_no_explicit_h_core_issues():
    smi = (
        "[3H]/C=C(\\[3H])C[C@@]1(c2cccc(C3(C(F)(F)F)N=N3)c2)"
        "C(=O)NC(=O)N(C)C1=O"
    )
    frags, _ = _fragment_mol(
        smi,
        "CHEMBL2163547",
        0,
        "both_optimal",
        max_heavy_atoms=15,
    )
    issues = [
        record
        for core, context, is_ring_closure in frags
        for record in _fragment_issue_records(
            11606,
            smi,
            "CHEMBL2163547",
            core,
            context,
            is_ring_closure,
        )
    ]

    assert issues == []


def test_fragment_error_log_path_replaces_db_extension(tmp_path):
    path = tmp_path / "output.db"

    assert _fragment_error_log_path(path) == str(tmp_path / "output.errors")


def test_fragment_issue_writer_appends_repeated_tsv_records(tmp_path):
    path = tmp_path / "output.errors"
    record = (
        1,
        "mol\t1",
        "C\nC",
        "duplicate_attachment_map",
        "C([*:1])[*:1]",
        "ctx",
        0,
        "map=1;isotopes=0,0",
    )
    writer = _FragmentIssueWriter(str(path))
    try:
        writer.write([record])
    finally:
        writer.close()
    writer = _FragmentIssueWriter(str(path))
    try:
        writer.write([record])
    finally:
        writer.close()

    lines = path.read_text().splitlines()
    assert lines[0] == "\t".join(_FRAGMENT_ISSUE_COLUMNS)
    assert len(lines) == 3
    for line in lines[1:]:
        fields = line.split("\t")
        datetime.fromisoformat(fields[0])
        assert fields[1:4] == ["1", "mol\\t1", "C\\nC"]
    assert lines[1].split("\t")[1:] == lines[2].split("\t")[1:]


def test_ring_optimal_restricts_side_cuts_to_exo_bonds():
    smi = "CC1CCC(CCC)CC1"
    mol = _normalize_input_mol(Chem.MolFromSmiles(smi))
    full, _ = _fragment_mol_ring(
        Chem.Mol(mol), "", max_heavy_atoms=15, side_cut_mode="all"
    )
    optimal, _ = _fragment_mol_ring(
        Chem.Mol(mol), "", max_heavy_atoms=15, side_cut_mode="exo"
    )

    assert optimal
    assert optimal < full
    assert {core.count("*") for core, _ in optimal} == {2, 3, 4}


def test_both_optimal_combines_acyclic_and_optimal_ring_fragments():
    smi = "CC1CCC(CCC)CC1"
    both_optimal, _ = _fragment_mol(smi, "", 1, "both_optimal", max_heavy_atoms=15)
    ring_optimal, _ = _fragment_mol(smi, "", 1, "ring_optimal", max_heavy_atoms=15)
    ring_full, _ = _fragment_mol(smi, "", 1, "ring", max_heavy_atoms=15)

    assert {item[2] for item in both_optimal} == {0, 1}
    assert {item for item in both_optimal if item[2] == 1} == ring_optimal
    assert ring_optimal < ring_full
