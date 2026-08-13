import functools
import sqlite3
from datetime import datetime

import pytest

from crem.db import add_fragment_props, create_db, merge_dbs
from crem.scripts.cremdb_create import DB_SCHEMA_VERSION

CORPUS_A = [
    "CCO mol1", "c1ccccc1 mol2", "CCCO mol3",
    "c1ccccc1N mol4", "CC(=O)O mol5",
]
CORPUS_B = [
    "CCNCC mol6", "c1ccc(O)cc1 mol7", "CC(C)O mol8",
    "c1ccncc1 mol9", "CCOCC mol10",
]
CORPUS_RING = ["CC1CCC(CCC)CC1 ring1"]


def _frag_count(db):
    with sqlite3.connect(db) as c:
        return c.execute("SELECT count(*) FROM frags").fetchone()[0]


def _tables(db):
    with sqlite3.connect(db) as c:
        return {r[0] for r in c.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}


def _indices(db):
    with sqlite3.connect(db) as c:
        return {r[0] for r in c.execute(
            "SELECT name FROM sqlite_master WHERE type='index'"
        )}


def _cols(db, table):
    with sqlite3.connect(db) as c:
        return {r[1] for r in c.execute(f"PRAGMA table_info({table})")}


def _null_count(db, table, col):
    with sqlite3.connect(db) as c:
        return c.execute(
            f"SELECT count(*) FROM {table} WHERE {col} IS NULL"
        ).fetchone()[0]


def _provenance_counts(db, radius=2):
    """Rows per provenance, keyed 0 (acyclic cut) / 1 (ring cut).

    On v2 there is no is_ring_closure column: every ring-cut attachment point carries
    isotope 1, so the label in core_smi *is* the provenance.
    """
    with sqlite3.connect(db) as c:
        return dict(c.execute(
            f"SELECT CASE WHEN f.core_smi LIKE '%[1*%' THEN 1 ELSE 0 END AS is_ring_closure, "
            f"       count(*) "
            f"FROM radius{radius} r JOIN frags f ON r.core_smi_id = f.core_smi_id "
            f"GROUP BY is_ring_closure"
        ).fetchall())


# ---------------------------------------------------------------------------
# create_db
# ---------------------------------------------------------------------------

def test_create_from_file_path(tmp_path):
    smi = tmp_path / "mols.smi"
    smi.write_text("\n".join(CORPUS_A))
    db = str(tmp_path / "test.db")
    create_db(str(smi), db, "chembl", radii=(1, 2, 3), verbose=False)
    assert _frag_count(db) > 0


def test_create_from_path_object(tmp_path):
    smi = tmp_path / "mols.smi"
    smi.write_text("\n".join(CORPUS_A))
    db = tmp_path / "test.db"
    create_db(smi, db, "chembl", radii=(1, 2, 3), verbose=False)
    assert _frag_count(db) > 0


def test_create_from_iterable(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, "chembl", radii=(1, 2, 3), verbose=False)
    assert _frag_count(db) > 0


def test_create_custom_separator_forwarded(tmp_path):
    smi = tmp_path / "mols.smi"
    smi.write_text("CCO|mol1\nCCCO|mol2\n")
    db = str(tmp_path / "test.db")
    create_db(smi, db, "chembl", radii=(1,), sep="|", verbose=False)
    assert _frag_count(db) > 0


def test_create_valid_schema(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, "chembl", radii=(1, 2, 3), verbose=False)
    with sqlite3.connect(db) as c:
        ver = c.execute("PRAGMA user_version").fetchone()[0]
    assert ver == DB_SCHEMA_VERSION
    assert {"envs", "frags_h", "frags", "radius1", "radius2", "radius3"}.issubset(_tables(db))
    assert "idx_radius3_lookup" in _indices(db)


def test_create_only_requested_radii(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, "chembl", radii=(1, 2), verbose=False)
    tbls = _tables(db)
    assert "radius1" in tbls and "radius2" in tbls
    assert "radius3" not in tbls


def test_create_is_additive(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, "chembl", radii=(1, 2, 3), verbose=False)
    n1 = _frag_count(db)
    create_db(CORPUS_B, db, "chembl", radii=(1, 2, 3), verbose=False)
    n2 = _frag_count(db)
    assert n2 > n1


def test_create_set_name_str_creates_column(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, "myset", radii=(1, 2, 3), verbose=False)
    assert "myset" in _cols(db, "radius3")


def test_create_set_name_dict_all_mols(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, {"myset": None}, radii=(1, 2, 3), verbose=False)
    assert "myset" in _cols(db, "radius3")
    assert _frag_count(db) > 0


def test_create_set_name_dict_filtered(tmp_path):
    db_all = str(tmp_path / "all.db")
    db_filtered = str(tmp_path / "filtered.db")
    # All molecules in one set
    create_db(CORPUS_A, db_all, "s", radii=(1, 2, 3), verbose=False)
    # Only the first two molecule IDs
    id_filter = {"mol1", "mol2"}
    create_db(CORPUS_A, db_filtered, {"s": id_filter}, radii=(1, 2, 3), verbose=False)
    # Filtered DB should have fewer or equal fragments
    assert _frag_count(db_filtered) <= _frag_count(db_all)


def test_create_max_heavy_atoms(tmp_path):
    db_strict = str(tmp_path / "strict.db")
    db_loose = str(tmp_path / "loose.db")
    create_db(CORPUS_A, db_strict, "s", radii=(1, 2, 3), max_heavy_atoms=3, verbose=False)
    create_db(CORPUS_A, db_loose, "s", radii=(1, 2, 3), max_heavy_atoms=15, verbose=False)
    assert _frag_count(db_strict) <= _frag_count(db_loose)


def test_create_frag_mode_default_matches_both_optimal(tmp_path):
    db_default = str(tmp_path / "default.db")
    db_explicit = str(tmp_path / "explicit.db")
    create_db(CORPUS_RING, db_default, "s", radii=(2,), verbose=False)
    create_db(CORPUS_RING, db_explicit, "s", radii=(2,), frag_mode="both_optimal", verbose=False)
    assert _provenance_counts(db_default) == _provenance_counts(db_explicit)


def test_create_frag_mode_argument_forwarded(tmp_path):
    db = str(tmp_path / "ring.db")
    create_db(CORPUS_RING, db, "s", radii=(2,), frag_mode="ring", verbose=False)
    counts = _provenance_counts(db)
    assert counts.get(0, 0) == 0
    assert counts.get(1, 0) > 0


def test_create_fragment_error_log_argument_forwarded(tmp_path):
    db = tmp_path / "test.db"
    create_db(CORPUS_A, db, "s", radii=(1,), fragment_error_log=True, verbose=False)
    err = tmp_path / "test.errors"

    assert err.exists()


def test_create_processed_chunks_forwarded_for_file_input(tmp_path):
    smi = tmp_path / "mols.smi"
    smi.write_text("\n".join(CORPUS_A))
    db = tmp_path / "test.db"
    processed = tmp_path / "processed.txt"
    create_db(
        smi,
        db,
        "s",
        radii=(1,),
        processed_chunks=processed,
        chunk_size=1,
        flush_every=1,
        verbose=False,
    )

    assert processed.exists()
    assert processed.read_text().strip()


def test_create_processed_chunks_ignored_for_iterable(tmp_path):
    db = tmp_path / "test.db"
    processed = tmp_path / "processed.txt"
    processed.write_text("0\n")
    create_db(
        CORPUS_A,
        db,
        "s",
        radii=(1,),
        processed_chunks=processed,
        chunk_size=1,
        verbose=False,
    )

    assert _frag_count(db) > 0
    assert processed.read_text() == "0\n"


def test_create_parallel_arguments_forwarded(tmp_path):
    db = tmp_path / "parallel.db"
    create_db(
        CORPUS_RING,
        db,
        "s",
        radii=(1,),
        parallel_shards=2,
        ncpu=2,
        merge_parallel=1,
        fragment_error_log=True,
        verbose=False,
    )

    assert _frag_count(db) > 0
    assert (tmp_path / "parallel.errors").exists()


def test_create_invalid_set_name_type(tmp_path):
    db = str(tmp_path / "test.db")
    with pytest.raises(TypeError):
        create_db(CORPUS_A, db, 123, radii=(1, 2, 3), verbose=False)


# ---------------------------------------------------------------------------
# merge_dbs
# ---------------------------------------------------------------------------

def test_merge_combines_fragments(tmp_path):
    db_a = str(tmp_path / "a.db")
    db_b = str(tmp_path / "b.db")
    create_db(CORPUS_A, db_a, "s", radii=(1, 2, 3), verbose=False)
    create_db(CORPUS_B, db_b, "s", radii=(1, 2, 3), verbose=False)
    n_a = _frag_count(db_a)
    n_b = _frag_count(db_b)
    merge_dbs(db_a, [db_b], verbose=False)
    n_merged = _frag_count(db_a)
    assert n_merged > max(n_a, n_b)


def test_merge_no_duplicates(tmp_path):
    # Same corpus split into two DBs — merged total <= sum (dedup on duplicates)
    db_a = str(tmp_path / "a.db")
    db_b = str(tmp_path / "b.db")
    create_db(CORPUS_A, db_a, "s", radii=(1, 2, 3), verbose=False)
    create_db(CORPUS_A, db_b, "s", radii=(1, 2, 3), verbose=False)
    n_a = _frag_count(db_a)
    n_b = _frag_count(db_b)
    merge_dbs(db_a, [db_b], verbose=False)
    assert _frag_count(db_a) <= n_a + n_b


def test_merge_index_rebuilt(tmp_path):
    db_a = str(tmp_path / "a.db")
    db_b = str(tmp_path / "b.db")
    create_db(CORPUS_A, db_a, "s", radii=(1, 2, 3), verbose=False)
    create_db(CORPUS_B, db_b, "s", radii=(1, 2, 3), verbose=False)
    merge_dbs(db_a, [db_b], rebuild_index=True, verbose=False)
    assert "idx_radius3_lookup" in _indices(db_a)


def test_merge_rebuild_index_false(tmp_path):
    db_a = str(tmp_path / "a.db")
    db_b = str(tmp_path / "b.db")
    create_db(CORPUS_A, db_a, "s", radii=(1, 2, 3), verbose=False)
    create_db(CORPUS_B, db_b, "s", radii=(1, 2, 3), verbose=False)
    merge_dbs(db_a, [db_b], rebuild_index=False, verbose=False)
    # Covering index was dropped before merge and not recreated
    assert "idx_radius3_lookup" not in _indices(db_a)


def test_merge_history_timestamp_is_human_readable(tmp_path):
    db_a = str(tmp_path / "a.db")
    db_b = str(tmp_path / "b.db")
    create_db(CORPUS_A, db_a, "s", radii=(1, 2, 3), verbose=False)
    create_db(CORPUS_B, db_b, "s", radii=(1, 2, 3), verbose=False)
    merge_dbs(db_a, [db_b], verbose=False)

    with sqlite3.connect(db_a) as conn:
        columns = {
            row[1]: row[2]
            for row in conn.execute("PRAGMA table_info(_merge_history)")
        }
        source_id, merged_at, storage_type = conn.execute(
            "SELECT source_id, merged_at, typeof(merged_at) FROM _merge_history"
        ).fetchone()

    assert columns["merged_at"].upper() == "TEXT"
    assert source_id == "b.db"
    assert storage_type == "text"
    datetime.strptime(merged_at, "%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# add_fragment_props — built-in
# ---------------------------------------------------------------------------

@pytest.fixture
def small_db(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A + CORPUS_B, db, "s", radii=(1, 2, 3), verbose=False)
    return db


def test_props_builtin_columns_added(small_db):
    add_fragment_props(small_db, ["mw", "logp"])
    cols = _cols(small_db, "frags")
    assert "mw" in cols and "logp" in cols


def test_props_no_nulls_after_builtin(small_db):
    add_fragment_props(small_db, ["mw", "logp", "rtb"])
    assert _null_count(small_db, "frags", "mw") == 0
    assert _null_count(small_db, "frags", "logp") == 0
    assert _null_count(small_db, "frags", "rtb") == 0


def test_props_none_computes_all(small_db):
    add_fragment_props(small_db, None)
    for col in ("mw", "logp", "rtb", "tpsa", "fcsp3"):
        assert _null_count(small_db, "frags", col) == 0


def test_props_subset_only_requested_columns(small_db):
    add_fragment_props(small_db, ["mw"])
    cols = _cols(small_db, "frags")
    assert "mw" in cols
    assert "logp" not in cols


def test_props_idempotent(small_db):
    add_fragment_props(small_db, ["mw"])
    with sqlite3.connect(small_db) as c:
        values_first = c.execute("SELECT core_smi_id, mw FROM frags ORDER BY core_smi_id").fetchall()
    add_fragment_props(small_db, ["mw"])
    with sqlite3.connect(small_db) as c:
        values_second = c.execute("SELECT core_smi_id, mw FROM frags ORDER BY core_smi_id").fetchall()
    assert values_first == values_second


def test_props_only_fills_null_rows(tmp_path):
    db = str(tmp_path / "test.db")
    create_db(CORPUS_A, db, "s", radii=(1, 2, 3), verbose=False)
    add_fragment_props(db, ["mw"])
    # Manually overwrite one row with a sentinel value
    with sqlite3.connect(db) as c:
        row_id = c.execute("SELECT core_smi_id FROM frags LIMIT 1").fetchone()[0]
        c.execute("UPDATE frags SET mw = 999.0 WHERE core_smi_id = ?", (row_id,))
        c.commit()
    # Extend DB with new molecules (their mw will be NULL)
    create_db(CORPUS_B, db, "s", radii=(1, 2, 3), verbose=False)
    add_fragment_props(db, ["mw"])
    # Sentinel must be preserved — re-run does not touch non-NULL rows
    with sqlite3.connect(db) as c:
        mw_val = c.execute("SELECT mw FROM frags WHERE core_smi_id = ?", (row_id,)).fetchone()[0]
    assert mw_val == 999.0
    # All rows must now have mw
    assert _null_count(db, "frags", "mw") == 0


# ---------------------------------------------------------------------------
# add_fragment_props — custom
# ---------------------------------------------------------------------------

def test_props_custom_on_frags(small_db):
    add_fragment_props(small_db, custom_props={"n_c": lambda smi: smi.count("C")})
    assert "n_c" in _cols(small_db, "frags")
    assert _null_count(small_db, "frags", "n_c") == 0


def test_props_custom_on_frags_h(small_db):
    add_fragment_props(small_db, custom_props={"n_c_h": lambda smi: smi.count("C")},
                      table="frags_h")
    assert "n_c_h" in _cols(small_db, "frags_h")
    assert _null_count(small_db, "frags_h", "n_c_h") == 0


def test_props_custom_multiple_at_once(small_db):
    add_fragment_props(small_db, custom_props={
        "n_c": lambda smi: smi.count("C"),
        "n_n": lambda smi: smi.count("N"),
    })
    assert "n_c" in _cols(small_db, "frags")
    assert "n_n" in _cols(small_db, "frags")
    assert _null_count(small_db, "frags", "n_c") == 0
    assert _null_count(small_db, "frags", "n_n") == 0


def test_props_custom_picklable_function(small_db):
    def count_c(smi):
        return smi.count("C")
    add_fragment_props(small_db, custom_props={"n_c2": count_c})
    assert _null_count(small_db, "frags", "n_c2") == 0


def test_props_custom_partial_function(small_db):
    import functools
    count_char = functools.partial(str.count, "C")
    add_fragment_props(small_db, custom_props={"n_c3": lambda smi: count_char(smi)})
    assert _null_count(small_db, "frags", "n_c3") == 0


def test_props_custom_idempotent(small_db):
    add_fragment_props(small_db, custom_props={"n_c": lambda smi: smi.count("C")})
    with sqlite3.connect(small_db) as c:
        vals_first = c.execute("SELECT core_smi_id, n_c FROM frags ORDER BY core_smi_id").fetchall()
    # Second call: no NULL rows remain, so nothing changes
    add_fragment_props(small_db, custom_props={"n_c": lambda smi: -1})
    with sqlite3.connect(small_db) as c:
        vals_second = c.execute("SELECT core_smi_id, n_c FROM frags ORDER BY core_smi_id").fetchall()
    assert vals_first == vals_second


def test_props_invalid_table(small_db):
    with pytest.raises(ValueError, match="table must be one of"):
        add_fragment_props(small_db, custom_props={"x": lambda s: 0}, table="nonexistent")


def test_props_empty_properties_skips_builtin_computes_custom(small_db):
    # Explicit `[]` continues to skip built-ins (back-compat).
    add_fragment_props(small_db, [], custom_props={"n_c": lambda smi: smi.count("C")})
    cols = _cols(small_db, "frags")
    assert "n_c" in cols
    assert "mw" not in cols


def test_props_default_with_custom_skips_builtins(small_db):
    # Omitting `properties` while passing `custom_props` should NOT compute
    # built-ins — the user only asked for custom columns.
    add_fragment_props(small_db, custom_props={"n_c": lambda smi: smi.count("C")})
    cols = _cols(small_db, "frags")
    assert "n_c" in cols
    for col in ("mw", "logp", "rtb", "tpsa", "fcsp3"):
        assert col not in cols


def test_props_explicit_none_with_custom_computes_both(small_db):
    # Explicit `properties=None` (or 'all') opts back in to all built-ins.
    add_fragment_props(small_db, None, custom_props={"n_c": lambda smi: smi.count("C")})
    cols = _cols(small_db, "frags")
    assert "n_c" in cols
    for col in ("mw", "logp", "rtb", "tpsa", "fcsp3"):
        assert _null_count(small_db, "frags", col) == 0


def test_props_all_string_alias(small_db):
    add_fragment_props(small_db, "all")
    for col in ("mw", "logp", "rtb", "tpsa", "fcsp3"):
        assert _null_count(small_db, "frags", col) == 0
