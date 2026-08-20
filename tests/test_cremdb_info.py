import json
import shutil
import sqlite3
import subprocess
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from crem.db import add_fragment_props, get_db_info
from crem.scripts.cremdb_create import DB_SCHEMA_VERSION, _STRIDE_SHARD_SENTINEL


def run_info(*args):
    return subprocess.run([sys.executable, "-m", "crem.scripts.cremdb_info", *args],
                          capture_output=True, text=True)


def heavy_atoms(smi):
    """Module-level (picklable) custom property function for add_fragment_props."""
    return Descriptors.HeavyAtomCount(Chem.MolFromSmiles(smi))


@pytest.fixture
def v0_db(tmp_path):
    """v0 database: one table per radius, no sets. radius1 carries a property column."""
    path = tmp_path / "v0.db"
    with sqlite3.connect(path) as con:
        con.execute("CREATE TABLE radius1(env TEXT, core_smi TEXT, core_num_atoms INTEGER, "
                    "core_sma TEXT, dist2 INTEGER, freq INTEGER, mw REAL)")
        con.execute("CREATE TABLE radius3(env TEXT, core_smi TEXT, core_num_atoms INTEGER, "
                    "core_sma TEXT, dist2 INTEGER, freq INTEGER)")
    return str(path)


@pytest.fixture
def old_v1_db(tmp_path):
    """v1 database as older builds wrote it, with core_num_atoms on frags as well as on
    the radius table. Both are schema columns and must not be reported as properties."""
    path = tmp_path / "old_v1.db"
    with sqlite3.connect(path) as con:
        con.execute("PRAGMA user_version = 1")
        con.execute("CREATE TABLE envs(env_id INTEGER PRIMARY KEY, env TEXT)")
        con.execute("CREATE TABLE frags(core_smi_id INTEGER PRIMARY KEY, core_smi TEXT, "
                    "core_num_atoms INTEGER, core_smi_h_id INTEGER, mw REAL)")
        con.execute("CREATE TABLE frags_h(core_smi_h_id INTEGER PRIMARY KEY, smi TEXT)")
        con.execute("CREATE TABLE radius1(env_id INTEGER, core_smi_id INTEGER, "
                    "core_num_atoms INTEGER, dist2 INTEGER, is_ring_closure INTEGER, "
                    "chembl INTEGER)")
    return str(path)


def test_current_schema(db):
    info = get_db_info(db)
    assert info['version'] == DB_SCHEMA_VERSION
    assert info['has_sets'] is True
    assert info['radius_tables'] == {'radius1': ['test'], 'radius2': ['test'], 'radius3': ['test']}
    assert info['properties'] == {'frags': [], 'frags_h': []}
    assert info['is_shard'] is False


def test_v0_has_no_sets(v0_db):
    info = get_db_info(v0_db)
    assert info['version'] == 0
    assert info['has_sets'] is False
    # the pre-cremdb_info script reported env / core_smi / core_sma / freq as set names here
    assert info['radius_tables'] == {'radius1': [], 'radius3': []}
    assert info['properties'] == {'radius1': ['mw'], 'radius3': []}


def test_v1_schema_columns_are_not_properties(old_v1_db):
    info = get_db_info(old_v1_db)
    assert info['version'] == 1
    assert info['radius_tables'] == {'radius1': ['chembl']}
    assert info['properties'] == {'frags': ['mw'], 'frags_h': []}


def test_properties_on_frags_h(db, tmp_path):
    target = tmp_path / "props.db"
    shutil.copy(db, target)
    add_fragment_props(target, custom_props={'heavy': heavy_atoms}, table='frags_h')
    info = get_db_info(target)
    assert info['properties']['frags_h'] == ['heavy']
    assert info['properties']['frags'] == []


def test_unmerged_shard_is_flagged(db, tmp_path):
    shard = tmp_path / "shard.db"
    shutil.copy(db, shard)
    with sqlite3.connect(shard) as con:
        con.execute(f"PRAGMA application_id = {_STRIDE_SHARD_SENTINEL}")
    assert get_db_info(shard)['is_shard'] is True
    assert "unmerged stride shard" in run_info("-i", str(shard)).stdout


def test_missing_file_is_not_created(tmp_path):
    missing = tmp_path / "nope.db"
    result = run_info("-i", str(missing))
    assert result.returncode == 1
    assert "file not found" in result.stdout
    assert not missing.exists()


def test_not_a_crem_database(tmp_path):
    path = tmp_path / "other.db"
    with sqlite3.connect(path) as con:
        con.execute("CREATE TABLE foo(a INTEGER)")
    result = run_info("-i", str(path))
    assert result.returncode == 1
    assert "not a CReM database" in result.stdout


def test_several_databases_in_argv_order(db, tmp_path):
    second = tmp_path / "second.db"
    shutil.copy(db, second)
    result = run_info("-i", db, str(second))
    assert result.returncode == 0
    assert result.stdout.count("schema version") == 2
    assert result.stdout.index(db) < result.stdout.index(str(second))


def test_json_output(db, tmp_path):
    missing = tmp_path / "nope.db"
    result = run_info("-i", db, str(missing), "--json")
    assert result.returncode == 1
    records = json.loads(result.stdout)
    assert records[0]['version'] == DB_SCHEMA_VERSION
    assert records[0]['radius_tables']['radius1'] == ['test']
    assert 'error' in records[1]


def test_deprecated_script_keeps_its_output_format(db):
    result = subprocess.run(
        [sys.executable, "-m", "crem.scripts.cremdb_get_set_names", "-i", db],
        capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout == ("radius1: ['test']\n"
                             "radius2: ['test']\n"
                             "radius3: ['test']\n")
    assert "deprecated" in result.stderr
