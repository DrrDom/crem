# Command-line tools

Installing CReM adds the following console commands to your `PATH`.

## Command overview

| Command | Purpose | Guide |
|---|---|---|
| `cremdb_create` | Build a v2 database from SMILES (recommended) | [Build (v2)](../fragment-databases/build-v2.md) |
| `cremdb_convert` | Convert a v0 database to v2 (or v1) | [Convert](../fragment-databases/convert.md) |
| `cremdb_add_prop` | Add property columns to a database | [Properties](../fragment-databases/properties.md) |
| `cremdb_info` | Show the schema version, fragment sets and property columns of one or more databases | [Fragment sets](../fragment-databases/fragment-sets.md) |
| `cremdb_get_set_names` | *Deprecated* — use `cremdb_info` | [Fragment sets](../fragment-databases/fragment-sets.md) |
| `cremdb_radius0` | Derive a radius0 table from an existing v2 database | [Radius 0](../fragment-databases/radius0.md) |
| `cremdb_merge` | Merge shard databases into one | [Build (v2)](../fragment-databases/build-v2.md#merge-shards-manually-with-cremdb_merge) |
| `fragmentation` | v0 pipeline — fragment molecules | [Build (v0)](../fragment-databases/build-v0.md) |
| `frag_to_env` | v0 pipeline — standardize context/core | [Build (v0)](../fragment-databases/build-v0.md) |
| `env_to_db` | v0 pipeline — import into a v0 database | [Build (v0)](../fragment-databases/build-v0.md) |
| `crem_create_frag_db.sh` | v0 pipeline — all-in-one shell wrapper | [Build (v0)](../fragment-databases/build-v0.md) |
| `guacamol_test` | Run the GuacaMol benchmark | [Benchmarks](../about/benchmarks.md) |

---

## `cremdb_info`

Print the schema version, the fragment-set columns of each radius table, and the
property columns of one or more databases. The file is opened read-only and only
PRAGMAs are read, so the command is instant whatever the database size.

```bash
cremdb_info -i fragments.db
cremdb_info -i *.db --json
```

```text
fragments.db
  schema version : 2 (current)
  radius1        : chembl, drugs
  radius2        : chembl, drugs
  radius3        : chembl, drugs
  properties     : frags: mw, logp
```

Each database is reported as a block headed by the path as it was typed. A database
that cannot be read gets an `error:` line instead of a block, and the command exits
with status 1 while still reporting the others.

| Option | Description |
|---|---|
| `-i`, `--input` | One or more paths to databases *(required)* |
| `--json` | Print the same information as a JSON list |

The same information is available from Python as
[`crem.db.get_db_info`](db.md).

---

## `cremdb_get_set_names`

!!! warning "Deprecated"

    Superseded by [`cremdb_info`](#cremdb_info), which also reports the schema version
    and property columns and accepts several databases at once. This command still
    works and prints the same output as before, preceded by a deprecation notice on
    stderr, but it will be removed in a future release.

List the fragment-set columns in each radius table.

```bash
cremdb_get_set_names -i fragments.db
```

| Option | Description |
|---|---|
| `-i`, `--input` | Path to the database *(required)* |

---

## `guacamol_test`

Run the GuacaMol goal-directed benchmark with a CReM-based generator. Requires
the optional `guacamol` package (and `pandas`, `joblib`, `numpy`).

```bash
guacamol_test --smiles_file train.smi --db_fname fragments.db --output_dir results/
```

Key options: `--smiles_file`, `--db_fname`, `--selection_size` (`10`),
`--radius` (`3`), `--replacements` (`1000`), `--min_size` (`0`), `--max_size`
(`10`), `--min_inc` (`-7`), `--max_inc` (`7`), `--generations` (`1000`),
`--ncpu` (`1`), `--seed` (`42`), `--suite` (`v2`), `--output_dir`. Results are
written as `goal_directed_results.json` plus per-task `.smi` files. See
[Benchmarks](../about/benchmarks.md).
