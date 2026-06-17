# Command-line tools

Installing CReM adds the following console commands to your `PATH`.

## Command overview

| Command | Purpose | Guide |
|---|---|---|
| `cremdb_create` | Build a v1 database from SMILES (recommended) | [Build (v1)](../fragment-databases/build-v1.md) |
| `cremdb_convert` | Convert a v0 database to v1 | [Convert](../fragment-databases/convert.md) |
| `cremdb_add_prop` | Add property columns to a database | [Properties](../fragment-databases/properties.md) |
| `cremdb_get_set_names` | List the fragment sets in a database | [Fragment sets](../fragment-databases/fragment-sets.md) |
| `cremdb_merge` | Merge shard databases into one | [Build (v1)](../fragment-databases/build-v1.md#merge-shards-manually-with-cremdb_merge) |
| `fragmentation` | v0 pipeline — fragment molecules | [Build (v0)](../fragment-databases/build-v0.md) |
| `frag_to_env` | v0 pipeline — standardize context/core | [Build (v0)](../fragment-databases/build-v0.md) |
| `env_to_db` | v0 pipeline — import into a v0 database | [Build (v0)](../fragment-databases/build-v0.md) |
| `crem_create_frag_db.sh` | v0 pipeline — all-in-one shell wrapper | [Build (v0)](../fragment-databases/build-v0.md) |
| `guacamol_test` | Run the GuacaMol benchmark | [Benchmarks](../about/benchmarks.md) |

---

## `cremdb_get_set_names`

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
