# Build a database (v2)

`cremdb_create` builds a v2 fragment database directly from a SMILES file in one
step. It is the recommended way to create new databases. The same functionality
is available programmatically through [`crem.db.create_db`](#python-api).

## Basic command

```bash
cremdb_create -i input.smi -o fragments.db -s chembl
```

- `-i/--input` — input SMILES file (plain text, or `.zst`-compressed).
- `-o/--output` — output SQLite database.
- `-s/--set-name` — name of the [fragment set](fragment-sets.md) to create. May
  also be one or more membership files (see [Fragment sets](fragment-sets.md)).

## Input format

One molecule per line: a SMILES, optionally followed by an ID.

```text
CCO         mol_0001
c1ccccc1    mol_0002
```

The ID (second column) is optional, but it is **required** when assigning
molecules to sets via membership files. Columns are whitespace-separated by
default; use `--sep` for a custom delimiter.

## Common options

```bash
cremdb_create -i input.smi -o fragments.db -s chembl \
  --radii 1 2 3 4 5 \
  --ncpu 16 \
  --frag-mode both_optimal \
  --max-heavy-atoms 15
```

| Option | Default | Description |
|---|---|---|
| `-r`, `--radii` | `1 2 3 4 5` | Context radii to build |
| `-c`, `--ncpu` | `1` | Worker processes (capped at available CPUs) |
| `--max-heavy-atoms` | `15` | Maximum heavy atoms in a core fragment |
| `--mode` | `0` | Acyclic cut mode: `0` all atoms, `1` heavy only, `2` H only |
| `--frag-mode` | `both_optimal` | Fragmentation source: `acyclic`, `ring`, `both`, `ring_optimal`, `both_optimal` (see [frag modes](../concepts.md#fragmentation-modes-frag_mode)) |
| `--keep-stereo` | off | Retain stereochemistry in env/core |
| `--sep` | whitespace | Input column delimiter |
| `--chunk-size` | `100` | Input lines per worker task |
| `--flush-every` | `100` | Chunks accumulated in memory before each DB flush |
| `--prefetch` | `4` | In-flight task batches per worker |
| `--zstd` | off | Force zstd decompression of the input |
| `--log-every` | off | Print a progress line every N chunks |
| `--timings` | off | Print per-flush timing breakdown to stderr |
| `--fragment-error-log` | off | Write fragment validation issues to `<output>.errors` instead of stderr |

!!! note "`--frag-mode` and ring-based generation"
    Ring-derived rows (`is_ring_closure = 1`) are needed by
    [`make_cycle(..., ring_closures=True)`](../operations/make-cycle.md) **and**
    by [`mutate_mol(..., replace_cycles="partial_all"`/`"partial_exo")`](../operations/mutate-grow-link.md#replacing-cyclic-source-fragments).
    The *optimal* modes store only *exo* side cuts — a subset of the exhaustive
    `ring`/`both` cuts — so they are smaller. Therefore:

    * `partial_exo` and `make_cycle` work with **any** ring-capable mode
      (`ring`, `both`, `ring_optimal`, `both_optimal`).
    * `partial_all` needs the **exhaustive** `ring` or `both` (on an *optimal*
      DB it under-matches — the non-exo cuts are absent).

    Ordinary mutate/grow/link and `replace_cycles="forced"` use only acyclic
    rows and work with any mode. The default `both_optimal` covers ordinary
    generation, ring closure, and `partial_exo`.

Indices are created automatically at the end of the run.

## Building on an existing database

Running `cremdb_create` against an existing database is **additive**: new radii
and new set names are added, and existing data is preserved. This lets you
extend a database with more molecules or more sets later.

## Large datasets: sharded and parallel builds

For very large inputs, build in shards and merge:

- `--shard-size N` — write at most `N` input structures per shard 
  database created sequentially, then merge the shards into the output at the end. 
  Shard 0 is the output file; later shards get a `_NNN` suffix.
- `--parallel-shards N` — build `N` shards concurrently, each on a stride of the
  input, splitting `--ncpu` evenly across them. Intermediate parts live in
  `<output>.parts/` and are merged with a parallel binary-tree reduction.

```bash
# 8 concurrent shard builders across 32 CPUs
cremdb_create -i big_input.smi -o fragments.db -s chembl \
  --parallel-shards 8 --ncpu 32
```

`--shard-size` and `--parallel-shards > 1` are mutually exclusive.

### Merge shards manually with `cremdb_merge`

Shards or individual databases of the same schema version can also be merged by hand — for example 
to combine the per-shard databases from `--shard-size`, or to merge shards built 
on different machines. `cremdb_merge` merges source databases into an existing 
target; it is idempotent and resumable, so already-absorbed sources are skipped.

```bash
cremdb_merge -t base.db -i shard_001.db shard_002.db shard_003.db
```

| Option | Default | Description |
|---|---|---|
| `-t`, `--target` | — *(required)* | Target database; must already exist with schema and data |
| `-i`, `--input` | — *(required)* | One or more source shard databases to merge in |
| `--no-index` | off | Skip index creation after merge (useful when more shards will follow) |
| `--parallel N` | `1` | Merge with a binary-tree reduction using up to N concurrent pair-merges per round |

The same operation is available in Python as
[`crem.db.merge_dbs`](../reference/db.md).

## Resumable runs

```bash
cremdb_create -i input.smi -o fragments.db -s chembl \
  --processed-chunks chunks.done --log-every 100
```

If the run is interrupted, rerun the **same** command: chunks recorded in
`chunks.done` are skipped. Parallel and sharded builds manage their own
resume markers internally, so simply rerunning the command resumes them.

## Naming rules of fragment sets

- A set name must be a valid SQLite identifier: `[A-Za-z_][A-Za-z0-9_]*`.
- The reserved names `env_id` and `core_smi_id` are not allowed.

## Python API

`crem.db.create_db` exposes the same build process and accepts either a file
path or an iterable of `"SMILES [ID]"` strings:

```python
from crem.db import create_db

create_db(
    "input.smi",
    "fragments.db",
    set_name="chembl",
    radii=(1, 2, 3, 4, 5),
    ncpu=16,
    frag_mode="both_optimal",
    max_heavy_atoms=15,
)
```

To assign molecules to multiple sets, pass a dict mapping each set name to a set
of IDs (or `None` for "all molecules"):

```python
create_db(
    "input.smi",
    "fragments.db",
    set_name={"chembl": None, "focused": {"mol_0001", "mol_0007"}},
)
```

See [`crem.db`](../reference/db.md) for the full signature, and
[Fragment sets](fragment-sets.md) for the set model.
