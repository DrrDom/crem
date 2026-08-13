# Convert a v0 database

`cremdb_convert` migrates a legacy [v0 database](schema.md#v0-schema-legacy) to
the deduplicated [v2 schema](schema.md#v2-schema-current) (or to
[v1](schema.md#v1-schema-deprecated) with `--target-version 1`). The old `freq`
column becomes a named [fragment set](fragment-sets.md) in the new database.

!!! warning "v0 sources only"

    `cremdb_convert` rejects a v1 or v2 input. There is no v1 → v2 upgrade path:
    relabelling a v1 database's ring arcs would require re-standardising its stored
    `env` strings, and a radius-truncated `env` is not a valid molecule — feeding one
    back through canonicalisation fails or silently alters its aromaticity. To move a
    v1 database to v2, rebuild it from the source SMILES with
    [`cremdb_create`](build-v2.md).

## Basic conversion

```bash
cremdb_convert -i old.db -o fragments.db --set-name chembl
```

This creates `fragments.db` with the v2 tables (`envs`, `frags`, `frags_h`,
`radiusN`) and one fragment-set column (`chembl`) populated from the v0 `freq`
values.

| Option | Default | Description                                                        |
|---|---|--------------------------------------------------------------------|
| `-i`, `--input` | — *(required)* | Existing v0 database                                               |
| `-o`, `--output` | — *(required)* | New database (created)                                             |
| `--target-version` | `2` | Schema version to write: `2` (current) or `1` (will be deprecated) |
| `--set-name` | `undefined` | Name of the set column to create and fill from the old `freq`      |
| `--radii` | `1 2 3 4 5` | Radii to convert                                                   |
| `--batch-size` | `10000` | Rows processed per batch                                           |
| `--verify` | off | Verify the conversion after completion                             |
| `--quiet` | off | Suppress progress output                                           |

## Convert selected radii

```bash
cremdb_convert -i old.db -o fragments.db --radii 1 2 3 --set-name chembl
```

## Verify the conversion

```bash
cremdb_convert -i old.db -o fragments.db --set-name chembl --verify
```

`--verify` samples rows from each radius and checks that the reconstructed
environment, core, and SMARTS match the source database.

## What is preserved

- Environment, core SMILES, and `dist2` information.
- The old `freq` values, mapped into the new `--set-name` column.
- Additional non-standard fragment columns from the old database are carried
  over into the `frags` table where possible.

## What a v0 source cannot provide

v0 predates ring-bond fragmentation entirely: it has no provenance column and every
row it holds is an acyclic cut. A converted database therefore contains **no
ring-closure fragments**, whichever target version you choose, and
`replace_cycles='partial_*'` / `make_cycle(ring_closures=True)` will return nothing
against it. To obtain ring-closure fragments, build from SMILES with
[`cremdb_create --frag-mode`](build-v2.md).

This also means a v0 → v2 conversion equals a from-scratch v2 build only when that
build used `--frag-mode acyclic`.

## Naming and safety

- `--set-name` must be a valid SQLite identifier and cannot be `env_id` or
  `core_smi_id`.
- If the output file already exists, `cremdb_convert` asks for confirmation
  before overwriting it.
