# Convert v0 to v1

`cremdb_convert` migrates a legacy [v0 database](schema.md#v0-schema-legacy) to
the deduplicated [v1 schema](schema.md#v1-schema-current). The old `freq` column
becomes a named [fragment set](fragment-sets.md) in the new database.

## Basic conversion

```bash
cremdb_convert -i old.db -o fragments.db --set-name chembl
```

This creates `fragments.db` with the v1 tables (`envs`, `frags`, `frags_h`,
`radiusN`) and one fragment-set column (`chembl`) populated from the v0 `freq`
values.

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | — *(required)* | Existing v0 database |
| `-o`, `--output` | — *(required)* | New v1 database (created) |
| `--set-name` | `undefined` | Name of the set column to create and fill from the old `freq` |
| `--radii` | `1 2 3 4 5` | Radii to convert |
| `--batch-size` | `10000` | Rows processed per batch |
| `--verify` | off | Verify the conversion after completion |
| `--quiet` | off | Suppress progress output |

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

New `radiusN` tables also gain an `is_ring_closure` column (default `0`); v0
databases contain only acyclic-cut fragments, so converted rows are all
acyclic. To obtain ring-closure rows, build with
[`cremdb_create --frag-mode`](build-v1.md).

## Naming and safety

- `--set-name` must be a valid SQLite identifier and cannot be `env_id` or
  `core_smi_id`.
- If the output file already exists, `cremdb_convert` asks for confirmation
  before overwriting it.
