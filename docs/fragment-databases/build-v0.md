# Build a database (v0)

The original CReM database builder is a three-stage text pipeline:
`fragmentation` → `frag_to_env` → `env_to_db`. It produces a
[v0 database](schema.md#v0-schema-legacy).

!!! tip "Use `cremdb_create` for new databases"
    The one-step [`cremdb_create`](build-v1.md) builds the richer v1 format
    (fragment sets, ring-closure provenance, smaller files). Use the v0 pipeline
    only when you specifically need it; a v0 database can later be
    [converted to v1](convert.md).

## All-in-one wrapper

The shipped shell script runs the whole pipeline. It takes the input SMILES
file, an output directory for intermediate files and the final database, and an
optional CPU count (default `1`):

```bash
crem_create_frag_db.sh input.smi fragdb_dir 32
```

## Step by step

The manual steps below give you control over each stage and over which radii are
built.

### 1. Fragment the input structures

```bash
fragmentation -i input.smi -o frags.txt -c 32 -v
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | — | Input SMILES, optionally with an ID column |
| `-o`, `--out` | — | Output fragments file |
| `-s`, `--sep` | whitespace | Input column delimiter |
| `-d`, `--sep_out` | `,` | Output delimiter |
| `-m`, `--mode` | `0` | `0` all atoms, `1` heavy only, `2` H only |
| `-c`, `--ncpu` | `1` | Number of CPUs |
| `-v`, `--verbose` | off | Print progress |

### 2. Convert fragments to standardized context/core at a radius

```bash
frag_to_env -i frags.txt -o r3.txt -r 3 -c 32 -v
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | — | Fragmented molecules from step 1 |
| `-o`, `--out` | — | Output text file |
| `-d`, `--sep` | `,` | Input delimiter |
| `-r`, `--radius` | `1` | Context radius (in bonds) |
| `-a`, `--max_heavy_atoms` | `20` | Maximum heavy atoms in a core; larger fragments are discarded |
| `-k`, `--keep_mols` | — | File of molecule names to keep (others ignored) |
| `-s`, `--keep_stereo` | off | Keep stereochemistry in env/core |
| `-c`, `--ncpu` | `1` | Number of CPUs |
| `-v`, `--verbose` | off | Print progress |

The output may contain duplicate lines by design — they are counted in the next
step.

### 3. Count occurrences

`sort` and `uniq` are standard shell utilities:

```bash
sort r3.txt | uniq -c > r3_c.txt
```

This prepends an occurrence count to each unique line, which becomes the `freq`
column.

### 4. Import into the database

```bash
env_to_db -i r3_c.txt -o fragments.db -r 3 -c -v
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | — | Counted env/core file from step 3 |
| `-o`, `--out` | — | Output SQLite database |
| `-r`, `--radius` | — *(required)* | Radius of this table; an existing table for the radius is dropped |
| `-c`, `--counts` | off | Input has a leading occurrence count → adds a `freq` column |
| `-n`, `--ncpu` | `1` | Number of CPUs |
| `-v`, `--verbose` | off | Print progress |

## Building multiple radii

Repeat steps 2–4 for each radius, writing into the **same** database file. Each
radius becomes its own `radiusN` table.

```bash
for r in 1 2 3 4 5; do
  frag_to_env -i frags.txt -o r${r}.txt -r ${r} -c 32 -v
  sort r${r}.txt | uniq -c > r${r}_c.txt
  env_to_db -i r${r}_c.txt -o fragments.db -r ${r} -c -v
done
```

## Next step

To use fragment sets, ring closures, and the smaller v1 layout, convert the
result: see [Convert v0 to v1](convert.md).
