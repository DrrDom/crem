# Fragment properties

You can annotate fragments with molecular properties and then use those
properties to constrain which fragments are selected during generation. This is
a two-part workflow: **add** property columns to the database, then **filter**
on them at generation time.

## Add properties with `cremdb_add_prop`

```bash
cremdb_add_prop -i fragments.db -p mw logp rtb tpsa fcsp3 -c 8
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | — *(required)* | CReM fragment database |
| `-p`, `--properties` | all five | Properties to compute (see below) |
| `-c`, `--ncpu` | `1` | Number of CPUs |
| `--fetch-batch` | `50000` | Rows fetched per batch |
| `--write-batch` | `10000` | Rows per write |
| `--imap-chunk` | `500` | Multiprocessing chunk size |
| `-v`, `--verbose` | off | Print progress to stderr |

Built-in properties:

| Name | Property |
|---|---|
| `mw` | molecular weight |
| `logp` | Crippen logP |
| `rtb` | rotatable bonds |
| `tpsa` | topological polar surface area |
| `fcsp3` | fraction of sp³ carbons |

The command is schema-aware: in a v1 or v2 database it adds columns to the `frags`
table; in a v0 database it adds them to each `radiusN` table. It is also
**incremental** — only rows whose property values are still `NULL` are
computed, so rerunning after adding new fragments fills just the new rows.

## Filter on properties at generation time

Pass property ranges as extra keyword arguments to any generation function. A
value may be a single number (exact match) or a `(low, high)` tuple (inclusive
bounds):

```python
from rdkit import Chem
from crem.crem import mutate_mol

mols = list(mutate_mol(
    Chem.MolFromSmiles("c1ccccc1N"),
    db_name="fragments.db",
    set_names="chembl",
    min_freq=5,
    mw=(50, 180),
    logp=(0.0, 3.5),
))
```

!!! note "Where property columns live"
    The keyword names must correspond to existing columns. In a **v1/v2** database
    these are columns of `frags` (or `frags_h`); in a **v0** database they are
    columns of the `radiusN` tables. This is the one generation-time difference
    between the two formats. Run `cremdb_info -i fragments.db` to list the property
    columns a database actually has.

## Add custom properties (Python API)

`crem.db.add_fragment_props` adds the built-in properties and, optionally,
arbitrary custom columns computed by your own functions:

```python
from crem.db import add_fragment_props
from rdkit import Chem
from rdkit.Chem import Descriptors

def num_rings(smi):
    return Chem.MolFromSmiles(smi).GetRingInfo().NumRings()

# Built-ins only
add_fragment_props("fragments.db")

# Built-ins plus a custom column
add_fragment_props(
    "fragments.db",
    properties="all",
    custom_props={"nrings": num_rings},
    ncpu=8,
)

# Only a custom column, targeting the H-collapsed table
add_fragment_props(
    "fragments.db",
    custom_props={"heavy": lambda s: Descriptors.HeavyAtomCount(Chem.MolFromSmiles(s))},
    table="frags_h",
)
```

Picklable functions (named functions, `functools.partial`) run on `ncpu`
workers; lambdas and closures run serially. Custom properties can target either
`frags` (default, keyed on `core_smi`) or `frags_h` (keyed on the H-capped
SMILES). See [`crem.db`](../reference/db.md) for the full signature.
