# Radius 0 — no environment constraint

At radius 0 a replacement is not required to match its surroundings at all. The only thing
a candidate fragment has to agree on is **how many attachment points it has and which of
them close a ring** — everything else about the context is ignored.

That makes it the most permissive mode CReM offers. It is a deliberately exotic switch:
useful for exploring chemical space far from what the database has actually seen, and for
generating structures that are formally valid but unusual. It is not a good default.

!!! warning "Result sets are enormous"

    A radius-0 query is matched only by attachment-point class, so a single four-point cut
    site can match millions of fragments. Always pass `max_replacements`. With
    `max_replacements=None` CReM will retrieve every matching row, which on a
    ChEMBL-scale database is not something your machine will enjoy.

## The `RxAy` environment classes

Since no context survives, the radius-0 `env` is not a SMILES string but the attachment
point class, written `R<x>A<y>` for x ring-cut points and y acyclic ones:

| env | meaning | typical source |
|---|---|---|
| `R0A1` | one acyclic point | `grow`, single-cut `mutate` |
| `R0A2` | two acyclic points | `link`, broad `make_cycle` |
| `R0A3`, `R0A4` | three / four acyclic points | multi-cut `mutate` |
| `R2A0` | a ring arc, nothing else | strict `make_cycle`, `replace_cycles` |
| `R2A1`, `R2A2` | a ring arc with one / two exo substituents | `replace_cycles='partial_*'` |

Ring cuts always come in a pair, so x is only ever 0 or 2 — verified across every row of
the chembl36 fragment table. That leaves exactly these seven classes.

The class string is deliberately not valid SMILES: it cannot be mistaken for a real
[env](schema.md#v2-schema-current), and it states the counts without implying any
particular attachment-point numbering.

## The `radius0` table

Identical in shape to the other [radius tables](schema.md#radiusn-context-to-fragment-mapping-one-per-radius):

```sql
CREATE TABLE radius0(
    env_id          INTEGER NOT NULL,   -- -> the RxAy row in envs
    core_smi_id     INTEGER NOT NULL,
    core_num_atoms  INTEGER NOT NULL,
    dist2           INTEGER NOT NULL,
    chembl          INTEGER NOT NULL DEFAULT 0,
    UNIQUE (env_id, core_smi_id)
);
```

`dist2` follows the usual convention — it is populated for the two-point classes `R0A2`
and `R2A0` and left at 0 otherwise.

### Every orientation is a separate row

A fragment with several attachment points can be spliced in more than one way, so each
distinct **orientation** is stored as its own row. Uniform sampling over rows is then
uniform over (fragment, orientation), which is what a mode designed for unbiased
exploration needs.

Orientations that differ only by a symmetry of the fragment itself are **not** stored
twice, because they would splice to the identical product:

```
C([1*:1])[1*:2]  and  C([1*:2])[1*:1]   ->  one row, not two
O(C[*:1])[*:2]   and  O(C[*:2])[*:1]    ->  two rows: the points are not equivalent
```

Without that collapse a fully symmetric four-point fragment would be drawn up to 24 times
more often than an asymmetric one.

## Building it

### With real occurrence counts — during the build

Include 0 in `--radii`:

```bash
cremdb_create -i input.smi -o fragments.db -s chembl --radii 0 1 2 3
```

Counts are then true occurrence counts: one increment per fragmentation event, so
[`min_freq`](fragment-sets.md) works exactly as it does at other radii.

### From an existing database — `cremdb_radius0`

```bash
cremdb_radius0 -i fragments.db -c 8
```

This derives the table from `frags` alone, which is possible because under
[v2](schema.md#v2-schema-current) the ring-cut isotope label puts the provenance in the
fragment string. It requires a v2 database and refuses to overwrite an existing `radius0`
unless `--force` is given.

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | — *(required)* | v2 database to extend |
| `-c`, `--ncpu` | `1` | Worker processes |
| `--batch-size` | `50000` | Rows per batch |
| `--force` | off | Overwrite an existing `radius0` |
| `--quiet` | off | Suppress progress output |

!!! warning "A derived table has no occurrence counts"

    Counts in a radius table are occurrences multiplied by the environment's symmetry
    orbit size, and that factor depends on the radius — so no sum over them recovers a
    radius-independent occurrence count. A derived table therefore stores **1** where a
    fragment belongs to a [fragment set](fragment-sets.md) and **0** where it does not.
    Per-set membership *is* exact, because membership is a boolean OR rather than a sum.

    Consequences: `min_freq=0` and `min_freq=1` behave as "no filter", `min_freq >= 2`
    matches nothing, and `return_rxn_freq` reports 1. For real counts, rebuild with
    `cremdb_create --radii 0`.

Both routes produce the **same rows** — they differ only in the counts. The script also
adds any fragment orientations that no radius ≥ 1 row happens to reference, so `frags`
grows.

## Generating at radius 0

Nothing special is required beyond `radius=0`:

```python
from rdkit import Chem
from crem.crem import mutate_mol

m = Chem.MolFromSmiles('CCC1CCOCC1')
for smi in mutate_mol(m, 'fragments.db', radius=0, max_size=4, max_replacements=100):
    print(smi)
```

All four operations work — `mutate` (including every `replace_cycles` mode), `grow`,
`link` and `make_cycle` in both ring modes. Because radius 0 constrains nothing beyond the
attachment-point class, its output is a **superset** of what the same query returns at any
higher radius.

## What to expect

Radius 0 does not produce valence-broken structures: CReM joins fragments at matching
attachment points, so the result is always a chemically valid molecule. What it does
produce is a large, cheap supply of molecules that are valid but unusual — higher
synthetic-accessibility scores, more structural alerts, and in the ring-forming modes a
tendency toward larger rings than a context-aware radius would have closed.
