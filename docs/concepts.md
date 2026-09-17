# Concepts

This page defines the vocabulary used throughout the documentation and in the
function signatures. Understanding these terms makes the parameters of
`mutate_mol`, `grow_mol`, `link_mols`, and `make_cycle` straightforward.

## Fragment, core, and context

CReM cuts a molecule into two complementary parts:

- The **core** (also called the *fragment*) — the piece that will be replaced.
- The **context** (also called the *environment*, `env`) — the surrounding
  atoms that stay in place.

The bond(s) that were cut leave **attachment points**, written as `[*:1]`,
`[*:2]`, … in SMILES. A core and a context fit together at matching attachment
points.

```text
   molecule:        Cc1ccc(OC)cc1
   one cut gives →  context:  Cc1ccc([*:1])cc1
                    core:     [*:1]OC
```

A **replacement** is the act of swapping one core for another core that was
seen attached to the *same* context somewhere in the fragment database.

## Context radius

The full context of a fragment can be large, so CReM only keeps the atoms
within a fixed number of bonds from each attachment point. That number is the
**radius**.

- `radius=1` keeps only the directly bonded atoms — a permissive context that
  matches many fragments.
- `radius=3` (the default) keeps three bonds of context — a good balance
  between chemical specificity and database coverage.

A larger radius means a stricter match (fewer, more context-appropriate
replacements); a smaller radius means a looser match (more, but less
context-aware replacements). The radius used at generation time must exist as a
`radiusN` table in the database.

`radius=0` is the limiting case: no context is matched at all, only the number of
attachment points and which of them close a ring. It has to be built explicitly and
returns very large result sets — see [Radius 0](fragment-databases/radius0.md).

## Fragment size and size windows

Several parameters constrain the size (number of heavy atoms) of the fragments
involved:

- `min_size` / `max_size` — size of the **core being replaced** (`mutate_mol`).
  `min_size=0` allows replacing a hydrogen.
- `min_inc` / `max_inc` — how much the **replacing** fragment may change the
  heavy-atom count relative to the replaced one. `min_inc=-2, max_inc=2` allows
  the new fragment to be up to two atoms smaller or larger.
- `min_atoms` / `max_atoms` — size of the **incoming** fragment for `grow_mol`,
  `link_mols`, and `make_cycle`.
- `min_rel_size` / `max_rel_size` — size of the replaced fragment relative to
  the whole molecule.

## Fragment frequency and sets

Each context–core relationship in the database carries a **frequency**: how many
times that pair was observed in the molecules used to build the database.

- `min_freq` keeps only replacements seen at least that many times — a simple
  way to favour common, well-precedented chemistry.

A database can hold several **fragment sets** side by side (for example
`chembl`, `natural_products`, `my_library`). Each set is a separate frequency
column, so the same database can describe how common a fragment is in different
collections.

- `set_names` selects which set column(s) the `min_freq` threshold applies to.
  When several sets are named, a fragment passes if **any** of them meets the
  threshold (OR logic). `set_names=None` (default) uses all available sets.

See [Fragment sets](fragment-databases/fragment-sets.md) for how to build and
use them.

## Attachment-point count and `dist2`

- A **mutate** fragment (`mutate_mol`) can have **one to four** attachment
  points.
- A **grow** fragment (`grow_mol`) has **one** attachment point (it replaces a
  single hydrogen).
- A **linker** (`link_mols`) or a **ring-closing** fragment (`make_cycle`) has
  **two** attachment points.
- Partial-ring fragments may have **two to four** attachment points.

For two-attachment-point fragments, `dist2` is the topological distance (in
bonds) between the two attachment points, which is stored in the database. 
It lets CReM control linker geometry and ring size:

- `dist` (in `link_mols`) filters linkers by `dist2`.
- `ring_size` (in `make_cycle`) is translated per anchor pair into a `dist2`
  filter.

## Fragmentation modes (`frag_mode`)

When **building** a database you choose how molecules are fragmented. This is
the `--frag-mode` option of `cremdb_create` (and the `frag_mode` argument of
`crem.db.create_db`):

| `frag_mode` | What it cuts                                                                           | `is_ring_closure` |
|---|----------------------------------------------------------------------------------------|---|
| `acyclic` | acyclic (non-ring) single bonds — classic MMPA cuts                                    | `0` |
| `ring` | pairs of single bonds inside one ring (ring arcs) + exhaustive acyclic side cuts | `1` |
| `ring_optimal` | ring arcs + only *exo* side cuts adjacent to the arc                                   | `1` |
| `both` | `acyclic` + `ring`                                                                     | `0` and `1` |
| `both_optimal` *(default)* | `acyclic` + `ring_optimal`                                                             | `0` and `1` |

The *optimal* modes (`ring_optimal`, `both_optimal`) emit only the *exo* side
cuts — a **subset** of the exhaustive cuts emitted by `ring`/`both` — so they
build a smaller database. Because the exo rows are nested inside the exhaustive
ones, an exhaustive `ring`/`both` database also satisfies the exo queries, but
not the other way round.

Rows produced by ring cutting are flagged with `is_ring_closure = 1`. They are
queried by [`make_cycle`](operations/make-cycle.md) when `ring_closures=True`
and by [`mutate_mol`](operations/mutate-grow-link.md#replacing-cyclic-source-fragments)
when `replace_cycles` is `"partial_all"` or `"partial_exo"`. Acyclic rows
(`is_ring_closure = 0`) drive ordinary mutate/grow/link and
`replace_cycles="forced"`.

The separate, integer **`mode`** option (`0` all atoms, `1` heavy atoms only,
`2` hydrogen atoms only) controls whether hydrogen cuts are produced within the
acyclic fragmenter. `min_size=0` / `grow_mol` need hydrogen-cut rows, which mode
`0` produces.

## Database format: v2, v1 and v0

CReM reads three database formats:

- **v2** — the current, deduplicated schema (`PRAGMA user_version = 2`) built by
  [`cremdb_create`](fragment-databases/build-v2.md). It stores fragment sets and
  shared environment/fragment tables, and marks ring-cut attachment points with an
  isotope label so ring provenance lives in the fragment SMILES.
- **v1** — the deprecated predecessor of v2, which recorded ring provenance in an
  `is_ring_closure` column instead. Still readable, but it raises a `FutureWarning`
  and there is no converter to v2 — the upgrade is a rebuild from SMILES.
- **v0** — the legacy single-table-per-radius layout produced by the
  [pipeline](fragment-databases/build-v0.md)
  (`fragmentation` → `frag_to_env` → `env_to_db`).

All three formats work with every generation function. The
[schema page](fragment-databases/schema.md) describes them, and
[Convert a v0 database](fragment-databases/convert.md) shows how to upgrade a v0 one.

The only place the format is visible at generation time is property filtering
via `**kwargs`: those filters read columns from `radiusN` in a v0 database and
from `frags` / `frags_h` in v1 and v2.
