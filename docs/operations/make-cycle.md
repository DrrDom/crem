# Make cycle

`make_cycle` forms a **new ring** by connecting two atoms of the same molecule
with a two-attachment-point fragment from the database. It covers both
macrocyclization and the closure of smaller native rings.

```python
from rdkit import Chem
from crem.crem import make_cycle

m = Chem.MolFromSmiles("c1ccccc1N")
res = list(make_cycle(
    m,
    db_name="fragments.db",
    radius=3,
    ring_size=(5, 7),
    min_atoms=1,
    max_atoms=10,
    replace_ids=[1, 2],
))
```

`make_cycle` works on an H-expanded copy of the molecule internally.

## Two modes

The `ring_closures` argument selects which database rows are queried:

### `ring_closures=True` (default)

**Strict.** Only ring-closure (arc-cut) fragments are queried — rows with
`is_ring_closure = 1`. Useful for closing native, typically aliphatic rings.
Requires a database built so that ring-closure rows exist:
[`cremdb_create --frag-mode`](../fragment-databases/build-v2.md) with `ring`,
`both`, `ring_optimal`, or `both_optimal`. If the database has no such rows, the
call raises a clear error.

### `ring_closures=False`

**Broad.** Any linker fragment may be used; the `is_ring_closure` provenance is
not filtered. This may be useful for transplanting *intact* aromatic or scaffold
rings, since the donor fragment's own ring closures are preserved. This
can also create macrocycles.

```python
# Broad mode: any linker fragment
res = list(make_cycle(m, db_name="fragments.db", ring_closures=False,
                      ring_size=(5, 7), max_atoms=10))
```

## Key parameters

| Parameter | Meaning |
|---|---|
| `ring_size` | Size of the **new** ring (atoms = bonds). `int` for a single size, `(min, max)` for a window, `None` for no constraint. Translated per anchor pair into a `dist2` filter as `ring_size − d_in`, where `d_in` is the topological distance between the two anchor atoms in the input molecule. |
| `ring_closures` | Which provenance to query (see above). Default `True`. |
| `min_atoms` / `max_atoms` | Heavy-atom size window of the linking fragment. Defaults `1` / `10`. |
| `replace_ids` / `protected_ids` | Restrict which atoms may serve as ring-closure anchors. |
| `set_names` / `min_freq` | Fragment set and frequency threshold. |

`symmetry_fixes` is accepted for API compatibility with the other generation
functions but is not used by `make_cycle`.

## Prefer acyclic attachment points

When forming rings you often want the linker to attach to acyclic atoms rather
than to atoms already inside a ring. The built-in
`crem.utils.filter_acyclic_attachment_points` filter enforces this:

```python
from crem.utils import filter_acyclic_attachment_points

res = list(make_cycle(
    m,
    db_name="fragments.db",
    radius=3,
    ring_size=(5, 7),
    max_atoms=10,
    replace_ids=[1, 2],
    filter_func=filter_acyclic_attachment_points,
))
```

## Parallel use

For `multiprocessing`, use the list-returning wrapper `make_cycle2` — see
[Multiprocessing](multiprocessing.md).
