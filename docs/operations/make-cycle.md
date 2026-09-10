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

## Spiro closures

Both ends of the linker can also land on one and the same atom, closing a new
ring through it. That atom must be in a ring, so the new ring meets the existing
one there and nowhere else — a spiro centre. The `spiro` argument is off by
default:

| `spiro` | Behaviour |
|---|---|
| `False` (default) | Two-atom closures only. |
| `True` | Same-atom closures **in addition to** the two-atom ones. |
| `'only'` | Same-atom closures alone. |

```python
# spiro[2.4]heptane from cyclopentane: both linker ends on one CH2
res = list(make_cycle(Chem.MolFromSmiles("C1CCCC1"), db_name="fragments.db",
                      ring_size=3, spiro="only"))
```

Eligible atoms are ring atoms with two replaceable hydrogens. A ring atom carries
at least two ring bonds and therefore never more than two hydrogens, so "has two"
already names both of them and the pair is unambiguous; the same rule excludes
aromatic atoms, which is right — a spiro centre is sp3. An **acyclic** atom is not
eligible: a closure through one of those is an ordinary core replacement, already
produced by [`mutate_mol`](mutate-grow-link.md).

A same-atom context can only be matched by a fragment that was itself cut out of a
ring — the two cut bonds and the fragment close a cycle through the atom — so
`spiro='only'` needs a database with ring-cut rows (`--frag-mode ring|both`) and
returns the same rows whichever way `ring_closures` is set. `'only'` also skips the
two-cut MMPA pass that the two-atom enumeration needs, which makes a focused spiro
scan considerably cheaper than filtering the output of a full run.

Under `replace_ids` each named eligible atom is used as a spiro centre, while pairs
of named atoms go on closing ordinary rings between them; naming a single atom
therefore leaves its spiro closure as the only possibility.

## Key parameters

| Parameter | Meaning |
|---|---|
| `ring_size` | Size of the **new** ring (atoms = bonds). `int` for a single size, `(min, max)` for a window, `None` for no constraint. Translated per anchor pair into a `dist2` filter as `ring_size − d_in`, where `d_in` is the topological distance between the two anchor atoms in the input molecule. |
| `ring_closures` | Which provenance to query (see above). Default `True`. |
| `spiro` | Close a ring through a single ring atom: `False` (default), `True` (in addition), `'only'` (see above). |
| `min_atoms` / `max_atoms` | Heavy-atom size window of the linking fragment. Defaults `1` / `10`. |
| `replace_ids` / `protected_ids` | Restrict which atoms may serve as ring-closure anchors. |
| `set_names` / `min_freq` | Fragment set and frequency threshold. |
| `discard_ring_geometry` | Discard products whose new ring is geometrically impossible. Default `True` (see [Ring geometry filter](ring-geometry.md)). |

`symmetry_fixes` is accepted for API compatibility with the other generation
functions but is not used by `make_cycle`.

## Geometry of the new ring

Products whose new ring cannot exist in 3D — a six-membered ring bridging the
*meta* positions of a benzene ring, a short bridge across a naphthalene or a
biaryl — are discarded automatically. Set `discard_ring_geometry=False` to keep
them; see [Ring geometry filter](ring-geometry.md).

## Parallel use

For `multiprocessing`, use the list-returning wrapper `make_cycle2` — see
[Multiprocessing](multiprocessing.md).
