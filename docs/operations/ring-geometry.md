# Ring geometry filter

Ring-forming operations can join two positions of a molecule that are held too
far apart, or at the wrong mutual orientation, for the new ring to exist — a
six-membered ring bridging the *meta* positions of a benzene ring being the
canonical example. Such products are removed by a geometry filter that runs on
every assembled product of a ring-forming transformation. It is **on by
default**; `discard_ring_geometry=False` restores the unfiltered output.

```python
# keep everything, including geometrically impossible rings
res = list(make_cycle(m, db_name="fragments.db", ring_size=(5, 7),
                      discard_ring_geometry=False))
```

## Scope

Only ring-forming transformations are examined. They are detected structurally
rather than declared by the caller: the inserted fragment is tracked by atom
labels, and a transformation is a ring closure when its two attachment bonds
land on atoms of the same connected parent fragment. Substitution, growing and
linking are therefore never filtered.

| Function | `discard_ring_geometry` |
|---|---|
| [`make_cycle`](make-cycle.md) / `make_cycle2` | accepted, default `True` |
| [`mutate_mol`](mutate-grow-link.md#replacing-cyclic-source-fragments) / `mutate_mol2` | accepted, default `True` — acts on the `replace_cycles="partial_*"` modes |
| [`get_mols_from_replacements`](two-step-replacements.md) | accepted, default `True` |
| `grow_mol` / `grow_mol2` | not accepted — growing replaces a hydrogen and forms no ring |
| `link_mols` | not accepted — linking joins two separate molecules |

Closures with three or more attachment points are not analysed and are always
kept.

## Design constraint

The filter is a **sound rejecter**: it discards only structures whose
non-embeddability is established, and accepts anything uncertain. Distance-geometry
embedding is incomplete — a failure to embed is not proof of impossibility — so
every threshold sits a deliberate margin inside the empirically measured
boundary, and any internal error in the analysis falls through to acceptance.

## Rule 1 — rigid-arc bands

For every new ring of 4–9 atoms the filter identifies the contiguous arc of an
aromatic ring that the new ring contains. The rule is keyed on the arc, i.e. on
the ring's contents, and not on the atoms the fragment is bonded to, so it
applies equally when the closure starts on a side chain.

| Arc inside the new ring | New ring rejected at |
|---|---|
| 3 atoms of a six-membered aromatic ring (*meta* span) | 4–8 |
| 4 atoms of a six-membered aromatic ring (*para* span) | 5–9 |
| 3 atoms of a five-membered aromatic ring (1,3-span across its middle atom) | 4–8 |

The five-membered 1,3-span is the exact analogue of benzene-*meta*: two
divergent exocyclic vectors. Adjacent (*ortho*) arcs, sp³ spans and larger rings
pass. The upper band edges stop one step below the smallest bridge known to
close in practice, so isolable strained systems such as [6]metacyclophane and
[6]paracyclophane are retained.

## Rule 2 — reach test

The complementary failure is one of pure distance: a new ring spanning a fused
aromatic system (naphthalene 2,6 ≈ 5.0 Å) or a biaryl (4,4′ ≈ 7–9 Å) with a
bridge too short to get there. The span is measured on a conformer of the
*isolated* ring system, which is rigid and therefore fully determined by one
embedding (for a biaryl the single torsion is scanned in 15° steps and the
minimum distance kept); results are cached per scaffold. It is compared with the
remainder of the cycle stretched perfectly straight — the sum of its maximum
bond lengths (covalent radii + 0.1 Å), an upper bound no conformation can
exceed. The product is discarded when the span exceeds this reach by more than a
10 % tolerance. No bond angle enters the comparison, so strained but real rings
cannot be rejected.

## Calibration and validation

Band edges were determined on bare templates in which a (CH₂)*ₙ* bridge of
increasing length spans the relevant ring positions, using multi-conformer
ETKDGv3 with random-coordinate restarts followed by MMFF optimisation. Because
embedding alone accepts plainly impossible molecules — puckering an "aromatic"
ring by 50–100° or stretching bonds to 1.9 Å — realisability was judged from the
worst aromatic-ring torsion and the strain energy per heavy atom (reference
values: naphthalene 0.0°, indane 0.6°, [6]metacyclophane 14.6°).

On a benchmark of 104 945 unique ring closures generated from ten drug-like
parents against a ChEMBL 36 fragment database, the filter removes 26.5 % of ring
closures and reduces the fraction of surviving closures with no valid 3D
structure from 28.8 % to 12.0 %; the residual lies almost entirely in the
cyclophane-like grey zone where the conformational assay itself is unreliable.
Of 1 200 randomly sampled rejected products none was a realistic structure, and
no product independently verified as realisable was discarded. The check costs
0.2–0.3 ms per ring-forming product and nothing for other transformations.
