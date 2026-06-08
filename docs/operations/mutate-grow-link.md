# Mutate, grow, and link

The three core generation functions share most of their parameters. All are
**generators** of distinct SMILES (wrap in `list(...)`), and all can optionally
return the reaction, its frequency, and/or the RDKit `Mol`.

| Function | What it does | Incoming fragment |
|---|---|---|
| `mutate_mol` | replaces an existing fragment | 1–4 attachment points |
| `grow_mol` | replaces a hydrogen | 1 attachment point |
| `link_mols` | joins two molecules | 2 attachment points (linker) |

Common parameters: `radius`, `min_freq` / `set_names`, `replace_ids` /
`protected_ids`, `max_replacements`, `filter_func` / `sample_func`, `ncores`,
property filters via `**kwargs`. See [`crem.crem`](../reference/crem.md) for the
complete, authoritative parameter list.

## Mutate

```python
from rdkit import Chem
from crem.crem import mutate_mol

m = Chem.MolFromSmiles("c1cc(OC)ccc1C")
res = list(mutate_mol(m, db_name="fragments.db", radius=3, max_size=1, max_inc=3))
```

To also replace hydrogens, pass an H-expanded molecule:

```python
res = list(mutate_mol(Chem.AddHs(m), db_name="fragments.db", max_size=1))
```

### Replacing cyclic source fragments

By default `mutate_mol` only cuts acyclic bonds, so ring systems are left
intact. The `replace_cycles` argument changes this:

| `replace_cycles` | Behaviour |
|---|---|
| `"no"` *(default)* | ordinary acyclic-cut mutation only |
| `"forced"` | allow cyclic cores from ordinary fragmentation to be replaced, ignoring the size filters |
| `"partial_all"` | additionally replace partial ring arcs using exhaustive side cuts |
| `"partial_exo"` | additionally replace partial ring arcs using only *exo* side cuts adjacent to the arc |

`partial_all` enumerates exhaustive side cuts, so it needs a database built with
the exhaustive `--frag-mode ring` or `both`; on an *optimal* database it
under-matches because the non-exo cuts are absent. `partial_exo` enumerates only
the *exo* side cuts — a subset — so it works with **any** ring-capable database
(`ring_optimal`, `both_optimal`, `ring`, or `both`); it is faster and narrower
and may return fewer products. `"no"` and `"forced"` use ordinary acyclic rows
and work with any database.

```python
res = list(mutate_mol(
    m,
    db_name="fragments.db",
    radius=3,
    max_size=8,
    replace_cycles="partial_exo",
))
```

For forming *new* rings rather than swapping existing ones, see
[Make cycle](make-cycle.md).

## Grow

`grow_mol` adds hydrogens internally and replaces them — do **not** call
`Chem.AddHs` yourself.

```python
from rdkit import Chem
from crem.crem import grow_mol

m = Chem.MolFromSmiles("c1cc(OC)ccc1C")
res = list(grow_mol(m, db_name="fragments.db", radius=3, min_atoms=1, max_atoms=2))
```

## Link

```python
from rdkit import Chem
from crem.crem import link_mols

m1 = Chem.MolFromSmiles("c1cc(OC)ccc1C")
m2 = Chem.MolFromSmiles("NCC(=O)O")
res = list(link_mols(m1, m2, db_name="fragments.db", radius=3, min_atoms=1, max_atoms=3))
```

Constrain the linker geometry with `dist` — the topological distance between the
two attachment points (a single value or a `(low, high)` tuple):

```python
res = list(link_mols(m1, m2, db_name="fragments.db", radius=3,
                     dist=(2, 6), min_atoms=1, max_atoms=4))
```

## Restricting where changes happen

- `replace_ids` — only these atoms (and their hydrogens) may be modified.
- `protected_ids` — these atoms are never modified. `protected_ids` has higher
  priority than `replace_ids`.

When protecting positions that have symmetry-equivalent atoms, supply the ids of
**all** equivalent atoms (for example both meta carbons in toluene). For
hydrogen replacement, supply hydrogen ids only when the molecule was created
with explicit hydrogens.

`link_mols` takes per-molecule variants: `replace_ids_1` / `replace_ids_2` and
`protected_ids_1` / `protected_ids_2`.

## Returning transformations and frequencies

```python
res = list(mutate_mol(
    m,
    db_name="fragments.db",
    set_names="chembl",
    min_freq=10,
    return_rxn=True,
    return_rxn_freq=True,
))
# each item: [smiles, rxn, freq]
```

The optional return values are appended in this order: SMILES, then `rxn`
(`return_rxn`), then `freq` (`return_rxn_freq`, only alongside `return_rxn`),
then the `Mol` (`return_mol`). With no optional returns, the generator yields
plain SMILES strings.

## Limiting and reproducing output

- `max_replacements=N` returns at most `N` products, sampled uniformly at random
  from the available replacements (use `sample_func` to bias the sampling — see
  [Advanced fragment selection](advanced-selection.md)).
- `seed=...` makes that random selection reproducible.
