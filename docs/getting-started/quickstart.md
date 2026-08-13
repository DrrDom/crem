# Quick start

This page assumes you already have a fragment database. You can
[build one](../fragment-databases/build-v2.md) or download a precompiled ChEMBL
database from <http://www.qsar4u.com/pages/crem.php>. All examples use a file
named `fragments.db`.

## Import the API

```python
from rdkit import Chem
from crem.crem import mutate_mol, grow_mol, link_mols, make_cycle
```

## Mutate a molecule

Replace one fragment of the molecule with alternatives from the database. The
functions are **generators**, so wrap them in `list(...)` to materialize the
results.

```python
m = Chem.MolFromSmiles("c1cc(OC)ccc1C")  # methoxytoluene
mols = list(mutate_mol(m, db_name="fragments.db", max_size=1))
print(mols[:10])
```

To also replace hydrogens, pass an H-expanded molecule:

```python
mols = list(mutate_mol(Chem.AddHs(m), db_name="fragments.db", max_size=1))
```

## Grow a molecule

Replace hydrogens with fragments (scaffold decoration). Do **not** add
hydrogens explicitly — `grow_mol` does it internally.

```python
m = Chem.MolFromSmiles("c1cc(OC)ccc1C")
mols = list(grow_mol(m, db_name="fragments.db", radius=3, max_atoms=2))
print(mols[:10])
```

## Link two molecules

Connect two molecules with a linker fragment.

```python
m1 = Chem.MolFromSmiles("c1cc(OC)ccc1C")
m2 = Chem.MolFromSmiles("NCC(=O)O")  # glycine
mols = list(link_mols(m1, m2, db_name="fragments.db", radius=3, max_atoms=3))
print(mols[:10])
```

## Form a new ring

Close a ring between two atoms of the same molecule.

```python
m = Chem.MolFromSmiles("c1ccccc1N")
mols = list(make_cycle(m, db_name="fragments.db", ring_size=(5, 7), max_atoms=10))
print(mols[:10])
```

See [Make cycle](../operations/make-cycle.md) for the database requirements of
the two cyclization modes.

## Select fragments by set and frequency

In a database with named [fragment sets](../fragment-databases/fragment-sets.md),
restrict replacements to fragments frequent in a given set:

```python
mols = list(mutate_mol(m, db_name="fragments.db", set_names="chembl", min_freq=5))
```

## Next steps

- [Concepts](../concepts.md) — what radius, context, and sets mean.
- [Mutate, grow, link](../operations/mutate-grow-link.md) — all parameters.
- [Advanced fragment selection](../operations/advanced-selection.md) — custom
  filters and biased sampling.
