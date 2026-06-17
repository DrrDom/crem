# Two-step replacements

`mutate_mol`, `grow_mol`, and `link_mols` do two things at once: they query the
database for matching fragments and then assemble the products. CReM also
exposes these two steps separately:

- `get_replacements` — query the database and return the matching fragments.
- `get_mols_from_replacements` — assemble products from a set of replacements.

This is useful when you want to **inspect, cache, or filter** the candidate
fragments before building molecules, or build products from the same
replacement set more than once.

## Inspect matching fragments

By default `get_replacements` yields the SMILES of the candidate replacement
fragments:

```python
from rdkit import Chem
from crem.crem import get_replacements

m = Chem.MolFromSmiles("c1cc(OC)ccc1C")
frags = list(get_replacements(m, db_name="fragments.db", radius=3, max_size=1))
print(frags[:10])
```

## Get replacements, then build molecules

Set `return_frag_smi_only=False` to get full replacement tuples
`(source_core_smi, replacement_core_smi, freq, context_mol)`, then pass them to
`get_mols_from_replacements`:

```python
from rdkit import Chem
from crem.crem import get_replacements, get_mols_from_replacements

m = Chem.MolFromSmiles("c1cc(OC)ccc1C")
radius = 3

replacements = list(get_replacements(
    m,
    db_name="fragments.db",
    radius=radius,
    max_size=1,
    return_frag_smi_only=False,
))

# (optionally inspect / filter `replacements` here)

mols = list(get_mols_from_replacements(
    m,
    radius=radius,
    replacements=replacements,
    return_rxn=True,
))
# each item: [smiles, rxn]
```

The `radius` passed to `get_mols_from_replacements` must match the one used to
fetch the replacements.

## Linking

To search for linking fragments, pass a second molecule as `mol2` to
`get_replacements`, and the same `mol2` to `get_mols_from_replacements`:

```python
m1 = Chem.MolFromSmiles("c1cc(OC)ccc1C")
m2 = Chem.MolFromSmiles("NCC(=O)O")

replacements = list(get_replacements(
    m1, db_name="fragments.db", radius=3, mol2=m2,
    min_atoms=1, max_atoms=3, return_frag_smi_only=False,
))
mols = list(get_mols_from_replacements(m1, radius=3, replacements=replacements, mol2=m2))
```

See [`crem.crem`](../reference/crem.md) for the full signatures.
