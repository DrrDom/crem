# Iterative enumeration

`crem.utils.enumerate_compounds` applies generation **repeatedly**: molecules
produced in one iteration become the inputs to the next. It has two modes:

- `mode="scaffold"` — scaffold decoration, built on `grow_mol` (hydrogens are
  replaced with fragments).
- `mode="analogs"` — analog enumeration, built on `mutate_mol`.

It returns the distinct molecules generated across all iterations (as `Mol`
objects, or SMILES with `return_smi=True`).

!!! warning "Combinatorial explosion"
    The number of products grows quickly with `n_iterations`. Restrict the
    editable region (`replace_ids` / `protected_ids`) and/or cap output with
    `max_replacements`.

## Scaffold decoration

Decorate specific positions of a scaffold over several rounds. In `scaffold`
mode, `protect_added_frag` is forced to `True`, so fragments are never attached
to fragments added on a previous iteration — growth happens only at the chosen
positions.

```python
from rdkit import Chem
from crem.utils import enumerate_compounds

mol = Chem.MolFromSmiles("Cc1cccc(Cl)c1")  # 1-chloro-3-methylbenzene

mols = enumerate_compounds(
    mol,
    db_fname="fragments.db",
    mode="scaffold",
    n_iterations=3,
    radius=3,
    max_replacements=2,
    replace_ids=[2, 4, 6],     # only these atoms are decorated
    return_smi=True,
)
print(mols)
```

## Analog enumeration

In `analogs` mode the molecule is mutated. Supply explicit hydrogens
(`Chem.AddHs`) if you also want hydrogens replaced. Set `protect_added_frag=True`
to keep newly added fragments from being mutated further.

```python
mol = Chem.MolFromSmiles("c1ccccc1C")
analogs = enumerate_compounds(
    mol,
    db_fname="fragments.db",
    mode="analogs",
    n_iterations=2,
    max_replacements=50,
    return_smi=True,
)
```

## Key parameters

| Parameter | Meaning |
|---|---|
| `mode` | `"scaffold"` (decorate) or `"analogs"` (mutate). Default `"scaffold"`. |
| `n_iterations` | Number of generation rounds. Default `1`. |
| `max_replacements` | Max replacements per molecule per iteration (random subset). Default `None` = all. |
| `replace_ids` / `protected_ids` | Restrict the editable atoms. |
| `protect_added_frag` | Prevent edits to previously added fragments. Forced `True` in `scaffold` mode. Default `False`. |
| `min_freq` | Minimum fragment frequency. Default `0`. |
| `return_smi` | Return SMILES instead of `Mol` objects. Default `False`. |
| `ncpu` | Number of cores. `None` = all CPUs. |

Extra keyword arguments are forwarded to `grow_mol` (scaffold mode) or
`mutate_mol` (analogs mode) — for example `min_atoms` / `max_atoms` for
decoration, or `min_size` / `max_size` / `min_inc` / `max_inc` /
`replace_cycles` for analogs. See [`crem.utils`](../reference/utils.md) for the
full signature.
