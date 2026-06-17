# Multiprocessing

There are two independent ways to parallelize CReM.

## Parallelism within one molecule: `ncores`

Every generation function takes an `ncores` argument. It parallelizes the
replacements made *within a single molecule* across that many processes:

```python
from rdkit import Chem
from crem.crem import mutate_mol

m = Chem.MolFromSmiles("c1cc(OC)ccc1C")
mols = list(mutate_mol(m, db_name="fragments.db", max_size=1, ncores=8))
```

## Parallelism across many molecules: the `*2` wrappers

`mutate_mol`, `grow_mol`, `link_mols`, and `make_cycle` are **generators**.
Generators cannot be pickled, so they can't be passed directly to
`multiprocessing.Pool`. For that case CReM provides list-returning wrappers that
take exactly the same arguments and are picklable:

| Generator | List-returning wrapper |
|---|---|
| `mutate_mol` | `mutate_mol2` |
| `grow_mol` | `grow_mol2` |
| `link_mols` | `link_mols2` |
| `make_cycle` | `make_cycle2` |

Use them with `multiprocessing.Pool` to process several molecules at once:

```python
from multiprocessing import Pool
from functools import partial
from rdkit import Chem
from crem.crem import mutate_mol2

input_smi = ["c1ccccc1N", "NCC(=O)OC", "NCCCO"]
input_mols = [Chem.MolFromSmiles(s) for s in input_smi]

with Pool(2) as p:
    res = list(p.imap(
        partial(mutate_mol2, db_name="fragments.db", max_size=1),
        input_mols,
    ))
# res is a list (one per input molecule) of lists of product SMILES
```

When parallelizing across molecules this way, keep `ncores=1` (the default) in
each wrapper call so you don't oversubscribe the CPUs — let the outer `Pool`
own the parallelism.

!!! tip
    `crem.utils.enumerate_compounds` already parallelizes across molecules
    internally (via `joblib`), so you don't need to wrap it yourself — see
    [Iterative enumeration](enumeration.md).
