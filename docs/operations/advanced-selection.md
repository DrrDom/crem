# Advanced fragment selection

Beyond `radius`, size windows, `min_freq`, and `set_names`, three mechanisms
give fine-grained control over which fragments are used:

1. `filter_func` — arbitrary on-the-fly filtering of candidate fragments.
2. `sample_func` — biased on-the-fly sampling when `max_replacements` is set.
3. `**kwargs` — filtering on [property columns](../fragment-databases/properties.md)
   in the database.

All three are accepted by `mutate_mol`, `grow_mol`, `link_mols`, and
`make_cycle`.

## 1. Custom filtering with `filter_func`

A `filter_func` receives the candidate database row ids and returns the subset
to keep. Its first three arguments are fixed — `row_ids`, the database cursor
`cur`, and the context `radius` — and any further arguments are your own (bind
them with `functools.partial`). Use the helper `crem.crem._get_replacements` to
turn row ids into fragment SMILES.

```python
from collections import defaultdict
from functools import partial
from rdkit import Chem
from crem.crem import mutate_mol, _get_replacements

def filter_by_atom(row_ids, cur, radius, atom_number):
    """Keep only fragments containing an atom of the given atomic number."""
    if not row_ids:
        return []
    by_smi = defaultdict(list)
    for rowid, core_smi, _, _ in _get_replacements(cur, radius, row_ids):
        by_smi[core_smi].append(rowid)
    out = []
    for smi, ids in by_smi.items():
        mol = Chem.MolFromSmiles(smi)
        if mol and any(a.GetAtomicNum() == atom_number for a in mol.GetAtoms()):
            out.extend(ids)
    return out

# only fluorine-containing fragments will be used
mols = list(mutate_mol(
    Chem.MolFromSmiles("c1ccccc1C"),
    db_name="fragments.db",
    filter_func=partial(filter_by_atom, atom_number=9),
    max_size=1,
    max_inc=3,
))
```

### Built-in filters

`crem.utils` ships ready-made filters:

- `filter_max_ring_size(row_ids, cur, radius, max_size=6)` — drop fragments with
  a ring larger than `max_size`.
- `filter_acyclic_attachment_points(row_ids, cur, radius)` — keep only fragments
  whose attachment points sit on acyclic atoms (handy for
  [`make_cycle`](make-cycle.md)).

```python
from functools import partial
from crem.utils import filter_max_ring_size

mols = list(mutate_mol(
    Chem.MolFromSmiles("c1ccccc1C"),
    db_name="fragments.db",
    filter_func=partial(filter_max_ring_size, max_size=6),
))
```

## 2. Biased sampling with `sample_func`

When `max_replacements` is set, CReM samples that many replacements uniformly at
random. A `sample_func` replaces the uniform draw with your own selection. Its
first four arguments are fixed — `row_ids`, `cur`, `radius`, and `n` (the number
to return) — and it returns the selected row ids.

The built-in `crem.utils.sample_csp3` biases selection toward fragments with a
higher fraction of sp³ carbons:

```python
from crem.utils import sample_csp3

mols = list(mutate_mol(
    Chem.MolFromSmiles("c1ccccc1F"),
    db_name="fragments.db",
    max_replacements=10,
    sample_func=sample_csp3,
))
```

!!! note
    Complex sampling functions run for every replacement site and can slow down
    generation noticeably.

## 3. Property filters with `**kwargs`

After [adding property columns](../fragment-databases/properties.md), filter on
them by passing ranges as keyword arguments. A value is either an exact number
or an inclusive `(low, high)` tuple:

```python
mols = list(mutate_mol(
    Chem.MolFromSmiles("c1ccccc1N"),
    db_name="fragments.db",
    set_names="chembl",
    min_freq=5,
    mw=(50, 180),
    logp=(0.0, 3.5),
))
```

Standard properties (`mw`, `logp`, `rtb`, `tpsa`, `fcsp3`) are added with:

```bash
cremdb_add_prop -i fragments.db -p mw logp rtb tpsa fcsp3 -c 8
```

The keyword names must match existing columns: `frags` / `frags_h` columns in a
v1 database, or `radiusN` columns in a v0 database.
