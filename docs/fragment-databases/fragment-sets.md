# Fragment sets

A v1 or v2 database can hold several **fragment sets** in one file. Each set is a
separate frequency column on every `radiusN` table, so the same deduplicated
`envs` and `frags` tables are shared while each set records how often a fragment
occurs *within that set*.

This lets one database describe, for example, how common a fragment is in
ChEMBL versus in a focused in-house library, and lets you switch between those
views at generation time.

!!! note "Not available in v0"
    Fragment sets need the v1/v2 schema. v0 databases have a single `freq` column and
    ignore the `set_names` argument.

## Build with a single set

```bash
cremdb_create -i input.smi -o fragments.db -s chembl
```

Every fragment occurrence is counted in the `chembl` column.

## Build with set-membership files

Pass one or more membership files to `-s/--set-name`. Each existing file is
treated as a membership list, and its basename (without extension) becomes the
set column name:

```bash
cremdb_create -i input_with_ids.smi -o fragments.db \
  -s set1_ids.txt set2_ids.txt all_set
```

- `set1_ids.txt` → column `set1_ids`
- `set2_ids.txt` → column `set2_ids`
- `all_set` is not an existing file, so it becomes a default set containing
  **all** molecules.

A fragment occurrence contributes to a set's count only when the source
molecule's ID belongs to that set.

### Membership file format

One molecule ID per line:

```text
mol_0001
mol_0007
mol_0042
```

For membership to work, the input SMILES file must carry IDs in the second
column:

```text
CCO         mol_0001
c1ccccc1    mol_0002
```

## Inspect the sets in a database

```bash
cremdb_info -i fragments.db
```

prints the schema version and the set columns per radius table, e.g.:

```text
fragments.db
  schema version : 2 (current)
  radius1        : chembl
  radius2        : chembl
  radius3        : chembl
  properties     : (none)
```

Several databases may be given at once, and `--json` prints the same information in a
machine-readable form. (The older `cremdb_get_set_names` still works but is deprecated.)

Equivalently, with SQLite:

```sql
PRAGMA table_info(radius3);
```

Any column other than `env_id`, `core_smi_id`, `core_num_atoms`, `dist2`, and
`is_ring_closure` is a fragment-set frequency column.

## Use a set at generation time

Choose the set(s) with `set_names` and the threshold with `min_freq`:

```python
from rdkit import Chem
from crem.crem import mutate_mol

m = Chem.MolFromSmiles("c1ccccc1N")

# Fragments seen at least 5 times in set1_ids
res = list(mutate_mol(m, db_name="fragments.db", set_names="set1_ids", min_freq=5))
```

`set_names` accepts a single column name or a list. When several sets are
named, a fragment is included if **any** of them meets `min_freq` (OR logic).
With `set_names=None` (default), all set columns are considered. Naming a column
that does not exist raises a `ValueError` listing the available set names.

```python
# Frequent in either set
res = list(mutate_mol(m, db_name="fragments.db",
                      set_names=["set1_ids", "set2_ids"], min_freq=3))
```
