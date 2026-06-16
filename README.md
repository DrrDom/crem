<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/crem-logo-white.png">
    <img src="docs/assets/crem-logo.png" alt="CReM" width="300">
  </picture>
</div>

# CReM — chemically reasonable mutations

[![PyPI version](https://img.shields.io/pypi/v/crem.svg)](https://pypi.org/project/crem/)
[![Documentation](https://img.shields.io/badge/docs-crem.readthedocs.io-0f766e.svg)](https://crem.readthedocs.io/)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE.txt)

**CReM** is an open-source Python framework to generate chemical structures using
a fragment-based approach.

The idea is similar to matched molecular pairs: fragments that occur in the same
context are considered interchangeable. CReM stores such context–fragment
relationships in a database and uses them to generate chemically valid
structures.

## Features

- **Four generation modes** — `mutate`, `grow`, `link`, and `make_cycle`
  (ring closure / macrocyclization).
- **Custom fragment databases** built in one step with `cremdb_create`, or
  downloaded as precompiled ChEMBL databases.
- **Multiple fragment sets per database** — switch between them at generation
  time with `set_names` and a frequency threshold (`min_freq`).
- **Fine control** — context radius, fragment-size windows, replaceable/protected
  atoms, and `replace_cycles` for partial-ring replacement.
- **Custom selection** — bias or restrict fragments with `filter_func` /
  `sample_func`, or with molecular-property columns.
- **Reproducible and parallel** — `seed` for deterministic sampling; `ncores`
  and picklable `*_mol2` wrappers for multiprocessing.

## Links

- 📖 Documentation: <https://crem.readthedocs.io/>
- 🧪 Web app: <https://crem.imtm.cz/>
- ⬇️ Precompiled ChEMBL databases: <http://www.qsar4u.com/pages/crem.php>
- 📝 Changelog: [changelog](changelog)

## Installation

```bash
pip install crem
```

From source:

```bash
git clone https://github.com/DrrDom/crem
cd crem
pip install .
```

CReM requires `rdkit>=2025.3.5`. Optional extras: `guacamol` (to run the
benchmark) and `zstandard` (to read `.zst`-compressed input when building
databases).

## Quick start

All examples assume a fragment database `fragments.db` — [build one](#build-a-fragment-database)
or download a precompiled ChEMBL database.

```python
from rdkit import Chem
from crem.crem import mutate_mol, grow_mol, link_mols, make_cycle

m = Chem.MolFromSmiles('c1cc(OC)ccc1C')          # methoxytoluene

# replace an existing fragment
mutants = list(mutate_mol(m, db_name='fragments.db', max_size=1))

# decorate by replacing a hydrogen
grown = list(grow_mol(m, db_name='fragments.db'))

# link two molecules with a linker
m2 = Chem.MolFromSmiles('NCC(=O)O')              # glycine
linked = list(link_mols(m, m2, db_name='fragments.db'))

# form a new ring
cyclic = list(make_cycle(m, db_name='fragments.db', ring_size=(5, 7)))
```

All four are generators (wrap in `list(...)`) and share many options — `radius`,
size windows, `min_freq` / `set_names`, `replace_ids` / `protected_ids`,
`filter_func` / `sample_func`, `max_replacements`, `seed`, and `ncores`. See
[Mutate, grow, link](https://crem.readthedocs.io/en/latest/operations/mutate-grow-link/),
[Advanced fragment selection](https://crem.readthedocs.io/en/latest/operations/advanced-selection/),
and the [API reference](https://crem.readthedocs.io/en/latest/reference/crem/).

## Build a fragment database

Build a database directly from a SMILES file in one step:

```bash
cremdb_create -i input.smi -o fragments.db -s chembl
```

This produces the current database format with fragment-set support and
ring-closure fragments. For multiple sets, property columns, sharded/parallel
builds, conversion of older databases, and the programmatic `crem.db` API, see
[Fragment databases](https://crem.readthedocs.io/en/latest/fragment-databases/build-v1/).

## Benchmarks

GuacaMol goal-directed benchmark (scores marked `*` are from the original
GuacaMol publication):

|task|SMILES LSTM*|SMILES GA*|Graph GA*|Graph MCTS*|CReM
|---|:---:|:---:|:---:|:---:|:---:|
|Celecoxib rediscovery|**1.000**|0.732|**1.000**|0.355|**1.000**
|Troglitazone rediscovery|**1.000**|0.515|**1.000**|0.311|**1.000**
|Thiothixene rediscovery|**1.000**|0.598|**1.000**|0.311|**1.000**
|Aripiprazole similarity|**1.000**|0.834|**1.000**|0.380|**1.000**
|Albuterol similarity|**1.000**|0.907|**1.000**|0.749|**1.000**
|Mestranol similarity|**1.000**|0.79|**1.000**|0.402|**1.000**
|C11H24|**0.993**|0.829|0.971|0.410|0.966
|C9H10N2O2PF2Cl|0.879|0.889|**0.982**|0.631|0.940
|Median molecules 1|**0.438**|0.334|0.406|0.225|0.371
|Median molecules 2|0.422|0.38|0.432|0.170|**0.434**
|Osimertinib MPO|0.907|0.886|0.953|0.784|**0.995**
|Fexofenadine MPO|0.959|0.931|0.998|0.695|**1.000**
|Ranolazine MPO|0.855|0.881|0.92|0.616|**0.969**
|Perindopril MPO|0.808|0.661|0.792|0.385|**0.815**
|Amlodipine MPO|0.894|0.722|0.894|0.533|**0.902**
|Sitagliptin MPO|0.545|0.689|**0.891**|0.458|0.763
|Zaleplon MPO|0.669|0.413|0.754|0.488|**0.770**
|Valsartan SMARTS|0.978|0.552|0.990|0.04|**0.994**
|Deco Hop|0.996|0.970|**1.000**|0.590|**1.000**
|Scaffold Hop|0.998|0.885|**1.000**|0.478|**1.000**
|total score|17.341|14.398|17.983|9.011|17.919

## Limitations

- CReM builds structures only from fragments present in the database, so the
  ring systems that can appear depend on the database. `make_cycle` and
  `replace_cycles` form or replace rings using fragments observed in the
  database rather than inventing entirely new ring systems.
- Very large molecules are skipped in some workflows: a molecule with more than
  30 non-ring single bonds is not mutated, and one with more than 100 hydrogen
  atoms is not grown or linked.
- Context canonicalization relies on RDKit's SMILES output. A database is best
  used with the RDKit version it was built with (no incompatibilities observed
  so far); pin RDKit when sharing databases across machines.

## License

BSD-3-Clause. See [LICENSE.txt](LICENSE.txt).

## Citation

CReM: chemically reasonable mutations framework for structure generation
Pavel Polishchuk
*Journal of Cheminformatics* **2020**, 12, (1), 28
<https://doi.org/10.1186/s13321-020-00431-w>

Control of Synthetic Feasibility of Compounds Generated with CReM
Pavel Polishchuk
*Journal of Chemical Information and Modeling* **2020**, 60, 6074-6080
<https://dx.doi.org/10.1021/acs.jcim.0c00792>