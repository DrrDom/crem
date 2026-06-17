# Installation

## Requirements

- Python ≥ 3.8
- [RDKit](https://www.rdkit.org/) ≥ 2025.3.5
- `numpy`, `joblib`, `tqdm` (installed automatically)

!!! note "RDKit version and database compatibility"
    Context canonicalization relies on RDKit's SMILES output. A fragment
    database built with one RDKit version is only guaranteed to be readable by
    the same (or a compatible) RDKit version. Up to now we did not observe 
    incompatibilities. However, we recommend to pin RDKit if you share databases
    across machines.

!!! note "RDKit version"
    The previous RDKit versions may also be used. However, this may result in 
    a loss of stereoconfigurations in generated structures.   

Optional extras:

- `guacamol` — only needed to run the [GuacaMol benchmark](../about/benchmarks.md)
  (`guacamol_test`).
- `zstandard` — only needed to read `.zst`-compressed input when building
  databases.

## Install from PyPI

```bash
pip install crem
```

Installing the package makes the `crem` module importable in Python and adds
the [command-line tools](../reference/cli.md) (`cremdb_create`, `cremdb_convert`,
…) to your `PATH`.

## Install from source

```bash
git clone https://github.com/DrrDom/crem
cd crem
pip install .
```

To build distribution artifacts (PEP 517):

```bash
pip install build
python -m build          # writes wheel + sdist to dist/
```

## Uninstall

```bash
pip uninstall crem
```

## Build the documentation locally

The documentation is built with [MkDocs](https://www.mkdocs.org/) and the
Material theme. The API reference is generated from docstrings by
[mkdocstrings](https://mkdocstrings.github.io/), which imports `crem` — so the
package (and RDKit) must be installed in the same environment.

```bash
pip install -r docs/requirements.txt
pip install .                     # so the API pages can import crem

mkdocs serve            # live preview at http://127.0.0.1:8000
mkdocs build            # static site in site/
```
