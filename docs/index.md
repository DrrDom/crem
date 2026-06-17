# CReM — chemically reasonable mutations

[![GitHub repo](https://img.shields.io/badge/GitHub-DrrDom%2Fcrem-0f766e?logo=github&logoColor=white)](https://github.com/DrrDom/crem)
[![GitHub stars](https://img.shields.io/github/stars/DrrDom/crem?style=social)](https://github.com/DrrDom/crem/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/DrrDom/crem?style=social)](https://github.com/DrrDom/crem/network/members)

**CReM** is an open-source Python framework for generating chemical structures
using a fragment-based approach.

The core idea is borrowed from matched molecular pairs: fragments that occur in
the **same chemical context** are considered interchangeable. CReM stores such
context–fragment relationships in a SQLite database and uses them to grow,
mutate, link, and cyclize molecules so that every generated structure is built
only from fragment substitutions that have been observed in real molecules.

## What you can do

- **Build fragment databases** from your own molecular collections, or download
  precompiled ChEMBL databases.
- **Generate molecules** in four modes:
  [`mutate`](operations/mutate-grow-link.md),
  [`grow`](operations/mutate-grow-link.md),
  [`link`](operations/mutate-grow-link.md), and
  [`make_cycle`](operations/make-cycle.md).
- **Control the chemistry** through context radius, fragment-size windows, and
  per-set frequency thresholds.
- **Restrict where changes happen** with `replace_ids` / `protected_ids`.
- **Bias fragment selection** with custom `filter_func` and `sample_func`
  callbacks, or with property filters stored in the database.
- **Store several fragment sets in one database** and switch between them at
  generation time with `set_names` and `min_freq`.

## Where to start

| If you want to… | Read |
|---|---|
| Understand the vocabulary (context, radius, core, sets) | [Concepts](concepts.md) |
| Install CReM | [Installation](getting-started/installation.md) |
| Run your first generation | [Quick start](getting-started/quickstart.md) |
| Build a fragment database | [Build a database (v1)](fragment-databases/build-v1.md) |
| Look up a function or CLI flag | [Reference](reference/crem.md) |

## How CReM works in one picture

```text
            context (environment)            core fragment
        ┌───────────────────────────┐      ┌─────────────┐
   …—C—C—[*:1]                       and    [*:1]—CH3      ← stored together in the DB
        └───────────────────────────┘      └─────────────┘

   To modify a molecule, CReM:
     1. fragments it into (context, core) pairs;
     2. looks up the context in the database;
     3. retrieves every core seen in that same context;
     4. reassembles the molecule with each alternative core.
```

See [Concepts](concepts.md) for the full model and
[Database schema](fragment-databases/schema.md) for how this is stored.

## Online resources

- Documentation: <https://crem.readthedocs.io/>
- Web application: <https://crem.imtm.cz/>
- Precompiled ChEMBL databases: <http://www.qsar4u.com/pages/crem.php>
- Source code: <https://github.com/DrrDom/crem>

## License

BSD-3-Clause. See [Citation](about/citation.md) for the papers to cite.
