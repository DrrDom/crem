# Benchmarks

## GuacaMol

CReM was evaluated on the [GuacaMol](https://github.com/BenevolentAI/guacamol)
goal-directed benchmark suite. Scores marked `*` are taken from the original
GuacaMol publication; the `CReM` column is CReM's score.

| Task | SMILES LSTM* | SMILES GA* | Graph GA* | Graph MCTS* | CReM |
|---|:---:|:---:|:---:|:---:|:---:|
| Celecoxib rediscovery | **1.000** | 0.732 | **1.000** | 0.355 | **1.000** |
| Troglitazone rediscovery | **1.000** | 0.515 | **1.000** | 0.311 | **1.000** |
| Thiothixene rediscovery | **1.000** | 0.598 | **1.000** | 0.311 | **1.000** |
| Aripiprazole similarity | **1.000** | 0.834 | **1.000** | 0.380 | **1.000** |
| Albuterol similarity | **1.000** | 0.907 | **1.000** | 0.749 | **1.000** |
| Mestranol similarity | **1.000** | 0.79 | **1.000** | 0.402 | **1.000** |
| C11H24 | **0.993** | 0.829 | 0.971 | 0.410 | 0.966 |
| C9H10N2O2PF2Cl | 0.879 | 0.889 | **0.982** | 0.631 | 0.940 |
| Median molecules 1 | **0.438** | 0.334 | 0.406 | 0.225 | 0.371 |
| Median molecules 2 | 0.422 | 0.38 | 0.432 | 0.170 | **0.434** |
| Osimertinib MPO | 0.907 | 0.886 | 0.953 | 0.784 | **0.995** |
| Fexofenadine MPO | 0.959 | 0.931 | 0.998 | 0.695 | **1.000** |
| Ranolazine MPO | 0.855 | 0.881 | 0.92 | 0.616 | **0.969** |
| Perindopril MPO | 0.808 | 0.661 | 0.792 | 0.385 | **0.815** |
| Amlodipine MPO | 0.894 | 0.722 | 0.894 | 0.533 | **0.902** |
| Sitagliptin MPO | 0.545 | 0.689 | **0.891** | 0.458 | 0.763 |
| Zaleplon MPO | 0.669 | 0.413 | 0.754 | 0.488 | **0.770** |
| Valsartan SMARTS | 0.978 | 0.552 | 0.990 | 0.04 | **0.994** |
| Deco Hop | 0.996 | 0.970 | **1.000** | 0.590 | **1.000** |
| Scaffold Hop | 0.998 | 0.885 | **1.000** | 0.478 | **1.000** |
| **Total score** | 17.341 | 14.398 | 17.983 | 9.011 | 17.919 |

## Reproducing

The benchmark is driven by the `guacamol_test` command, which requires the
optional `guacamol` package. See its options in the
[CLI reference](../reference/cli.md#guacamol_test).
