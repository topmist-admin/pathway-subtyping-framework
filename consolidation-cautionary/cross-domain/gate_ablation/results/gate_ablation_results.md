# Gate Ablation Study — results (R3.10)

n=60 samples, p=50 pathways, 15 reps/condition, n_ref=120, n_bootstrap=40. Elapsed 4623.9s.

**TPR** = fraction of genuine-subtype datasets certified (higher is better). **FPR** = fraction of no-structure datasets (single Gaussian + continuum) falsely certified as subtypes (lower is better).

| Gate subset | TPR | FPR | continuum_1d cert rate |
|---|---|---|---|
| `stability_only` | 1.00 | 0.37 | 0.73 |
| `discreteness_only` | 0.97 | 0.00 | 0.00 |
| `stability+discreteness` | 0.97 | 0.00 | 0.00 |

**Reading:** `stability_only` (the pre-v0.8.0 gate set) false-certifies the 1-D continuum — the erratum finding. `stability+discreteness` (v0.8.0) drives the continuum false-positive rate down while retaining true positives, which is the quantitative answer to R3.10: removing the discreteness gate reintroduces continuum false positives.