# Gate Ablation Study — results (R3.10)

n=60 samples, p=50 pathways, 15 reps/condition, n_ref=100, n_bootstrap=40. Elapsed 7295.2s.

**TPR** = fraction of genuine-subtype datasets certified (higher is better). **FPR** = fraction of no-structure datasets (single Gaussian + continuum) falsely certified as subtypes (lower is better).

| Gate subset | TPR | FPR | continuum_1d cert rate |
|---|---|---|---|
| `stability_only` | 1.00 | 0.43 | 0.87 |
| `discreteness_only` | 1.00 | 0.00 | 0.00 |
| `stability+discreteness` | 1.00 | 0.00 | 0.00 |

**Reading:** `stability_only` (the pre-v0.8.0 gate set) false-certifies the 1-D continuum — the erratum finding. `stability+discreteness` (v0.8.0) removes those false certifications, which is the quantitative answer to R3.10: removing the discreteness gate reintroduces continuum false positives.

**Two caveats that must travel with the FPR number.** (1) The gate reaches its low FPR mostly by ABSTAINING rather than rejecting — on the negative controls it returns "not-testable (no reproducible k)" on 28/30 (93%), leaving a testable-negative denominator of 2. An FPR quoted against n=30 is the over-generous convention; against the testable n=2 the interval is [0, 0.658] and is nearly uninformative. Do not quote FPR=0.000 unqualified. (2) The gate is not free: TPR moves 1.000 -> 0.967 (one `discrete_k3` replicate lost). Exact McNemar on positives is p=1.0, so that cost is not distinguishable from zero — but the study has no power to exclude one. Report the point change and the p-value together; write neither "TPR held" nor a settled "3% cost".