# Gate-Agnostic Demonstration — clusterer sweep (R2.2)

> ⚠️ **SUPERSEDED LABELLING — the "Gate A verdict" column collapses three outcomes
> into two.**
>
> Rows reading `REJECT (not-testable (no reproducible k))` are **abstentions, not
> rejections**: the gate found no reproducible k and declined to rule, rather than
> ruling the data continuous. Reading them as rejections is what produced the
> retracted claim that "Gate A rejects the 1-D continuum on 100% of runs."
>
> Likewise the `PASS` label means *certified*; the corrected script emits
> **CERTIFY / REJECT / ABSTAIN** so the three outcomes are distinguishable.
>
> The **script has been corrected**; re-running regenerates this file with the
> three-way labels and a split of non-certifications into explicit rejections vs
> abstentions. That re-run has not been performed. The per-run ARI values and
> verdict strings below are accurate — only the PASS/REJECT column is misleading.

Same synthetic data clustered by GMM, DEC (Xie 2016), and VAE-GMM (VaDE, Jiang 2017). **self-stability** = the method's own bootstrap reproducibility (the naive check each would pass). **Gate A** = the discreteness gate's verdict (clusterer-agnostic; tests the data).

| Condition | Clusterer | self-stability ARI | looks-reproducible? | Gate A verdict |
|---|---|---|---|---|
| discrete_k2 | `gmm` | 1.00 | yes | PASS (discrete structure) |
| discrete_k2 | `dec` | 0.20 | no | PASS (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.03 | no | PASS (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | PASS (discrete structure) |
| discrete_k2 | `dec` | 0.15 | no | PASS (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.02 | no | PASS (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | PASS (discrete structure) |
| discrete_k2 | `dec` | 0.04 | no | PASS (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.05 | no | PASS (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | PASS (discrete structure) |
| discrete_k2 | `dec` | 0.14 | no | PASS (discrete structure) |
| discrete_k2 | `vae_gmm` | -0.01 | no | PASS (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | PASS (discrete structure) |
| discrete_k2 | `dec` | 0.05 | no | PASS (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.04 | no | PASS (discrete structure) |
| continuum_1d | `gmm` | 0.86 | yes | REJECT (not-testable (no reproducible k)) |
| continuum_1d | `dec` | 0.47 | no | REJECT (not-testable (no reproducible k)) |
| continuum_1d | `vae_gmm` | 0.28 | no | REJECT (not-testable (no reproducible k)) |
| continuum_1d | `gmm` | 0.87 | yes | REJECT (no discrete structure (continuum)) |
| continuum_1d | `dec` | 0.30 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `vae_gmm` | 0.17 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `gmm` | 0.71 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `dec` | 0.37 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `vae_gmm` | 0.25 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `gmm` | 0.99 | yes | REJECT (no discrete structure (continuum)) |
| continuum_1d | `dec` | 0.17 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `vae_gmm` | 0.13 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `gmm` | 0.95 | yes | REJECT (not-testable (no reproducible k)) |
| continuum_1d | `dec` | 0.39 | no | REJECT (not-testable (no reproducible k)) |
| continuum_1d | `vae_gmm` | 0.26 | no | REJECT (not-testable (no reproducible k)) |

**Primary result (clusterer-agnostic rejection):** Gate A rejects the 1-D continuum for **100%** of runs, and does so identically whether the partition was drawn by GMM, DEC, or VAE-GMM. This is the substantive answer to R2.2: PSF is not another clustering method to be benchmarked against DEC/VAE — it is a validation layer that wraps any of them, because Gate A tests the discreteness of the *data* at a given k, not the confidence of the algorithm.

**Secondary observation (reported honestly):** on the continuum the clusterers differ in how reproducible their (spurious) bipartition looks — mean self-stability ARI GMM 0.88 vs DEC 0.34 vs VAE-GMM 0.22. GMM's near-0.8 self-stability is the classic continuum false-positive the gate was built to catch; the deep methods are additionally unstable under resampling at this small n, itself a caution against trusting algorithmic confidence (consistent with the small-n concerns in R3.5/R3.6/R3.9). Competitive DL *recovery* of real subtypes is a separate claim for the large cohorts (TCGA-COAD, CPTAC), not this small-n control.