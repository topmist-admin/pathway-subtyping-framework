# Clusterer sweep — ⛔ SUPERSEDED framing, retained for the record — clusterer sweep (R2.2)

> ⛔ **NON-REPRODUCIBLE ARTIFACT — do not cite any number below.** `clusterer_sweep_raw.csv`
> was produced before the salted-`hash()` seeding defect was fixed, and cannot be regenerated
> from the current code: recomputing `testable` for `continuum_1d` from the deposited
> generator disagrees with this file on 4 of 5 reps (the ablation CSV, by contrast, matches
> on 5 of 5). Its derived figures — the 60%/40% reject-vs-abstain split and the 0.88/0.34/0.22
> self-stability means — are therefore unverifiable, and the 40% abstention rate contradicts
> Result 1's 93% on the same generator. The manuscript withdraws this sweep as a result
> (Result 1, clusterer-invariance retraction); it is retained only as a record.


Same synthetic data clustered by GMM, DEC (Xie 2016), and VAE-GMM (VaDE, Jiang 2017). **self-stability** = the method's own bootstrap reproducibility (the naive check each would pass). **Gate A** = the discreteness gate's verdict (re-clusters internally; takes no partition argument, so the three rows are copies of one call).

| Condition | Clusterer | self-stability ARI | looks-reproducible? | Gate A verdict |
|---|---|---|---|---|
| discrete_k2 | `gmm` | 1.00 | yes | CERTIFY (discrete structure) |
| discrete_k2 | `dec` | 0.20 | no | CERTIFY (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.03 | no | CERTIFY (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | CERTIFY (discrete structure) |
| discrete_k2 | `dec` | 0.15 | no | CERTIFY (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.02 | no | CERTIFY (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | CERTIFY (discrete structure) |
| discrete_k2 | `dec` | 0.04 | no | CERTIFY (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.05 | no | CERTIFY (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | CERTIFY (discrete structure) |
| discrete_k2 | `dec` | 0.14 | no | CERTIFY (discrete structure) |
| discrete_k2 | `vae_gmm` | -0.01 | no | CERTIFY (discrete structure) |
| discrete_k2 | `gmm` | 1.00 | yes | CERTIFY (discrete structure) |
| discrete_k2 | `dec` | 0.05 | no | CERTIFY (discrete structure) |
| discrete_k2 | `vae_gmm` | 0.04 | no | CERTIFY (discrete structure) |
| continuum_1d | `gmm` | 0.86 | yes | ABSTAIN (not-testable (no reproducible k)) |
| continuum_1d | `dec` | 0.47 | no | ABSTAIN (not-testable (no reproducible k)) |
| continuum_1d | `vae_gmm` | 0.28 | no | ABSTAIN (not-testable (no reproducible k)) |
| continuum_1d | `gmm` | 0.87 | yes | REJECT (no discrete structure (continuum)) |
| continuum_1d | `dec` | 0.30 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `vae_gmm` | 0.17 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `gmm` | 0.71 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `dec` | 0.37 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `vae_gmm` | 0.25 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `gmm` | 0.99 | yes | REJECT (no discrete structure (continuum)) |
| continuum_1d | `dec` | 0.17 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `vae_gmm` | 0.13 | no | REJECT (no discrete structure (continuum)) |
| continuum_1d | `gmm` | 0.95 | yes | ABSTAIN (not-testable (no reproducible k)) |
| continuum_1d | `dec` | 0.39 | no | ABSTAIN (not-testable (no reproducible k)) |
| continuum_1d | `vae_gmm` | 0.26 | no | ABSTAIN (not-testable (no reproducible k)) |

**⛔ SUPERSEDED FRAMING — do not cite.** Gate A declines to certify the 1-D continuum on **100%** of runs. The arms below are IDENTICAL BY CONSTRUCTION, not by agreement: the gate is called once per condition and its verdict copied into all three clusterer rows. This shows only that the gate's INTERFACE accepts any clusterer. That total splits into **60% explicit rejection** and **40% abstention** ("not-testable — no reproducible k"). The abstentions must not be reported as rejections: in those runs the gate declined to rule rather than ruling against, and counting them as rejections is what produced the retracted "rejects on 100% of runs" claim. This is the substantive answer to R2.2: PSF is not another clustering method to be benchmarked against DEC/VAE — it is a validation layer that wraps any of them, because Gate A tests the discreteness of the *data* at a given k, not the confidence of the algorithm.

**Secondary observation (reported honestly):** on the continuum the clusterers differ in how reproducible their (spurious) bipartition looks — mean self-stability ARI GMM 0.88 vs DEC 0.34 vs VAE-GMM 0.22. GMM's near-0.8 self-stability is the classic continuum false-positive the gate was built to catch; the deep methods are additionally unstable under resampling at this small n, itself a caution against trusting algorithmic confidence (consistent with the small-n concerns in R3.5/R3.6/R3.9). Competitive DL *recovery* of real subtypes is a separate claim for the large cohorts (TCGA-COAD, CPTAC), not this small-n control.