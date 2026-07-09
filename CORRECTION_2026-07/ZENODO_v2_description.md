# Zenodo new-version description (paste into the deposit)

**Version 2 — Corrected (2026-07-08).** This version corrects errors in the original release.
The adaptive bootstrap-threshold model previously reported at R^2=0.889 does not reproduce and
is retracted (refit R^2~0.11 on released data; ~0.015 on valid rows). Fourteen benchmark rows
have degenerate ground truth (n_true_clusters=1) and are flagged invalid — the empty-input
adjusted-Rand-index artifact gave 13 of them a spurious ARI=1.0 (the 14th, GSE136196, returned 0.0); one row had impossible labels; the claim that the benchmark excludes the
manuscript analysis dataset (TCGA-COAD) was incorrect and is withdrawn; dataset counts and sample
totals are restated (41 unique datasets, not 47; 40,778 samples in passing rows, not 36,551). See ERRATUM_2026-07-08.md.
The corrected benchmark shows that discrete-subtype reproducibility is dataset-dependent and low on
most independent cohorts — a negative methodological result; it does not support a silhouette-
calibrated threshold. This version (v2.0) supersedes v1.0 (10.5281/zenodo.19323754) and v1.1 (10.5281/zenodo.19324360), both retained per Zenodo DOI permanence. Cite the concept DOI going forward.
