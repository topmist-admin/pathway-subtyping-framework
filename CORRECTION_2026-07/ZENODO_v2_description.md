# Zenodo new-version description (paste into the deposit)

**Version 2 — Corrected (2026-07-08).** This version corrects errors in the original release.
The adaptive bootstrap-threshold model previously reported at R^2=0.889 does not reproduce and
is retracted (refit R^2~0.11 on released data; ~0.015 on valid rows). Fourteen benchmark rows
carried an empty-input adjusted-Rand-index artifact (ARI=1.0 on degenerate ground truth) and are
flagged invalid; one row had impossible labels; the claim that the benchmark excludes the
manuscript analysis dataset (TCGA-COAD) was incorrect and is withdrawn; dataset counts and sample
totals are restated (41 unique datasets, not 47; 40,778 samples in passing rows, not 36,551). See ERRATUM_2026-07-08.md.
The corrected benchmark shows that discrete-subtype reproducibility is dataset-dependent and low on
most independent cohorts — a negative methodological result; it does not support a silhouette-
calibrated threshold. Superseded version retained per Zenodo DOI permanence.
