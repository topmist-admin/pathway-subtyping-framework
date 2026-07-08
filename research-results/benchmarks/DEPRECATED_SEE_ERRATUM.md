# ⚠️ DEPRECATED — see erratum (2026-07-08)

The benchmark file in this directory,
`bootstrap_threshold_calibration_47datasets_final_zenodo.csv` (and its variants), is the
**uncorrected** release. It contains:

- 14 rows with an empty-input `adjusted_rand_score` artifact (ARI=1.0 on degenerate ground
  truth, `n_true_clusters=1`) — spurious, not measured stability;
- 1 row with impossible labels (GSE92332: 1,957 "true clusters" for 533 samples);
- **TCGA-COAD present** despite the model file / manuscript stating the benchmark excludes it
  (the independence claim is withdrawn);
- count discrepancies (41 unique datasets, not 47; 40,778 samples in passing rows, not 36,551).

**Use the corrected file instead:** [`../../CORRECTION_2026-07/corrected_benchmark_47datasets_v2.csv`](../../CORRECTION_2026-07/corrected_benchmark_47datasets_v2.csv)
(all original rows retained, with `valid` / `invalid_reason` columns).

The adaptive-threshold model calibrated from this benchmark (R²=0.889) **does not reproduce**
and is **retracted**. Full notice: [`../../CORRECTION_2026-07/ERRATUM_2026-07-08.md`](../../CORRECTION_2026-07/ERRATUM_2026-07-08.md).
