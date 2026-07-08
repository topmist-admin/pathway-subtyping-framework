# ⚠️ RETRACTED — `threshold_model_real47.json` (2026-07-08)

The adaptive bootstrap-threshold model in this directory
(`threshold_model_real47.json`, reported at **R²=0.889**, slope +0.914) **does not reproduce**
and is **retracted**. Refit on the released benchmark gives R²≈0.11 (slope reversed); on the
*valid* rows only, R²≈0.015. Silhouette does **not** predict bootstrap-ARI reproducibility in
these data, and the reported R²=0.889 does not reproduce from the released benchmark and was not
verified against it before release. The pipeline does **not** consume this file (the shipped stability threshold is
the sample-size × cluster-count grid in `threshold_calibration.py`).

Do not use this model. Corrected/retraction record:
[`RETRACTED replacement`](../../CORRECTION_2026-07/threshold_model_real47_CORRECTED.json) ·
[`ERRATUM`](../../CORRECTION_2026-07/ERRATUM_2026-07-08.md).
