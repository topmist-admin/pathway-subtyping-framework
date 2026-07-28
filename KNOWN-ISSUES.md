# Known Issues

> **⚠️ 2026-07-08 — Correction issued.** The empty-input `adjusted_rand_score` artifact documented below (NB17) was **also present, uncorrected, in the 47-dataset calibration benchmark** — it produced 14 degenerate rows (`n_true_clusters=1`), 13 of which carried a spurious ARI=1.0, that drove the retracted R²=0.889 adaptive-threshold model. The calibration model is **retracted** and the benchmark **corrected**. See [`CORRECTION_2026-07/ERRATUM_2026-07-08.md`](CORRECTION_2026-07/ERRATUM_2026-07-08.md).
>
> **⚠️ 2026-07 — Scope of this notice.** The correction above covers the ARI artifact
> and the retracted R²=0.889 model. It does **not** cover the **CMS4 recovery
> 75.9%/76%** figure quoted later in this file, which is also superseded as a headline
> validation claim — see the "Framing" section at the end, whose guidance has been
> **reversed**.
>
> **✅ 2026-07-09 — Fixed at the source.** The empty-input / single-true-cluster ARI artifact is now guarded framework-wide by `pathway_subtyping.utils.metrics.safe_adjusted_rand_score` / `ari_with_validity`, which return NaN (not a spurious score) on degenerate ground truth and are wired into `benchmark.py` and `pipeline.py`. The guard keys on ground-truth *structure* (`n_true_clusters < 2`), so it also catches the one degenerate row that returned ARI=0.0 (GSE136196), which a value-based (`== 1.0`) guard would miss. Regression tests: `tests/test_metrics_ari_guard.py`.

## NB17: MSI ARI reports 1.0 when no MSI data available

**Status:** FIXED — `valid_msi.sum() >= 2` guard added (Feb 25, 2026); results_summary.json was manually corrected to `null` earlier
**File:** `examples/notebooks/17_tcga_cancer_validation.ipynb`
**Severity:** Low (cosmetic — only affects NB17 output JSON and heatmap)

### Problem

The `msi_status` column exists in the TCGA-COAD clinical metadata (created during GDC API fetch) but contains no valid values — all are `None`. The MSI ARI computation block guards on `if msi_col_final:` which is truthy because the column name string `'msi_status'` exists, even though there are 0 valid samples.

`sklearn.metrics.adjusted_rand_score([], [])` returns `1.0` for empty inputs, producing a misleading perfect ARI.

### Location (approximate cell in Section 11)

```python
# MSI ARI computation
if msi_col_final:
    valid_msi = (
        clinical[msi_col_final].notna() &
        (clinical[msi_col_final].astype(str) != 'nan') &
        (clinical[msi_col_final].astype(str) != '')
    )
    our_labels_msi = labels[valid_msi.values]
    msi_labels_raw = clinical.loc[valid_msi, msi_col_final].values
    msi_unique = sorted(set(msi_labels_raw))
    msi_int = np.array([msi_unique.index(m) for m in msi_labels_raw])
    ari_vs_msi = adjusted_rand_score(msi_int, our_labels_msi)  # returns 1.0 on empty
```

### Fix

Add a sample-count guard after the validity filter:

```python
if msi_col_final:
    valid_msi = (
        clinical[msi_col_final].notna() &
        (clinical[msi_col_final].astype(str) != 'nan') &
        (clinical[msi_col_final].astype(str) != '')
    )
    if valid_msi.sum() >= 2:  # need at least 2 samples for meaningful ARI
        our_labels_msi = labels[valid_msi.values]
        msi_labels_raw = clinical.loc[valid_msi, msi_col_final].values
        msi_unique = sorted(set(msi_labels_raw))
        msi_int = np.array([msi_unique.index(m) for m in msi_labels_raw])
        ari_vs_msi = adjusted_rand_score(msi_int, our_labels_msi)
    else:
        ari_vs_msi = None
        print(f'[Note] MSI column "{msi_col_final}" has {valid_msi.sum()} valid samples — skipping ARI.')
```

The same guard should be applied to the `results_summary` serialization and the visualization panel.

### Manual correction applied

`results_summary.json`: `"msi_ari": 1.0` changed to `"msi_ari": null` (2026-02-25).

---

## NB17: Low global CMS ARI (0.10) due to k=3 vs k=4 granularity mismatch

**Status:** RESOLVED — k=4 sensitivity analysis executed (Feb 25, 2026). k=3 confirmed as primary.
**File:** `examples/notebooks/17_tcga_cancer_validation.ipynb` (cells 48–51)
**Severity:** Low (documented, not a defect — granularity mismatch is expected)

### Problem

The k=3 TCGA-COAD analysis recovers CMS4 at 75.9% (Fisher OR=16.71, p=1.37e-25), but the global ARI is only 0.1016. This is because 3 subtypes cannot biject to 4 CMS classes — Subtype 2 (n=235, 52%) is a merger subtype absorbing CMS1+CMS2+CMS3.

### k=4 Sensitivity Analysis Results

Refitted GMM at k=4 on the same 452 samples, reusing existing CMS NTP predictions (k-independent).

| Metric | k=3 | k=4 |
|--------|-----|-----|
| Silhouette | 0.1644 | 0.1625 |
| Global CMS ARI | 0.1016 | 0.0824 |
| Cluster sizes | [86, 131, 235] | [124, 220, 26, 82] |

**k=4 Subtype-CMS Mapping (Fisher exact):**

| Subtype | n | Dominant CMS | Concordance | OR | p-value |
|---------|---|-------------|-------------|------|---------|
| 0 | 124 (27%) | CMS2 | 54.5% | 3.35 | 7.67e-08 *** |
| 1 | 220 (49%) | CMS2 | 32.1% | 0.83 | 0.419 ns |
| 2 | 26 (6%) | CMS2 | 52.2% | 2.20 | 0.072 ns |
| 3 | 82 (18%) | CMS4 | 76.2% | 16.50 | 8.84e-25 *** |

### Conclusion

k=4 **worsened** global ARI (0.0824 vs 0.1016, delta=-0.019) and did not produce a clean 1:1 CMS mapping. Three of four subtypes map to CMS2 (no CMS1 or CMS3 recovery). CMS4 recovery is maintained (76.2%) but the additional cluster (Subtype 2, n=26) is a small near-degenerate split. **k=3 remains the primary result.**

### Framing — ⚠️ REVERSED 2026-07 (this section previously said the opposite)

**Superseded guidance.** This section used to read: *"Use CMS4 recovery at 76%
(Fisher p=1.4e-25) as the validation line, not global ARI."* **Do not follow that.**

Preferring a single-subtype enrichment figure over the k-way ARI is exactly the
metric substitution the framework's own cautionary analysis identifies as the way
subtyping results get flattered. The same partition can look strong or weak purely
by metric choice — on TCGA-BRCA, one partition gives k-way ARI 0.218 while LumA
single-subtype enrichment reads 87.6%. Quoting only the second is not a validation
line, it is a selection effect.

**Current guidance.** Report the k-way ARI as the primary recovery metric, and
report single-subtype enrichment *alongside* it, never instead of it. If the two
disagree, that disagreement is the finding and should be stated. The underlying
observation in this section is still valid and worth keeping: 3 subtypes cannot
biject onto 4 CMS classes, so some ARI loss here is genuine granularity mismatch
rather than framework failure. That explains part of the gap — it does not license
swapping in the more flattering number.

See [`consolidation-cautionary/cross-domain/cancer_r38/README.md`](consolidation-cautionary/cross-domain/cancer_r38/README.md)
for the worked example and [`CORRECTION_2026-07/ERRATUM_2026-07-08.md`](CORRECTION_2026-07/ERRATUM_2026-07-08.md).
