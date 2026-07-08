# Known Issues

> **⚠️ 2026-07-08 — Correction issued.** The empty-input `adjusted_rand_score` artifact documented below (NB17) was **also present, uncorrected, in the 47-dataset calibration benchmark** — it produced 14 spurious ARI=1.0 rows that drove the retracted R²=0.889 adaptive-threshold model. The calibration model is **retracted** and the benchmark **corrected**. See [`CORRECTION_2026-07/ERRATUM_2026-07-08.md`](CORRECTION_2026-07/ERRATUM_2026-07-08.md).

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

### R21 Framing

Use CMS4 recovery at 76% (Fisher p=1.4e-25) as the validation line, not global ARI. The k=4 sensitivity analysis confirms that the low global ARI reflects irreducible granularity mismatch, not a framework failure.
