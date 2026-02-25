# Known Issues

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
