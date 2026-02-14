# Understanding Validation Gates

The framework includes three validation gates to prevent overfitting and ensure discovered subtypes are meaningful.

---

## Why Validation Gates?

Unsupervised clustering can always find "clusters" in data — the question is whether they're real or artifacts. Validation gates answer:

1. **Are clusters better than random?** (Label Shuffle)
2. **Do pathways add signal?** (Random Genes)
3. **Are clusters reproducible?** (Bootstrap Stability)

All three must pass for subtypes to be considered robust.

> **New in v0.3:** Validation thresholds can be **automatically calibrated** based on your dataset's sample size and cluster count, replacing the fixed defaults. See [Threshold Calibration](#threshold-calibration) below.

---

## Gate 1: Label Shuffle Test

### What It Does
Randomly permutes sample labels and re-runs clustering. If real subtypes exist, shuffled data should NOT cluster as well.

### How It Works
```
1. Shuffle pathway scores across samples (break sample-pathway associations)
2. Run GMM clustering on shuffled data
3. Compare shuffled clusters to original clusters using ARI
4. Repeat N times (default: 100)
5. Report mean ARI across iterations
```

### Interpretation

| ARI | Meaning | Action |
|-----|---------|--------|
| < 0.10 | Strong pass | Clusters are not random |
| 0.10 - 0.15 | Pass | Acceptable |
| 0.15 - 0.20 | Borderline | Review cluster count |
| > 0.20 | Fail | Likely overfitting |

### If It Fails
- Reduce number of clusters
- Increase cohort size
- Check for batch effects in data

---

## Gate 2: Random Genes Test

### What It Does
Replaces disease pathways with random gene sets and re-runs the pipeline. If your pathways capture real biology, random genes should NOT work as well.

### How It Works
```
1. Generate random "pathways" with same sizes as real pathways
2. Score samples using random gene sets
3. Run GMM clustering
4. Compare random clusters to real pathway clusters using ARI
5. Repeat N times (default: 100)
6. Report mean ARI across iterations
```

### Interpretation

| ARI | Meaning | Action |
|-----|---------|--------|
| < 0.10 | Strong pass | Pathways capture disease biology |
| 0.10 - 0.15 | Pass | Pathways add signal |
| 0.15 - 0.20 | Borderline | Review pathway definitions |
| > 0.20 | Fail | Pathways may not be informative |

### If It Fails
- Pathways may be too generic
- Curate more specific gene sets
- Review gene evidence quality
- See [Pathway Curation Guide](pathway-curation-guide.md)

---

## Gate 3: Bootstrap Stability Test

### What It Does
Resamples data with replacement and re-clusters. Robust subtypes should be stable across resamples.

### How It Works
```
1. Sample N patients with replacement (bootstrap sample)
2. Run GMM clustering on bootstrap sample
3. Compare bootstrap clusters to original using ARI
4. Repeat N times (default: 100)
5. Report mean and std of ARI across iterations
```

### Interpretation

| Mean ARI | Std ARI | Meaning | Action |
|----------|---------|---------|--------|
| > 0.85 | < 0.05 | Strong pass | Highly stable subtypes |
| 0.80 - 0.85 | < 0.10 | Pass | Acceptable stability |
| 0.70 - 0.80 | any | Borderline | Consider merging clusters |
| < 0.70 | any | Fail | Subtypes are unstable |

### If It Fails
- Reduce number of clusters (k)
- Clusters may be overlapping — try merging
- Increase cohort size if possible
- Some clusters may be "real" while others are noise

---

## Threshold Calibration

### Why Calibrate?

The default thresholds (0.15 for null ARI, 0.8 for stability) don't account for sample size or number of clusters:

- **Small samples** (n=30) produce noisier ARI distributions — 0.15 may be too permissive
- **Large samples** (n=500) have tighter distributions — 0.15 may be too strict
- **More clusters** inflate chance ARI — the same threshold is wrong for k=2 vs k=8

### How It Works

The framework uses **pre-computed lookup tables** derived from empirical simulations across a grid of (n_samples, n_clusters) configurations. At runtime:

1. **Lookup**: If your (n, k) matches a grid point, use the pre-computed threshold
2. **Interpolate**: If between grid points, bilinearly interpolate
3. **Simulate**: If outside the grid range, run fresh simulations on-the-fly

### Auto-Calibration (Recommended)

Set thresholds to `null` in your config to enable automatic calibration:

```yaml
validation:
  run_gates: true
  calibrate: true          # Enable auto-calibration (default: true)
  stability_threshold: null  # null = auto-calibrate based on n_samples and n_clusters
  null_ari_max: null         # null = auto-calibrate
  alpha: 0.05               # Significance level for calibration
  n_permutations: 100
  n_bootstrap: 50
```

### Manual Override

To use specific thresholds (disables auto-calibration for that threshold):

```yaml
validation:
  stability_threshold: 0.85  # Explicit value overrides auto-calibration
  null_ari_max: 0.10          # Stricter null threshold
```

### Programmatic Calibration

```python
from pathway_subtyping import calibrate_thresholds

# Auto-calibrate for your data
ct = calibrate_thresholds(n_samples=150, n_clusters=4)
print(f"Null ARI threshold: {ct.null_ari_threshold:.4f}")
print(f"Stability threshold: {ct.stability_threshold:.4f}")
print(ct.format_report())
```

---

## Running Validation

### Enable All Gates with Auto-Calibration (Recommended)
```yaml
validation:
  run_gates: true
  calibrate: true
  n_permutations: 100
  n_bootstrap: 50
```

### Enable All Gates with Fixed Thresholds
```yaml
validation:
  run_gates: true
  stability_threshold: 0.8
  null_ari_max: 0.15
  n_permutations: 100
  n_bootstrap: 50
```

### Quick Validation (Development)
```yaml
validation:
  run_gates: true
  n_permutations: 20
  n_bootstrap: 20
```

### Skip Validation (Not Recommended)
```yaml
validation:
  run_gates: false
```

---

## Understanding ARI (Adjusted Rand Index)

ARI measures similarity between two clusterings:

| ARI Value | Meaning |
|-----------|---------|
| 1.0 | Perfect agreement |
| 0.5 - 1.0 | High similarity |
| 0.0 | No better than random |
| < 0.0 | Worse than random |

For negative controls (shuffle, random genes), we WANT low ARI.
For stability (bootstrap), we WANT high ARI.

---

## Validation Output

The framework generates a validation report:

```json
{
  "validation_gates": {
    "all_passed": true,
    "tests": [
      {
        "name": "label_shuffle",
        "status": "PASS",
        "threshold": 0.15,
        "observed": 0.03,
        "iterations": 100
      },
      {
        "name": "random_genes",
        "status": "PASS",
        "threshold": 0.15,
        "observed": 0.08,
        "iterations": 100
      },
      {
        "name": "bootstrap_stability",
        "status": "PASS",
        "threshold": 0.80,
        "observed": 0.89,
        "std": 0.04,
        "iterations": 100
      }
    ]
  }
}
```

---

## Decision Tree

```
                    Run Pipeline
                         │
                         ▼
              ┌─────────────────────┐
              │  Label Shuffle Pass? │
              └──────────┬──────────┘
                    │
           ┌────────┴────────┐
           │                 │
          Yes               No
           │                 │
           ▼                 ▼
   ┌───────────────┐   Reduce clusters
   │ Random Genes  │   or increase N
   │    Pass?      │
   └───────┬───────┘
           │
    ┌──────┴──────┐
    │             │
   Yes           No
    │             │
    ▼             ▼
┌──────────┐   Refine pathways
│ Bootstrap │   (curation guide)
│   Pass?   │
└─────┬────┘
      │
  ┌───┴───┐
  │       │
 Yes     No
  │       │
  ▼       ▼
PUBLISH  Merge clusters
         or increase N
```

---

## Common Questions

### Q: Can I publish if one gate fails?
**A:** Not recommended. Report the failure transparently and investigate the cause before claiming subtypes.

### Q: How many iterations are enough?
**A:** 100 is standard. For publication, consider 1000 for tighter confidence intervals.

### Q: My bootstrap is 0.78 — is that okay?
**A:** Borderline. Try reducing k by 1 and re-running. If bootstrap improves significantly, the original k was too high.

### Q: Random genes ARI is 0.12 — should I worry?
**A:** That's acceptable (<0.15), but close to threshold. Review pathway quality to see if you can improve specificity.

---

## References

- Hubert & Arabie (1985). Comparing partitions. Journal of Classification.
- Lange et al. (2004). Stability-based validation of clustering solutions. Neural Computation.
- Hennig (2007). Cluster-wise assessment of cluster stability. Computational Statistics & Data Analysis.
- Vinh NX, Epps J, Bailey J (2010). Information Theoretic Measures for Clusterings Comparison. J Mach Learn Res.
