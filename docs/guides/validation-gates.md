# Understanding Validation Gates

The framework includes three validation gates to prevent overfitting and ensure discovered subtypes are meaningful.

---

## Why Validation Gates?

Unsupervised clustering can always find "clusters" in data — the question is whether they're real or artifacts. Validation gates answer:

1. **Are clusters better than random?** (Label Shuffle)
2. **Do pathways add signal?** (Random Genes)
3. **Are clusters reproducible?** (Bootstrap Stability)

All three must pass for subtypes to be considered robust.

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

## Running Validation

### Enable All Gates (Recommended)
```yaml
validation:
  run_gates: true
  label_shuffle_iterations: 100
  random_genes_iterations: 100
  bootstrap_iterations: 100
  stability_threshold: 0.8
  negative_control_threshold: 0.15
```

### Quick Validation (Development)
```yaml
validation:
  run_gates: true
  label_shuffle_iterations: 20
  random_genes_iterations: 20
  bootstrap_iterations: 20
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
