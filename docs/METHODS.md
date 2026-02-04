# Statistical Methods Documentation

This document provides detailed statistical methodology for the Pathway Subtyping Framework, suitable for methods sections in peer-reviewed publications.

## Overview

The Pathway Subtyping Framework identifies molecular subtypes in genetically heterogeneous diseases by:
1. Computing gene-level burden scores from rare variant data
2. Aggregating burdens into pathway-level scores
3. Clustering samples based on pathway profiles
4. Validating discoveries through rigorous statistical tests

## Gene Burden Scoring

### Variant Selection

Qualifying variants are selected based on:
- **Minor Allele Frequency (MAF)**: < 1% in gnomAD (configurable threshold)
- **Consequence**: Loss-of-function (LoF) or missense variants
- **Quality**: Standard VCF quality filters (PASS status)

### Burden Weighting

Gene burdens are computed as weighted sums of variant counts:

```
burden_g = Σ_v (w_consequence × w_CADD × allele_count_v)
```

Where:
- `w_consequence`: LoF = 2.0, missense = 1.0 (default)
- `w_CADD`: CADD score normalized by percentile rank
- `allele_count`: Number of alternate alleles (0, 1, or 2)

#### Literature-Based Weight Schemes

The framework supports multiple weighting schemes (see `BurdenWeightScheme`):

| Scheme | LoF Weight | High-Impact Missense | Moderate Missense | VUS |
|--------|------------|---------------------|-------------------|-----|
| DEFAULT | 2.0 | 1.5 | 1.0 | 0.5 |
| GNOMAD_CONSTRAINT | pLI-scaled | CADD-scaled | CADD-scaled | 0.1 |
| ACMG_INSPIRED | 3.0 | 2.0 | 1.5 | 0.5 |
| UNIFORM | 1.0 | 1.0 | 1.0 | 1.0 |

### CADD Score Handling

For variants without CADD annotations, consequence-based defaults are applied:
- LoF variants: CADD = 35 (highly deleterious)
- Missense variants: CADD = 20 (moderate impact)
- Other variants: CADD = 10 (low impact)

## Pathway Score Computation

### Gene-to-Pathway Aggregation

Pathway scores are computed from gene burdens using GMT file definitions:

```
pathway_score_p = normalize(aggregate(burden_g for g in pathway_p))
```

### Aggregation Methods

The framework supports multiple aggregation approaches:

| Method | Formula | Use Case |
|--------|---------|----------|
| MEAN | Σburden / n | Default; simple average |
| MEDIAN | median(burdens) | Robust to outliers |
| SIZE_NORMALIZED | Σburden / √n | Corrects for pathway size bias |
| PCA | PC1(burdens) | Captures dominant signal |

### Normalization

After aggregation, pathway scores are Z-score normalized across samples:

```
z_p,i = (x_p,i - μ_p) / σ_p
```

This ensures equal contribution from each pathway regardless of scale.

## Clustering Methodology

### Gaussian Mixture Model (Default)

The primary clustering method uses Gaussian Mixture Models (GMM):

```python
GMM(n_components=k, covariance_type='full', n_init=10, reg_covar=1e-6)
```

Key parameters:
- **n_components**: Number of clusters (determined by model selection)
- **covariance_type**: Full covariance matrices (captures correlations)
- **n_init**: Multiple random initializations for robustness
- **reg_covar**: Regularization for numerical stability

### Alternative Algorithms

For validation, multiple clustering algorithms are supported:

| Algorithm | Characteristics | Best For |
|-----------|-----------------|----------|
| GMM | Soft assignments, BIC selection | Primary analysis |
| K-means | Fast, spherical clusters | Large datasets |
| Hierarchical | Dendogram, no K required | Exploratory |
| Spectral | Nonlinear boundaries | Complex structure |

### Model Selection (Number of Clusters)

The optimal number of clusters is determined by:

**BIC Method** (default):
```
BIC = -2 × log(L) + k × log(n)
```
Select k that minimizes BIC.

**Silhouette Method** (alternative):
```
silhouette = (b - a) / max(a, b)
```
Select k that maximizes mean silhouette score.

## Validation Framework

### Negative Control 1: Label Shuffle Test

Tests whether discovered subtypes are spurious:
1. Shuffle sample labels randomly
2. Recompute pathway scores
3. Recluster and compute ARI vs original labels
4. **Pass criterion**: ARI < 0.15

### Negative Control 2: Random Gene Sets Test

Tests whether pathway structure is meaningful:
1. Replace curated pathways with random gene sets (same sizes)
2. Recompute pathway scores with random genes
3. Recluster and compute ARI vs original labels
4. **Pass criterion**: ARI < 0.15

### Stability Test: Bootstrap Resampling

Tests cluster stability across samples:
1. Generate 100 bootstrap samples (with replacement)
2. Recluster each bootstrap sample
3. Compare to original clustering via ARI
4. **Pass criterion**: Mean ARI ≥ 0.80

### Cross-Validation

K-fold cross-validation for clustering stability:
1. Split data into K folds
2. Train clustering on K-1 folds
3. Predict labels for held-out fold
4. Compare held-out predictions to training labels
5. Report mean and standard deviation of fold ARIs

## Statistical Corrections

### Multiple Testing Correction

All pathway-level comparisons use Benjamini-Hochberg FDR correction:

```
p_adjusted_i = min(1, p_i × m / rank_i)
```

Where m = total tests, rank = ascending p-value rank.

### Permutation-Based P-Values

For effect size significance:
```
p_perm = (# permutations with |effect| ≥ |observed|) / n_permutations
```

Using ≥1000 permutations for robust estimation.

### Confidence Intervals

95% bootstrap confidence intervals for effect sizes:
1. Generate 1000 bootstrap samples
2. Compute effect size for each
3. Report 2.5th and 97.5th percentiles

## Effect Size Metrics

### Cohen's d

For continuous pathway scores between clusters:
```
d = (μ_1 - μ_2) / σ_pooled
```

Interpretation:
- |d| < 0.2: Negligible
- 0.2 ≤ |d| < 0.5: Small
- 0.5 ≤ |d| < 0.8: Medium
- |d| ≥ 0.8: Large

### Adjusted Rand Index (ARI)

For comparing cluster assignments:
```
ARI = (RI - E[RI]) / (max(RI) - E[RI])
```

Range: -1 to 1 (1 = perfect agreement, 0 = random)

### Normalized Mutual Information (NMI)

Alternative clustering comparison metric:
```
NMI = 2 × I(U,V) / (H(U) + H(V))
```

Range: 0 to 1 (1 = perfect agreement)

## Power Analysis

### Type I Error Estimation

Estimate false positive rate under null:
1. Generate random data with no true structure
2. Apply clustering pipeline
3. Measure rate of ARI > threshold

Reported as Type I error rate at threshold = 0.15.

### Power Curves

Power as a function of effect size:
1. Generate data with planted subtypes at various effect sizes
2. Apply clustering pipeline
3. Measure recovery rate (ARI ≥ 0.80)
4. Report power (proportion successful) at each effect size

### Sample Size Recommendations

Based on power analysis:
- **Minimum**: n = 50 per subtype (effect size ≥ 1.5)
- **Recommended**: n = 100 per subtype (effect size ≥ 1.0)
- **Optimal**: n = 200 per subtype (effect size ≥ 0.5)

## Reproducibility

### Random Seed Control

All stochastic operations use controlled random seeds:
```python
np.random.seed(config.seed)
sklearn_model(..., random_state=config.seed)
```

### Metadata Logging

Each analysis logs:
- Software versions (package, Python, sklearn)
- Configuration parameters
- Input file checksums
- Timestamp and execution environment

## Simulation Validation

### Synthetic Data Generation

For validation, synthetic data with known structure is generated:
1. Define n_subtypes with specified proportions
2. Assign subset of pathways to each subtype
3. Elevate pathway scores by effect_size × noise_level
4. Add Gaussian noise

### Ground Truth Recovery

Clustering performance is evaluated against planted labels:
- **ARI**: Primary metric (accounts for label permutation)
- **NMI**: Secondary metric (information-theoretic)
- **Correct K**: Whether optimal K matches true K

## Reporting Standards

### Required Outputs

All analyses report:
1. Number of samples and pathways
2. Optimal K and selection method
3. Per-cluster sample sizes
4. Top contributing pathways per cluster (FDR < 0.05)
5. Validation gate results (pass/fail)
6. Effect sizes with 95% CIs

### Visualization Standards

Recommended plots:
- Heatmap of pathway scores by cluster
- PCA/UMAP projection colored by cluster
- BIC/silhouette curves for K selection
- Bootstrap ARI distribution

## References

- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. J R Stat Soc B.
- Cohen J (1988). Statistical Power Analysis for the Behavioral Sciences.
- Hubert L, Arabie P (1985). Comparing partitions. J Classif.
- McLachlan GJ, Peel D (2000). Finite Mixture Models.
- Strehl A, Ghosh J (2002). Cluster ensembles. J Mach Learn Res.

---

*Document version: 0.2.0 | Last updated: February 2026*
