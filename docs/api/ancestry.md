# Ancestry Correction API

> **Module**: `pathway_subtyping.ancestry`

Detects and corrects population stratification in pathway-based subtype analysis. Provides PCA-based ancestry inference from genotype data, regression-based score correction, independence testing for cluster validation, and stratified analysis within ancestry groups.

---

## Quick Example

```python
from pathway_subtyping import (
    compute_ancestry_pcs,
    adjust_pathway_scores,
    check_ancestry_independence,
    AncestryMethod,
)

# Step 1: Compute ancestry PCs from genotype data
pcs = compute_ancestry_pcs(genotype_matrix, n_components=10, seed=42)
print(f"Top 10 PCs explain {sum(pcs.explained_variance_ratio):.1%} of variance")

# Step 2: Correct pathway scores
result = adjust_pathway_scores(pathway_scores, pcs, method=AncestryMethod.REGRESS_OUT)
print(result.format_report())
# Use result.adjusted_scores for downstream clustering

# Step 3: Verify clusters are ancestry-independent
report = check_ancestry_independence(cluster_labels, pcs)
print(report.format_report())
print(f"Independent of ancestry: {report.overall_independence_passed}")
```

---

## Enums

### `AncestryMethod`

| Value | Description |
|-------|-------------|
| `REGRESS_OUT` | Regress pathway scores on ancestry PCs and take residuals. Simple, widely used (Price et al., 2006). |
| `COVARIATE_AWARE` | Include ancestry PCs as covariates. Currently uses same residualization as REGRESS_OUT; reserved for future mixed-model approaches. |
| `STRATIFIED` | Cluster independently within each ancestry group, then compare concordance. Use with `stratified_analysis()`. |

---

## Functions

### `compute_ancestry_pcs()`

```python
def compute_ancestry_pcs(
    genotype_matrix: pd.DataFrame,
    n_components: int = 10,
    seed: Optional[int] = None,
) -> AncestryPCs
```

Compute principal components from genotype data for ancestry inference. Standardizes variants (mean-center, unit variance per variant), handles missing genotypes, and runs PCA.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `genotype_matrix` | pd.DataFrame | required | Genotype dosages (samples x variants), values 0/1/2 |
| `n_components` | int | 10 | Number of PCs to compute (Price et al. recommend 10) |
| `seed` | int or None | None | Random seed for PCA |

**Returns:** `AncestryPCs` with components DataFrame and variance explained.

---

### `adjust_pathway_scores()`

```python
def adjust_pathway_scores(
    pathway_scores: pd.DataFrame,
    ancestry_pcs: AncestryPCs,
    method: AncestryMethod = AncestryMethod.REGRESS_OUT,
    n_pcs: Optional[int] = None,
    confounding_threshold: float = 0.1,
) -> AncestryAdjustmentResult
```

Remove ancestry-correlated variance from pathway scores to prevent spurious subtypes driven by population structure.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `ancestry_pcs` | AncestryPCs | required | From `compute_ancestry_pcs()` |
| `method` | AncestryMethod | REGRESS_OUT | Correction method |
| `n_pcs` | int or None | All | Number of PCs to use |
| `confounding_threshold` | float | 0.1 | R^2 threshold for flagging confounded pathways |

**Returns:** `AncestryAdjustmentResult` with corrected scores and per-pathway R^2 diagnostics.

**Raises:** `ValueError` if method is `STRATIFIED` (use `stratified_analysis()` instead).

---

### `check_ancestry_independence()`

```python
def check_ancestry_independence(
    cluster_labels: np.ndarray,
    ancestry_pcs: AncestryPCs,
    significance_threshold: float = 0.05,
    n_pcs_to_test: Optional[int] = None,
) -> AncestryStratificationReport
```

Test whether discovered clusters are independent of ancestry PCs using Kruskal-Wallis H-test with Bonferroni correction.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_labels` | np.ndarray | required | Cluster assignments |
| `ancestry_pcs` | AncestryPCs | required | From `compute_ancestry_pcs()` |
| `significance_threshold` | float | 0.05 | P-value threshold (before Bonferroni) |
| `n_pcs_to_test` | int or None | min(5, n_components) | Number of top PCs to test |

**Returns:** `AncestryStratificationReport` with per-PC p-values and overall pass/fail.

---

### `stratified_analysis()`

```python
def stratified_analysis(
    pathway_scores: pd.DataFrame,
    ancestry_groups: np.ndarray,
    n_clusters: int,
    seed: Optional[int] = None,
) -> Dict[str, Any]
```

Cluster independently within each ancestry group, then evaluate cross-group concordance. Groups with fewer than `n_clusters * 3` samples are skipped.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `ancestry_groups` | np.ndarray | required | Ancestry group labels per sample |
| `n_clusters` | int | required | Number of clusters per group |
| `seed` | int or None | None | Random seed |

**Returns:** Dictionary with `per_group_results`, `cross_group_concordance`, and metadata.

---

### `compute_ancestry_correlation()`

```python
def compute_ancestry_correlation(
    pathway_scores: pd.DataFrame,
    ancestry_pcs: AncestryPCs,
    n_pcs: Optional[int] = None,
) -> pd.DataFrame
```

Compute Pearson correlation between each pathway score and each ancestry PC. Useful for identifying confounded pathways before and after correction.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `ancestry_pcs` | AncestryPCs | required | From `compute_ancestry_pcs()` |
| `n_pcs` | int or None | min(5, n_components) | Number of PCs to correlate |

**Returns:** DataFrame of Pearson correlations (pathways x PCs).

---

## Dataclasses

### `AncestryPCs`

| Attribute | Type | Description |
|-----------|------|-------------|
| `components` | pd.DataFrame | Ancestry PCs (samples x n_components) |
| `explained_variance_ratio` | np.ndarray | Variance explained by each PC |
| `n_components` | int | Number of PCs computed |
| `sample_ids` | List[str] | Sample identifiers |

**Methods:** `to_dict()`

### `AncestryAdjustmentResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `adjusted_scores` | pd.DataFrame | Corrected pathway scores (samples x pathways) |
| `method` | AncestryMethod | Method used |
| `r_squared_per_pathway` | Dict[str, float] | R^2 of ancestry PCs for each pathway |
| `n_pcs_used` | int | Number of ancestry PCs used |
| `highly_confounded_pathways` | List[str] | Pathways with R^2 > threshold |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `AncestryStratificationReport`

| Attribute | Type | Description |
|-----------|------|-------------|
| `ancestry_cluster_correlation` | Dict[str, Dict[str, float]] | Per-PC per-cluster mean values |
| `independence_test_pvalues` | Dict[str, float] | Per-PC Kruskal-Wallis p-values |
| `overall_independence_passed` | bool | Whether clusters are ancestry-independent |
| `significance_threshold` | float | Alpha level used (default 0.05) |

**Methods:** `to_dict()`, `format_report()`

---

## Interpretation

| Metric | Good | Concerning |
|--------|------|------------|
| R^2 per pathway (before correction) | < 0.05 | > 0.10 |
| R^2 per pathway (after correction) | < 0.02 | > 0.05 |
| Independence test | All PCs p > 0.05/n | Any PC p < 0.05/n |
| Cross-group concordance | > 0.7 | < 0.5 |

Highly confounded pathways (R^2 > 0.1) have substantial ancestry-driven variance that could create spurious subtypes. After correction, R^2 should decrease substantially.

---

## References

- Price AL, et al. Principal components analysis corrects for stratification in genome-wide association studies. *Nat Genet*. 2006;38(8):904-909.
- Patterson N, et al. Population structure and eigenanalysis. *PLoS Genet*. 2006;2(12):e190.

---

## See Also

- [Batch Correction](batch_correction.md) — technical batch effect removal
- [Validation Gates](validation.md) — cluster quality validation
- [Sensitivity](sensitivity.md) — robustness testing
