# Sensitivity Analysis API

> **Module**: `pathway_subtyping.sensitivity`

Assesses robustness of subtype discovery to analytical parameter choices. Systematically varies clustering algorithm, cluster count, feature set, and normalization method, measuring concordance (pairwise ARI) to identify which parameters most affect results.

---

## Quick Example

```python
from pathway_subtyping import run_sensitivity_analysis

result = run_sensitivity_analysis(
    pathway_scores=scores_df,
    n_clusters=3,
    seed=42,
    robustness_threshold=0.7,
)

print(result.format_report())
print(f"Overall stability: {result.overall_stability:.3f}")
print(f"Robust: {result.is_robust}")
print(f"Most sensitive to: {result.most_sensitive_parameter}")
```

---

## Enums

### `SensitivityParameter`

| Value | Description | Variations |
|-------|-------------|------------|
| `CLUSTERING_ALGORITHM` | Test different algorithms | GMM, K-means, Hierarchical, Spectral |
| `N_CLUSTERS` | Test different cluster counts | k = n-1 to n+2 |
| `NORMALIZATION` | Test normalization methods | Z-score, min-max, robust, rank |
| `FEATURE_SUBSET` | Leave-one-out pathway | Each pathway removed once |

---

## Functions

### `run_sensitivity_analysis()`

```python
def run_sensitivity_analysis(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
    robustness_threshold: float = 0.7,
    parameters: Optional[List[SensitivityParameter]] = None,
    show_progress: bool = True,
) -> SensitivityAnalysisResult
```

Run comprehensive sensitivity analysis across multiple parameters.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `n_clusters` | int | required | Reference number of clusters |
| `seed` | int or None | None | Random seed |
| `robustness_threshold` | float | 0.7 | Minimum mean ARI for robustness |
| `parameters` | List or None | All four | Which parameters to test |
| `show_progress` | bool | True | Show tqdm progress bar |

**Returns:** `SensitivityAnalysisResult` with per-parameter and overall stability metrics.

---

### `vary_clustering_algorithm()`

```python
def vary_clustering_algorithm(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> ParameterVariationResult
```

Test sensitivity to clustering algorithm choice. Runs GMM, K-means, Hierarchical, and Spectral on the same data.

---

### `vary_n_clusters()`

```python
def vary_n_clusters(
    pathway_scores: pd.DataFrame,
    cluster_range: Tuple[int, int] = (2, 6),
    seed: Optional[int] = None,
) -> ParameterVariationResult
```

Test sensitivity to number of clusters. Runs GMM with each k in the specified range.

---

### `vary_feature_subset()`

```python
def vary_feature_subset(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> ParameterVariationResult
```

Leave-one-out feature sensitivity. Removes each pathway once and re-clusters on the remaining features.

---

### `vary_normalization()`

```python
def vary_normalization(
    pathway_scores_raw: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> ParameterVariationResult
```

Test sensitivity to normalization method. Applies Z-score, min-max, robust (IQR), and rank normalization.

---

## Dataclasses

### `SensitivityAnalysisResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `parameter_results` | Dict[str, ParameterVariationResult] | Per-parameter results |
| `overall_stability` | float | Mean ARI across all parameters |
| `most_sensitive_parameter` | str | Lowest mean ARI parameter |
| `least_sensitive_parameter` | str | Highest mean ARI parameter |
| `is_robust` | bool | overall_stability >= threshold |
| `robustness_threshold` | float | Threshold (default 0.7) |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `ParameterVariationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `parameter` | SensitivityParameter | Which parameter was varied |
| `configurations` | List[str] | Configuration values tested |
| `labels_per_config` | Dict[str, np.ndarray] | Labels for each config |
| `pairwise_ari` | Dict[str, float] | ARI between each pair |
| `mean_ari` | float | Mean pairwise ARI |
| `min_ari` | float | Worst-case pairwise ARI |
| `reference_ari` | Dict[str, float] | ARI of each config vs reference |

**Methods:** `to_dict()`

---

## Interpretation

| Metric | Robust | Concerning |
|--------|--------|------------|
| `overall_stability` | > 0.7 | < 0.5 |
| Per-parameter `mean_ari` | > 0.8 | < 0.5 |
| Per-parameter `min_ari` | > 0.5 | < 0.2 |

A result is **robust** if `overall_stability >= robustness_threshold` (default 0.7), meaning the subtypes are stable regardless of analytical choices.

The **most sensitive parameter** identifies which analytical decision matters most — useful for focusing validation effort.

---

## See Also

- [Clustering](clustering.md) — clustering algorithms
- [Benchmark](benchmark.md) — method comparison with ground truth
- [Validation Gates](validation.md) — cluster quality validation
