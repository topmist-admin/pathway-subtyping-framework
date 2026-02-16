# Clustering API

> **Module**: `pathway_subtyping.clustering`

Implements multiple clustering algorithms for molecular subtype discovery, model selection (optimal K), cross-validation, and algorithm comparison.

---

## Quick Example

```python
from pathway_subtyping import run_clustering, ClusteringAlgorithm

# Run GMM clustering
result = run_clustering(
    data=pathway_scores.values,
    n_clusters=3,
    algorithm=ClusteringAlgorithm.GMM,
    seed=42,
)
print(f"Silhouette: {result.silhouette:.3f}")
print(f"BIC: {result.bic:.1f}")

# Select optimal K
from pathway_subtyping import select_n_clusters

selection = select_n_clusters(pathway_scores.values, k_range=[2, 3, 4, 5, 6], seed=42)
print(f"Optimal K: {selection.optimal_k}")
```

---

## Enums

### `ClusteringAlgorithm`

| Value | Description | Key Characteristics |
|-------|-------------|---------------------|
| `GMM` | Gaussian Mixture Model | Soft assignments, BIC for K selection, full covariance |
| `KMEANS` | K-means | Fast, spherical clusters, hard assignments |
| `HIERARCHICAL` | Agglomerative (Ward linkage) | Dendrogram-based, deterministic |
| `SPECTRAL` | Spectral clustering (RBF kernel) | Nonlinear cluster boundaries |

---

## Functions

### `run_clustering()`

```python
def run_clustering(
    data: np.ndarray,
    n_clusters: int,
    algorithm: ClusteringAlgorithm = ClusteringAlgorithm.GMM,
    seed: Optional[int] = None,
    **kwargs,
) -> ClusteringResult
```

Run clustering with specified algorithm.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | np.ndarray | required | Feature matrix (samples x features) |
| `n_clusters` | int | required | Number of clusters |
| `algorithm` | ClusteringAlgorithm | `GMM` | Clustering algorithm |
| `seed` | int or None | None | Random seed |

**Algorithm-specific kwargs:**
- GMM: `covariance_type` (default: `"full"`), `n_init` (default: 10), `reg_covar` (default: 1e-6)
- K-means: `n_init` (default: 10), `max_iter` (default: 300)
- Hierarchical: `linkage_method` (default: `"ward"`)
- Spectral: `affinity` (default: `"rbf"`)

**Returns:** `ClusteringResult` with labels and quality metrics.

---

### `select_n_clusters()`

```python
def select_n_clusters(
    data: np.ndarray,
    k_range: List[int] = None,
    method: str = "bic",
    seed: Optional[int] = None,
) -> ModelSelectionResult
```

Select optimal number of clusters.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | np.ndarray | required | Feature matrix (samples x features) |
| `k_range` | List[int] or None | `[2..8]` | K values to test |
| `method` | str | `"bic"` | Selection method (`"bic"` or `"silhouette"`) |
| `seed` | int or None | None | Random seed |

**Returns:** `ModelSelectionResult` with optimal K and per-K scores.

---

### `cross_validate_clustering()`

```python
def cross_validate_clustering(
    data: np.ndarray,
    n_clusters: int,
    algorithm: ClusteringAlgorithm = ClusteringAlgorithm.GMM,
    n_folds: int = 5,
    seed: Optional[int] = None,
) -> CrossValidationResult
```

Cross-validate clustering stability using K-fold splits.

**Returns:** `CrossValidationResult` with mean/std ARI across folds.

---

### `compare_algorithms()`

```python
def compare_algorithms(
    data: np.ndarray,
    n_clusters: int,
    algorithms: Optional[List[ClusteringAlgorithm]] = None,
    seed: Optional[int] = None,
) -> AlgorithmComparisonResult
```

Compare multiple clustering algorithms on the same data.

**Returns:** `AlgorithmComparisonResult` with pairwise ARI, per-algorithm metrics, and most stable algorithm.

---

## Dataclasses

### `ClusteringResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `labels` | np.ndarray | Cluster assignments per sample |
| `n_clusters` | int | Number of clusters |
| `algorithm` | ClusteringAlgorithm | Algorithm used |
| `silhouette` | float | Silhouette score |
| `calinski_harabasz` | float | Calinski-Harabasz index |
| `davies_bouldin` | float | Davies-Bouldin index |
| `bic` | float or None | BIC (GMM only) |
| `converged` | bool | Whether algorithm converged |

**Methods:** `to_dict()`

### `ModelSelectionResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `optimal_k` | int | Recommended number of clusters |
| `k_range` | List[int] | K values tested |
| `bic_values` | Dict[int, float] | BIC per K (GMM) |
| `silhouette_values` | Dict[int, float] | Silhouette per K |
| `method` | str | Selection method used |

**Methods:** `to_dict()`

### `CrossValidationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `mean_ari` | float | Mean ARI across folds |
| `std_ari` | float | Standard deviation of ARI |
| `fold_aris` | List[float] | ARI per fold |
| `n_folds` | int | Number of folds |

**Methods:** `to_dict()`

### `AlgorithmComparisonResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `results` | Dict[str, ClusteringResult] | Per-algorithm results |
| `pairwise_ari` | Dict[str, float] | ARI between algorithm pairs |
| `consensus_labels` | np.ndarray | Labels from most stable algorithm |
| `most_stable_algorithm` | str | Algorithm with highest average ARI |

**Methods:** `to_dict()`

---

## Clustering Quality Metrics

| Metric | Range | Interpretation |
|--------|-------|----------------|
| Silhouette | -1 to 1 | Higher = better separated clusters |
| Calinski-Harabasz | 0 to inf | Higher = denser, well-separated clusters |
| Davies-Bouldin | 0 to inf | Lower = better separation |
| BIC | -inf to inf | Lower = better model fit (GMM) |

---

## See Also

- [Sensitivity Analysis](sensitivity.md) — robustness to algorithm choice
- [Validation Gates](validation.md) — cluster quality validation
- [Benchmark](benchmark.md) — method comparison
