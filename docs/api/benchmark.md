# Benchmark Comparison API

> **Module**: `pathway_subtyping.benchmark`

Compares the framework's pathway-based GMM approach against established methods: NMF + K-means, PCA + K-means, direct gene-level K-means, and a random baseline. Supports single comparisons, multi-condition sweeps, and publication-ready reports.

---

## Quick Example

```python
from pathway_subtyping import run_benchmark_comparison

result = run_benchmark_comparison(
    gene_burdens=burden_df,
    pathway_scores=scores_df,
    pathways=pathway_dict,
    true_labels=true_labels,  # optional ground truth
    n_clusters=3,
    seed=42,
)

print(result.format_report())
print(f"Best method: {result.best_method}")
```

---

## Enums

### `BenchmarkMethod`

| Value | Description | Input Data |
|-------|-------------|------------|
| `PATHWAY_GMM` | Framework's pathway-based GMM (reference) | Pathway scores |
| `NMF_CLUSTERING` | NMF on gene burdens + K-means | Gene burdens |
| `PCA_KMEANS` | PCA on pathway scores + K-means | Pathway scores |
| `GENE_KMEANS` | Direct K-means on gene burdens | Gene burdens |
| `RANDOM_BASELINE` | Random label assignment (null model) | N/A |

---

## Functions

### `run_single_benchmark()`

```python
def run_single_benchmark(
    method: BenchmarkMethod,
    gene_burdens: pd.DataFrame,
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    true_labels: Optional[np.ndarray] = None,
    seed: Optional[int] = None,
) -> BenchmarkResult
```

Run a single benchmark method and evaluate performance.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | BenchmarkMethod | required | Which method to run |
| `gene_burdens` | pd.DataFrame | required | Gene burden matrix (samples x genes) |
| `pathway_scores` | pd.DataFrame | required | Pathway score matrix (samples x pathways) |
| `n_clusters` | int | required | Number of clusters |
| `true_labels` | np.ndarray or None | None | Ground truth for ARI/NMI |
| `seed` | int or None | None | Random seed |

**Returns:** `BenchmarkResult` with labels, metrics, and timing.

---

### `run_benchmark_comparison()`

```python
def run_benchmark_comparison(
    gene_burdens: pd.DataFrame,
    pathway_scores: pd.DataFrame,
    pathways: Dict[str, List[str]],
    true_labels: Optional[np.ndarray] = None,
    n_clusters: Optional[int] = None,
    methods: Optional[List[BenchmarkMethod]] = None,
    seed: Optional[int] = None,
) -> BenchmarkComparisonResult
```

Compare multiple benchmark methods on the same data.

If `n_clusters` is None, uses the number of unique `true_labels` (or defaults to 3). Methods are ranked by ARI (if ground truth provided) or silhouette score.

**Returns:** `BenchmarkComparisonResult` with ranked results.

---

### `run_benchmark_sweep()`

```python
def run_benchmark_sweep(
    effect_sizes: Optional[List[float]] = None,
    sample_sizes: Optional[List[int]] = None,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    methods: Optional[List[BenchmarkMethod]] = None,
    seed: Optional[int] = None,
) -> BenchmarkSweepResult
```

Sweep benchmark across multiple (effect_size, sample_size) conditions using synthetic data. Produces publication-ready results showing which method performs best under different signal strengths.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `effect_sizes` | List[float] or None | `[0.25, 0.5, 1.0, 1.5, 2.0]` | Effect sizes to test |
| `sample_sizes` | List[int] or None | `[100]` | Sample sizes to test |
| `n_pathways` | int | 15 | Pathways per simulation |
| `n_subtypes` | int | 3 | Planted subtypes |
| `methods` | List or None | All methods | Methods to benchmark |
| `seed` | int or None | None | Random seed |

**Returns:** `BenchmarkSweepResult` with per-condition results and summary.

---

## Dataclasses

### `BenchmarkResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `method` | BenchmarkMethod | Method used |
| `predicted_labels` | np.ndarray | Cluster assignments |
| `n_clusters_found` | int | Unique clusters in output |
| `ari` | float or None | ARI vs ground truth |
| `nmi` | float or None | NMI vs ground truth |
| `silhouette` | float | Silhouette score |
| `runtime_seconds` | float | Wall-clock time |
| `converged` | bool | Whether method converged |

**Methods:** `to_dict()`

### `BenchmarkComparisonResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `method_results` | Dict[str, BenchmarkResult] | Per-method results |
| `best_method` | str | Top-ranked method |
| `ranking` | List[str] | Methods sorted by primary metric |
| `n_samples` | int | Number of samples |
| `n_clusters` | int | Number of clusters |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `BenchmarkSweepResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `conditions` | List[Dict] | (effect_size, n_samples) tested |
| `results_per_condition` | List[BenchmarkComparisonResult] | Results per condition |
| `method_mean_ari` | Dict[str, float] | Mean ARI per method |
| `method_wins` | Dict[str, int] | Conditions each method ranked first |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

---

## See Also

- [Clustering](clustering.md) — clustering algorithms
- [Simulation](simulation.md) — synthetic data generation for benchmarks
- [Sensitivity Analysis](sensitivity.md) — parameter robustness
