# Simulation & Power Analysis API

> **Module**: `pathway_subtyping.simulation`

Generates synthetic data with known ground truth for framework validation: planted subtype structures, configurable effect sizes, Type I error estimation, power analysis, and sample size recommendations. Also supports expression data simulation.

---

## Quick Example

```python
from pathway_subtyping import SimulationConfig, generate_synthetic_data, evaluate_recovery

# Generate synthetic data with 3 planted subtypes
config = SimulationConfig(
    n_samples=200,
    n_pathways=15,
    n_subtypes=3,
    effect_size=1.0,
    seed=42,
)
sim_data = generate_synthetic_data(config)

# Cluster and evaluate recovery
from pathway_subtyping import run_clustering
result = run_clustering(sim_data.pathway_scores.values, n_clusters=3, seed=42)
recovery = evaluate_recovery(result.labels, sim_data.true_labels)
print(f"ARI: {recovery.ari:.3f}, Correct K: {recovery.correct_k}")
```

---

## Functions

### `generate_synthetic_data()`

```python
def generate_synthetic_data(
    config: SimulationConfig,
) -> SimulatedData
```

Generate synthetic variant-based data with planted subtype structure. Creates pathway scores where a subset of pathways have different means between subtypes. Optionally adds simulated ancestry/population structure.

**Returns:** `SimulatedData` with gene burdens, pathway scores, pathways, and ground truth labels.

---

### `generate_synthetic_expression_data()`

```python
def generate_synthetic_expression_data(
    config: ExpressionSimulationConfig,
) -> SimulatedExpressionData
```

Generate synthetic bulk RNA-seq expression data with planted subtypes. Simulates log2-scale expression with zero-inflation (dropout), subtype-specific pathway effects, and configurable noise.

**Returns:** `SimulatedExpressionData` with gene expression matrix, pathways, and labels.

---

### `evaluate_recovery()`

```python
def evaluate_recovery(
    predicted_labels: np.ndarray,
    true_labels: np.ndarray,
) -> RecoveryResult
```

Evaluate how well clustering recovered the true subtype structure.

**Returns:** `RecoveryResult` with ARI, NMI, and whether correct K was found.

---

### `estimate_type_i_error()`

```python
def estimate_type_i_error(
    n_samples: int = 100,
    n_pathways: int = 15,
    n_clusters: int = 3,
    ari_threshold: float = 0.15,
    n_simulations: int = 100,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> TypeIErrorResult
```

Estimate Type I error rate (false positive clustering). Generates random data with no true structure and measures how often clustering finds spurious patterns.

**Returns:** `TypeIErrorResult` with null ARI distribution and false positive rate.

---

### `run_power_analysis()`

```python
def run_power_analysis(
    n_samples: int = 100,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    effect_sizes: Optional[List[float]] = None,
    ari_threshold: float = 0.8,
    n_simulations_per_effect: int = 50,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> PowerAnalysisResult
```

Run power analysis across effect sizes. Determines what effect size is needed for reliable subtype recovery.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_samples` | int | 100 | Samples per simulation |
| `n_pathways` | int | 15 | Number of pathways |
| `n_subtypes` | int | 3 | Planted subtypes |
| `effect_sizes` | List or None | `[0.1, 0.25, ..., 2.0]` | Effect sizes to test |
| `ari_threshold` | float | 0.8 | ARI for "successful" recovery |
| `n_simulations_per_effect` | int | 50 | Simulations per effect size |
| `seed` | int or None | None | Random seed |

**Returns:** `PowerAnalysisResult` with power at each effect size and recommended minimum.

---

### `run_sample_size_analysis()`

```python
def run_sample_size_analysis(
    effect_size: float = 1.0,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    sample_sizes: Optional[List[int]] = None,
    ari_threshold: float = 0.8,
    n_simulations: int = 50,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> SampleSizeAnalysisResult
```

Analyze power as a function of sample size at a fixed effect size.

**Returns:** `SampleSizeAnalysisResult` with power by sample size and recommended N.

---

### `validate_framework()`

```python
def validate_framework(
    n_samples: int = 200,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    effect_size: float = 1.0,
    n_runs: int = 10,
    seed: Optional[int] = None,
) -> Dict[str, Any]
```

Run comprehensive framework validation: multiple recovery runs + Type I error estimation.

---

## Dataclasses

### `SimulationConfig`

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_samples` | int | 200 | Number of samples |
| `n_pathways` | int | 15 | Number of pathways |
| `n_genes_per_pathway` | int | 20 | Average genes per pathway |
| `n_subtypes` | int | 3 | Planted subtypes |
| `effect_size` | float | 1.0 | Cohen's d for subtype differences |
| `noise_level` | float | 1.0 | Background noise std |
| `subtype_proportions` | List[float] or None | Equal | Relative sizes |
| `seed` | int or None | 42 | Random seed |
| `n_ancestry_groups` | int | 0 | Ancestry groups (0 = none) |
| `ancestry_effect_size` | float | 0.5 | Ancestry effect on scores |
| `ancestry_confounding` | float | 0.0 | Correlation between subtypes and ancestry |

### `SimulatedData`

| Attribute | Type | Description |
|-----------|------|-------------|
| `gene_burdens` | pd.DataFrame | Gene burden matrix (samples x genes) |
| `pathway_scores` | pd.DataFrame | Z-normalized pathway scores |
| `pathways` | Dict[str, List[str]] | Pathway definitions |
| `true_labels` | np.ndarray | Ground truth subtype labels |
| `config` | SimulationConfig | Config used |
| `subtype_pathway_effects` | Dict[int, List[str]] | Which pathways differ per subtype |
| `ancestry_labels` | np.ndarray or None | Ancestry group assignments |
| `ancestry_pcs` | pd.DataFrame or None | Simulated ancestry PCs |

**Methods:** `to_dict()`

### `ExpressionSimulationConfig`

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_samples` | int | 200 | Number of samples |
| `n_genes` | int | 500 | Total genes |
| `n_pathways` | int | 15 | Number of pathways |
| `n_subtypes` | int | 3 | Planted subtypes |
| `effect_size` | float | 1.5 | Log2 fold-change |
| `noise_level` | float | 1.0 | Background noise |
| `dropout_rate` | float | 0.1 | Zero-inflation rate |
| `seed` | int or None | 42 | Random seed |

### `RecoveryResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `ari` | float | Adjusted Rand Index |
| `nmi` | float | Normalized Mutual Information |
| `n_clusters_predicted` | int | Clusters found |
| `n_clusters_true` | int | True number of subtypes |
| `correct_k` | bool | Whether correct K was found |

**Methods:** `to_dict()`

### `TypeIErrorResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `null_ari_mean` | float | Mean ARI under null |
| `null_ari_std` | float | Std of null ARI |
| `type_i_rate` | float | Proportion exceeding threshold |
| `threshold` | float | ARI threshold used |
| `n_simulations` | int | Simulations run |

**Methods:** `to_dict()`

### `PowerAnalysisResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `effect_sizes` | List[float] | Effect sizes tested |
| `recovery_rates` | Dict[float, List[float]] | ARIs per effect size |
| `power_at_threshold` | Dict[float, float] | Power per effect size |
| `threshold` | float | ARI threshold for success |
| `recommended_effect_size` | float or None | Min effect for 80% power |

**Methods:** `to_dict()`

### `SampleSizeAnalysisResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `sample_sizes` | List[int] | Sample sizes tested |
| `power_by_size` | Dict[int, float] | Power per sample size |
| `effect_size` | float | Effect size used |
| `threshold` | float | ARI threshold |
| `recommended_n` | int or None | Min N for 80% power |

**Methods:** `to_dict()`

---

## See Also

- [Benchmark](benchmark.md) — comparing methods using synthetic data
- [METHODS.md](../METHODS.md) — statistical methodology
- [Clustering](clustering.md) — clustering algorithms
