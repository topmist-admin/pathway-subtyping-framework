# Statistical Rigor API

> **Module**: `pathway_subtyping.statistical_rigor`

Publication-quality statistical methods: FDR correction (Benjamini-Hochberg), permutation-based p-values, Cohen's d effect sizes with bootstrap confidence intervals, literature-based burden weight schemes, and pathway aggregation normalization.

---

## Quick Example

```python
from pathway_subtyping import run_statistical_analysis, BurdenWeightScheme

result = run_statistical_analysis(
    pathway_scores=scores_df,
    cluster_labels=labels,
    weight_scheme=BurdenWeightScheme.GNOMAD_CONSTRAINT,
    fdr_alpha=0.05,
    n_permutations=1000,
    n_bootstrap=1000,
    seed=42,
)

print(result.format_report())
print(f"Significant pathways: {result.get_significant_pathways()}")
```

---

## Enums

### `BurdenWeightScheme`

| Value | Description | Reference |
|-------|-------------|-----------|
| `DEFAULT` | Ad-hoc weights (backwards compatibility) | — |
| `GNOMAD_CONSTRAINT` | gnomAD constraint metrics | Lek et al., Nature 2016 |
| `ACMG_INSPIRED` | ACMG variant classification | Richards et al., Genet Med 2015 |
| `CADD_PERCENTILE` | CADD score percentile-based | Kircher et al., Nat Genet 2014 |
| `UNIFORM` | All damaging variants equal weight | Sensitivity analysis |

### `PathwayNormalization`

| Value | Description |
|-------|-------------|
| `MEAN` | Simple mean of gene burdens |
| `MEDIAN` | Median (robust to outliers) |
| `SIZE_NORMALIZED` | Mean / sqrt(n_genes) — corrects for pathway size bias |
| `PCA` | First principal component of gene burdens |

---

## Functions

### `run_statistical_analysis()`

```python
def run_statistical_analysis(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    weight_scheme: BurdenWeightScheme = BurdenWeightScheme.GNOMAD_CONSTRAINT,
    normalization: PathwayNormalization = PathwayNormalization.SIZE_NORMALIZED,
    fdr_alpha: float = 0.05,
    n_permutations: int = 1000,
    n_bootstrap: int = 1000,
    seed: Optional[int] = None,
) -> StatisticalRigorResult
```

Run comprehensive statistical analysis on clustering results.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `cluster_labels` | np.ndarray | required | Cluster assignments |
| `weight_scheme` | BurdenWeightScheme | `GNOMAD_CONSTRAINT` | Burden weighting scheme used |
| `normalization` | PathwayNormalization | `SIZE_NORMALIZED` | Pathway aggregation method |
| `fdr_alpha` | float | 0.05 | FDR significance threshold |
| `n_permutations` | int | 1000 | Permutations for p-value calculation |
| `n_bootstrap` | int | 1000 | Bootstrap iterations for confidence intervals |
| `seed` | int or None | None | Random seed |

**Returns:** `StatisticalRigorResult` with FDR-corrected results, effect sizes, and CIs.

---

### `benjamini_hochberg()`

```python
def benjamini_hochberg(
    p_values: np.ndarray,
    alpha: float = 0.05,
) -> np.ndarray
```

Apply Benjamini-Hochberg FDR correction.

**Returns:** Array of q-values (FDR-adjusted p-values).

---

### `compute_pathway_pvalues()`

```python
def compute_pathway_pvalues(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    n_permutations: int = 1000,
    seed: Optional[int] = None,
) -> Dict[str, float]
```

Compute permutation-based p-values for pathway associations using ANOVA F-statistic.

---

### `compute_cohens_d()`

```python
def compute_cohens_d(
    group1: np.ndarray,
    group2: np.ndarray,
) -> float
```

Compute Cohen's d effect size with pooled standard deviation.

---

### `compute_pathway_effect_sizes()`

```python
def compute_pathway_effect_sizes(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
) -> Dict[str, float]
```

Compute maximum pairwise Cohen's d for each pathway across all cluster pairs.

---

### `bootstrap_effect_size_ci()`

```python
def bootstrap_effect_size_ci(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathway: str,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    seed: Optional[int] = None,
) -> Tuple[float, float]
```

Compute bootstrap confidence interval for effect size.

**Returns:** `(lower_bound, upper_bound)` tuple.

---

### `compute_variant_weight()`

```python
def compute_variant_weight(
    consequence: str,
    cadd_score: float,
    weights: BurdenWeights,
) -> float
```

Compute burden weight for a variant using a literature-based scheme.

---

### `aggregate_pathway_scores()`

```python
def aggregate_pathway_scores(
    gene_burdens: pd.DataFrame,
    pathway_genes: List[str],
    method: PathwayNormalization = PathwayNormalization.MEAN,
) -> pd.Series
```

Aggregate gene burdens into a pathway score with normalization.

---

## Dataclasses

### `StatisticalRigorResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `pathway_results` | List[FDRResult] | FDR-corrected results per pathway |
| `n_pathways_tested` | int | Total pathways tested |
| `n_significant_pathways` | int | Pathways with q < alpha |
| `fdr_alpha` | float | Significance threshold |
| `weight_scheme` | str | Burden weight scheme used |
| `weight_citations` | List[str] | Literature citations for weights |
| `normalization_method` | str | Aggregation method used |
| `n_permutations` | int | Permutations for p-values |
| `n_bootstrap` | int | Bootstrap iterations for CIs |

**Methods:** `to_dict()`, `format_report()`, `get_significant_pathways()`

### `FDRResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `pathway` | str | Pathway name |
| `p_value` | float | Raw permutation p-value |
| `q_value` | float | FDR-adjusted q-value |
| `significant` | bool | Whether q < alpha |
| `effect_size` | float | Cohen's d |
| `ci_lower` | float | 95% CI lower bound |
| `ci_upper` | float | 95% CI upper bound |

**Methods:** `to_dict()`

### `BurdenWeights`

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `scheme` | BurdenWeightScheme | `DEFAULT` | Weighting scheme |
| `lof_weight` | float | 1.0 | Loss-of-function weight |
| `missense_high_weight` | float | 0.8 | High-impact missense |
| `missense_moderate_weight` | float | 0.5 | Moderate-impact missense |
| `missense_low_weight` | float | 0.2 | Low-impact missense |
| `other_weight` | float | 0.1 | Other variant types |
| `cadd_high_threshold` | float | 25.0 | CADD "high impact" cutoff |
| `cadd_moderate_threshold` | float | 20.0 | CADD "moderate" cutoff |

**Methods:** `get_citations()`, `from_scheme(scheme)`

---

## Weight Scheme Comparison

| Scheme | LoF | High Missense | Moderate Missense | Low Missense | Other |
|--------|-----|---------------|-------------------|-------------|-------|
| DEFAULT | 1.0 | 0.5 | 0.1 | 0.1 | 0.1 |
| GNOMAD_CONSTRAINT | 1.0 | 0.8 | 0.4 | 0.1 | 0.05 |
| ACMG_INSPIRED | 1.0 | 0.9 | 0.5 | 0.1 | 0.05 |
| CADD_PERCENTILE | 1.0 | 0.7 | 0.4 | 0.2 | 0.1 |
| UNIFORM | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |

---

## See Also

- [Characterization](characterization.md) — subtype profiling using these metrics
- [Benchmark](benchmark.md) — method comparison
- [METHODS.md](../METHODS.md) — full statistical methodology
