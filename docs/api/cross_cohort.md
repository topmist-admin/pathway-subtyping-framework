# cross_cohort — Cross-Cohort Validation

Compare subtype definitions across independent cohorts to validate reproducibility.

## Classes

### `CohortResult`

Results from a single cohort analysis.

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Cohort identifier |
| `pathway_scores` | `pd.DataFrame` | Pathway scores matrix (samples x pathways) |
| `cluster_labels` | `np.ndarray` | Cluster assignments per sample |
| `cluster_names` | `Dict[int, str]` | Mapping of cluster IDs to names |
| `n_samples` | `int` | Number of samples |
| `n_clusters` | `int` | Number of clusters |

**Methods:**

- `to_dict() -> Dict[str, Any]` — Serialize to dictionary (pathway_scores converted to dict, cluster_labels to list)

### `CrossCohortResult`

Results from cross-cohort validation.

| Attribute | Type | Description |
|-----------|------|-------------|
| `cohort_a` | `str` | Name of reference cohort |
| `cohort_b` | `str` | Name of validation cohort |
| `transfer_ari` | `float` | ARI from transfer learning validation |
| `projection_ari` | `float` | ARI from projection validation |
| `pathway_correlation` | `float` | Pearson correlation of pathway variances |
| `shared_subtypes` | `List[str]` | Subtype labels found in both cohorts |
| `details` | `Dict[str, Any]` | Additional metadata (common pathways, sample counts) |

**Methods:**

- `to_dict() -> Dict[str, Any]` — Serialize all fields to dictionary
- `format_report() -> str` — Human-readable summary with interpretation guidance
- `get_citations() -> List[str]` — Return relevant citations (Hubert & Arabie 1985)

## Functions

### `compare_cohorts(cohort_a, cohort_b, seed=42) -> CrossCohortResult`

Compare subtype definitions between two cohorts.

Performs:
1. **Transfer Learning**: Train GMM on cohort A, predict on cohort B
2. **Projection**: Project both cohorts into shared PCA space
3. **Pathway Correlation**: Correlate pathway variances

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cohort_a` | `CohortResult` | required | Reference cohort |
| `cohort_b` | `CohortResult` | required | Validation cohort |
| `seed` | `int` | `42` | Random seed |

**Returns:** `CrossCohortResult`

**Raises:** `ValueError` if fewer than 2 common pathways exist.

```python
result = compare_cohorts(cohort_a, cohort_b)
print(f"Transfer ARI: {result.transfer_ari:.3f}")
```

### `batch_compare_cohorts(output_dirs, report_path=None, seed=42) -> List[CrossCohortResult]`

Compare all pairs of cohorts from pipeline output directories.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_dirs` | `List[str]` | required | Pipeline output directories |
| `report_path` | `Optional[str]` | `None` | Path to save markdown report |
| `seed` | `int` | `42` | Random seed |

**Returns:** List of `CrossCohortResult` for all pairwise comparisons.

**Raises:** `ValueError` if fewer than 2 cohorts load successfully.

```python
results = batch_compare_cohorts(
    ["outputs/cohort_a/", "outputs/cohort_b/", "outputs/cohort_c/"],
    report_path="cross_cohort_report.md",
)
```

### `load_cohort_result(output_dir) -> CohortResult`

Load results from a completed pipeline run.

Expects the directory to contain:
- `pathway_scores.csv` — pathway score matrix
- `subtype_assignments.csv` — cluster assignments with `cluster_id` and `cluster_label` columns
- `report.json` (optional) — pipeline metadata

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `output_dir` | `str` | Path to pipeline output directory |

**Returns:** `CohortResult`

**Raises:** `FileNotFoundError` if required files are missing.

### `generate_cross_cohort_report(results, output_path) -> None`

Generate a markdown cross-cohort validation report.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `results` | `List[CrossCohortResult]` | Comparison results to include |
| `output_path` | `str` | Path to save the markdown report |

### `generate_synthetic_cohort_pair(...) -> Tuple[CohortResult, CohortResult]`

Generate a pair of synthetic cohorts for testing and demos.

Uses the simulation framework to create two cohorts with matching subtype structure.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_samples_a` | `int` | `100` | Samples in cohort A |
| `n_samples_b` | `int` | `80` | Samples in cohort B |
| `n_subtypes` | `int` | `2` | Number of planted subtypes |
| `n_pathways` | `int` | `10` | Number of pathways |
| `effect_size` | `float` | `1.5` | Cohen's d for subtype differences |
| `seed` | `Optional[int]` | `None` | Random seed |

**Returns:** Tuple of two `CohortResult` objects.

```python
from pathway_subtyping import generate_synthetic_cohort_pair, compare_cohorts

a, b = generate_synthetic_cohort_pair(seed=42)
result = compare_cohorts(a, b)
print(result.format_report())
```
