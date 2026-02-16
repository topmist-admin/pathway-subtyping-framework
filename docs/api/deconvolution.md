# Bulk Deconvolution API

> **Module**: `pathway_subtyping.deconvolution`

Estimates cell-type proportions from bulk RNA-seq data using a single-cell reference expression profile. Uses non-negative least squares (NNLS) to solve: `bulk ≈ reference.T @ proportions`. The estimated proportions can be combined with pathway scores for cell-type-aware subtype discovery.

---

## Quick Example

```python
from pathway_subtyping import (
    build_reference_profile,
    deconvolve_bulk,
    combine_features,
)

# Build reference from single-cell data
reference = build_reference_profile(sc_expression, cell_type_labels)

# Estimate cell-type proportions from bulk RNA-seq
result = deconvolve_bulk(bulk_expression, reference, seed=42)
print(result.format_report())

# Combine with pathway scores for clustering
combined = combine_features(
    pathway_scores=pathway_scores_df,
    cell_type_proportions=result.cell_type_proportions,
    proportion_weight=0.5,
)
# combined: samples x (pathways + CELLTYPE:* columns)
```

---

## Enums

### `DeconvolutionMethod`

| Value | Description |
|-------|-------------|
| `NNLS` | Non-negative least squares (scipy.optimize.nnls). Default and recommended. |

---

## Functions

### `build_reference_profile()`

```python
def build_reference_profile(
    reference_expr: pd.DataFrame,
    cell_type_labels: np.ndarray,
    min_cells_per_type: int = 5,
) -> pd.DataFrame
```

Aggregate single-cell expression into cell-type mean profiles.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reference_expr` | pd.DataFrame | required | Single-cell expression matrix (cells x genes) |
| `cell_type_labels` | np.ndarray | required | Cell type label per cell |
| `min_cells_per_type` | int | 5 | Exclude cell types with fewer cells |

**Returns:** `pd.DataFrame` — cell_types x genes mean expression profile.

**Raises:** `ValueError` if no cell types pass the minimum cell count filter.

---

### `deconvolve_bulk()`

```python
def deconvolve_bulk(
    bulk_expr: pd.DataFrame,
    reference_profile: pd.DataFrame,
    method: DeconvolutionMethod = DeconvolutionMethod.NNLS,
    min_genes: int = 50,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> DeconvolutionResult
```

Estimate cell-type proportions from bulk RNA-seq using a reference profile.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bulk_expr` | pd.DataFrame | required | Bulk expression matrix (samples x genes) |
| `reference_profile` | pd.DataFrame | required | Reference profile from `build_reference_profile()` |
| `method` | DeconvolutionMethod | `NNLS` | Deconvolution algorithm |
| `min_genes` | int | 50 | Minimum shared genes required |
| `seed` | int or None | None | Random seed |
| `show_progress` | bool | True | Show tqdm progress bar |

**Returns:** `DeconvolutionResult` with cell-type proportion matrix.

**Raises:** `ValueError` if shared genes < `min_genes` or inputs are empty.

---

### `combine_features()`

```python
def combine_features(
    pathway_scores: pd.DataFrame,
    cell_type_proportions: pd.DataFrame,
    proportion_weight: float = 0.5,
) -> pd.DataFrame
```

Merge pathway scores and cell-type proportions into a unified feature matrix for clustering.

**Key behaviors:**
- Z-normalizes both inputs independently
- Scales proportions by `proportion_weight` and pathway scores by `(1 - proportion_weight)`
- Concatenates column-wise with `CELLTYPE:{name}_proportion` prefix
- Aligns on shared sample IDs

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `cell_type_proportions` | pd.DataFrame | required | Proportions from `deconvolve_bulk()` |
| `proportion_weight` | float | 0.5 | Weight for cell-type features (0 = pathway only, 1 = proportions only) |

**Returns:** `pd.DataFrame` — combined features ready for clustering.

**Raises:** `ValueError` if no shared samples or invalid weight.

---

### `generate_synthetic_bulk()`

```python
def generate_synthetic_bulk(
    reference_profile: pd.DataFrame,
    n_samples: int = 100,
    n_subtypes: int = 3,
    noise_level: float = 0.1,
    seed: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]
```

Generate synthetic bulk expression from known cell-type proportions for testing.

**Returns:** `(bulk_expression, true_proportions, subtype_labels)` tuple.

Each subtype has a distinct cell-type composition pattern (e.g., subtype 0 = high Neuron, subtype 1 = high Astrocyte).

---

## Dataclasses

### `DeconvolutionResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `cell_type_proportions` | pd.DataFrame | Samples x cell_types (values in [0,1], rows sum to 1) |
| `method` | DeconvolutionMethod | Algorithm used |
| `quality_report` | DeconvolutionQualityReport | Quality metrics |
| `reference_cell_types` | List[str] | Cell type names |
| `n_samples` | int | Number of deconvolved samples |
| `n_cell_types` | int | Number of cell types |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `DeconvolutionQualityReport`

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_samples` | int | Number of bulk samples |
| `n_genes_bulk` | int | Genes in bulk data |
| `n_genes_reference` | int | Genes in reference |
| `n_genes_shared` | int | Genes in common |
| `n_cell_types` | int | Number of cell types |
| `gene_coverage` | float | shared / reference genes |
| `proportion_stats` | Dict | Per-cell-type min/max/mean proportions |
| `warnings` | List[str] | Quality warnings |
| `is_usable` | bool | Whether results are reliable |

**Methods:** `to_dict()`

---

## Algorithm Details

### NNLS Deconvolution

For each bulk sample, solve:

```
minimize ||reference.T @ proportions - bulk_vector||^2
subject to: proportions >= 0
```

Post-processing: normalize proportions to sum to 1 (if sum > 0; else uniform).

NNLS is biologically appropriate because cell-type proportions cannot be negative. This approach is used in CIBERSORT (Newman et al., 2015) and MuSiC (Wang et al., 2019).

---

## End-to-End Workflow

```python
from pathway_subtyping import (
    load_expression_matrix, score_pathways_from_expression,
    build_reference_profile, deconvolve_bulk, combine_features,
    run_clustering,
)

# 1. Load bulk expression and score pathways
bulk_expr, _ = load_expression_matrix("bulk_rnaseq.csv")
pathway_result = score_pathways_from_expression(bulk_expr, pathways, seed=42)

# 2. Build reference from single-cell data
reference = build_reference_profile(sc_expression, cell_type_labels)

# 3. Deconvolve bulk to get cell-type proportions
deconv = deconvolve_bulk(bulk_expr, reference, seed=42)

# 4. Combine pathway scores + proportions
combined = combine_features(
    pathway_result.pathway_scores,
    deconv.cell_type_proportions,
    proportion_weight=0.3,  # 30% cell-type, 70% pathway
)

# 5. Cluster on combined features
clustering = run_clustering(combined.values, n_clusters=3, seed=42)
```

---

## See Also

- [Expression Scoring](expression.md) — bulk RNA-seq pathway scoring
- [Single-Cell Scoring](single_cell.md) — building the reference
- [Multi-Omic Fusion](multi_omic.md) — integrating deconvolution as a modality
- [Cross-Modal Validation](cross_modal_validation.md) — validating multi-omic subtypes
