# Single-Cell Pathway Scoring API

> **Requires:** `pip install pathway-subtyping[sc]` (installs `anndata>=0.9.0`)

The single-cell module extends the framework to scRNA-seq data, enabling pathway-based subtyping at the cell-type level. Pseudobulk methods reuse `expression.py` internals (ssGSEA, GSVA, mean-Z), so scoring is methodologically identical to bulk RNA-seq.

---

## Enums

### `SingleCellScoringMethod`

```python
from pathway_subtyping import SingleCellScoringMethod

SingleCellScoringMethod.MEAN_Z              # Per-cell mean-Z (fast, chunked)
SingleCellScoringMethod.PSEUDOBULK_MEAN_Z   # Pseudobulk mean-Z
SingleCellScoringMethod.PSEUDOBULK_SSGSEA   # Pseudobulk ssGSEA (recommended)
SingleCellScoringMethod.PSEUDOBULK_GSVA     # Pseudobulk GSVA
```

**String conversion:**
```python
method = SingleCellScoringMethod.from_string("pseudobulk_ssgsea")
```

### `SingleCellInputType`

```python
from pathway_subtyping import SingleCellInputType

SingleCellInputType.RAW_COUNTS       # Raw UMI counts (will be log-normalized)
SingleCellInputType.LOG_NORMALIZED   # Already log-normalized
SingleCellInputType.H5AD             # AnnData h5ad file (auto-detected)
```

---

## Dataclasses

### `SingleCellQualityReport`

Quality control report returned by `load_single_cell_data()`.

| Field | Type | Description |
|-------|------|-------------|
| `n_cells` | int | Total cells loaded |
| `n_genes` | int | Total genes |
| `n_cell_types` | int | Number of unique cell types |
| `cell_type_counts` | Dict[str, int] | Cells per cell type |
| `sparsity` | float | Fraction of zeros in the matrix |
| `median_genes_per_cell` | float | Median detected genes per cell |
| `median_umi_per_cell` | float | Median UMI counts per cell |
| `n_pathways_with_coverage` | int | Pathways with at least 1 gene in data |
| `n_pathways_total` | int | Total pathways provided |
| `excluded_cell_types` | List[str] | Types excluded due to min_cells filter |
| `warnings` | List[str] | QC warnings |
| `is_usable` | bool | Whether data passes minimum QC |

**Methods:** `.to_dict()`

### `SingleCellScoringResult`

Result returned by `score_single_cell_pathways()`.

| Field | Type | Description |
|-------|------|-------------|
| `pathway_scores` | pd.DataFrame | Cell types x pathways (Z-normalized) |
| `per_cell_scores` | pd.DataFrame or None | Cells x pathways (only for MEAN_Z) |
| `pseudobulk_expression` | pd.DataFrame or None | Cell types x genes (mean expression) |
| `method` | SingleCellScoringMethod | Method used |
| `quality_report` | SingleCellQualityReport | QC report |
| `n_cells_per_type` | Dict[str, int] | Cell counts per type |

**Methods:** `.to_dict()`, `.format_report()`, `.get_citations()`

**Pipeline compatibility:**
```python
# Directly compatible with run_clustering, ValidationGates, characterize_subtypes
from pathway_subtyping import run_clustering, characterize_subtypes

labels = run_clustering(result.pathway_scores.values, n_clusters=3, seed=42).labels
characterize_subtypes(result.pathway_scores, labels)
```

---

## Functions

### `load_single_cell_data()`

Load and validate single-cell data from h5ad or CSV/TSV.

```python
from pathway_subtyping import load_single_cell_data, SingleCellInputType

adata, qc_report = load_single_cell_data(
    path="data/pbmc_3k.h5ad",
    cell_type_column="cell_type",      # Column in adata.obs
    input_type=SingleCellInputType.H5AD,  # Auto-detected for .h5ad
    min_genes_per_cell=200,            # Filter cells with few genes
    min_cells_per_gene=3,              # Filter rarely-detected genes
    pathways=pathway_dict,             # Optional: for coverage QC
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | str | required | Path to h5ad or CSV/TSV file |
| `cell_type_column` | str | required | Column name for cell type labels |
| `input_type` | SingleCellInputType | `H5AD` | Input data format |
| `min_genes_per_cell` | int | 200 | Minimum genes for a cell to pass QC |
| `min_cells_per_gene` | int | 3 | Minimum cells for a gene to be retained |
| `pathways` | Dict[str, List[str]] | None | Pathway dict for coverage reporting |

**Returns:** `Tuple[anndata.AnnData, SingleCellQualityReport]`

**Notes:**
- Raw counts are auto-normalized via `log1p(x / total * 10000)`
- Sparse matrices (CSR/CSC) are preserved throughout
- CSV/TSV files are loaded with genes as columns, cells as rows

### `score_single_cell_pathways()`

Main entry point for single-cell pathway scoring.

```python
from pathway_subtyping import score_single_cell_pathways, SingleCellScoringMethod

result = score_single_cell_pathways(
    adata=adata,
    pathways=pathway_dict,
    cell_type_column="cell_type",
    method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
    min_cells_per_type=10,    # Exclude rare cell types
    chunk_size=5000,          # For per-cell scoring memory efficiency
    seed=42,
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | anndata.AnnData | required | Loaded AnnData object |
| `pathways` | Dict[str, List[str]] | required | Pathway gene sets |
| `cell_type_column` | str | required | Column in adata.obs |
| `method` | SingleCellScoringMethod | `PSEUDOBULK_SSGSEA` | Scoring method |
| `min_cells_per_type` | int | 10 | Minimum cells per type |
| `chunk_size` | int | 5000 | Cells per chunk (MEAN_Z only) |
| `seed` | int or None | None | Random seed for reproducibility |

**Returns:** `SingleCellScoringResult`

---

## Scoring Methods in Detail

### Pseudobulk Methods (Recommended)

Pseudobulk methods aggregate single-cell expression per cell type, then score the resulting cell_types x genes matrix using bulk expression scoring methods:

1. Compute mean expression per cell type (sparse-aware)
2. Pass pseudobulk matrix to `expression._score_mean_z()`, `_score_ssgsea()`, or `_score_gsva()`
3. Z-normalize output

This approach is statistically sound because the pseudobulk matrix is structurally identical to a bulk expression matrix.

### Per-Cell Mean-Z

The `MEAN_Z` method scores each individual cell:

1. For each cell, extract expression for pathway genes
2. Z-normalize across genes
3. Mean of Z-scores = pathway score for that cell
4. Aggregate per-cell scores to cell-type level (mean)
5. Z-normalize cell-type scores

Per-cell scoring is memory-efficient via chunked processing (default 5000 cells per chunk) and works with sparse matrices.

**Note:** Per-cell ssGSEA is not offered because it would be prohibitively slow at 50K+ cells.

---

## End-to-End Example

```python
from pathway_subtyping import (
    load_single_cell_data,
    score_single_cell_pathways,
    SingleCellScoringMethod,
    run_clustering,
    ValidationGates,
    characterize_subtypes,
)

# 1. Load data
adata, qc = load_single_cell_data("pbmc.h5ad", cell_type_column="cell_type")

# 2. Score pathways
result = score_single_cell_pathways(
    adata, pathways, cell_type_column="cell_type",
    method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA, seed=42,
)

# 3. Cluster cell types into subtypes
clustering = run_clustering(result.pathway_scores.values, n_clusters=3, seed=42)

# 4. Validate
gates = ValidationGates(seed=42)
validation = gates.run_all(
    pathway_scores=result.pathway_scores,
    cluster_labels=clustering.labels,
    pathways=pathways,
    gene_burdens=result.pseudobulk_expression,
    n_clusters=3, gmm_seed=42,
)

# 5. Characterize
profiles = characterize_subtypes(result.pathway_scores, clustering.labels)
```

---

## See Also

- [Expression Scoring](../framework_overview.md) - Bulk RNA-seq scoring (reused by pseudobulk)
- [Clustering](../api/index.md) - Clustering API
- [Data Formats](../data_formats.md) - h5ad and CSV format specifications

---

*Last updated: February 2026*
