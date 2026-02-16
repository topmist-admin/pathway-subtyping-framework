# Expression Pathway Scoring API

> **Module**: `pathway_subtyping.expression`

Computes pathway-level scores from bulk RNA-seq gene expression matrices using ssGSEA, GSVA, or mean-Z methods. Produces the same `pathway_scores` DataFrame format as variant-based scoring, enabling unified downstream clustering and validation.

---

## Quick Example

```python
from pathway_subtyping import (
    load_expression_matrix,
    score_pathways_from_expression,
    ExpressionScoringMethod,
    ExpressionInputType,
)

# Load expression matrix (auto-detects orientation)
expr, qc = load_expression_matrix(
    "data/expression.csv",
    input_type=ExpressionInputType.TPM,
)

# Load pathways
import yaml
pathways = {}  # your GMT-loaded dict: pathway_name -> [gene1, gene2, ...]

# Score pathways (ssGSEA recommended)
result = score_pathways_from_expression(
    expr, pathways,
    method=ExpressionScoringMethod.SSGSEA,
    seed=42,
)

# result.pathway_scores is samples x pathways, Z-normalized
print(result.format_report())
```

---

## Enums

### `ExpressionScoringMethod`

| Value | Description | Speed | Best For |
|-------|-------------|-------|----------|
| `MEAN_Z` | Z-normalize genes, average per pathway | Fast | Quick exploration, large datasets |
| `SSGSEA` | Single-sample GSEA with rank-based enrichment | Medium | **Recommended default** |
| `GSVA` | Gene Set Variation Analysis (empirical CDF + KS) | Medium | Alternative to ssGSEA |

### `ExpressionInputType`

| Value | Description | Transformation |
|-------|-------------|---------------|
| `COUNTS` | Raw read counts | `log2(x + 1)` applied automatically |
| `TPM` | Transcripts per million | `log2(x + 1)` if max > 20, else no-op |
| `FPKM` | Fragments per kilobase per million | Same as TPM |
| `LOG2` | Already log2-transformed | No transformation |

---

## Functions

### `load_expression_matrix()`

```python
def load_expression_matrix(
    path: str,
    input_type: ExpressionInputType = ExpressionInputType.TPM,
    gene_column: Optional[str] = None,
    sample_column: Optional[str] = None,
    min_genes_per_sample: int = 100,
    min_samples_per_gene: int = 3,
) -> Tuple[pd.DataFrame, ExpressionDataQualityReport]
```

Load and validate an expression matrix from CSV/TSV.

**Key behaviors:**
- Auto-detects orientation (genes as rows vs columns) using gene symbol pattern matching
- Applies appropriate log transformation based on `input_type`
- Removes all-zero, low-occurrence, and zero-variance genes
- Returns samples x genes DataFrame

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | str | required | Path to CSV or TSV file |
| `input_type` | ExpressionInputType | `TPM` | Type of expression values |
| `gene_column` | str or None | None | Column containing gene symbols (if genes are rows) |
| `sample_column` | str or None | None | Column containing sample IDs |
| `min_genes_per_sample` | int | 100 | Minimum genes for a usable dataset |
| `min_samples_per_gene` | int | 3 | Minimum non-zero samples per gene |

**Returns:** `(gene_expression, quality_report)` tuple.

**Raises:**
- `FileNotFoundError` if file does not exist
- `ValueError` if data is empty or unusable

---

### `score_pathways_from_expression()`

```python
def score_pathways_from_expression(
    gene_expression: pd.DataFrame,
    pathways: Dict[str, List[str]],
    method: ExpressionScoringMethod = ExpressionScoringMethod.SSGSEA,
    min_genes_per_pathway: int = 2,
    alpha: float = 0.25,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> ExpressionScoringResult
```

Main entry point for expression-based pathway scoring.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_expression` | pd.DataFrame | required | Preprocessed samples x genes matrix |
| `pathways` | Dict[str, List[str]] | required | Pathway name to gene list mapping |
| `method` | ExpressionScoringMethod | `SSGSEA` | Scoring algorithm |
| `min_genes_per_pathway` | int | 2 | Skip pathways with fewer overlapping genes |
| `alpha` | float | 0.25 | ssGSEA rank weight exponent |
| `seed` | int or None | None | Random seed for reproducibility |
| `show_progress` | bool | True | Show tqdm progress bar |

**Returns:** `ExpressionScoringResult` with Z-normalized pathway scores.

---

## Dataclasses

### `ExpressionScoringResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `pathway_scores` | pd.DataFrame | Z-normalized scores (samples x pathways) |
| `gene_expression` | pd.DataFrame | Preprocessed expression matrix |
| `method` | ExpressionScoringMethod | Method used |
| `quality_report` | ExpressionDataQualityReport | Quality metrics |
| `n_pathways_scored` | int | Number of pathways successfully scored |
| `n_pathways_skipped` | int | Number of pathways skipped |
| `skipped_pathways` | List[str] | Names of skipped pathways |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `ExpressionDataQualityReport`

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_samples` | int | Number of samples |
| `n_genes` | int | Number of genes after filtering |
| `n_genes_before_filter` | int | Number of genes before filtering |
| `n_zero_genes` | int | Genes with all-zero expression |
| `n_low_variance_genes` | int | Zero-variance genes removed |
| `n_pathways_covered` | int | Pathways with sufficient gene overlap |
| `n_pathways_total` | int | Total pathways provided |
| `mean_pathway_gene_coverage` | float | Average fraction of pathway genes found |
| `input_type` | str | Input type used |
| `orientation_detected` | str | `"genes_as_rows"` or `"genes_as_columns"` |
| `was_transposed` | bool | Whether matrix was transposed |
| `warnings` | List[str] | Quality warnings |
| `is_usable` | bool | Whether data passes minimum thresholds |

**Methods:** `to_dict()`

---

## Scoring Method Details

### Mean-Z

1. Z-score normalize each gene across samples
2. For each pathway: compute mean of member gene Z-scores
3. Z-normalize resulting pathway scores

Fast and simple. Best for quick exploration.

### ssGSEA (Recommended)

1. Rank genes by expression within each sample
2. Walk along ranked list; step up by |rank|^alpha for pathway genes, step down for non-pathway genes
3. Enrichment score = sum of running statistic

Based on Barbie et al. (2009). Captures pathway enrichment without assuming normality.

### GSVA (Simplified)

1. Compute empirical CDF per gene across samples
2. Rank CDF values within each sample
3. KS-like statistic: max deviation between pathway and non-pathway gene distributions

Based on Hanzelmann et al. (2013). For publication-grade GSVA, precompute scores using the R GSVA package and import with `input_type=LOG2`.

---

## Integration with Pipeline

Expression scoring produces the same `pathway_scores` DataFrame format as VCF-based scoring. This enables:

```python
from pathway_subtyping import run_clustering, ValidationGates

# Cluster expression-based pathway scores
clustering = run_clustering(result.pathway_scores.values, n_clusters=3, seed=42)

# Validate
gates = ValidationGates(seed=42)
validation = gates.run_all(
    pathway_scores=result.pathway_scores,
    cluster_labels=clustering.labels,
    pathways=pathways,
    gene_burdens=result.gene_expression,
    n_clusters=3,
)
```

---

## See Also

- [Multi-Omic Fusion](multi_omic.md) — combine expression scores with VCF or single-cell
- [Single-Cell Scoring](single_cell.md) — scRNA-seq pathway scoring
- [Framework Overview](../framework_overview.md) — architecture and pipeline flow
