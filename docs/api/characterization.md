# Subtype Characterization API

> **Module**: `pathway_subtyping.characterization`

Interprets and describes discovered subtypes through pathway enrichment analysis (Kruskal-Wallis + FDR correction), gene-level contribution scores (Cohen's d effect sizes), publication-quality heatmaps, and export to CSV/Excel.

---

## Quick Example

```python
from pathway_subtyping import characterize_subtypes

result = characterize_subtypes(
    pathway_scores=scores_df,
    cluster_labels=labels,
    gene_burdens=burden_df,
    pathways=pathway_dict,
    confidence_scores=gmm_proba,
    fdr_alpha=0.05,
    top_n_genes=20,
    seed=42,
)

print(result.format_report())
```

---

## Functions

### `characterize_subtypes()`

```python
def characterize_subtypes(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    gene_burdens: Optional[pd.DataFrame] = None,
    pathways: Optional[Dict[str, List[str]]] = None,
    cluster_names: Optional[Dict[int, str]] = None,
    confidence_scores: Optional[np.ndarray] = None,
    fdr_alpha: float = 0.05,
    top_n_genes: int = 20,
    seed: Optional[int] = None,
) -> CharacterizationResult
```

Main entry point for subtype characterization.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `cluster_labels` | np.ndarray | required | Cluster assignments (one per sample) |
| `gene_burdens` | pd.DataFrame or None | None | Gene burdens (samples x genes) for gene-level analysis |
| `pathways` | Dict or None | None | Pathway name to gene list mapping (required for gene analysis) |
| `cluster_names` | Dict[int, str] or None | None | Human-readable names per cluster ID |
| `confidence_scores` | np.ndarray or None | None | Assignment confidence (e.g., GMM posterior probability) |
| `fdr_alpha` | float | 0.05 | FDR significance threshold |
| `top_n_genes` | int | 20 | Top genes per subtype |
| `seed` | int or None | None | Random seed (reserved for future use) |

**Returns:** `CharacterizationResult` with complete subtype profiles.

**Raises:** `ValueError` if `pathway_scores` and `cluster_labels` have mismatched lengths.

---

### `pathway_enrichment_analysis()`

```python
def pathway_enrichment_analysis(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    fdr_alpha: float = 0.05,
) -> Dict[int, List[PathwayEnrichment]]
```

Run pathway enrichment analysis across all subtypes.

For each pathway, performs a Kruskal-Wallis test across all subtypes, then computes Cohen's d effect size for each subtype vs the rest. P-values are FDR-corrected per subtype using Benjamini-Hochberg.

**Returns:** Dict mapping subtype ID to list of `PathwayEnrichment` results (sorted by absolute effect size descending).

---

### `gene_contribution_scores()`

```python
def gene_contribution_scores(
    gene_burdens: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathways: Dict[str, List[str]],
    top_n: int = 20,
) -> Dict[int, List[GeneContribution]]
```

Identify top genes driving each subtype's pathway profile. For each gene in each pathway, computes Cohen's d effect size comparing the subtype to the rest.

**Returns:** Dict mapping subtype ID to top N `GeneContribution` results.

---

### `generate_subtype_heatmap()`

```python
def generate_subtype_heatmap(
    result: CharacterizationResult,
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 8),
) -> Optional[matplotlib.Figure]
```

Generate a subtypes x pathways heatmap of mean pathway Z-scores. Uses a diverging RdBu_r colormap centered at 0.

**Returns:** matplotlib Figure, or None if matplotlib is unavailable.

---

### `generate_gene_heatmap()`

```python
def generate_gene_heatmap(
    result: CharacterizationResult,
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (14, 10),
    top_n: int = 15,
) -> Optional[matplotlib.Figure]
```

Generate a subtypes x top-genes heatmap of effect sizes (Cohen's d).

---

### `export_characterization()`

```python
def export_characterization(
    result: CharacterizationResult,
    output_dir: str,
    formats: Optional[List[str]] = None,
) -> List[str]
```

Export characterization results to CSV and/or Excel.

**Generated files:**
- `subtype_summary.{csv,xlsx}` — one row per subtype
- `pathway_enrichment.{csv,xlsx}` — all enrichment results
- `gene_contributions.{csv,xlsx}` — top gene contributions
- `pathway_scores_matrix.{csv,xlsx}` — subtypes x pathways mean scores

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `result` | CharacterizationResult | required | From `characterize_subtypes()` |
| `output_dir` | str | required | Directory to write output files |
| `formats` | List[str] or None | `["csv"]` | Formats to export (`"csv"`, `"excel"`) |

**Returns:** List of file paths written.

---

## Dataclasses

### `CharacterizationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `subtype_profiles` | List[SubtypeProfile] | One profile per discovered subtype |
| `n_subtypes` | int | Number of subtypes |
| `n_samples` | int | Total number of samples |
| `n_pathways` | int | Number of pathways analyzed |
| `n_genes` | int | Number of genes analyzed (0 if gene burdens not provided) |
| `fdr_alpha` | float | FDR significance threshold used |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `SubtypeProfile`

| Attribute | Type | Description |
|-----------|------|-------------|
| `subtype_id` | int | Cluster ID |
| `subtype_label` | str | Human-readable label |
| `n_samples` | int | Samples in this subtype |
| `fraction` | float | Fraction of total samples |
| `mean_confidence` | float | Mean GMM assignment confidence |
| `enriched_pathways` | List[PathwayEnrichment] | Sorted by effect size |
| `top_genes` | List[GeneContribution] | Top contributing genes |
| `pathway_score_means` | Dict[str, float] | Mean score per pathway |

**Methods:** `to_dict()`

### `PathwayEnrichment`

| Attribute | Type | Description |
|-----------|------|-------------|
| `pathway` | str | Pathway name |
| `mean_score` | float | Mean score in this subtype |
| `overall_mean` | float | Mean score across all samples |
| `fold_change` | float | Subtype mean / overall mean |
| `effect_size` | float | Cohen's d (subtype vs rest) |
| `p_value` | float | Kruskal-Wallis p-value |
| `q_value` | float | FDR-corrected p-value |
| `significant` | bool | Whether q_value < alpha |

**Methods:** `to_dict()`

### `GeneContribution`

| Attribute | Type | Description |
|-----------|------|-------------|
| `gene` | str | Gene name |
| `pathway` | str | Pathway the gene belongs to |
| `mean_burden` | float | Mean gene burden in this subtype |
| `overall_mean` | float | Mean gene burden across all samples |
| `fold_change` | float | Subtype mean / overall mean |
| `effect_size` | float | Cohen's d (subtype vs rest) |

**Methods:** `to_dict()`

---

## Interpretation

| Metric | Good | Concerning |
|--------|------|------------|
| `effect_size` (Cohen's d) | > 0.8 (large) | < 0.2 (negligible) |
| `q_value` | < 0.05 | > 0.1 |
| `fold_change` | > 1.5 or < 0.5 | ~1.0 (no difference) |
| `mean_confidence` | > 0.8 | < 0.5 |

---

## See Also

- [Clustering](clustering.md) — subtype discovery
- [Statistical Rigor](statistical_rigor.md) — FDR correction, effect sizes
- [Visualization](visualization.md) — interactive plots
