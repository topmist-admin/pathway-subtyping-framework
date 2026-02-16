# API Reference

This document provides comprehensive API documentation for the Pathway Subtyping Framework.

## Modules Overview

| Module | Description |
|--------|-------------|
| [`pipeline`](pipeline.md) | Main pipeline orchestrator and configuration |
| [`expression`](expression.md) | Bulk RNA-seq pathway scoring (ssGSEA, GSVA, mean-Z) |
| [`single_cell`](single_cell.md) | Single-cell scRNA-seq pathway scoring (pseudobulk + per-cell) |
| [`multi_omic`](multi_omic.md) | Multi-omic pathway score fusion (concatenate, weighted, intersection) |
| [`deconvolution`](deconvolution.md) | Bulk deconvolution for cell-type proportion estimation (NNLS) |
| [`cross_modal_validation`](cross_modal_validation.md) | Cross-modal validation gate (Gate 5) for multi-omic analyses |
| [`visualization`](visualization.md) | Interactive Plotly reports, UMAP/t-SNE, multi-format export |
| [`variant_qc`](variant_qc.md) | Variant quality control filters (QUAL, HWE, MAF, call rate) |
| [`clustering`](clustering.md) | Clustering algorithms (GMM, K-means, Hierarchical, Spectral) and model selection |
| [`characterization`](characterization.md) | Subtype characterization, pathway enrichment, gene contributions |
| [`statistical_rigor`](statistical_rigor.md) | FDR correction, permutation tests, effect sizes, burden weighting |
| [`benchmark`](benchmark.md) | Method comparison (Pathway GMM vs NMF, PCA+K-means, gene-level) |
| [`simulation`](simulation.md) | Synthetic data generation, power analysis, sample size estimation |
| [`sensitivity`](sensitivity.md) | Parameter sensitivity analysis and robustness testing |
| [`ancestry`](ancestry.md) | Ancestry PCA, population stratification correction, independence testing |
| [`batch_correction`](batch_correction.md) | Batch effect detection, ComBat/mean-center/standardize correction |
| [`validation`](validation.md) | Validation gates and stability testing |
| [`threshold_calibration`](threshold_calibration.md) | Data-driven validation threshold calibration |
| [`cross_cohort`](cross_cohort.md) | Cross-cohort validation and replication |
| [`config`](config.md) | Configuration loading and validation utilities |
| [`cli`](cli.md) | Command-line interface |

## Quick Start

### Running the Pipeline Programmatically

```python
from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

# Load configuration
config = PipelineConfig.from_yaml("configs/my_config.yaml")

# Run pipeline
pipeline = DemoPipeline(config)
pipeline.run()

# Access results
print(f"Found {pipeline.n_clusters} subtypes")
print(pipeline.cluster_assignments.head())
```

### Running Validation Only

```python
from pathway_subtyping.validation import ValidationGates

# Initialize validator
validator = ValidationGates(
    seed=42,
    n_permutations=100,
    n_bootstrap=50,
    stability_threshold=0.8,
    null_ari_max=0.15
)

# Run validation
result = validator.run_all(
    pathway_scores=pathway_scores_df,
    cluster_labels=labels,
    pathways=pathway_dict,
    gene_burdens=burden_df,
    n_clusters=4,
    gmm_seed=42
)

print(f"All passed: {result.all_passed}")
print(result.summary)
```

### Custom Configuration

```python
from pathway_subtyping.pipeline import PipelineConfig

# Create config programmatically
config = PipelineConfig(
    name="my_analysis",
    output_dir="outputs/my_analysis",
    seed=42,
    vcf_path="data/my_cohort.vcf",
    phenotype_path="data/my_phenotypes.csv",
    pathway_db="data/pathways/my_pathways.gmt",
    n_clusters_range=[2, 6],  # Test 2-6 clusters
)
```

## Guides

| Guide | Description |
|-------|-------------|
| [Performance & Hardware](../guides/performance-and-hardware.md) | Hardware recommendations, memory estimation, chunked processing, benchmarking |
| [Cross-Cohort Validation](../guides/cross-cohort-validation.md) | Comparing subtypes across independent cohorts |
| [Validation Gates](../guides/validation-gates.md) | Understanding and configuring validation gates |

## Expression Pathway Scoring

The expression module computes pathway scores from bulk RNA-seq gene expression matrices, producing the same `pathway_scores` DataFrame format as variant-based scoring.

```python
from pathway_subtyping import (
    load_expression_matrix, score_pathways_from_expression,
    ExpressionScoringMethod, ExpressionInputType,
)

# Load and score
expr, qc = load_expression_matrix("expression.csv", input_type=ExpressionInputType.TPM)
result = score_pathways_from_expression(expr, pathways, method=ExpressionScoringMethod.SSGSEA, seed=42)
# result.pathway_scores: samples x pathways, Z-normalized
```

| Method | Speed | Best For |
|--------|-------|----------|
| `MEAN_Z` | Fast | Quick exploration |
| `SSGSEA` | Medium | **Recommended default** |
| `GSVA` | Medium | Alternative to ssGSEA |

See [expression.md](expression.md) for full API reference.

---

## Single-Cell Pathway Scoring

The single-cell module (`pip install pathway-subtyping[sc]`) provides per-cell and pseudobulk pathway scoring from scRNA-seq data. Pseudobulk methods reuse expression.py internals, so ssGSEA/GSVA scoring is identical to bulk RNA-seq.

### Quick Example

```python
from pathway_subtyping import (
    load_single_cell_data,
    score_single_cell_pathways,
    SingleCellScoringMethod,
    run_clustering,
)

# Load h5ad file with cell type annotations
adata, qc_report = load_single_cell_data(
    "data/pbmc_3k.h5ad",
    cell_type_column="cell_type",
)
print(qc_report.format_report())  # QC summary

# Score pathways (pseudobulk ssGSEA — recommended)
result = score_single_cell_pathways(
    adata,
    pathways=pathway_dict,
    cell_type_column="cell_type",
    method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
    seed=42,
)

# Result is cell_types × pathways — directly compatible with clustering
clustering = run_clustering(result.pathway_scores.values, n_clusters=3, seed=42)
```

### Scoring Methods

| Method | Level | Speed | Best For |
|--------|-------|-------|----------|
| `MEAN_Z` | Per-cell | Fast | Large datasets, cell-level resolution |
| `PSEUDOBULK_MEAN_Z` | Cell-type | Fast | Quick exploration |
| `PSEUDOBULK_SSGSEA` | Cell-type | Medium | Recommended default |
| `PSEUDOBULK_GSVA` | Cell-type | Medium | Alternative to ssGSEA |

See [single_cell.md](single_cell.md) for full API reference.

---

## Multi-Omic Fusion

Fuse pathway scores from multiple modalities into a unified feature matrix:

```python
from pathway_subtyping import (
    ModalityType, FusionStrategy, prepare_modality, fuse_modalities,
)

vcf_mod = prepare_modality(ModalityType.VCF, vcf_scores, label="WES")
expr_mod = prepare_modality(ModalityType.EXPRESSION, expr_scores, label="RNA-seq")
result = fuse_modalities([vcf_mod, expr_mod], strategy=FusionStrategy.CONCATENATE, seed=42)
# result.fused_pathway_scores: unified feature matrix
# result.per_modality_scores: for cross-modal validation (Gate 5)
```

| Strategy | Description |
|----------|-------------|
| `CONCATENATE` | Column-bind with prefixes; preserves all info |
| `WEIGHTED_AVERAGE` | Weighted mean of shared pathways |
| `INTERSECTION_ONLY` | Restrict to shared samples and pathways |

See [multi_omic.md](multi_omic.md) for full API reference.

---

## Bulk Deconvolution

Estimate cell-type proportions from bulk RNA-seq using a single-cell reference:

```python
from pathway_subtyping import build_reference_profile, deconvolve_bulk, combine_features

reference = build_reference_profile(sc_expression, cell_type_labels)
result = deconvolve_bulk(bulk_expression, reference, seed=42)
combined = combine_features(pathway_scores, result.cell_type_proportions, proportion_weight=0.5)
```

See [deconvolution.md](deconvolution.md) for full API reference.

---

## Cross-Modal Validation (Gate 5)

Validates that subtypes are consistent across data modalities. Automatically included by `ValidationGates.run_all()` when `per_modality_scores` is provided:

```python
from pathway_subtyping import ValidationGates

gates = ValidationGates(seed=42, n_permutations=100)
result = gates.run_all(
    pathway_scores=fused_scores,
    cluster_labels=labels,
    pathways=pathways,
    gene_burdens=gene_data,
    n_clusters=3,
    per_modality_scores={"WES": vcf_scores, "RNA-seq": expr_scores},
)
# Gate 5 included automatically when >= 2 modalities
```

See [cross_modal_validation.md](cross_modal_validation.md) for full API reference.

---

## Visualization

The visualization module (`pip install pathway-subtyping[viz]`) provides interactive Plotly charts and multi-format figure export. Static matplotlib fallbacks work with the base install.

### Interactive HTML Report

```python
from pathway_subtyping import create_interactive_report, ReportConfig, DimReductionMethod

# Generate a self-contained HTML report
config = ReportConfig(
    title="My Cohort Analysis",
    dim_reduction=DimReductionMethod.UMAP,
    disclaimer="Research use only.",
)

result = create_interactive_report(
    pathway_scores=pathway_scores_df,
    labels=cluster_labels,
    output_path="outputs/report.html",
    config=config,
    seed=42,
)
# Open outputs/report.html in any browser
```

### Individual Charts

```python
from pathway_subtyping import (
    plot_interactive_scatter,
    plot_interactive_heatmap,
    plot_subtype_trajectories,
    plot_cluster_distribution,
    export_figure,
    FigureFormat,
    DimReductionMethod,
)

# UMAP scatter plot
fig = plot_interactive_scatter(
    pathway_scores_df, labels,
    method=DimReductionMethod.UMAP, seed=42,
)

# Export in multiple formats
export_figure(fig, "outputs/scatter", [FigureFormat.PNG, FigureFormat.SVG, FigureFormat.HTML])

# Subtype radar chart
fig = plot_subtype_trajectories(pathway_scores_df, labels, top_n=8)
export_figure(fig, "outputs/radar", [FigureFormat.HTML])
```

### Static Fallback (No Plotly Required)

```python
from pathway_subtyping import plot_static_scatter, DimReductionMethod

fig = plot_static_scatter(
    pathway_scores_df, labels,
    method=DimReductionMethod.PCA,
    output_path="outputs/scatter.png",
    seed=42,
)
```

See [visualization.md](visualization.md) for full API reference.

## Data Structures

### Input Formats

| Data Type | Format | Required Columns/Fields |
|-----------|--------|------------------------|
| VCF | `.vcf` | INFO fields: `GENE`, `CONSEQUENCE`, `CADD` |
| Phenotypes | `.csv` | `sample_id` (required), other columns optional |
| Pathways | `.gmt` | Standard GMT: `NAME\tDESCRIPTION\tGENE1\tGENE2...` |

### Output Formats

| File | Format | Description |
|------|--------|-------------|
| `pathway_scores.csv` | CSV | Z-normalized pathway burden scores per sample |
| `subtype_assignments.csv` | CSV | Cluster assignments with confidence scores |
| `report.json` | JSON | Machine-readable analysis report |
| `report.md` | Markdown | Human-readable analysis report |
| `figures/summary.png` | PNG | Static visualization of clustering results |
| `figures/interactive_report.html` | HTML | Interactive Plotly report (when `generate_interactive_report: true`) |
| `run_metadata.yaml` | YAML | Reproducibility metadata |

## Error Handling

All major operations raise descriptive exceptions:

```python
from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

try:
    config = PipelineConfig.from_yaml("missing.yaml")
except FileNotFoundError as e:
    print(f"Config file not found: {e}")

try:
    pipeline = DemoPipeline(config)
    pipeline.run()
except ValueError as e:
    print(f"Invalid configuration: {e}")
```

## Thread Safety

The pipeline is designed for single-threaded execution. For parallel processing of multiple cohorts, run separate pipeline instances:

```python
from concurrent.futures import ProcessPoolExecutor
from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

def run_cohort(config_path):
    config = PipelineConfig.from_yaml(config_path)
    pipeline = DemoPipeline(config)
    pipeline.run()
    return pipeline.cluster_assignments

# Run cohorts in parallel processes
config_files = ["cohort_a.yaml", "cohort_b.yaml", "cohort_c.yaml"]
with ProcessPoolExecutor(max_workers=3) as executor:
    results = list(executor.map(run_cohort, config_files))
```
