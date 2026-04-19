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
| [`signaling_databases`](signaling_databases.md) | Cell-cell signaling pathway databases (CellPhoneDB, CellChatDB) |
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

### Optional Layers (v0.5)

| Module | Extra | Description |
|--------|-------|-------------|
| [`qc`](qc.md) | `[qc]` | 12-feature molecular QC for cell manufacturing (cascade, dosage, drift, off-target, stress) |
| [`gnn`](gnn.md) | `[gnn]` | TransE/RotatE KG embeddings, OntologyAwareGNN, biological attention *(experimental)* |
| [`autism`](autism.md) | `[autism]` | Neuro-symbolic rules (R1-R7), therapeutic hypothesis ranking *(autism-only)* |

### v0.6 Layers — Rigor + Foundation-Model Interface

All v0.6 modules ship a deterministic fallback alongside the production backend, so CI and test suites run without heavyweight checkpoints. The production backend is gated on the extra + a locally-cached checkpoint (see each module's in-package docstring for checkpoint discovery). User-facing walkthroughs in the `Guide` column; notebook numbers in the `Notebook` column mirror [notebook-guide.md](../notebook-guide.md).

| Module | Feature | Extra | Guide | Notebook | What it does |
|--------|---------|-------|-------|----------|--------------|
| [`uncertainty`](../guides/uncertainty.md) | F1 | (core) | [uncertainty](../guides/uncertainty.md) | 21 | `ConformalPathwayPredictor`, `BootstrapMSV`, `BayesianPathwayGMM`, `CalibrationReport` — calibrated intervals on pathway scores with a distribution-free coverage guarantee |
| [`harmonize`](../guides/cross-platform.md) | F2 | `[harmonize]` | [cross-platform](../guides/cross-platform.md) | 22 | `UCEEmbedder`, `CrossPlatformAligner`, `HarmonizationReport`, `CrossPlatformBenchmark` — align pathway scores across 10x / Smart-seq2 / bulk / spatial on an embedding anchor |
| `knowledge_graph` (refresh) | F3 | (core) | [v0.5→v0.6 KG migration](../migration/v05-to-v06-kg.md) | 16 | `sources.py` (pinned OmniPath/SIGNOR/Reactome manifest), `diff.py`, `regression.py` — KG provenance + diff |
| `qc.alphamissense` | F4 | `[qc]` | *(in QC guide)* | — | AlphaMissense-aware variant cascade extension; `gene_weights=None` is bit-identical to v0.5 |
| [`perturb`](../guides/perturbation.md) | F5 | `[perturb]` | [perturbation](../guides/perturbation.md) | 23 | `GeneformerPerturber`, `MSVFromEmbedding`, `PerturbationScreen`, `PerturbationReport`, `FallbackPerturber`, `OfficialBackend` (Geneformer V2 104M; content-hashed embedding cache via `cache_dir=`) |
| [`embed`](../guides/embeddings.md) | F6, F8 | `[embed]` | [embeddings](../guides/embeddings.md) | 24, 26 | `Embedder` ABC, `scGPTEmbedder`, `NicheformerEmbedder`, `EmbeddingCache`, `embed_joint` for dissociated + spatial |
| [`genesets`](../guides/genesets.md) | F7 | `[genesets]` | *(see module docstring)* | 25 | `RegulatoryGeneSetExpander` with `BorzoiBackend` (opt-in) + `CoexpressionBackend` (fallback) — expand curated gene sets by regulatory co-occurrence |
| `qc.offtarget_sequence` | F9 | `[qc-sequence]` | *(in QC guide)* | 27 | `Evo2OffTargetScorer`, `SimulatedEvo2Backend`, `SimilarityBackend` — sequence-level CRISPR off-target scoring |
| [`omics`](../guides/multi-omics.md) | F10 | (core) | *(see module docstring)* | 28 | `ATACScorer`, `ProteomicsScorer`, `MultiOmicsFusion`, `FusionWeights`, `flag_discordant_pathways` — weighted fusion of per-modality pathway scores |
| [`causal`](../guides/causal.md) | F11 | (core) | *(see module docstring)* | 29 | `InvariantPathwayPredictor` — Invariant Causal Prediction (ICP) with combined mean+variance invariance test |
| [`active`](../guides/active.md) | F12 | (core) | *(see module docstring)* | 30 | `ActiveSampleSelector` — uncertainty / diversity / hybrid sample selection strategies |

### Real-Data Validation Scripts (v0.6.3)

Each script fetches a public cohort, computes the roadmap acceptance metric, and writes a skip-on-absent JSON artefact consumed by the corresponding `tests/test_*_real_data.py` suite. Dataset provenance + SHA-256 in [docs/provenance.md](../provenance.md).

| Script | Cohort | Asserts |
|--------|--------|---------|
| `scripts/validate_f1_real_data.py` | TCGA-COAD (n=57) + GSE28521 (n=79) | Conformal oracle deviation within ±2% at 80/90/95% target |
| `scripts/validate_f2_real_data.py` | GSE28521 × GSE80655 | Pathway-mean Spearman rho uplift ≥ 0.10 post-alignment |
| `scripts/validate_f5_real_data.py` | TCGA-COAD + GSE123753 (MECP2-KO, via `scripts/fetch_gse123753.py`) | Directional agreement ≥ 70% on curated edges; 50/50 on real WT-vs-KO cohort with Geneformer |
| `scripts/validate_f10_real_data.py` | 10x `pbmc_1k_protein_v3` CITE-seq | Fused-vs-RNA-only 1-NN classification uplift ≥ 3% with CI_low > 0 |

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
| [How It Works](../how-it-works.md) | Plain-language conceptual guide — pathway scoring, validation gates, 5-layer architecture |
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
