# Framework Overview

> **RESEARCH USE ONLY** — This framework is for research purposes only. Not for clinical decision-making. See [DISCLAIMER.md](../DISCLAIMER.md).

This document provides a high-level overview of the Pathway Subtyping Framework architecture, design principles, and workflow.

---

## Purpose

The Pathway Subtyping Framework enables **pathway-based molecular subtype discovery** in genetically heterogeneous diseases. It addresses a fundamental challenge: how to find meaningful subgroups when hundreds of genes contribute to disease risk.

### The Problem

Genetically complex diseases (autism, schizophrenia, epilepsy, etc.) involve:
- Thousands of associated genetic variants
- Hundreds of implicated genes
- No single "disease gene"

Traditional approaches struggle because:
- Individual variants are too rare for statistical power
- Gene-level analysis loses biological context
- Unsupervised clustering on raw variants is noisy

### The Solution

Aggregate rare variants into **pathway-level scores**, then cluster to discover reproducible molecular subtypes.

```
Variants → Genes → Pathways → Subtypes
```

---

## Architecture

### Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                     INPUT DATA                                   │
├─────────────────────────────────────────────────────────────────┤
│  VCF File     Expression CSV    scRNA-seq h5ad    Pathways GMT  │
│  (variants)   (bulk RNA-seq)    (single-cell)     (gene sets)   │
└────────┬────────────┬────────────────┬──────────────┬───────────┘
         │                │                    │                │
         ▼                ▼                    ▼                │
┌─────────────────────────────────────────────────────────────────┐
│              PATHWAY SCORING (modality-specific)                 │
├─────────────────────────────────────────────────────────────────┤
│  VCF path:     Variant QC → Gene Burden → Pathway Aggregation   │
│  Expression:   Load matrix → ssGSEA / GSVA / mean-Z scoring    │
│  Single-cell:  Load h5ad → Pseudobulk or per-cell scoring       │
│  Output:       samples/cell_types × pathways (Z-normalized)     │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              VARIANT QUALITY CONTROL (VCF path only)             │
├─────────────────────────────────────────────────────────────────┤
│  • Filter by QUAL score (default ≥ 30)                          │
│  • Remove low call rate variants (default ≥ 90%)                │
│  • Test Hardy-Weinberg equilibrium (chi-squared)                │
│  • Filter by minor allele frequency (default ≤ 1%)              │
│  • Apply per-genotype GQ/DP filters (optional)                  │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                   GENE BURDEN SCORING                            │
├─────────────────────────────────────────────────────────────────┤
│  • Parse VCF variants                                            │
│  • Weight by consequence (LoF > missense > other)               │
│  • Weight by pathogenicity (CADD score)                         │
│  • Aggregate to gene-level burden scores                         │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                   PATHWAY SCORING                                │
├─────────────────────────────────────────────────────────────────┤
│  • Map genes to pathways                                         │
│  • Aggregate gene burdens to pathway scores                     │
│  • Normalize by pathway size (√n)                               │
│  • Z-score normalize across samples                              │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              ANCESTRY CORRECTION (optional)                      │
├─────────────────────────────────────────────────────────────────┤
│  • Load or compute ancestry PCs from genotype data              │
│  • Regress out ancestry-correlated variance per pathway         │
│  • Flag confounded pathways (R² > 0.1)                          │
│  • Re-normalize adjusted scores                                  │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              BATCH CORRECTION (optional)                          │
├─────────────────────────────────────────────────────────────────┤
│  • Detect batch effects via ANOVA (eta-squared)                  │
│  • ComBat empirical Bayes correction (default)                   │
│  • Validate variance reduction & signal preservation             │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              MULTI-OMIC FUSION (optional)                        │
├─────────────────────────────────────────────────────────────────┤
│  • Combine VCF + expression + single-cell pathway scores        │
│  • Strategies: concatenate, weighted average, intersection       │
│  • Handle partial sample overlap (impute_zero/mean or drop)     │
│  Output: unified samples × pathways matrix                       │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              BULK DECONVOLUTION (optional)                        │
├─────────────────────────────────────────────────────────────────┤
│  • Build reference profile from single-cell data                │
│  • NNLS deconvolution: estimate cell-type proportions           │
│  • Combine proportions + pathway scores for cell-type-aware     │
│    subtype discovery                                             │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                   CLUSTERING                                     │
├─────────────────────────────────────────────────────────────────┤
│  • GMM (default), K-means, Hierarchical, Spectral               │
│  • Select optimal K by BIC (GMM) or cross-validation            │
│  • Assign samples to clusters with confidence scores            │
│  • Algorithm comparison with pairwise ARI                        │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                   VALIDATION GATES                               │
├─────────────────────────────────────────────────────────────────┤
│  Gate 1: Label Shuffle     → Verify no spurious patterns        │
│  Gate 2: Random Gene Sets  → Verify pathways drive clustering   │
│  Gate 3: Bootstrap         → Verify cluster stability           │
│  Gate 4: Ancestry Indep.   → Verify no ancestry confounding     │
│  Gate 5: Cross-Modal       → Verify subtypes replicate across   │
│           Concordance        data modalities (multi-omic only)  │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              SENSITIVITY ANALYSIS (optional)                      │
├─────────────────────────────────────────────────────────────────┤
│  • Vary clustering algorithm (GMM, K-means, Hierarchical)       │
│  • Vary number of clusters (k sweep)                             │
│  • Leave-one-out feature sensitivity                             │
│  • Vary normalization method                                     │
│  • Report overall robustness score                               │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              VISUALIZATION (optional)                             │
├─────────────────────────────────────────────────────────────────┤
│  • Interactive HTML report (Plotly) — scatter, heatmap, radar   │
│  • UMAP / t-SNE / PCA dimensionality reduction                  │
│  • Multi-format export (PNG, SVG, PDF, HTML)                    │
│  • Static matplotlib fallback when Plotly unavailable            │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                     OUTPUTS                                      │
├─────────────────────────────────────────────────────────────────┤
│  pathway_scores.csv       subtype_assignments.csv               │
│  report.json              report.md                              │
│  interactive_report.html  figures/ (PNG, SVG, PDF)              │
└─────────────────────────────────────────────────────────────────┘
```

---

## Core Components

### 1. Data Loaders (`io/`)

**Purpose**: Parse input files into internal data structures.

| Module | Input | Output |
|--------|-------|--------|
| `load_vcf()` | VCF file | Variant DataFrame |
| `load_pathways()` | GMT file | Dict[pathway → genes] |
| `load_phenotypes()` | CSV file | Phenotype DataFrame |

### 2. Variant QC (`variant_qc.py`)

**Purpose**: Remove technical artifacts and common variants before burden computation.

**Filters** (applied sequentially):
- QUAL score filtering
- Per-variant call rate
- Hardy-Weinberg equilibrium test (chi-squared, 1 df)
- Minor allele frequency (MAF) filtering
- Per-genotype GQ/DP masking (optional)

```python
from pathway_subtyping import VariantQCConfig, filter_variants

config = VariantQCConfig(min_qual=30, min_call_rate=0.95, max_maf=0.01)
filtered_variants, filtered_genotypes, result = filter_variants(
    variants_df, genotypes_df, config
)
print(result.format_report())
```

### 3. Pipeline (`pipeline.py`)

**Purpose**: Orchestrate the analysis workflow.

```python
from pathway_subtyping import DemoPipeline, PipelineConfig

config = PipelineConfig(
    run_name="my_analysis",
    vcf_path="variants.vcf",
    pathways_path="pathways.gmt",
    seed=42
)

pipeline = DemoPipeline(config)
results = pipeline.run()
```

**Key Classes**:
- `PipelineConfig`: Configuration dataclass
- `DemoPipeline`: Main pipeline orchestrator

### 4. Clustering (`clustering.py`)

**Purpose**: Discover molecular subtypes via unsupervised clustering.

**Algorithms**:
- **GMM** (default): Soft assignments, automatic K selection via BIC
- **K-means**: Fast, spherical clusters
- **Hierarchical**: Dendrogram-based, no K required
- **Spectral**: Nonlinear cluster boundaries

**Workflow**:
1. Fit models for each K in range [min_k, max_k]
2. Select optimal K (BIC for GMM, silhouette or cross-validation for others)
3. Return hard assignments, soft probabilities (GMM), and confidence scores
4. Compare algorithms with pairwise ARI for robustness

```python
from pathway_subtyping import run_clustering, ClusteringMethod

result = run_clustering(
    pathway_scores.values,
    n_clusters=3,
    method=ClusteringMethod.GMM,
    seed=42,
)
print(f"Labels: {result.labels}")
print(f"BIC: {result.bic}")
```

### 5. Validation (`validation.py`)

**Purpose**: Ensure clustering quality via statistical tests.

**Validation Gates**:

| Gate | Test | Pass Condition |
|------|------|----------------|
| Label Shuffle | Shuffle labels, re-cluster | ARI < 0.1 |
| Random Genes | Random gene sets, re-cluster | ARI < 0.1 |
| Bootstrap | Resample, re-cluster | ARI >= 0.7 |
| Ancestry Independence | Kruskal-Wallis per PC | No significant association (Bonferroni) |
| Cross-Modal Concordance | Cluster each modality independently | ARI > null 95th percentile |

Gate 5 (Cross-Modal Concordance) is automatically included when `per_modality_scores` is provided to `run_all()`.

```python
from pathway_subtyping import ValidationGates

gates = ValidationGates(seed=42)
results = gates.run_all(
    pathway_scores=scores,
    cluster_labels=labels,
    pathways=pathway_dict,
    gene_burdens=gene_data,
    n_clusters=3,
    per_modality_scores={"WES": vcf_scores, "RNA-seq": expr_scores},  # enables Gate 5
)
```

### 6. Cross-Cohort (`cross_cohort.py`)

**Purpose**: Compare subtypes across independent cohorts.

**Methods**:
- Direct comparison: Cluster both cohorts, compare with ARI
- Transfer: Project cohort B onto cohort A's model

### 7. Ancestry Correction (`ancestry.py`)

**Purpose**: Detect and correct for population stratification confounding.

**Methods**:
- `compute_ancestry_pcs()`: PCA on genotype matrix for ancestry inference
- `adjust_pathway_scores()`: Regress ancestry PCs out of pathway scores
- `check_ancestry_independence()`: Kruskal-Wallis test for cluster-ancestry association
- `stratified_analysis()`: Per-ancestry-group clustering with concordance

```python
from pathway_subtyping import (
    compute_ancestry_pcs,
    adjust_pathway_scores,
    check_ancestry_independence,
)

# Compute ancestry PCs
pcs = compute_ancestry_pcs(genotype_matrix, n_components=10, seed=42)

# Adjust pathway scores
result = adjust_pathway_scores(pathway_scores, pcs)
adjusted_scores = result.adjusted_scores

# Verify independence after clustering
report = check_ancestry_independence(cluster_labels, pcs)
print(f"Independent of ancestry: {report.passed}")
```

### 8. Batch Correction (`batch_correction.py`)

**Purpose**: Detect and correct batch effects from multi-site or multi-run data.

**Methods**:
- `detect_batch_effects()`: ANOVA-based detection with eta-squared variance explained
- `correct_batch_effects()`: ComBat empirical Bayes, mean centering, or standardization
- `validate_batch_correction()`: Post-correction validation of variance reduction

```python
from pathway_subtyping import detect_batch_effects, correct_batch_effects

# Detect batch effects
report = detect_batch_effects(pathway_scores, batch_labels)
print(f"Batch effect detected: {report.overall_batch_effect}")

# Correct batch effects
result = correct_batch_effects(pathway_scores, batch_labels)
corrected_scores = result.corrected_scores
```

### 9. Sensitivity Analysis (`sensitivity.py`)

**Purpose**: Evaluate robustness of clustering results to parameter choices.

**Methods**:
- `vary_clustering_algorithm()`: Compare results across GMM, K-means, Hierarchical
- `vary_n_clusters()`: Sweep k range with pairwise ARI concordance
- `vary_feature_subset()`: Leave-one-out pathway sensitivity
- `vary_normalization()`: Compare z-score, min-max, robust, rank normalization
- `run_sensitivity_analysis()`: Full sensitivity analysis across all axes

```python
from pathway_subtyping import run_sensitivity_analysis

result = run_sensitivity_analysis(pathway_scores, n_clusters=3, seed=42)
print(f"Overall stability: {result.overall_stability:.3f}")
print(f"Robust: {result.is_robust}")
print(f"Most sensitive: {result.most_sensitive_parameter}")
```

### 10. Single-Cell Scoring (`single_cell.py`)

**Purpose**: Score pathways from scRNA-seq data at per-cell or cell-type level.

**Requires**: `pip install pathway-subtyping[sc]` (anndata).

**Methods**:
- `load_single_cell_data()`: Load h5ad or CSV, auto-normalize raw counts, QC
- `score_single_cell_pathways()`: Pseudobulk (ssGSEA, GSVA, mean-Z) or per-cell (mean-Z) scoring

```python
from pathway_subtyping import (
    load_single_cell_data,
    score_single_cell_pathways,
    SingleCellScoringMethod,
)

adata, qc = load_single_cell_data("pbmc.h5ad", cell_type_column="cell_type")
result = score_single_cell_pathways(
    adata, pathways, cell_type_column="cell_type",
    method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA, seed=42,
)
# result.pathway_scores: cell_types × pathways (Z-normalized)
```

### 11. Visualization (`visualization.py`)

**Purpose**: Generate interactive and publication-quality visualizations.

**Requires**: `pip install pathway-subtyping[viz]` for interactive features. Static matplotlib fallbacks work with base install.

**Functions**:
- `compute_dim_reduction()`: PCA, t-SNE, or UMAP dimensionality reduction
- `plot_interactive_scatter()`: Plotly scatter plot of samples in reduced space
- `plot_interactive_heatmap()`: Plotly heatmap of mean pathway Z-scores per subtype
- `plot_cluster_distribution()`: Plotly bar chart of cluster sizes
- `plot_subtype_trajectories()`: Plotly radar chart of subtype pathway profiles
- `plot_static_scatter()`: Static matplotlib fallback
- `export_figure()`: Multi-format export (PNG, SVG, PDF, HTML)
- `create_interactive_report()`: Self-contained interactive HTML report
- `generate_all_figures()`: Generate all visualizations at once

```python
from pathway_subtyping import create_interactive_report, ReportConfig, DimReductionMethod

result = create_interactive_report(
    pathway_scores=scores_df,
    labels=labels,
    output_path="report.html",
    config=ReportConfig(dim_reduction=DimReductionMethod.UMAP),
    seed=42,
)
# Open report.html in any browser — no server needed
```

### 12. Expression Scoring (`expression.py`)

**Purpose**: Score pathways from bulk RNA-seq gene expression data.

**Methods**:
- **ssGSEA** (recommended): Rank-based single-sample gene set enrichment
- **GSVA**: Gene Set Variation Analysis via empirical CDF + KS statistic
- **Mean-Z**: Fast Z-score averaging per pathway (best for large datasets)

**Key functions**:
- `load_expression_matrix()`: Load CSV/TSV with auto-orientation detection and log transformation
- `score_pathways_from_expression()`: Main scoring entry point — produces same Z-normalized pathway scores as VCF-based scoring

```python
from pathway_subtyping import (
    load_expression_matrix, score_pathways_from_expression,
    ExpressionScoringMethod, ExpressionInputType,
)

expr, qc = load_expression_matrix("expression.csv", input_type=ExpressionInputType.TPM)
result = score_pathways_from_expression(expr, pathways, method=ExpressionScoringMethod.SSGSEA, seed=42)
# result.pathway_scores: samples × pathways (Z-normalized)
```

### 13. Multi-Omic Fusion (`multi_omic.py`)

**Purpose**: Fuse pathway scores from multiple data modalities (VCF, expression, single-cell, deconvolution) into a unified feature matrix.

**Fusion strategies**:
- **Concatenate** (default): Column-bind with modality prefixes (e.g., `WES:PATHWAY_1`)
- **Weighted Average**: Weighted mean of shared pathway scores
- **Intersection Only**: Restrict to shared samples and pathways

**Key functions**:
- `prepare_modality()`: Wrap a modality's scores into `ModalityInput`
- `fuse_modalities()`: Fuse 2+ modalities with configurable strategy and missing data handling
- `correlation_analysis()`: Pairwise pathway-level correlations between modalities

```python
from pathway_subtyping import (
    ModalityType, FusionStrategy, prepare_modality, fuse_modalities,
)

vcf_mod = prepare_modality(ModalityType.VCF, vcf_scores, label="WES")
expr_mod = prepare_modality(ModalityType.EXPRESSION, expr_scores, label="RNA-seq")
fusion = fuse_modalities([vcf_mod, expr_mod], strategy=FusionStrategy.CONCATENATE, seed=42)
# fusion.fused_pathway_scores: unified samples × pathways
# fusion.per_modality_scores: dict for cross-modal validation (Gate 5)
```

### 14. Bulk Deconvolution (`deconvolution.py`)

**Purpose**: Estimate cell-type proportions from bulk RNA-seq using a single-cell reference profile, then combine with pathway scores for cell-type-aware subtype discovery.

**Algorithm**: Non-negative least squares (NNLS) — `scipy.optimize.nnls`. Biologically appropriate because cell-type proportions cannot be negative (used in CIBERSORT, MuSiC).

**Key functions**:
- `build_reference_profile()`: Aggregate single-cell expression to cell-type mean profiles
- `deconvolve_bulk()`: NNLS deconvolution with quality checks
- `combine_features()`: Merge pathway scores + cell-type proportions into unified feature matrix
- `generate_synthetic_bulk()`: Create synthetic bulk expression from known proportions for testing

```python
from pathway_subtyping import build_reference_profile, deconvolve_bulk, combine_features

reference = build_reference_profile(sc_expression, cell_type_labels)
result = deconvolve_bulk(bulk_expression, reference, seed=42)
combined = combine_features(pathway_scores, result.cell_type_proportions, proportion_weight=0.3)
# combined: samples × (pathways + CELLTYPE:* columns), ready for clustering
```

### 15. Cross-Modal Validation (`cross_modal_validation.py`)

**Purpose**: Test whether molecular subtypes are consistent across data modalities (Validation Gate 5).

**How it works**:
1. For each pair of modalities, cluster each independently using GMM
2. Measure agreement via ARI and NMI (concordance)
3. Test transfer: train on modality A, predict on modality B (and vice versa)
4. Build null distribution by permuting labels across `n_permutations` runs
5. Gate passes if observed ARI > 95th percentile of null distribution

**Key functions**:
- `cross_modal_concordance()`: Main entry point for cross-modal validation
- `single_cell_composition_test()`: Test whether subtypes have distinct cell-type compositions
- `generate_synthetic_multimodal_data()`: Generate synthetic multi-modal data with planted subtypes

```python
from pathway_subtyping import cross_modal_concordance

result = cross_modal_concordance(
    per_modality_scores={"WES": vcf_scores, "RNA-seq": expr_scores},
    cluster_labels=fused_labels,
    fused_sample_ids=list(fused_scores.index),
    n_clusters=3,
    seed=42,
)
print(f"Gate passed: {result.gate_passed}")
print(result.format_report())
```

---

## Design Principles

### 1. Disease-Agnostic

The framework makes no assumptions about specific diseases:
- User provides pathways relevant to their condition
- Works with any VCF + pathway combination
- No hardcoded gene lists or disease logic

### 2. Validated by Default

Every run includes mandatory validation gates:
- Prevents reporting spurious clusters
- Ensures reproducibility
- Builds user confidence

### 3. Reproducible

All randomness is controlled by a single seed:
- Same inputs + same seed = same outputs
- Facilitates debugging and verification
- Enables exact replication

### 4. Interpretable

Outputs are designed for understanding:
- Pathway-level scores (not opaque embeddings)
- Clear cluster labels
- Confidence scores for each assignment

### 5. Extensible

Modular design allows customization:
- Swap clustering algorithms
- Add validation gates
- Customize variant weighting

---

## Data Flow

### Input Requirements

| Data | Required | Format | Purpose |
|------|----------|--------|---------|
| Variants | Yes* | VCF | Raw genetic data |
| Expression | Yes* | CSV/TSV | Bulk RNA-seq gene expression |
| Single-Cell | Yes* | h5ad/CSV | scRNA-seq data with cell types |
| Pathways | Yes | GMT | Biological groupings |
| Phenotypes | No | CSV | Sample metadata |

*One input modality required (VCF, expression, or single-cell).

### Internal Representations

| Stage | Shape | Description |
|-------|-------|-------------|
| Variant matrix | (n_variants, n_samples) | Genotype calls |
| Gene burden | (n_samples, n_genes) | Per-gene scores |
| Pathway scores | (n_samples, n_pathways) | Per-pathway scores |
| Cluster labels | (n_samples,) | Subtype assignments |

### Output Files

| File | Format | Content |
|------|--------|---------|
| `pathway_scores.csv` | CSV | Sample × Pathway matrix |
| `subtype_assignments.csv` | CSV | Cluster assignments |
| `report.json` | JSON | Machine-readable summary |
| `report.md` | Markdown | Human-readable report |
| `figures/summary.png` | PNG | Static clustering visualization |
| `figures/interactive_report.html` | HTML | Interactive Plotly report (optional) |

---

## Configuration

### Minimal Config

```yaml
run_name: my_analysis
input:
  vcf_path: data/variants.vcf
  pathways_path: data/pathways.gmt
pipeline:
  seed: 42
```

### Full Config

```yaml
run_name: my_analysis

input:
  vcf_path: data/variants.vcf
  phenotypes_path: data/phenotypes.csv
  pathways_path: data/pathways.gmt

output:
  output_dir: outputs/my_analysis

pipeline:
  seed: 42
  min_samples_per_cluster: 10

clustering:
  n_clusters_range: [2, 8]
  covariance_type: full

variant_qc:
  enabled: true
  min_qual: 30
  min_call_rate: 0.95
  hwe_p_threshold: 1e-6
  max_maf: 0.01

validation:
  run_validation: true
  n_permutations: 100
  n_bootstrap: 50
  ari_threshold: 0.7
```

---

## Extension Points

### Custom Variant Weighting

```python
def custom_weight(consequence, cadd_score):
    if consequence in LOF_CONSEQUENCES:
        return 3.0
    elif consequence == 'missense_variant':
        return min(cadd_score / 20.0, 2.0)
    else:
        return 0.5
```

### Custom Clustering

```python
from sklearn.cluster import SpectralClustering

def custom_cluster(pathway_scores, n_clusters):
    model = SpectralClustering(n_clusters=n_clusters)
    return model.fit_predict(pathway_scores)
```

### Additional Validation Gates

```python
def custom_validation_gate(pathway_scores, labels):
    # Your custom validation logic
    metric = calculate_custom_metric(pathway_scores, labels)
    passed = metric > threshold
    return ValidationResult(
        name="Custom Gate",
        status="PASS" if passed else "FAIL",
        metric_value=metric
    )
```

---

## Performance Considerations

### Memory Usage

| Cohort Size | Pathways | Approx. Memory |
|-------------|----------|----------------|
| 100 samples | 50 | ~100 MB |
| 1,000 samples | 50 | ~500 MB |
| 10,000 samples | 50 | ~4 GB |

### Optimization Utilities

For large cohorts, use chunked processing:

```python
from pathway_subtyping.utils import (
    chunked_vcf_reader,
    compute_gene_burdens_chunked,
    estimate_memory_usage
)

# Check memory requirements first
estimate_memory_usage(n_samples=10000, n_genes=20000)

# Process in chunks
for chunk in chunked_vcf_reader(vcf_path, chunk_size=1000):
    process(chunk)
```

---

## Data Provenance

**This framework contains zero proprietary, commercial, or third-party patient data.** All data shipped with the repository is either:

1. **Computationally generated** — Synthetic VCF and phenotype files created by `SyntheticDataGenerator` using random number generators with fixed seeds. They contain no real patient or clinical data.
2. **Curated from public scientific literature** — Pathway GMT files contain standard HGNC gene symbols from publicly available databases (SFARI Gene, KEGG, Reactome, MSigDB, Gene Ontology).
3. **Open-source code only** — All algorithms are original implementations or use standard open-source libraries (scikit-learn, scipy, numpy, pandas).

**No data from any employer, client, institution, or commercial entity was used at any stage of this project.** The framework is designed so that users supply their own data; it does not depend on any private or restricted datasets.

For full details, see [DISCLAIMER.md](../DISCLAIMER.md).

---

## See Also

- [Quickstart Guide](quickstart.md) - Getting started
- [API Reference](api/index.md) - Detailed API docs
- [Data Formats](data_formats.md) - Input/output specs
- [Validation Gates](guides/validation-gates.md) - Understanding validation
