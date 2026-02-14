# API Reference

This document provides comprehensive API documentation for the Pathway Subtyping Framework.

## Modules Overview

| Module | Description |
|--------|-------------|
| [`pipeline`](pipeline.md) | Main pipeline orchestrator and configuration |
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
| `figures/summary.png` | PNG | Visualization of clustering results |
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
