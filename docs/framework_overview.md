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
│  VCF File          Pathways GMT       Phenotypes CSV            │
│  (variants)        (gene sets)        (metadata)                │
└────────┬────────────────┬────────────────────┬──────────────────┘
         │                │                    │
         ▼                ▼                    │
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
│                   GMM CLUSTERING                                 │
├─────────────────────────────────────────────────────────────────┤
│  • Fit Gaussian Mixture Models for K = 2..max_k                 │
│  • Select optimal K by BIC                                       │
│  • Assign samples to clusters                                    │
│  • Calculate confidence (posterior probability)                  │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                   VALIDATION GATES                               │
├─────────────────────────────────────────────────────────────────┤
│  Gate 1: Label Shuffle     → Verify no spurious patterns        │
│  Gate 2: Random Gene Sets  → Verify pathways drive clustering   │
│  Gate 3: Bootstrap         → Verify cluster stability           │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                     OUTPUTS                                      │
├─────────────────────────────────────────────────────────────────┤
│  pathway_scores.csv       subtype_assignments.csv               │
│  report.json              report.md                              │
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

### 2. Pipeline (`pipeline.py`)

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

### 3. Clustering (`clustering.py`)

**Purpose**: Discover molecular subtypes via GMM.

**Algorithm**:
1. Fit GMM for each K in range [min_k, max_k]
2. Calculate BIC for each model
3. Select K with lowest BIC
4. Return hard assignments and soft probabilities

```python
from sklearn.mixture import GaussianMixture

gmm = GaussianMixture(
    n_components=k,
    covariance_type='full',
    random_state=seed
)
gmm.fit(pathway_scores)
labels = gmm.predict(pathway_scores)
probabilities = gmm.predict_proba(pathway_scores)
```

### 4. Validation (`validation.py`)

**Purpose**: Ensure clustering quality via statistical tests.

**Validation Gates**:

| Gate | Test | Pass Condition |
|------|------|----------------|
| Label Shuffle | Shuffle labels, re-cluster | ARI < 0.1 |
| Random Genes | Random gene sets, re-cluster | ARI < 0.1 |
| Bootstrap | Resample, re-cluster | ARI >= 0.7 |

```python
from pathway_subtyping import ValidationGates

gates = ValidationGates(seed=42)
results = gates.run_all_gates(
    pathway_scores=scores,
    cluster_labels=labels,
    pathways=pathway_dict
)
```

### 5. Cross-Cohort (`cross_cohort.py`)

**Purpose**: Compare subtypes across independent cohorts.

**Methods**:
- Direct comparison: Cluster both cohorts, compare with ARI
- Transfer: Project cohort B onto cohort A's model

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
| Variants | Yes | VCF | Raw genetic data |
| Pathways | Yes | GMT | Biological groupings |
| Phenotypes | No | CSV | Sample metadata |

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

## See Also

- [Quickstart Guide](quickstart.md) - Getting started
- [API Reference](api/index.md) - Detailed API docs
- [Data Formats](data_formats.md) - Input/output specs
- [Validation Gates](guides/validation-gates.md) - Understanding validation
