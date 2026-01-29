# Pipeline API Reference

The `pipeline` module provides the main orchestration for running the pathway subtyping analysis.

## Classes

### PipelineConfig

Configuration dataclass for the pipeline.

```python
from pathway_subtyping.pipeline import PipelineConfig
```

#### Attributes

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `name` | `str` | `"demo_run"` | Pipeline run name (used in outputs) |
| `output_dir` | `str` | `"outputs/demo_run"` | Output directory path |
| `seed` | `Optional[int]` | `42` | Random seed for reproducibility |
| `verbose` | `bool` | `True` | Enable verbose logging |
| `vcf_path` | `str` | `""` | Path to annotated VCF file |
| `phenotype_path` | `str` | `""` | Path to phenotype CSV |
| `pathway_db` | `str` | `""` | Path to pathway GMT file |
| `n_clusters` | `Optional[int]` | `None` | Fixed number of clusters (None = auto-select) |
| `n_clusters_range` | `List[int]` | `[2, 8]` | Range for automatic cluster selection |
| `disclaimer` | `str` | `"Research use only..."` | Disclaimer text for reports |

#### Methods

##### `from_yaml(yaml_path: str) -> PipelineConfig`

Load configuration from a YAML file.

```python
config = PipelineConfig.from_yaml("configs/my_config.yaml")
```

**Parameters:**
- `yaml_path`: Path to YAML configuration file

**Returns:** `PipelineConfig` instance

**Raises:** `FileNotFoundError` if file doesn't exist

**YAML Format:**
```yaml
pipeline:
  name: my_analysis
  output_dir: outputs/my_analysis
  seed: 42
  verbose: true

data:
  vcf_path: data/my_cohort.vcf
  phenotype_path: data/my_phenotypes.csv
  pathway_db: data/pathways/autism_pathways.gmt

clustering:
  n_clusters: null  # null = auto-select via BIC
  n_clusters_range: [2, 8]

output:
  disclaimer: "Research use only. Not for clinical decisions."
```

---

### DemoPipeline

Main pipeline orchestrator class.

```python
from pathway_subtyping.pipeline import DemoPipeline
```

#### Constructor

```python
pipeline = DemoPipeline(config: PipelineConfig)
```

**Parameters:**
- `config`: `PipelineConfig` instance

#### Attributes (after `run()`)

| Attribute | Type | Description |
|-----------|------|-------------|
| `variants_df` | `pd.DataFrame` | Parsed variant data |
| `phenotypes_df` | `pd.DataFrame` | Phenotype data |
| `pathways` | `Dict[str, List[str]]` | Pathway gene sets |
| `gene_burdens` | `pd.DataFrame` | Gene-level burden scores |
| `pathway_scores` | `pd.DataFrame` | Pathway-level scores (z-normalized) |
| `cluster_assignments` | `pd.DataFrame` | Final cluster assignments |
| `n_clusters` | `int` | Number of clusters identified |
| `validation_result` | `ValidationGatesResult` | Validation test results |

#### Methods

##### `run() -> None`

Execute the complete pipeline.

```python
pipeline = DemoPipeline(config)
pipeline.run()
```

This runs the following steps in order:
1. `setup()` - Create output directories, configure logging
2. `load_data()` - Load VCF, phenotypes, pathways
3. `compute_gene_burdens()` - Calculate per-gene burden scores
4. `compute_pathway_scores()` - Aggregate to pathway level
5. `cluster_samples()` - GMM clustering with BIC selection
6. `run_validation_gates()` - Execute validation tests
7. `generate_outputs()` - Save results and reports

**Raises:**
- `FileNotFoundError` if input files don't exist
- `ValueError` for invalid data formats

---

##### `setup() -> None`

Initialize output directories and logging.

```python
pipeline.setup()
```

Creates:
- Output directory (`output_dir/`)
- Figures subdirectory (`output_dir/figures/`)
- Log file (`output_dir/pipeline.log`)

---

##### `load_data() -> None`

Load all input data files.

```python
pipeline.load_data()
```

After calling:
- `pipeline.variants_df` - DataFrame of variants
- `pipeline.phenotypes_df` - DataFrame of phenotypes
- `pipeline.pathways` - Dict of pathway gene lists
- `pipeline.samples` - List of sample IDs

---

##### `compute_gene_burdens() -> None`

Calculate gene-level burden scores.

```python
pipeline.compute_gene_burdens()
```

Burden scoring:
- Loss-of-function variants: weight = 1.0
- Missense (CADD > 25): weight = 0.5
- Other missense: weight = 0.1
- Final score: `sum(genotype * weight * CADD/40)`

After calling:
- `pipeline.gene_burdens` - DataFrame (samples × genes)

---

##### `compute_pathway_scores() -> None`

Aggregate gene burdens to pathway level.

```python
pipeline.compute_pathway_scores()
```

For each pathway:
1. Find genes present in both pathway and burden data
2. Compute mean burden across pathway genes
3. Z-score normalize across samples

After calling:
- `pipeline.pathway_scores` - DataFrame (samples × pathways)

---

##### `cluster_samples() -> None`

Cluster samples into molecular subtypes.

```python
pipeline.cluster_samples()
```

Uses Gaussian Mixture Model (GMM) clustering:
- If `config.n_clusters` is set, use that value
- Otherwise, test range and select by lowest BIC

After calling:
- `pipeline.cluster_assignments` - DataFrame with columns:
  - `sample_id`: Sample identifier
  - `cluster_id`: Numeric cluster (0, 1, 2, ...)
  - `cluster_label`: Biological label (synaptic, chromatin, etc.)
  - `confidence`: GMM posterior probability
- `pipeline.n_clusters` - Number of clusters

---

##### `run_validation_gates() -> None`

Execute validation tests.

```python
pipeline.run_validation_gates()
```

Runs three tests:
1. **Label Shuffle**: Ensure clustering doesn't find spurious patterns
2. **Random Gene Sets**: Ensure biological pathways matter
3. **Bootstrap Stability**: Ensure clusters are robust

After calling:
- `pipeline.validation_result` - `ValidationGatesResult` object

---

##### `generate_outputs() -> None`

Save all output artifacts.

```python
pipeline.generate_outputs()
```

Creates:
- `pathway_scores.csv` - Normalized pathway scores
- `subtype_assignments.csv` - Cluster assignments
- `report.json` - Machine-readable report
- `report.md` - Human-readable report
- `figures/summary.png` - Visualization
- `run_metadata.yaml` - Reproducibility info

---

## Example: Full Custom Pipeline

```python
from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

# Create configuration
config = PipelineConfig(
    name="schizophrenia_cohort_2026",
    output_dir="outputs/scz_2026",
    seed=123,
    vcf_path="data/scz_exomes.vcf",
    phenotype_path="data/scz_phenotypes.csv",
    pathway_db="data/pathways/schizophrenia_pathways.gmt",
    n_clusters_range=[3, 7],  # Expect 3-7 subtypes
)

# Run pipeline
pipeline = DemoPipeline(config)
pipeline.run()

# Analyze results
print(f"Identified {pipeline.n_clusters} subtypes")
print(f"Validation: {'PASSED' if pipeline.validation_result.all_passed else 'FAILED'}")

# Get subtype distribution
subtype_counts = pipeline.cluster_assignments['cluster_label'].value_counts()
print(subtype_counts)

# Export for downstream analysis
pipeline.cluster_assignments.to_csv("my_subtypes.csv", index=False)
```

## Example: Step-by-Step Execution

```python
from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

config = PipelineConfig.from_yaml("my_config.yaml")
pipeline = DemoPipeline(config)

# Run steps individually for debugging
pipeline.setup()
print("Setup complete")

pipeline.load_data()
print(f"Loaded {len(pipeline.samples)} samples")

pipeline.compute_gene_burdens()
print(f"Computed burdens for {len(pipeline.gene_burdens.columns)} genes")

pipeline.compute_pathway_scores()
print(f"Scored {len(pipeline.pathway_scores.columns)} pathways")

# Inspect pathway scores before clustering
import matplotlib.pyplot as plt
pipeline.pathway_scores.boxplot(figsize=(12, 6), rot=45)
plt.title("Pathway Score Distributions")
plt.tight_layout()
plt.show()

# Continue with clustering
pipeline.cluster_samples()
pipeline.run_validation_gates()
pipeline.generate_outputs()
```
