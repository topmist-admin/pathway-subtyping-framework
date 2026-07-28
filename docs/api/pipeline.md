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
| `input_type` | `str` | `"vcf"` | Input modality: `"vcf"`, `"expression"`, or `"multi_omic"` |
| `expression_path` | `str` | `""` | Path to expression matrix (used when `input_type="expression"`) |
| `expression_input_type` | `str` | `"tpm"` | Expression units: `"counts"`, `"tpm"`, `"fpkm"`, `"log2"` |
| `expression_scoring_method` | `str` | `"ssgsea"` | Pathway scoring method: `"mean_z"`, `"ssgsea"`, `"gsva"` |
| `ssgsea_alpha` | `float` | `0.25` | Rank-weighting exponent passed to ssGSEA scoring |
| `ancestry_pcs_path` | `Optional[str]` | `None` | Path to ancestry PC CSV (sample IDs as index, PC columns) |
| `ancestry_correction` | `Optional[str]` | `None` | Ancestry correction method (`AncestryMethod` value: `"regress_out"`, `"covariate_aware"`, `"stratified"`); `None` = disabled |
| `ancestry_n_pcs` | `int` | `10` | Number of ancestry PCs used for correction |
| `disclaimer` | `str` | `"Research use only. Not medical advice."` | Disclaimer text for reports |
| `validation_run_gates` | `bool` | `True` | Whether to run validation gates |
| `validation_n_permutations` | `int` | `100` | Permutations for null tests |
| `validation_n_bootstrap` | `int` | `50` | Bootstrap iterations for stability |
| `validation_stability_threshold` | `Optional[float]` | `None` | Stability threshold (None = auto-calibrate) |
| `validation_null_ari_max` | `Optional[float]` | `None` | Null ARI threshold (None = auto-calibrate) |
| `validation_calibrate` | `bool` | `True` | Enable auto-calibration when thresholds are None |
| `validation_alpha` | `float` | `0.05` | Significance level for threshold calibration |
| `variant_qc_enabled` | `bool` | `False` | Enable variant QC before burden computation |
| `variant_qc_min_qual` | `float` | `30.0` | Minimum QUAL score |
| `variant_qc_min_call_rate` | `float` | `0.9` | Minimum genotype call rate |
| `variant_qc_hwe_p_threshold` | `float` | `1e-6` | HWE p-value threshold |
| `variant_qc_max_maf` | `float` | `0.01` | Maximum minor allele frequency |
| `multi_omic_modalities` | `List[Dict[str, Any]]` | `[]` | Modality entries (`type`, `path`, `label`, …) for `input_type="multi_omic"` |
| `multi_omic_fusion_strategy` | `str` | `"concatenate"` | Fusion strategy: `"concatenate"`, `"weighted_average"`, `"intersection_only"` |
| `multi_omic_missing_strategy` | `str` | `"impute_zero"` | Missing-sample handling: `"impute_zero"`, `"impute_mean"`, `"drop"` |
| `multi_omic_weights` | `Optional[Dict[str, float]]` | `None` | Per-modality weights (label → weight, summing to 1.0) |
| `multi_omic_renormalize` | `bool` | `True` | Z-normalize the fused pathway-score matrix |
| `use_chunked_processing` | `bool` | `False` | Compute gene burdens by streaming the VCF (`compute_gene_burdens_chunked`) instead of loading it in memory |
| `chunk_size` | `int` | `1000` | Chunk size used when `use_chunked_processing` is enabled |
| `generate_interactive_report` | `bool` | `False` | Emit an interactive HTML report alongside the static outputs |
| `interactive_dim_reduction` | `str` | `"pca"` | Embedding for the interactive report: `"pca"`, `"tsne"`, `"umap"` |

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

variant_qc:
  enabled: true
  min_qual: 30
  min_call_rate: 0.95
  hwe_p_threshold: 1e-6
  max_maf: 0.01

validation:
  run_gates: true
  calibrate: true             # Auto-calibrate thresholds
  stability_threshold: null   # null = auto-calibrate
  null_ari_max: null           # null = auto-calibrate
  alpha: 0.05
  n_permutations: 100
  n_bootstrap: 50

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
| `variant_qc_result` | `Optional[VariantQCResult]` | Variant QC results (if enabled) |
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
2. `load_data()` - Load VCF/expression/multi-omic inputs, phenotypes, pathways
3. Scoring, branching on `config.input_type`:
   - `"multi_omic"` - pathway scores were already computed during `load_data()`
   - `"expression"` - `compute_expression_pathway_scores()`
   - otherwise (`"vcf"`) - `run_variant_qc()` (if `variant_qc_enabled`), `compute_gene_burdens()`, `compute_pathway_scores()`
4. Ancestry correction (optional) - load PCs if `ancestry_pcs_path` is set, then adjust scores if `ancestry_correction` is set
5. `cluster_samples()` - GMM clustering with BIC selection
6. `run_validation_gates()` - Execute validation tests
7. `characterize()` - Subtype characterization (best-effort; a failure is logged as a warning and does not stop the run)
8. `generate_outputs()` - Save results and reports

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

##### `compute_expression_pathway_scores() -> None`

Score pathways directly from an expression matrix (used when `config.input_type == "expression"`).

```python
pipeline.compute_expression_pathway_scores()
```

Calls `score_pathways_from_expression()` with the method named by
`config.expression_scoring_method` (`mean_z`, `ssgsea`, or `gsva`),
`alpha=config.ssgsea_alpha`, and `seed=config.seed`.

After calling:
- `pipeline.pathway_scores` - DataFrame (samples × pathways)
- `pipeline.expression_scoring_result` - full `score_pathways_from_expression` result
- `pipeline.gene_burdens` - set to the expression matrix, so validation gates and
  characterization can consume it the same way as burden-based runs

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

Execute validation tests with auto-calibrated or explicit thresholds.

```python
pipeline.run_validation_gates()
```

Threshold determination (in priority order):
1. **Explicit config values**: If `validation_stability_threshold` or `validation_null_ari_max` are set, use those
2. **Auto-calibration**: If thresholds are `None` and `validation_calibrate=True`, call `calibrate_thresholds()` based on n_samples and n_clusters
3. **Defaults**: Fall back to 0.15 (null ARI) and 0.8 (stability)

Always runs three tests:
1. **Label Shuffle**: Ensure clustering doesn't find spurious patterns
2. **Random Gene Sets**: Ensure biological pathways matter
3. **Bootstrap Stability**: Ensure clusters are robust

Two further gates run conditionally, because the pipeline forwards the data they need
to [`ValidationGates.run_all()`](validation.md):
- **Ancestry Independence** — when `pipeline.ancestry_pcs` has been loaded
- **Cross-Modal Concordance** — when a multi-omic run produced per-modality scores

The remaining `run_all()` gates (confound association, genetic/somatic anchoring) take
annotations the pipeline does not collect; call `ValidationGates` directly to use them.

After calling:
- `pipeline.validation_result` - `ValidationGatesResult` object
- `pipeline.calibrated_thresholds` - `CalibratedThresholds` object (if auto-calibrated)
- `pipeline.ancestry_report` - ancestry independence report (if ancestry PCs were loaded)

---

##### `characterize() -> None`

Run subtype characterization on the clustering result.

```python
pipeline.characterize()
```

Calls `characterize_subtypes()` with the pathway scores, cluster labels, gene burdens,
pathways, the `cluster_label` names from `cluster_assignments`, the `confidence` column
(if present), and `config.seed`.

After calling:
- `pipeline.characterization_result` - `CharacterizationResult` object
- `output_dir/figures/subtype_heatmap.png` - characterization heatmap
- `output_dir/characterization/` - exported characterization CSVs

Invoked automatically by `run()` between `run_validation_gates()` and
`generate_outputs()`, wrapped in a `try/except` so a characterization failure is logged
as a warning rather than failing the pipeline.

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

pipeline.run_variant_qc()
if pipeline.variant_qc_result:
    print(f"Variant QC: {pipeline.variant_qc_result.passed_variants}/{pipeline.variant_qc_result.total_variants} retained")

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
