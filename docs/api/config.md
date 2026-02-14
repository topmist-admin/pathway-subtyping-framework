# Configuration API Reference

The `config` module provides utilities for loading and validating pipeline configuration.

## Classes

### `ConfigValidationError`

Custom exception for configuration validation errors with actionable suggestions.

```python
from pathway_subtyping.config import ConfigValidationError

# Raised automatically by validation functions
try:
    validate_config(invalid_config)
except ConfigValidationError as e:
    print(e.field)        # "data.vcf_path"
    print(e.suggestions)  # ["Add 'vcf_path: /path/to/file' under the 'data:' section"]
```

**Attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `field` | `Optional[str]` | The config field that caused the error |
| `suggestions` | `List[str]` | Actionable suggestions to fix the error |

**Features:**
- Inherits from `ValueError` for backward compatibility
- Formats error messages with field location and numbered suggestions
- Provides context-aware fix recommendations

---

## Functions

### `load_config(config_path: str) -> Dict[str, Any]`

Load configuration from a YAML file.

```python
from pathway_subtyping.config import load_config

config = load_config("configs/my_config.yaml")
print(config["pipeline"]["name"])
```

**Parameters:**
- `config_path`: Path to YAML configuration file

**Returns:** Configuration dictionary

**Raises:**
- `FileNotFoundError`: If config file doesn't exist
- `ConfigValidationError`: If file is empty or invalid YAML
- `yaml.YAMLError`: If YAML syntax is invalid

---

### `validate_config(config: Dict[str, Any], check_files: bool = True) -> bool`

Validate configuration structure and required fields.

```python
from pathway_subtyping.config import load_config, validate_config

config = load_config("configs/my_config.yaml")

# Full validation including file existence checks
validate_config(config)

# Skip file checks (useful for testing)
validate_config(config, check_files=False)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config` | `Dict[str, Any]` | required | Configuration dictionary |
| `check_files` | `bool` | `True` | Whether to verify input files exist |

**Returns:** `True` if valid

**Raises:** `ConfigValidationError` with field and suggestions if invalid

**Validation checks:**
- Required sections: `pipeline`, `data`
- Required data fields: `vcf_path`, `phenotype_path`, `pathway_db`
- File existence (if `check_files=True`)
- Seed must be integer (if provided)
- Cluster range validation (min ≥ 2, max ≥ min)

---

### `validate_gmt_file(gmt_path: str) -> Dict[str, List[str]]`

Validate and load a GMT (Gene Matrix Transposed) pathway file.

```python
from pathway_subtyping.config import validate_gmt_file

# Returns validated pathways or raises ConfigValidationError
pathways = validate_gmt_file("data/pathways/autism.gmt")

print(f"Loaded {len(pathways)} pathways")
for name, genes in pathways.items():
    print(f"  {name}: {len(genes)} genes")
```

**Parameters:**
- `gmt_path`: Path to GMT file

**Returns:** Dictionary mapping pathway names to gene lists

**Raises:**
- `FileNotFoundError`: If GMT file doesn't exist
- `ConfigValidationError`: If GMT format is invalid

**Validation checks:**
- Each line must have at least 3 tab-separated fields
- Each pathway must have at least 2 genes
- No duplicate pathway names allowed
- Empty pathway names are rejected
- Lines starting with `#` are treated as comments

**Example error output:**
```
GMT file has 2 error(s):
Line 5: Expected at least 3 tab-separated fields, got 2
Line 12: Pathway 'SYNAPTIC' has fewer than 2 genes

Suggested fixes:
  1. GMT format: PATHWAY_NAME<TAB>DESCRIPTION<TAB>GENE1<TAB>GENE2<TAB>...
  2. Each line must have at least 3 tab-separated fields
  3. Each pathway must have at least 2 genes
```

---

## Configuration File Format

### Complete Example

```yaml
# Pipeline metadata
pipeline:
  name: autism_cohort_analysis
  output_dir: outputs/autism_2026
  seed: 42
  verbose: true

# Input data paths
data:
  vcf_path: data/cohorts/autism_exomes.vcf
  phenotype_path: data/cohorts/autism_phenotypes.csv
  pathway_db: data/pathways/autism_pathways.gmt

# Clustering parameters
clustering:
  n_clusters: null           # null = auto-select via BIC
  n_clusters_range: [2, 8]   # Range to test for BIC selection

# Variant quality control (optional)
variant_qc:
  enabled: true
  min_qual: 30
  min_call_rate: 0.95
  hwe_p_threshold: 1e-6
  max_maf: 0.01

# Output settings
output:
  disclaimer: "Research use only. Not for clinical diagnosis."
```

### Minimal Example

```yaml
pipeline:
  name: quick_test

data:
  vcf_path: data/test.vcf
  phenotype_path: data/test.csv
  pathway_db: data/pathways/test.gmt
```

---

## Configuration Sections

### `pipeline`

General pipeline settings.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | string | `"demo_run"` | Run identifier |
| `output_dir` | string | `"outputs/demo_run"` | Output directory |
| `seed` | integer/null | `42` | Random seed (null for random) |
| `verbose` | boolean | `true` | Enable verbose logging |

### `data`

Input file paths.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `vcf_path` | string | Yes | Path to VCF file |
| `phenotype_path` | string | Yes | Path to phenotype CSV |
| `pathway_db` | string | Yes | Path to pathway GMT file |

### `clustering`

Clustering algorithm parameters.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `n_clusters` | integer/null | `null` | Fixed cluster count (null = auto) |
| `n_clusters_range` | list[int] | `[2, 8]` | Range for BIC selection |

### `validation`

Validation gate configuration. When thresholds are `null`, the pipeline auto-calibrates them based on sample size and cluster count.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `run_gates` | boolean | `true` | Whether to run validation gates |
| `calibrate` | boolean | `true` | Auto-calibrate thresholds when null |
| `stability_threshold` | float/null | `null` | Stability ARI threshold (null = auto) |
| `null_ari_max` | float/null | `null` | Null ARI max threshold (null = auto) |
| `alpha` | float | `0.05` | Significance level for calibration |
| `n_permutations` | integer | `100` | Permutations for null tests |
| `n_bootstrap` | integer | `50` | Bootstrap iterations for stability |

**Validation checks (performed by `_validate_validation_section()`):**
- `stability_threshold` must be in [0, 1] (if provided)
- `null_ari_max` must be in [0, 1] (if provided)
- `alpha` must be in (0, 1)
- `n_permutations` must be a positive integer
- `n_bootstrap` must be a positive integer

### `variant_qc`

Variant quality control filters. Applied between data loading and burden computation.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | boolean | `false` | Whether to apply variant QC filters |
| `min_qual` | float | `30.0` | Minimum QUAL score (≥0) |
| `min_call_rate` | float | `0.9` | Minimum fraction of non-missing genotypes (0.0–1.0) |
| `hwe_p_threshold` | float | `1e-6` | HWE p-value threshold; variants below are removed (0 to disable) |
| `max_maf` | float | `0.01` | Maximum minor allele frequency (0.0–1.0; 1.0 to disable) |
| `min_gq` | integer/null | `null` | Minimum genotype quality; genotypes below are set to missing |
| `min_dp` | integer/null | `null` | Minimum read depth; genotypes below are set to missing |

**Validation checks (performed by `_validate_variant_qc_section()`):**
- `min_qual` must be ≥ 0
- `min_call_rate` must be in [0, 1]
- `hwe_p_threshold` must be in [0, 1]
- `max_maf` must be in [0, 1]
- `min_gq` must be a non-negative integer (if provided)
- `min_dp` must be a non-negative integer (if provided)

### `output`

Output formatting options.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `disclaimer` | string | `"Research use only..."` | Report disclaimer |

---

## Example: Programmatic Configuration

```python
from pathway_subtyping.pipeline import PipelineConfig

# Create config without YAML file
config = PipelineConfig(
    name="programmatic_run",
    output_dir="outputs/prog_run",
    seed=123,
    verbose=True,
    vcf_path="data/my_data.vcf",
    phenotype_path="data/my_pheno.csv",
    pathway_db="data/pathways/my_pathways.gmt",
    n_clusters=None,  # Auto-select
    n_clusters_range=[3, 6],
)
```

---

## Example: Config Validation

```python
from pathway_subtyping.config import (
    load_config,
    validate_config,
    validate_gmt_file,
    ConfigValidationError,
)

# Load and validate config
config = load_config("my_config.yaml")

try:
    validate_config(config)
    print("Configuration is valid")
except ConfigValidationError as e:
    print(f"Invalid configuration: {e}")
    if e.field:
        print(f"Problem field: {e.field}")
    if e.suggestions:
        print("Suggestions:")
        for suggestion in e.suggestions:
            print(f"  - {suggestion}")

# Validate GMT file separately
try:
    pathways = validate_gmt_file(config["data"]["pathway_db"])
    print(f"Loaded {len(pathways)} valid pathways")
except ConfigValidationError as e:
    print(f"Invalid GMT file: {e}")
```
