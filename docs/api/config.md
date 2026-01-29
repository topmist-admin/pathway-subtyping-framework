# Configuration API Reference

The `config` module provides utilities for loading and validating pipeline configuration.

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

**Raises:** `FileNotFoundError` if file doesn't exist

---

### `validate_config(config: Dict[str, Any]) -> bool`

Validate configuration structure and required fields.

```python
from pathway_subtyping.config import load_config, validate_config

config = load_config("configs/my_config.yaml")
validate_config(config)  # Raises ValueError if invalid
```

**Parameters:**
- `config`: Configuration dictionary

**Returns:** `True` if valid

**Raises:** `ValueError` with details if invalid

**Validation checks:**
- Required sections: `pipeline`, `data`
- Required data fields: `vcf_path`, `phenotype_path`, `pathway_db`

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
from pathway_subtyping.config import load_config, validate_config

# Load and validate
config = load_config("my_config.yaml")

try:
    validate_config(config)
    print("Configuration is valid")
except ValueError as e:
    print(f"Invalid configuration: {e}")
```
