# Quickstart Guide

> **RESEARCH USE ONLY** — This framework is for research purposes only. Not for clinical decision-making. See [DISCLAIMER.md](../DISCLAIMER.md).

Get the Pathway Subtyping Framework running in under 30 minutes.

---

## Prerequisites

- **Python 3.9+** (3.11 recommended)
- **Git**
- **8 GB RAM** minimum (16 GB recommended for large cohorts)
- **2 GB disk space**

---

## Step 1: Install the Package

### Option A: Install from PyPI (Recommended)

```bash
pip install pathway-subtyping
```

### Option B: Install from Source

```bash
# Clone the repository
git clone https://github.com/topmist-admin/pathway-subtyping-framework.git
cd pathway-subtyping-framework

# Create virtual environment
python3 -m venv pathwayenv
source pathwayenv/bin/activate  # Linux/macOS
# or: pathwayenv\Scripts\activate  # Windows

# Install in development mode
pip install -e ".[dev]"
```

---

## Step 2: Verify Installation

```bash
# Check CLI is available
pathway-subtyping --version

# Or use the short alias
psf --help
```

Expected output:
```
pathway-subtyping 0.1.0
```

---

## Step 3: Run the Demo Pipeline

The demo uses a synthetic 60-sample dataset with planted ground truth subtypes.

```bash
# Run with the test configuration
pathway-subtyping --config configs/test_synthetic.yaml
```

Or using Python:
```bash
python -m pathway_subtyping --config configs/test_synthetic.yaml
```

**Expected runtime:** 1-3 minutes on a standard laptop.

---

## Step 4: Check Outputs

After the pipeline completes, outputs are in `outputs/synthetic_test/`:

```
outputs/synthetic_test/
├── pathway_scores.csv        # Sample x Pathway disruption scores
├── subtype_assignments.csv   # Cluster assignments with confidence
├── report.json               # Machine-readable results
└── report.md                 # Human-readable summary
```

### Key Files to Review

1. **`report.md`** - Start here for a summary of results
2. **`subtype_assignments.csv`** - Sample cluster assignments
3. **`pathway_scores.csv`** - Pathway-level scores per sample

---

## Step 5: Understand the Results

### Pipeline Flow

```
VCF Variants → Gene Burdens → Pathway Scores → GMM Clustering → Validation Gates
```

### Validation Gates

The pipeline runs three validation checks to ensure clustering quality:

| Gate | Purpose | Pass Condition |
|------|---------|----------------|
| **Label Shuffle** | Verify clusters aren't spurious | Null ARI < 0.1 |
| **Random Gene Sets** | Verify pathways drive clustering | Random ARI < 0.1 |
| **Bootstrap Stability** | Verify clusters are robust | Stability ARI >= 0.7 |

### Output Interpretation

See [outputs_dictionary.md](outputs_dictionary.md) for detailed guidance on:
- What each output file contains
- How to interpret pathway scores
- What NOT to infer from results

---

## Step 6: Run on Your Own Data

### Prepare Your Data

1. **VCF file** with annotated variants (GENE, CONSEQUENCE, CADD fields)
2. **Phenotypes CSV** (optional) with sample metadata
3. **Pathways GMT file** (or use provided defaults)

### Create a Configuration File

Copy and modify the example config:

```bash
cp configs/test_synthetic.yaml configs/my_analysis.yaml
```

Edit `my_analysis.yaml`:

```yaml
run_name: my_cohort_analysis

input:
  vcf_path: path/to/your/variants.vcf
  phenotypes_path: path/to/phenotypes.csv  # optional
  pathways_path: data/pathways/your_pathways.gmt

output:
  output_dir: outputs/my_cohort

pipeline:
  seed: 42
  min_samples_per_cluster: 10

clustering:
  n_clusters_range: [2, 8]

validation:
  run_validation: true
```

### Run Your Analysis

```bash
pathway-subtyping --config configs/my_analysis.yaml
```

---

## Quick Command Reference

| Command | Description |
|---------|-------------|
| `pathway-subtyping --config FILE` | Run pipeline with config |
| `pathway-subtyping --version` | Show version |
| `psf --help` | Show help (short alias) |
| `pytest tests/` | Run test suite |

---

## Common Issues

### Pipeline fails to start

```bash
# Ensure you're in the virtual environment
source pathwayenv/bin/activate

# Verify installation
pathway-subtyping --version
```

### Module not found errors

```bash
# Reinstall the package
pip install -e .
```

### Memory errors

- Use a machine with more RAM (16+ GB recommended)
- For large cohorts, see [Performance Optimization](api/index.md#performance-utilities)

See [troubleshooting.md](troubleshooting.md) for more solutions.

---

## Next Steps

| Resource | Description |
|----------|-------------|
| [Adapting for Your Disease](guides/adapting-for-your-disease.md) | Customize pathways for your condition |
| [Validation Gates Guide](guides/validation-gates.md) | Understanding validation results |
| [API Reference](api/index.md) | Python API documentation |
| [Data Formats](data_formats.md) | Input/output specifications |

---

## Getting Help

- **Issues:** [GitHub Issues](https://github.com/topmist-admin/pathway-subtyping-framework/issues)
- **Discussions:** [GitHub Discussions](https://github.com/topmist-admin/pathway-subtyping-framework/discussions)

---

*Last updated: January 2026*
