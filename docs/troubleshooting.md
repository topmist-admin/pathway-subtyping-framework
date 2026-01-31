# Troubleshooting Guide

> **RESEARCH USE ONLY** — This framework is for research purposes only. Not for clinical decision-making. See [DISCLAIMER.md](../DISCLAIMER.md).

This guide covers common issues when running the Pathway Subtyping Framework and their solutions.

---

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Environment Setup](#environment-setup)
3. [Pipeline Errors](#pipeline-errors)
4. [Data Issues](#data-issues)
5. [Reproducibility Issues](#reproducibility-issues)
6. [Performance Issues](#performance-issues)
7. [Validation Gate Failures](#validation-gate-failures)
8. [Platform-Specific Issues](#platform-specific-issues)
9. [Getting Help](#getting-help)

---

## Installation Issues

### `pip install` fails with dependency conflicts

**Symptom:**
```
ERROR: Cannot install pathway-subtyping because these package versions have conflicting dependencies.
```

**Solution:**
1. Create a fresh virtual environment:
   ```bash
   python -m venv pathwayenv
   source pathwayenv/bin/activate  # Linux/macOS
   # or: pathwayenv\Scripts\activate  # Windows
   ```

2. Upgrade pip and install:
   ```bash
   pip install --upgrade pip
   pip install pathway-subtyping
   ```

### Missing system dependencies on Linux

**Symptom:**
```
error: command 'gcc' failed with exit status 1
```

**Solution:**
Install build dependencies:
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install build-essential python3-dev libbz2-dev liblzma-dev

# CentOS/RHEL
sudo yum groupinstall "Development Tools"
sudo yum install python3-devel bzip2-devel xz-devel
```

### scikit-learn installation issues

**Symptom:**
```
ImportError: cannot import name 'GaussianMixture' from 'sklearn.mixture'
```

**Solution:**
```bash
pip install --upgrade scikit-learn
```

---

## Environment Setup

### Python version mismatch

**Symptom:**
```
This package requires Python >=3.9
```

**Solution:**
1. Check your Python version:
   ```bash
   python --version
   ```

2. Install Python 3.9+ using pyenv:
   ```bash
   pyenv install 3.11.0
   pyenv local 3.11.0
   ```

### Module import errors after installation

**Symptom:**
```python
ModuleNotFoundError: No module named 'pathway_subtyping'
```

**Solution:**
Install in editable mode (for development):
```bash
pip install -e .
```

Or reinstall:
```bash
pip uninstall pathway-subtyping
pip install pathway-subtyping
```

### CLI not found

**Symptom:**
```
command not found: pathway-subtyping
```

**Solution:**
1. Ensure the virtual environment is activated
2. Check if scripts are in PATH:
   ```bash
   which pathway-subtyping
   # or
   python -m pathway_subtyping --help
   ```

---

## Pipeline Errors

### FileNotFoundError for input files

**Symptom:**
```
FileNotFoundError: VCF file not found: path/to/file.vcf
```

**Solution:**
1. Verify the file exists:
   ```bash
   ls -la path/to/file.vcf
   ```

2. Use absolute paths in config:
   ```yaml
   input:
     vcf_path: /full/path/to/variants.vcf
   ```

3. Check for typos in the config file

### Configuration file not found

**Symptom:**
```
FileNotFoundError: [Errno 2] No such file or directory: 'configs/my_config.yaml'
```

**Solution:**
Run from the repository root or use absolute path:
```bash
pathway-subtyping --config /full/path/to/config.yaml
```

### Invalid YAML syntax

**Symptom:**
```
yaml.scanner.ScannerError: mapping values are not allowed here
```

**Solution:**
1. Check indentation (use spaces, not tabs)
2. Validate YAML syntax:
   ```bash
   python -c "import yaml; yaml.safe_load(open('config.yaml'))"
   ```

### Memory error during clustering

**Symptom:**
```
MemoryError: Unable to allocate array
```

**Solution:**
1. Close other applications to free memory
2. Use chunked processing for large datasets:
   ```python
   from pathway_subtyping.utils import compute_gene_burdens_chunked
   ```
3. Reduce dataset size or use machine with more RAM

---

## Data Issues

### VCF parsing errors

**Symptom:**
```
ValueError: Could not parse VCF line
```

**Solution:**
1. Verify VCF format is valid:
   ```bash
   bcftools view -h your_file.vcf | head -20
   ```

2. Check for required INFO fields (GENE, CONSEQUENCE):
   ```bash
   bcftools query -f '%INFO/GENE\n' your_file.vcf | head
   ```

3. Use the annotation helper script if fields are missing:
   ```bash
   python scripts/annotate_vcf.py --help
   ```

### Missing gene annotations

**Symptom:**
```
WARNING: 80% of variants have no gene annotation
```

**Solution:**
1. Annotate your VCF with VEP or ANNOVAR first
2. Use the provided annotation helper:
   ```bash
   python scripts/annotate_vcf.py \
     --vcf input.vcf \
     --annotations vep_output.tsv \
     --output annotated.vcf
   ```

### GMT file format errors

**Symptom:**
```
ValueError: Invalid GMT format
```

**Solution:**
1. Check GMT format (tab-separated):
   ```
   PATHWAY_NAME<TAB>DESCRIPTION<TAB>GENE1<TAB>GENE2<TAB>...
   ```

2. Verify no empty lines or malformed entries:
   ```bash
   awk -F'\t' 'NF < 3 {print NR": "$0}' pathways.gmt
   ```

### No samples pass QC

**Symptom:**
```
ERROR: No samples remaining after QC filtering
```

**Solution:**
1. Check your VCF has variant calls (not all reference)
2. Verify sample IDs match between VCF and phenotypes
3. Lower QC thresholds if appropriate for your data

---

## Reproducibility Issues

### Different results on different machines

**Symptom:**
Outputs differ between runs or machines despite using the same seed.

**Diagnosis:**
1. Check Python version:
   ```bash
   python --version
   ```

2. Check package versions:
   ```bash
   pip freeze | grep -E "(numpy|pandas|scikit-learn)"
   ```

3. Verify seed is set in config:
   ```yaml
   pipeline:
     seed: 42
   ```

**Solution:**
1. Use identical Python and package versions
2. Set environment variables for full determinism:
   ```bash
   export PYTHONHASHSEED=42
   export OMP_NUM_THREADS=1
   ```

### Results differ across runs

**Symptom:**
Same config produces different results.

**Solution:**
1. Ensure seed is set (not null):
   ```yaml
   pipeline:
     seed: 42  # Must be an integer, not null
   ```

2. Clear cached data:
   ```bash
   rm -rf outputs/<run_name>
   ```

---

## Performance Issues

### Pipeline runs very slowly

**Symptom:**
Analysis takes much longer than expected.

**Solution:**
1. Check dataset size vs. available RAM
2. Use chunked processing for large VCFs:
   ```python
   from pathway_subtyping.utils import chunked_vcf_reader
   ```

3. Reduce validation iterations:
   ```yaml
   validation:
     n_permutations: 50  # Reduce from 100
     n_bootstrap: 25     # Reduce from 50
   ```

4. Narrow cluster search range:
   ```yaml
   clustering:
     n_clusters_range: [2, 5]  # Instead of [2, 10]
   ```

### High memory usage

**Symptom:**
System becomes unresponsive or swaps heavily.

**Solution:**
1. Monitor memory:
   ```bash
   watch -n 1 free -h
   ```

2. Use memory-efficient utilities:
   ```python
   from pathway_subtyping.utils import estimate_memory_usage
   estimate_memory_usage(n_samples=1000, n_genes=20000)
   ```

3. Consider downsampling for initial exploration:
   ```python
   from pathway_subtyping.utils import downsample_cohort
   ```

---

## Validation Gate Failures

### Label Shuffle fails

**Symptom:**
```
Label Shuffle: FAIL (ARI = 0.25)
```

**Meaning:**
Clustering is finding structure in random data, suggesting overfitting.

**Investigation:**
1. Check if cluster count is too high
2. Review if dataset is too small
3. Consider using fewer pathways

### Random Gene Sets fails

**Symptom:**
```
Random Gene Sets: FAIL (ARI = 0.18)
```

**Meaning:**
Random gene sets produce similar cluster structure.

**Investigation:**
1. Check pathway overlap with cohort genes
2. Verify pathway database integrity
3. This can fail on small datasets

### Bootstrap Stability fails

**Symptom:**
```
Bootstrap Stability: FAIL (ARI = 0.55)
```

**Meaning:**
Clusters are not robust to resampling.

**Investigation:**
1. Dataset may be too small for stable clustering
2. Try reducing the number of clusters
3. Check for outlier samples

---

## Platform-Specific Issues

### macOS

**Issue:** OpenMP library conflict
```
OMP: Error #15: Initializing libiomp5.dylib...
```

**Solution:**
```bash
export KMP_DUPLICATE_LIB_OK=TRUE
```

**Issue:** SSL certificate errors
```
ssl.SSLCertVerificationError: certificate verify failed
```

**Solution:**
```bash
/Applications/Python\ 3.11/Install\ Certificates.command
```

### Linux (Ubuntu)

**Issue:** Missing BLAS libraries
```
numpy.linalg.LinAlgError: LAPACK library not found
```

**Solution:**
```bash
sudo apt-get install libopenblas-dev liblapack-dev
```

### Windows

**Issue:** Path issues with backslashes
```
FileNotFoundError: path\to\file.vcf
```

**Solution:**
Use forward slashes or raw strings in config:
```yaml
input:
  vcf_path: "C:/Users/name/data/variants.vcf"
```

---

## Getting Help

### Before asking for help

1. **Check this guide** for your specific error
2. **Search existing issues** on GitHub
3. **Collect diagnostic information:**
   ```bash
   python --version
   pip freeze | grep -E "(pathway|numpy|pandas|sklearn)"
   cat outputs/<run>/pipeline.log | tail -50
   ```

### Reporting issues

When opening a GitHub issue, include:

1. **Environment:**
   - OS and version
   - Python version
   - Package versions (`pip freeze`)

2. **Steps to reproduce:**
   - Exact commands run
   - Configuration file used

3. **Error message:**
   - Full traceback
   - Relevant log output

4. **Expected vs actual behavior**

### Community resources

- **GitHub Issues:** [Report bugs](https://github.com/topmist-admin/pathway-subtyping-framework/issues)
- **Discussions:** [Ask questions](https://github.com/topmist-admin/pathway-subtyping-framework/discussions)

---

## Quick Reference

### Common Commands

```bash
# Verify installation
pathway-subtyping --version

# Run pipeline
pathway-subtyping --config configs/test_synthetic.yaml

# Run tests
pytest tests/ -v

# Check environment
python -c "import pathway_subtyping; print(pathway_subtyping.__version__)"
```

### Environment Variables

| Variable | Purpose | Example |
|----------|---------|---------|
| `PYTHONHASHSEED` | Hash randomization seed | `42` |
| `OMP_NUM_THREADS` | OpenMP threads | `1` |
| `KMP_DUPLICATE_LIB_OK` | macOS OpenMP fix | `TRUE` |

---

*Last updated: January 2026*
