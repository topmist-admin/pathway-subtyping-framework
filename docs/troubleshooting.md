# Troubleshooting Guide

> **RESEARCH USE ONLY** — This framework is for research purposes only. Not for clinical decision-making. See [DISCLAIMER.md](../DISCLAIMER.md).

This guide covers common issues when running the Pathway Subtyping Framework and their solutions.

---

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Environment Setup](#environment-setup)
3. [Pipeline Errors](#pipeline-errors)
4. [Data Issues](#data-issues)
5. [Ancestry Correction Issues](#ancestry-correction-issues)
6. [Batch Correction Issues](#batch-correction-issues)
7. [Sensitivity Analysis Issues](#sensitivity-analysis-issues)
8. [Reproducibility Issues](#reproducibility-issues)
7. [Performance Issues](#performance-issues)
8. [Validation Gate Failures](#validation-gate-failures)
9. [Platform-Specific Issues](#platform-specific-issues)
10. [Getting Help](#getting-help)

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

### Google Colab installation

**Symptom:**
```
ValueError: numpy.dtype size changed, may indicate binary incompatibility
```

**Cause:**
This happens when an older version of the package (< 0.2.2) is installed. Versions before 0.2.2 included `pysam` as a required dependency, which forced NumPy to downgrade from 2.x to 1.x, breaking Colab's pre-compiled scipy/scikit-learn binaries.

**Solution:**
Install version 0.2.2 or later:
```bash
pip install pathway-subtyping>=0.2.2
```

If you see `v0.2.0` in the install output, force the latest:
```bash
pip install pathway-subtyping==0.2.2
```

No runtime restart is needed with v0.2.2+.

### VCF processing with pysam

**Note:** As of v0.2.2, `pysam` is an optional dependency. Install it when you need to process VCF files:
```bash
pip install pathway-subtyping[vcf]
```

The core framework (simulation, clustering, validation, benchmarks) works without pysam.

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

### Validating Your VCF Before Running

Before running the pipeline, validate your VCF has the required annotations:

```bash
# Quick validation
python scripts/annotate_vcf.py your_file.vcf --validate-only

# Detailed validation with statistics
python scripts/annotate_vcf.py your_file.vcf --validate-only --detailed
```

Or programmatically:
```python
from pathway_subtyping import validate_vcf_for_pipeline
is_valid, report, suggestions = validate_vcf_for_pipeline("your_file.vcf")
```

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
or
```
VCFDataQualityError: VCF data quality is insufficient for analysis
```

**Solution:**
1. Check your annotation coverage:
   ```bash
   python scripts/annotate_vcf.py your_file.vcf --validate-only
   ```

2. Annotate with VEP:
   ```bash
   # Install VEP (if not already installed)
   # See: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html

   # Run VEP
   vep -i input.vcf -o output.tsv \
     --tab --symbol --pick \
     --fields "Uploaded_variation,Gene,SYMBOL,Consequence,CADD_PHRED"

   # Apply annotations
   python scripts/annotate_vcf.py input.vcf annotated.vcf --vep-tsv output.tsv
   ```

3. Or annotate with ANNOVAR:
   ```bash
   # Convert VCF to ANNOVAR input format
   convert2annovar.pl -format vcf4 input.vcf > input.avinput

   # Run ANNOVAR
   table_annovar.pl input.avinput humandb/ \
     -buildver hg38 \
     -out output \
     -protocol refGene,cadd \
     -operation g,f

   # Apply annotations
   python scripts/annotate_vcf.py input.vcf annotated.vcf --annovar-file output.hg38_multianno.txt
   ```

### Multi-allelic variants

**Symptom:**
```
INFO: Expanded 150 multi-allelic variants to 320 bi-allelic records
```

**Explanation:**
The framework automatically expands multi-allelic variants (variants with multiple alternate alleles) into separate bi-allelic records. This is logged for transparency.

**Example:**
```
# Original multi-allelic
chr1  100  .  A  G,T  .  .  .

# Becomes two bi-allelic records
chr1  100  .  A  G  .  .  .
chr1  100  .  A  T  .  .  .
```

This is expected behavior and requires no action.

### Low CADD score coverage

**Symptom:**
```
WARNING: Low CADD score coverage (25.0%). Consider adding CADD annotations.
```

**Solution:**
CADD scores improve variant weighting but are optional. To add them:

1. Download CADD scores: https://cadd.gs.washington.edu/download
2. Annotate with VEP using CADD plugin, or
3. Use ANNOVAR with CADD database

### GMT file format errors

**Symptom:**
```
ConfigValidationError: GMT file has 2 error(s):
Line 5: Expected at least 3 tab-separated fields, got 2
Line 12: Pathway 'SYNAPTIC' has fewer than 2 genes
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

3. Validate the GMT file programmatically:
   ```python
   from pathway_subtyping.config import validate_gmt_file, ConfigValidationError

   try:
       pathways = validate_gmt_file("pathways.gmt")
       print(f"Valid: {len(pathways)} pathways")
   except ConfigValidationError as e:
       print(f"Errors:\n{e}")
   ```

4. Validate gene symbols match your VCF annotations

### Configuration validation errors

**Symptom:**
```
ConfigValidationError: Missing required field: data.vcf_path
Field: data.vcf_path

Suggested fixes:
  1. Add 'vcf_path: /path/to/file' under the 'data:' section
```

**Solution:**
The error message includes the field that's missing and suggested fixes. Common issues:

1. **Missing required fields**: Add the specified field to your config
2. **Invalid seed value**: Seed must be an integer, not a string
3. **Invalid cluster range**: Ensure min ≥ 2 and max ≥ min

**Example fix:**
```yaml
# Before (missing required field)
data:
  phenotype_path: pheno.csv

# After (field added)
data:
  vcf_path: variants.vcf        # Added
  phenotype_path: pheno.csv
  pathway_db: pathways.gmt      # Also required
```

### Insufficient pathways after filtering

**Symptom:**
```
ValueError: Insufficient pathways after filtering: 1 remaining.
Need at least 2 pathways with non-zero variance for clustering.
```

**Cause:**
Pathways are filtered out during processing for these reasons:
- No genes from the pathway found in the burden data
- Pathway has fewer than 2 genes overlapping with burden data
- Pathway scores have zero variance (all samples have same score)

**Solution:**
1. Check pathway-gene overlap:
   ```python
   # See which genes from pathways are in your data
   pathway_genes = set(gene for genes in pathways.values() for gene in genes)
   burden_genes = set(gene_burdens.columns)
   overlap = pathway_genes & burden_genes
   print(f"Overlap: {len(overlap)} genes")
   ```

2. Use broader pathway definitions with more genes

3. Ensure pathways contain genes that vary across samples

4. Review the warning logs for details on filtered pathways:
   ```
   WARNING: Removing 3 zero-variance pathway(s): ['PATHWAY_A', 'PATHWAY_B', ...]
   ```

### No samples pass QC

**Symptom:**
```
ERROR: No samples remaining after QC filtering
```

**Solution:**
1. Check your VCF has variant calls (not all reference)
2. Verify sample IDs match between VCF and phenotypes:
   ```bash
   # List VCF samples
   bcftools query -l your_file.vcf

   # Compare with phenotypes
   head -1 phenotypes.csv
   ```
3. Lower QC thresholds if appropriate for your data

### Data Quality Report Interpretation

The pipeline generates a data quality report with these key metrics:

| Metric | Threshold | Action if below |
|--------|-----------|-----------------|
| Gene coverage | ≥50% | Annotate with VEP/ANNOVAR |
| Consequence coverage | ≥50% | Re-run annotation |
| CADD coverage | ≥30% | Optional - add CADD scores |

Example output:
```
Data Quality Report
========================================
Total variants: 1000
Parsed variants: 998
Skipped variants: 2

Annotation Coverage:
  GENE: 95.0% (950/998)
  CONSEQUENCE: 90.0% (898/998)
  CADD: 80.0% (798/998)

Multi-allelic variants: 50 (expanded to 110)

Data Quality Status: PASS
```

---

## Ancestry Correction Issues

### Ancestry independence gate fails

**Symptom:**
```
Ancestry Independence: FAIL (min_p_value = 0.001)
```

**Meaning:**
Cluster assignments are significantly associated with ancestry principal components, suggesting population stratification is confounding the subtypes.

**Solution:**
1. Enable ancestry correction in your config:
   ```yaml
   ancestry:
     pcs_path: data/ancestry_pcs.csv
     correction: regress_out
     n_pcs: 10
   ```

2. Or correct programmatically:
   ```python
   from pathway_subtyping import compute_ancestry_pcs, adjust_pathway_scores

   pcs = compute_ancestry_pcs(genotype_matrix, n_components=10, seed=42)
   result = adjust_pathway_scores(pathway_scores, pcs)
   # Use result.adjusted_scores for clustering
   ```

3. If correction is already enabled, try increasing `n_pcs` (e.g., 15 or 20)

### Many pathways flagged as confounded

**Symptom:**
```
WARNING: 8 of 10 pathways have R² > 0.1 with ancestry PCs
```

**Meaning:**
Most pathway scores have substantial ancestry-correlated variance. This is common in multi-ethnic cohorts.

**Investigation:**
1. Check ancestry correlation details:
   ```python
   from pathway_subtyping import compute_ancestry_correlation
   corr = compute_ancestry_correlation(pathway_scores, ancestry_pcs)
   print(corr)  # Shows per-pathway, per-PC Pearson r
   ```

2. Verify correction was applied (R² should decrease after adjustment)
3. This is expected when ancestry groups have different genetic burden profiles — the correction handles it

### Ancestry PCs file format errors

**Symptom:**
```
ValueError: Ancestry PCs file must have sample IDs as first column
```

**Solution:**
The ancestry PCs CSV must have:
- First column: sample IDs matching your VCF
- Remaining columns: PC values (PC1, PC2, ...)

```csv
sample_id,PC1,PC2,PC3,...
SAMPLE_001,0.012,-0.034,0.008,...
SAMPLE_002,-0.015,0.022,0.001,...
```

Generate PCs using `compute_ancestry_pcs()` or external tools (PLINK, EIGENSOFT).

### Too few samples for stratified analysis

**Symptom:**
```
WARNING: Skipping ancestry group 'AFR' with 8 samples (minimum 20 required)
```

**Meaning:**
Ancestry groups with too few samples are excluded from stratified analysis because GMM clustering requires sufficient samples.

**Solution:**
1. Merge small ancestry groups if biologically appropriate
2. Use `adjust_pathway_scores()` (regression) instead of `stratified_analysis()` for small groups
3. Collect more samples from underrepresented populations

### sklearn warnings during ancestry PCA

**Symptom:**
```
RuntimeWarning: invalid value encountered in matmul
```

**Meaning:**
Numerical precision issues in PCA computation, typically from near-monomorphic variants.

**Solution:**
These warnings are generally benign — the framework automatically filters monomorphic variants. The PCA results remain valid. To suppress:
```python
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
```

---

## Batch Correction Issues

### Batch effects persist after correction

**Symptom:** `validate_batch_correction()` shows `batch_variance_reduced: False`.

**Solutions:**
1. Try a different correction method:
   ```python
   # If ComBat didn't work, try standardize
   result = correct_batch_effects(scores, batch_labels, method=BatchCorrectionMethod.STANDARDIZE)
   ```
2. Check if batches are confounded with biology — if biological groups are entirely within one batch, correction will remove signal
3. Verify batch labels are correct (misassigned labels will worsen the problem)

### ComBat produces NaN values

**Symptom:** NaN values in `corrected_scores` after ComBat correction.

**Solutions:**
1. Check for zero-variance pathways — the framework handles these, but very small batches (< 3 samples) may cause instability
2. Use `MEAN_CENTER` instead for very small batch sizes
3. Ensure at least 2 batches with ≥ 3 samples each

### Not sure if batch effects exist

**Solution:** Run detection first before correcting:
```python
report = detect_batch_effects(pathway_scores, batch_labels)
print(report.format_report())
# Only correct if report.overall_batch_effect is True
```

---

## Sensitivity Analysis Issues

### Sensitivity analysis shows low stability

**Symptom:** `run_sensitivity_analysis()` returns `is_robust: False`.

**Solutions:**
1. This may indicate genuine instability — consider whether your data supports clear subtypes
2. Check which parameter is most sensitive:
   ```python
   result = run_sensitivity_analysis(scores, n_clusters=3, seed=42)
   print(f"Most sensitive: {result.most_sensitive_parameter}")
   ```
3. If `n_clusters` is most sensitive, your optimal K may be uncertain — use `select_n_clusters()` to verify
4. If `feature_subset` is most sensitive, some pathways may be dominating — review pathway definitions

### Sensitivity analysis is slow

**Symptom:** `run_sensitivity_analysis()` takes too long.

**Solutions:**
1. Run only specific parameter axes:
   ```python
   result = run_sensitivity_analysis(
       scores, n_clusters=3, seed=42,
       parameters=[SensitivityParameter.CLUSTERING_ALGORITHM, SensitivityParameter.N_CLUSTERS]
   )
   ```
2. Reduce data size for exploratory analysis, then run full sensitivity on final dataset

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

### Threshold too strict or too lenient for your dataset

**Symptom:**
Gates fail on reasonably structured data, or pass on clearly artifactual data.

**Meaning:**
The default thresholds (0.15/0.8) don't account for sample size and cluster count. Small datasets produce noisier ARI distributions, while more clusters inflate chance ARI.

**Solution:**
Enable auto-calibration in your config:
```yaml
validation:
  calibrate: true
  stability_threshold: null   # Auto-calibrate based on n_samples and n_clusters
  null_ari_max: null           # Auto-calibrate
```

Or calibrate programmatically:
```python
from pathway_subtyping import calibrate_thresholds

ct = calibrate_thresholds(n_samples=150, n_clusters=4)
print(f"Recommended null ARI: {ct.null_ari_threshold:.4f}")
print(f"Recommended stability: {ct.stability_threshold:.4f}")
```

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

### GMM convergence warnings

**Symptom:**
```
WARNING: GMM did not converge during transfer learning validation
```
or
```
WARNING: Label shuffle test: No GMM fits converged
```

**Meaning:**
The Gaussian Mixture Model optimization did not converge within the maximum iterations.

**Investigation:**
1. **Small datasets**: GMM may struggle with very small sample sizes (<20)
2. **Numerical instability**: Data may have numerical issues (NaN, Inf values)
3. **Too many clusters**: The cluster count may be too high for the data

**Solution:**
1. The framework uses `reg_covar=1e-6` for numerical stability
2. If warnings persist, check your data for:
   ```python
   import numpy as np
   print(f"NaN values: {np.isnan(pathway_scores).sum().sum()}")
   print(f"Inf values: {np.isinf(pathway_scores).sum().sum()}")
   ```
3. Try reducing the cluster count
4. Increase sample size if possible

**Note:** The validation framework handles non-converged fits gracefully by skipping them. If all fits fail to converge, the test returns a failing result with diagnostic information.

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

*Last updated: February 2026*
