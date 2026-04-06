# Benchmark Datasets

This directory contains benchmark datasets for validating the pathway subtyping framework.

## Available Benchmarks

### 1. 47-Dataset Bootstrap Threshold Calibration Benchmark

**Zenodo DOI:** 10.5281/zenodo.19324360
**URL:** https://zenodo.org/records/19324360

**Description:** 47 publicly available datasets (12 TCGA + 35 GEO) spanning oncology, psychiatry, immunology, single-cell RNA-seq, and other domains. Used to calibrate the adaptive bootstrap stability threshold model (Methods §1.5 in manuscript).

**Results:** 39/47 PASS (83.0%): 12/12 TCGA (100%), 27/35 GEO (77.1%). 8 GEO failures are data-quality issues (insufficient samples, missing values, format incompatibility), not framework deficiencies.

**Files in this directory:**
- `dataset_candidate_list_47.csv` — 47 datasets with accession IDs, domains, sources, selection justifications
- `benchmark_47datasets_manifest.csv` — Complete metadata manifest for automation

**Results files (in `research-results/benchmarks/`):**
- `bootstrap_threshold_calibration_47datasets_zenodo.csv` — Final validated results (silhouette + bootstrap ARI for all 47)
- `ZENODO_47DATASETS_README.md` — Comprehensive documentation

**Selection criteria:** See `archive/benchmark-calibration-2026-03/BENCHMARK_SELECTION_CRITERIA.md`

**Circularity prevention:** None of the 47 benchmarks include TCGA-COAD, GSE28521, GSE64018, or GSE80655 (the manuscript analysis datasets).

---

### 2. Synthetic 4-Subtype Benchmark

**Location:** `data/sample/synthetic_cohort.vcf`, `data/sample/synthetic_phenotypes.csv`

**Description:** 60 synthetic samples with 4 planted subtypes (synaptic, chromatin, ion_channel, mtor), each with 15 samples.

**Expected Results:**
- Optimal clusters: 4
- ARI vs planted subtypes: > 0.95
- Validation gates: All PASS

**Run:**
```bash
psf --config configs/test_synthetic.yaml
```

---

## Expected Results Summary

| Benchmark | Samples | Planted Subtypes | Expected ARI | Validation |
|-----------|---------|------------------|--------------|------------|
| synthetic_4subtype | 60 | 4 | > 0.95 | PASS |

---

## Creating New Benchmarks

To create a benchmark dataset:

1. **Generate synthetic VCF** with known subtype structure
2. **Create phenotype CSV** with `sample_id` and `planted_subtype` columns
3. **Run pipeline** and record expected metrics
4. **Document** expected results in this README

### Benchmark Validation Script

```python
from pathway_subtyping import DemoPipeline, PipelineConfig
from sklearn.metrics import adjusted_rand_score

# Run benchmark
config = PipelineConfig.from_yaml("configs/test_synthetic.yaml")
pipeline = DemoPipeline(config)
pipeline.run()

# Validate against expected
assignments = pipeline.cluster_assignments
ari = adjusted_rand_score(
    assignments["planted_subtype"],
    assignments["cluster_label"]
)

assert ari > 0.95, f"ARI {ari} below expected 0.95"
assert pipeline.validation_result.all_passed, "Validation gates failed"

print(f"Benchmark PASSED: ARI = {ari:.4f}")
```

---

## Benchmark Results Log

Document benchmark results here for tracking:

| Date | Version | Benchmark | ARI | Validation | Notes |
|------|---------|-----------|-----|------------|-------|
| 2026-01-29 | 0.1.0 | synthetic_4subtype | 1.000 | PASS | Initial release |

---

## Data Provenance

**All benchmark data is computationally generated using the `SyntheticDataGenerator` class with fixed random seeds.** No real patient, clinical, or proprietary data is used in any benchmark.

- The synthetic VCF files contain **randomly generated variant calls** with no connection to real individuals
- The synthetic phenotype CSVs contain **randomly assigned sample metadata** (no real clinical information)
- Gene symbols in the data (e.g., SHANK3, CHD8) are standard HGNC identifiers used as labels only
- All benchmark results are fully reproducible from the random seed alone

**No data from any employer, client, institution, or commercial entity was used to create these benchmarks.**

For full provenance details, see [DISCLAIMER.md](../../DISCLAIMER.md).
