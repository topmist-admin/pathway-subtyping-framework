# Benchmark Datasets

This directory contains benchmark datasets for validating the pathway subtyping framework.

## Available Benchmarks

### 1. Synthetic 4-Subtype Benchmark

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
