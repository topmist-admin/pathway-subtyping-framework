# Cross-Cohort Validation Guide

## Why Cross-Cohort Validation?

Molecular subtypes discovered in a single dataset may reflect real biological groupings — or they may be artifacts of batch effects, sampling bias, or noise. Cross-cohort validation tests whether subtypes **replicate** in an independent dataset, which is essential before publishing or acting on subtype definitions.

If subtypes are real, a model trained on cohort A should recover similar groupings in cohort B, even though the samples are different.

## How It Works

The `cross_cohort` module compares subtype definitions across two cohorts using three complementary metrics:

1. **Transfer ARI**: Train a GMM on cohort A's pathway scores, apply it to cohort B, and compare the predicted labels to cohort B's discovered labels using the Adjusted Rand Index.

2. **Projection ARI**: Project both cohorts into a shared PCA space, assign cohort B samples to the nearest cohort A centroid, and compare assignments to cohort B's labels.

3. **Pathway Correlation**: Compute the Pearson correlation of pathway variances between cohorts. High correlation means the same pathways drive variation in both datasets.

## Quick Start

```python
from pathway_subtyping import (
    generate_synthetic_cohort_pair,
    compare_cohorts,
)

# Generate two synthetic cohorts with matching subtype structure
cohort_a, cohort_b = generate_synthetic_cohort_pair(seed=42)

# Compare them
result = compare_cohorts(cohort_a, cohort_b, seed=42)

# View results
print(result.format_report())
```

## Interpreting Results

| Metric | Threshold | Interpretation |
|--------|-----------|----------------|
| Transfer ARI | > 0.5 | Good replication — subtypes are consistent |
| Transfer ARI | 0.3 – 0.5 | Moderate — partially shared subtypes |
| Transfer ARI | < 0.3 | Weak — subtypes may not replicate |
| Projection ARI | > 0.5 | Cluster geometry preserved in shared space |
| Projection ARI | < 0.3 | Cluster structure differs between cohorts |
| Pathway Correlation | > 0.7 | Similar pathway importance |
| Pathway Correlation | 0.3 – 0.7 | Partially shared pathway effects |
| Pathway Correlation | < 0.3 | Different pathways drive each cohort |

**Tip**: All three metrics should be reported together. A high Transfer ARI with low Pathway Correlation may indicate that the clusters transfer but for different biological reasons.

## Real-World Workflow

When comparing actual cohorts (not synthetic data):

### Step 1: Run the pipeline on each cohort separately

```python
from pathway_subtyping import DemoPipeline, PipelineConfig

config_a = PipelineConfig.from_yaml("configs/cohort_a.yaml")
pipeline_a = DemoPipeline(config_a)
pipeline_a.run()

config_b = PipelineConfig.from_yaml("configs/cohort_b.yaml")
pipeline_b = DemoPipeline(config_b)
pipeline_b.run()
```

### Step 2: Load results and compare

```python
from pathway_subtyping import load_cohort_result, compare_cohorts

cohort_a = load_cohort_result("outputs/cohort_a/")
cohort_b = load_cohort_result("outputs/cohort_b/")

result = compare_cohorts(cohort_a, cohort_b)
print(result.format_report())
```

### Step 3: Generate a report

```python
from pathway_subtyping import generate_cross_cohort_report

generate_cross_cohort_report([result], "outputs/cross_cohort_report.md")
```

### Step 4: Batch comparison (3+ cohorts)

```python
from pathway_subtyping import batch_compare_cohorts

results = batch_compare_cohorts(
    ["outputs/cohort_a/", "outputs/cohort_b/", "outputs/cohort_c/"],
    report_path="outputs/batch_cross_cohort_report.md",
)
```

## Common Pitfalls

- **Different pathway sets**: The comparison uses only pathways present in both cohorts. If cohorts have very different pathway annotations, the overlap may be too small for meaningful comparison. Check `result.details["common_pathways"]`.

- **Very different sample sizes**: Large imbalances (e.g., 1000 vs 20 samples) can affect GMM training quality. The Transfer ARI metric is most reliable when cohorts are within 5x of each other in size.

- **Batch effects**: If cohorts were processed differently (different sequencing platforms, pipelines, normalization), apparent replication failures may reflect technical rather than biological differences. Apply batch correction first using `correct_batch_effects()`.

- **Different numbers of clusters**: The transfer learning approach trains on cohort A's cluster count. If cohort B naturally has a different number of subtypes, Transfer ARI will be low even if some subtypes overlap.

## Example Report Output

```
Cross-Cohort Validation Report
========================================
Cohort A: cohort_a
Cohort B: cohort_b

Metrics:
  Transfer ARI:      0.650
  Projection ARI:    0.580
  Pathway Correlation: 0.820

Transfer ARI > 0.5: Good replication of subtypes.
Pathway correlation > 0.7: Similar pathway importance across cohorts.
Shared subtype labels: synaptic

Details:
  Common pathways:  10
  Cohort A samples: 100
  Cohort B samples: 80
```

## CLI Example

A standalone example script is included:

```bash
python scripts/cross_cohort_example.py --output outputs/cross_cohort/
```

## Citations

- Hubert L, Arabie P. Comparing partitions. J Classif. 1985;2(1):193-218.
