# Performance and Hardware Guide

## Recommended Hardware

| Sample Count | RAM   | CPU Cores | Est. Time (full validation) |
|-------------|-------|-----------|----------------------------|
| 100         | 1 GB  | 1         | < 1 minute                 |
| 1,000       | 2 GB  | 2         | 2–5 minutes                |
| 5,000       | 4 GB  | 2         | 10–20 minutes              |
| 10,000      | 8 GB  | 4         | 15–30 minutes              |
| 25,000      | 16 GB | 4         | 30–60 minutes              |
| 50,000      | 32 GB | 8         | 1–3 hours                  |

Validation time scales primarily with `n_permutations` and `n_bootstrap` (default: 100 and 50 respectively). Reducing these to 50/25 halves validation time with minimal statistical impact.

## Memory Estimation

Use the built-in estimator before running large analyses:

```python
from pathway_subtyping.utils.performance import estimate_memory_usage

est = estimate_memory_usage(
    n_samples=10000,
    n_variants=500000,
    n_genes=5000,
    n_pathways=50,
)
for component, mb in est.items():
    print(f"  {component}: {mb:.0f} MB")
```

This returns per-component estimates (genotype matrix, gene burdens, pathway scores, GMM working set) so you can identify which step will be the bottleneck.

## Progress Bars

Long-running operations display `tqdm` progress bars by default:

```
Label shuffle:    100%|██████████| 100/100 [00:45<00:00]
Random gene sets: 100%|██████████| 100/100 [00:38<00:00]
Bootstrap:        100%|██████████| 50/50   [00:22<00:00]
```

To suppress progress bars (e.g., in batch jobs or CI):

```python
# Validation
vg = ValidationGates(seed=42, show_progress=False)

# Simulation analysis
from pathway_subtyping.simulation import estimate_type_i_error
result = estimate_type_i_error(n_simulations=100, show_progress=False)

# Expression scoring
from pathway_subtyping.expression import score_pathways_from_expression
result = score_pathways_from_expression(expr, pathways, show_progress=False)

# Sensitivity analysis
from pathway_subtyping.sensitivity import vary_feature_subset
result = vary_feature_subset(scores, n_clusters=3, show_progress=False)
```

## Chunked Processing for Large VCFs

For large VCF files that exceed available memory, enable chunked processing:

```python
from pathway_subtyping.pipeline import PipelineConfig, DemoPipeline

config = PipelineConfig(
    name="large_cohort",
    vcf_path="data/large_cohort.vcf",
    use_chunked_processing=True,
    chunk_size=1000,  # samples per chunk
)
pipeline = DemoPipeline(config)
pipeline.run()
```

This processes the VCF in chunks of `chunk_size` samples, computing gene burdens incrementally. The final pathway scores matrix is identical to non-chunked mode.

You can also use the chunked utilities directly:

```python
from pathway_subtyping.utils.performance import (
    chunked_vcf_reader,
    compute_gene_burdens_chunked,
)

# Stream VCF in chunks
for chunk_df in chunked_vcf_reader("large.vcf", chunk_size=1000):
    process(chunk_df)

# Compute gene burdens without loading full VCF
burdens = compute_gene_burdens_chunked("large.vcf", chunk_size=1000)
```

## Google Colab Constraints

The free tier provides ~12 GB RAM and 2 CPU cores. Recommendations:

- **Max samples**: ~5,000 with full validation, ~10,000 with reduced permutations
- **Reduce permutations**: `n_permutations=50, n_bootstrap=25`
- **Use mean_z**: Faster than ssGSEA for expression scoring
- **Avoid pysam**: Use the base install (`pip install pathway-subtyping`), not `[vcf]`

For larger analyses, use Colab Pro (25 GB RAM) or a local machine.

## Performance Tips

1. **Reduce validation iterations**: The defaults (100 permutations, 50 bootstrap) provide strong statistical power. For exploratory work, 50/25 is often sufficient:
   ```python
   vg = ValidationGates(n_permutations=50, n_bootstrap=25, seed=42)
   ```

2. **Use mean_z over ssGSEA**: For expression data, `mean_z` is 10–50x faster than `ssgsea` with comparable results for well-curated pathways:
   ```python
   from pathway_subtyping.expression import ExpressionScoringMethod
   result = score_pathways_from_expression(
       expr, pathways, method=ExpressionScoringMethod.MEAN_Z
   )
   ```

3. **Subset pathways**: Fewer pathways = faster clustering and validation. Focus on biologically relevant pathways rather than running all of Reactome.

4. **Chunked mode for VCFs > 1 GB**: Enable `use_chunked_processing=True` in `PipelineConfig` to avoid loading the full VCF into memory.

5. **Disable progress bars in batch mode**: Set `show_progress=False` to avoid tqdm overhead in non-interactive environments.

## Benchmarking

Run the included benchmark script to test performance on your hardware:

```bash
# Default: 10K samples, 15 pathways, 3 subtypes
python scripts/benchmark_performance.py

# Custom size
python scripts/benchmark_performance.py --n-samples 5000 --n-pathways 20

# Skip validation (benchmark data gen + clustering only)
python scripts/benchmark_performance.py --n-samples 50000 --skip-validation

# Reduced validation for faster benchmarking
python scripts/benchmark_performance.py --n-permutations 25 --n-bootstrap 10
```

The script reports wall-clock time and peak memory per step, and checks against targets (30 minutes, 8 GB for 10K samples).
