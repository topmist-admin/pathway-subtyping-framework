# Variant QC API Reference

The `variant_qc` module provides standard genetic variant quality control filters applied before burden computation.

## Classes

### `VariantQCConfig`

Configuration for variant quality control filters.

```python
from pathway_subtyping.variant_qc import VariantQCConfig

config = VariantQCConfig(
    min_qual=30.0,
    min_call_rate=0.95,
    hwe_p_threshold=1e-6,
    max_maf=0.01,
    min_gq=20,
    min_dp=10,
)
```

#### Attributes

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_qual` | `float` | `30.0` | Minimum QUAL score (0 to disable) |
| `min_call_rate` | `float` | `0.9` | Minimum genotype call rate per variant (0.0–1.0) |
| `hwe_p_threshold` | `float` | `1e-6` | HWE p-value threshold; variants below are removed (0 to disable) |
| `max_maf` | `float` | `0.01` | Maximum minor allele frequency (1.0 to disable) |
| `min_gq` | `Optional[int]` | `None` | Minimum genotype quality; genotypes below are set to missing |
| `min_dp` | `Optional[int]` | `None` | Minimum read depth; genotypes below are set to missing |

#### Methods

- **`to_dict() -> Dict[str, Any]`** — Serialize to dictionary

---

### `VariantQCResult`

Result from variant quality control filtering.

```python
result = filter_variants(variants_df, genotypes_df, config)
print(result.format_report())
```

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_variants` | `int` | Variants before QC |
| `passed_variants` | `int` | Variants after QC |
| `removed_variants` | `int` | Number of variants removed |
| `removal_reasons` | `Dict[str, int]` | Count per filter (e.g., `{"low_qual": 5, "high_maf": 3}`) |
| `per_variant_metrics` | `Optional[pd.DataFrame]` | QC metrics per variant (call_rate, maf, hwe_p, passed) |
| `config` | `Optional[VariantQCConfig]` | The QC config used |

#### Methods

- **`to_dict() -> Dict[str, Any]`** — Serialize to dictionary (includes retention rate)
- **`format_report() -> str`** — Human-readable QC report
- **`get_citations() -> List[str]`** — Relevant citations (Wigginton 2005, Anderson 2010)

---

## Functions

### `filter_variants(variants_df, genotypes_df, config=None)`

Apply variant-level QC filters. This is the main entry point.

```python
from pathway_subtyping.variant_qc import VariantQCConfig, filter_variants

config = VariantQCConfig(min_qual=30, max_maf=0.01)
filtered_variants, filtered_genotypes, result = filter_variants(
    variants_df, genotypes_df, config
)

print(f"Retained {result.passed_variants}/{result.total_variants} variants")
print(result.format_report())
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `variants_df` | `pd.DataFrame` | Variant metadata (must have `qual` column for QUAL filtering) |
| `genotypes_df` | `pd.DataFrame` | Allele counts (variants x samples), values 0/1/2 |
| `config` | `Optional[VariantQCConfig]` | QC configuration (defaults used if None) |

**Returns:** `Tuple[pd.DataFrame, pd.DataFrame, VariantQCResult]`
- Filtered variants DataFrame
- Filtered genotypes DataFrame
- QC result with metrics

**Filters applied (in order):**
1. **QUAL** — Remove variants with QUAL < `min_qual`
2. **Call rate** — Remove variants with call rate < `min_call_rate`
3. **HWE** — Remove variants with HWE p-value < `hwe_p_threshold`
4. **MAF** — Remove variants with MAF > `max_maf`

---

### `compute_call_rate(genotypes_df)`

Compute per-variant genotype call rate (proportion of non-missing genotypes).

```python
from pathway_subtyping.variant_qc import compute_call_rate

rates = compute_call_rate(genotypes_df)
print(rates.describe())
```

**Parameters:**
- `genotypes_df` (`pd.DataFrame`): Allele counts (variants x samples)

**Returns:** `pd.Series` of call rates indexed by variant

---

### `compute_maf(genotypes_df)`

Compute minor allele frequency per variant (diploid assumption).

```python
from pathway_subtyping.variant_qc import compute_maf

maf = compute_maf(genotypes_df)
rare = maf[maf < 0.01]
print(f"{len(rare)} rare variants (MAF < 1%)")
```

**Parameters:**
- `genotypes_df` (`pd.DataFrame`): Allele counts (0, 1, 2) per variant x sample

**Returns:** `pd.Series` of MAF values (range 0.0–0.5)

---

### `check_hwe(genotypes_df)`

Test Hardy-Weinberg equilibrium per variant using chi-squared goodness of fit (1 df).

```python
from pathway_subtyping.variant_qc import check_hwe

pvals = check_hwe(genotypes_df)
violations = pvals[pvals < 1e-6]
print(f"{len(violations)} HWE violations at p < 1e-6")
```

**Parameters:**
- `genotypes_df` (`pd.DataFrame`): Allele counts (0, 1, 2) per variant x sample

**Returns:** `pd.Series` of HWE p-values. Returns 1.0 for monomorphic variants or insufficient data.

---

### `apply_genotype_filters(genotypes_df, genotype_fields=None, min_gq=None, min_dp=None)`

Apply per-genotype quality filters. Sets genotypes below thresholds to NaN.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `genotypes_df` | `pd.DataFrame` | required | Allele counts (variants x samples) |
| `genotype_fields` | `Optional[pd.DataFrame]` | `None` | DataFrame with `GQ` and/or `DP` columns |
| `min_gq` | `Optional[int]` | `None` | Minimum genotype quality |
| `min_dp` | `Optional[int]` | `None` | Minimum read depth |

**Returns:** `Tuple[pd.DataFrame, int]` — (filtered genotypes, count of masked genotypes)

---

## Pipeline Integration

Enable variant QC in your YAML config:

```yaml
variant_qc:
  enabled: true
  min_qual: 30
  min_call_rate: 0.95
  hwe_p_threshold: 1e-6
  max_maf: 0.01
```

When enabled, variant QC runs between `load_data()` and `compute_gene_burdens()` in the pipeline. Results are included in both JSON and Markdown reports.

---

## References

- Wigginton JE, Cutler DJ, Gravel A. A note on exact tests of Hardy-Weinberg equilibrium. *Am J Hum Genet*. 2005;76(5):887-93.
- Anderson CA et al. Data quality control in genetic case-control association studies. *Nat Protoc*. 2010;5(9):1564-73.
