# Batch Correction API

> **Module**: `pathway_subtyping.batch_correction`

Detects and corrects technical batch effects in pathway scores using empirical Bayes (ComBat), mean-centering, or standardization. Includes post-correction validation to verify batch variance was reduced without destroying biological signal.

---

## Quick Example

```python
from pathway_subtyping import (
    detect_batch_effects,
    correct_batch_effects,
    validate_batch_correction,
    BatchCorrectionMethod,
)

# Step 1: Detect batch effects
report = detect_batch_effects(pathway_scores, batch_labels, batch_variable="site")
print(report.format_report())

if report.overall_batch_effect:
    # Step 2: Correct
    result = correct_batch_effects(
        pathway_scores, batch_labels,
        method=BatchCorrectionMethod.COMBAT,
        batch_variable="site",
    )
    print(result.format_report())

    # Step 3: Validate
    val = validate_batch_correction(
        pathway_scores, result.corrected_scores,
        batch_labels, biological_labels=subtype_labels,
    )
    print(f"Batch variance reduced: {val['batch_variance_reduced']}")
    print(f"Biological signal preserved: {val.get('biological_signal_preserved', 'N/A')}")
```

---

## Enums

### `BatchCorrectionMethod`

| Value | Algorithm | Speed | Best For |
|-------|-----------|-------|----------|
| `COMBAT` | Empirical Bayes (Johnson et al., 2007) — shrinks batch-specific location and scale parameters toward a common prior | Medium | **Recommended default** — best preservation of biological variance |
| `MEAN_CENTER` | Subtract batch-specific mean per feature | Fast | Quick correction when variance differences are minimal |
| `STANDARDIZE` | Per-batch Z-score, then re-scale to grand distribution | Fast | When batches differ in both mean and variance |

---

## Functions

### `detect_batch_effects()`

```python
def detect_batch_effects(
    pathway_scores: pd.DataFrame,
    batch_labels: np.ndarray,
    batch_variable: str = "batch",
    significance_threshold: float = 0.05,
) -> BatchEffectReport
```

Detect batch effects using one-way ANOVA per feature with Benjamini-Hochberg FDR correction. Reports F-statistics, p-values, and eta-squared (variance explained by batch).

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `batch_labels` | np.ndarray | required | Batch assignment per sample |
| `batch_variable` | str | "batch" | Name for reporting (e.g., "sequencing_site") |
| `significance_threshold` | float | 0.05 | FDR threshold for significance |

**Returns:** `BatchEffectReport` with per-feature statistics and overall detection result.

**Raises:** `ValueError` if fewer than 2 batches.

---

### `correct_batch_effects()`

```python
def correct_batch_effects(
    pathway_scores: pd.DataFrame,
    batch_labels: np.ndarray,
    method: BatchCorrectionMethod = BatchCorrectionMethod.COMBAT,
    batch_variable: str = "batch",
    preserve_biological: Optional[np.ndarray] = None,
) -> BatchCorrectionResult
```

Correct batch effects in pathway scores. ComBat uses empirical Bayes shrinkage to adjust batch-specific mean and variance while preserving biological variation.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `batch_labels` | np.ndarray | required | Batch assignment per sample |
| `method` | BatchCorrectionMethod | COMBAT | Correction algorithm |
| `batch_variable` | str | "batch" | Name for reporting |
| `preserve_biological` | np.ndarray or None | None | Known biological group labels to preserve (ComBat) |

**Returns:** `BatchCorrectionResult` with corrected scores, pre/post variance, and reduction metrics.

---

### `validate_batch_correction()`

```python
def validate_batch_correction(
    original_scores: pd.DataFrame,
    corrected_scores: pd.DataFrame,
    batch_labels: np.ndarray,
    biological_labels: Optional[np.ndarray] = None,
) -> Dict[str, Any]
```

Validate that correction removed batch effects without destroying biological signal.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `original_scores` | pd.DataFrame | required | Pre-correction pathway scores |
| `corrected_scores` | pd.DataFrame | required | Post-correction pathway scores |
| `batch_labels` | np.ndarray | required | Batch assignments |
| `biological_labels` | np.ndarray or None | None | Known biological groups (for signal preservation check) |

**Returns:** Dictionary with:

| Key | Type | Description |
|-----|------|-------------|
| `batch_variance_before` | float | Mean eta-squared before correction |
| `batch_variance_after` | float | Mean eta-squared after correction |
| `batch_variance_reduced` | bool | Whether batch variance decreased |
| `biological_variance_before` | float | Biological eta-squared before (if labels provided) |
| `biological_variance_after` | float | Biological eta-squared after (if labels provided) |
| `biological_signal_preserved` | bool | Whether >=50% of biological variance retained |
| `mean_correlation_with_original` | float | Mean Pearson r between original and corrected per feature |

---

## Dataclasses

### `BatchEffectReport`

| Attribute | Type | Description |
|-----------|------|-------------|
| `batch_variable` | str | Name of the batch variable |
| `n_batches` | int | Number of unique batches |
| `batch_sizes` | Dict[str, int] | Samples per batch |
| `f_statistics` | Dict[str, float] | Per-feature F-statistic (one-way ANOVA) |
| `p_values` | Dict[str, float] | Per-feature p-value |
| `significant_features` | List[str] | Features with FDR < threshold |
| `variance_explained` | Dict[str, float] | Per-feature eta-squared |
| `overall_batch_effect` | bool | Whether any significant batch effects detected |
| `significance_threshold` | float | FDR threshold used |

**Methods:** `to_dict()`, `format_report()`

### `BatchCorrectionResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `corrected_scores` | pd.DataFrame | Batch-corrected pathway scores |
| `method` | BatchCorrectionMethod | Correction method used |
| `batch_variable` | str | Name of the batch variable |
| `n_batches` | int | Number of batches |
| `pre_correction_variance` | Dict[str, float] | Per-feature batch eta-squared before |
| `post_correction_variance` | Dict[str, float] | Per-feature batch eta-squared after |
| `variance_reduction` | Dict[str, float] | Per-feature variance reduction |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

---

## Method Comparison

| Criterion | ComBat | Mean-Center | Standardize |
|-----------|--------|-------------|-------------|
| Corrects mean shift | Yes | Yes | Yes |
| Corrects variance differences | Yes (EB shrinkage) | No | Yes (Z-score) |
| Preserves biological signal | Best | Good | Moderate |
| Minimum batch size | 3 samples | 2 samples | 3 samples |
| Can use biological covariates | Yes (`preserve_biological`) | No | No |

---

## Interpretation

| Metric | Good | Concerning |
|--------|------|------------|
| `batch_variance_before` (eta-squared) | < 0.05 | > 0.10 |
| `batch_variance_after` | < 0.02 | > 0.05 |
| `variance_reduction` | > 50% | < 20% |
| `biological_signal_preserved` | True | False |
| `mean_correlation_with_original` | > 0.8 | < 0.5 |

Always run `detect_batch_effects()` first. If `overall_batch_effect` is `False`, correction is unnecessary and could introduce noise.

---

## References

- Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*. 2007;8(1):118-127.
- Leek JT, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nat Rev Genet*. 2010;11(10):733-739.

---

## See Also

- [Ancestry Correction](ancestry.md) — population stratification correction
- [Validation Gates](validation.md) — cluster quality validation
- [Sensitivity](sensitivity.md) — robustness testing
