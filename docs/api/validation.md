# Validation API Reference

The `validation` module provides validation gates to verify clustering quality and prevent overfitting.

## Overview

The framework implements three mandatory validation tests:

| Test | Purpose | Pass Criteria |
|------|---------|---------------|
| **Label Shuffle** | Detect spurious clustering | ARI < 0.15 |
| **Random Gene Sets** | Verify pathway biology matters | ARI < 0.15 |
| **Bootstrap Stability** | Ensure robust clusters | ARI ≥ 0.80 |

## Classes

### ValidationResult

Result from a single validation test.

```python
from pathway_subtyping.validation import ValidationResult
```

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Test name |
| `passed` | `bool` | Whether test passed |
| `metric_name` | `str` | Name of metric (e.g., "mean_null_ARI") |
| `metric_value` | `float` | Computed metric value |
| `threshold` | `float` | Pass/fail threshold |
| `comparison` | `str` | Comparison operator ("<", ">", ">=") |
| `details` | `Dict[str, Any]` | Additional details |

#### Properties

##### `status -> str`

Returns `"PASS"` or `"FAIL"`.

```python
print(result.status)  # "PASS" or "FAIL"
```

#### Methods

##### `to_dict() -> Dict[str, Any]`

Convert to dictionary for JSON serialization.

```python
result_dict = validation_result.to_dict()
# {
#     "name": "Negative Control 1: Label Shuffle",
#     "status": "PASS",
#     "metric": "mean_null_ARI",
#     "value": 0.0234,
#     "threshold": 0.15,
#     "comparison": "<",
#     "details": {...}
# }
```

---

### ValidationGatesResult

Aggregated results from all validation gates.

```python
from pathway_subtyping.validation import ValidationGatesResult
```

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `results` | `List[ValidationResult]` | Individual test results |
| `all_passed` | `bool` | Whether all tests passed |
| `summary` | `str` | Human-readable summary |

#### Methods

##### `to_dict() -> Dict[str, Any]`

Convert to dictionary for JSON serialization.

```python
gates_dict = validation_gates_result.to_dict()
# {
#     "all_passed": true,
#     "summary": "All 3 validation gates PASSED",
#     "tests": [...]
# }
```

---

### ValidationGates

Main validation orchestrator class.

```python
from pathway_subtyping.validation import ValidationGates
```

#### Constructor

```python
validator = ValidationGates(
    seed: Optional[int] = 42,
    n_permutations: int = 100,
    n_bootstrap: int = 50,
    stability_threshold: float = 0.8,
    null_ari_max: float = 0.15
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seed` | `Optional[int]` | `42` | Random seed for reproducibility |
| `n_permutations` | `int` | `100` | Permutations for null tests |
| `n_bootstrap` | `int` | `50` | Bootstrap iterations for stability |
| `stability_threshold` | `float` | `0.8` | Minimum ARI for stability test |
| `null_ari_max` | `float` | `0.15` | Maximum ARI under null hypothesis |

#### Methods

##### `run_all(...) -> ValidationGatesResult`

Run all validation gates.

```python
result = validator.run_all(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathways: Dict[str, List[str]],
    gene_burdens: pd.DataFrame,
    n_clusters: int,
    gmm_seed: Optional[int] = None
)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pathway_scores` | `pd.DataFrame` | Pathway scores matrix (samples × pathways) |
| `cluster_labels` | `np.ndarray` | Cluster assignments (integer labels) |
| `pathways` | `Dict[str, List[str]]` | Pathway name → gene list mapping |
| `gene_burdens` | `pd.DataFrame` | Gene burden matrix (samples × genes) |
| `n_clusters` | `int` | Number of clusters |
| `gmm_seed` | `Optional[int]` | Seed for GMM in validation tests |

**Returns:** `ValidationGatesResult`

---

##### `negative_control_label_shuffle(...) -> ValidationResult`

Test if clustering can recover randomly shuffled labels.

```python
result = validator.negative_control_label_shuffle(
    pathway_scores: pd.DataFrame,
    original_labels: np.ndarray,
    n_clusters: int,
    gmm_seed: Optional[int] = None
)
```

**Logic:**
1. Randomly shuffle cluster labels
2. Re-run GMM clustering
3. Compute ARI between new clusters and shuffled labels
4. Repeat `n_permutations` times
5. **PASS** if mean ARI < `null_ari_max` (clustering doesn't find patterns in noise)

---

##### `negative_control_random_gene_sets(...) -> ValidationResult`

Test if clusters are driven by biological pathways vs. random genes.

```python
result = validator.negative_control_random_gene_sets(
    gene_burdens: pd.DataFrame,
    real_pathways: Dict[str, List[str]],
    original_labels: np.ndarray,
    n_clusters: int,
    gmm_seed: Optional[int] = None
)
```

**Logic:**
1. Replace each pathway with random genes (same size)
2. Compute pathway scores using random gene sets
3. Cluster on random pathway scores
4. Compute ARI with original clustering
5. Repeat `n_permutations` times
6. **PASS** if mean ARI < `null_ari_max` (random genes don't replicate structure)

---

##### `stability_test_bootstrap(...) -> ValidationResult`

Test if clusters are robust to resampling.

```python
result = validator.stability_test_bootstrap(
    pathway_scores: pd.DataFrame,
    original_labels: np.ndarray,
    n_clusters: int,
    gmm_seed: Optional[int] = None
)
```

**Logic:**
1. Bootstrap sample (with replacement)
2. Re-run GMM clustering on bootstrap sample
3. Compute ARI with original labels
4. Repeat `n_bootstrap` times
5. **PASS** if mean ARI ≥ `stability_threshold` (clusters are stable)

---

## Utility Functions

### `format_validation_report(result: ValidationGatesResult) -> str`

Format validation results as Markdown.

```python
from pathway_subtyping.validation import format_validation_report

markdown = format_validation_report(validation_result)
print(markdown)
```

**Output:**
```markdown
## Validation Gates

**Overall Status:** PASS

All 3 validation gates PASSED

| Test | Status | Metric | Value | Threshold |
|------|--------|--------|-------|-----------|
| Negative Control 1: Label Shuffle | ✓ PASS | mean_null_ARI | 0.023 | < 0.15 |
| Negative Control 2: Random Gene Sets | ✓ PASS | mean_random_ARI | 0.048 | < 0.15 |
| Stability Test: Bootstrap | ✓ PASS | mean_bootstrap_ARI | 0.921 | >= 0.8 |
```

---

## Example: Standalone Validation

```python
import pandas as pd
import numpy as np
from pathway_subtyping.validation import ValidationGates

# Load your data
pathway_scores = pd.read_csv("pathway_scores.csv", index_col=0)
gene_burdens = pd.read_csv("gene_burdens.csv", index_col=0)
assignments = pd.read_csv("assignments.csv")

# Load pathways
pathways = {}
with open("pathways.gmt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            pathways[parts[0]] = parts[2:]

# Run validation
validator = ValidationGates(
    seed=42,
    n_permutations=200,  # More permutations for stricter testing
    n_bootstrap=100,
    stability_threshold=0.85,  # Stricter stability requirement
)

result = validator.run_all(
    pathway_scores=pathway_scores,
    cluster_labels=assignments["cluster_id"].values,
    pathways=pathways,
    gene_burdens=gene_burdens,
    n_clusters=assignments["cluster_id"].nunique(),
    gmm_seed=42,
)

# Check results
if result.all_passed:
    print("Clustering validated!")
else:
    print("Validation FAILED")
    for test in result.results:
        if not test.passed:
            print(f"  - {test.name}: {test.metric_value:.3f} (threshold: {test.comparison} {test.threshold})")
```

---

## Example: Custom Thresholds

```python
from pathway_subtyping.validation import ValidationGates

# Stricter validation for publication
validator = ValidationGates(
    seed=42,
    n_permutations=500,      # More permutations
    n_bootstrap=200,          # More bootstrap samples
    stability_threshold=0.9,  # Require 90% stability
    null_ari_max=0.10,        # Stricter null threshold
)

# More lenient for exploratory analysis
validator_exploratory = ValidationGates(
    seed=42,
    n_permutations=50,
    n_bootstrap=30,
    stability_threshold=0.7,
    null_ari_max=0.20,
)
```

---

## Technical Implementation Details

### GMM Reliability Improvements

The validation module includes several reliability improvements for GMM clustering:

**Covariance Regularization:**
All GMM fits use `reg_covar=1e-6` to prevent numerical instability:
```python
gmm = GaussianMixture(
    n_components=n_clusters,
    covariance_type="full",
    n_init=5,
    random_state=seed,
    reg_covar=1e-6,  # Regularization for numerical stability
)
```

**Convergence Checking:**
All validation tests check GMM convergence and handle non-converged fits gracefully:
```python
gmm.fit(data)
if not gmm.converged_:
    logger.warning("GMM did not converge")
    continue  # Skip this iteration
```

**Empty Results Handling:**
When no GMM fits converge during permutation testing, the validation test returns a failing result with clear diagnostic information:
```python
if not ari_values:
    return ValidationResult(
        name="Test Name",
        passed=False,
        metric_value=1.0,  # Worst case
        details={"error": "No GMM fits converged", "n_attempted": n_permutations}
    )
```

---

## Understanding Validation Failures

### Label Shuffle Fails (ARI too high)

**Meaning:** Clustering finds structure even in random labels — suggests overfitting or method artifacts.

**Potential fixes:**
- Reduce number of pathways (too many features)
- Increase sample size
- Use different clustering method
- Check for data quality issues

### Random Gene Sets Fails (ARI too high)

**Meaning:** Random gene groupings work as well as biological pathways — suggests clusters aren't driven by biology.

**Potential fixes:**
- Review pathway definitions (may be too broad)
- Check variant annotation quality
- Ensure pathways are disease-relevant
- Consider gene set enrichment analysis

### Bootstrap Fails (ARI too low)

**Meaning:** Clusters aren't robust to resampling — suggests unstable or weak signal.

**Potential fixes:**
- Increase sample size
- Use fewer clusters (simpler model)
- Check for outlier samples
- Review variant filtering criteria
