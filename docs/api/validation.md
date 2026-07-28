# Validation API Reference

The `validation` module provides validation gates to verify clustering quality and prevent overfitting.

## Overview

The framework implements three mandatory validation tests:

| Test | Purpose | Pass Criteria (default) |
|------|---------|------------------------|
| **Label Shuffle** | Detect spurious clustering | ARI < 0.15 |
| **Random Gene Sets** | Verify pathway biology matters | ARI < 0.15 |
| **Bootstrap Stability** | Ensure robust clusters | ARI ≥ 0.80 |

`run_all()` adds five further gates, each of which runs only when the data it needs is
supplied:

| Gate | Runs when | Pass Criteria (default) |
|------|-----------|------------------------|
| **Ancestry Independence** | `ancestry_pcs` given | min ancestry p > 0.05 / n_PCs (Bonferroni) |
| **Cross-Modal Concordance** | `per_modality_scores` has ≥ 2 modalities | mean concordance ARI > permutation null 95th pct |
| **Confound Association** | `confounds` given | no nuisance confound both significant and Cramér's V ≥ 0.30 |
| **Genetic Anchoring** | `genetic_anchoring` given | ≥ 1 subtype significantly enriched, fold > 1.0 |
| **Genetic Anchoring (somatic)** | `somatic_anchoring` given | ≥ 1 stratum significant with Cramér's V ≥ 0.30 |

The first four gates above are *negative* controls and confound checks — they can only
show a partition is reproducible and not confounded. The two anchoring gates have
inverted polarity: they are **positive-evidence**, low-power tests, so a failure is weak
evidence of absence rather than proof a partition lacks genetic grounding.

> **Threshold Calibration:** These default thresholds (0.15/0.8) can be automatically calibrated based on sample size and cluster count using the [`threshold_calibration`](threshold_calibration.md) module. Set `validation.stability_threshold: null` and `validation.null_ari_max: null` in your config to enable auto-calibration.

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
#     "summary": "All 3 validation gates PASSED",  # count = gates actually run
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
    null_ari_max: float = 0.15,
    show_progress: bool = True
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
| `show_progress` | `bool` | `True` | Show tqdm progress bars for long-running loops |

#### Methods

##### `run_all(...) -> ValidationGatesResult`

Run all validation gates. The three mandatory tests always run; each optional argument
below switches on one additional gate.

```python
result = validator.run_all(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathways: Dict[str, List[str]],
    gene_burdens: pd.DataFrame,
    n_clusters: int,
    gmm_seed: Optional[int] = None,
    ancestry_pcs=None,
    per_modality_scores: Optional[Dict[str, pd.DataFrame]] = None,
    confounds: Optional[Dict[str, Any]] = None,
    genetic_anchoring: Optional[Dict[str, Any]] = None,
    somatic_anchoring: Optional[Dict[str, Any]] = None
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | `pd.DataFrame` | required | Pathway scores matrix (samples × pathways) |
| `cluster_labels` | `np.ndarray` | required | Cluster assignments (integer labels) |
| `pathways` | `Dict[str, List[str]]` | required | Pathway name → gene list mapping |
| `gene_burdens` | `pd.DataFrame` | required | Gene burden matrix (samples × genes) |
| `n_clusters` | `int` | required | Number of clusters |
| `gmm_seed` | `Optional[int]` | `None` | Seed for GMM in validation tests |
| `ancestry_pcs` | `AncestryPCs` (unannotated in source) | `None` | Ancestry PCs; enables the ancestry independence gate |
| `per_modality_scores` | `Optional[Dict[str, pd.DataFrame]]` | `None` | Modality label → pathway score DataFrame; enables the cross-modal concordance gate when ≥ 2 modalities are given |
| `confounds` | `Optional[Dict[str, Any]]` | `None` | Confound name (e.g. `region`, `batch`, `diagnosis`) → per-sample categorical values; enables the confound association gate. Diagnosis-like keys are treated as biology-of-interest and never fail the gate |
| `genetic_anchoring` | `Optional[Dict[str, Any]]` | `None` | Kwargs forwarded to `genetic_anchoring_gate()` (requires at least `subtype_gene_sets`, `risk_genes`, `background_universe`) |
| `somatic_anchoring` | `Optional[Dict[str, Any]]` | `None` | Kwargs forwarded to `somatic_anchoring_gate()` (requires `somatic_strata`); `cluster_labels` is supplied by `run_all()` |

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

##### `negative_control_ancestry_independence(...) -> ValidationResult`

Test whether clusters are confounded with ancestry (population structure rather than
biology). Runs from `run_all()` only when `ancestry_pcs` is supplied.

```python
result = validator.negative_control_ancestry_independence(
    cluster_labels: np.ndarray,
    ancestry_pcs
)
```

**Logic:**
1. Delegate to `ancestry.check_ancestry_independence()` with `significance_threshold=0.05`
2. Take the minimum p-value across the ancestry PCs as the test statistic
3. Compare against the Bonferroni threshold `0.05 / n_PCs`
4. **PASS** if no ancestry PC is significantly associated with the clusters
   (`report.overall_independence_passed`)

Reported as `min_ancestry_pvalue` (comparison `>`); `details` carries `n_pcs_tested` and
the per-PC p-values.

---

##### `cross_modal_concordance_gate(...) -> ValidationResult`

**Gate 5.** Test whether subtypes discovered on fused data are consistent across the
individual data modalities. Runs from `run_all()` only when `per_modality_scores`
contains at least two modalities.

```python
result = validator.cross_modal_concordance_gate(
    per_modality_scores: Dict[str, pd.DataFrame],
    cluster_labels: np.ndarray,
    fused_sample_ids: List[str],
    n_clusters: int
)
```

**Logic:**
1. Delegate to `cross_modal_validation.cross_modal_concordance()`, passing the
   validator's `n_permutations`, `seed`, and `show_progress`
2. Cluster each modality independently and measure pairwise agreement with the fused
   labels
3. Compare the mean concordance ARI against the permutation null's 95th percentile
4. **PASS** if mean concordance ARI exceeds the null 95th percentile
   (`result.gate_passed`)

Reported as `mean_concordance_ARI` (comparison `>`); `details` carries
`mean_concordance_nmi`, `mean_transfer_ari`, `null_ari_95th`, `n_modality_pairs`, and
the per-pair results.

---

##### `confound_association_gate(...) -> ValidationResult`

**Gate 6.** Test whether the partition is explained by a technical or anatomical
confound (brain region, sequencing batch, …) rather than the biology of interest.
Mandatory whenever confound annotations are available; runs from `run_all()` when
`confounds` is non-empty.

```python
result = validator.confound_association_gate(
    cluster_labels: np.ndarray,
    confounds: Dict[str, Any],
    cramers_v_max: float = 0.30,
    alpha: float = 0.05
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_labels` | `np.ndarray` | required | Cluster assignments (n_samples,) |
| `confounds` | `Dict[str, Any]` | required | Confound name → per-sample categorical values (length n_samples). Confounds of the wrong length are skipped with a warning |
| `cramers_v_max` | `float` | `0.30` | Effect-size threshold above which an association counts as a confound (0.30 ≈ medium; 0.10 small, 0.50 large) |
| `alpha` | `float` | `0.05` | Significance level for the BH-adjusted chi-square p-values |

**Logic:**
1. For each confound, cross-tabulate cluster label × confound level
2. Chi-square test of independence (no continuity correction) plus Cramér's V effect size
3. Benjamini–Hochberg-adjust the p-values across the confounds actually tested
4. Evaluate over *nuisance* confounds only — keys in
   `{diagnosis, dx, disease, condition, phenotype}` are treated as biology-of-interest,
   reported for context but never able to fail the gate (a partition *should* track
   diagnosis)
5. **PASS** if no nuisance confound is both significant (adjusted p < `alpha`) **and**
   non-trivially associated (Cramér's V ≥ `cramers_v_max`)

Reported as `max_nuisance_cramers_v` (comparison `<`); `details` carries
`worst_confound`, `failing_confounds`, and the full `per_confound` breakdown.

**Why this gate exists.** Its absence is what let an anatomy artifact pass the rest of
the battery. On GSE80655 a k=3 partition passed every stability control (bootstrap ARI
~0.92) yet corresponded almost exactly to brain region (Cramér's V ~0.67, p ~4e-26)
while being independent of diagnosis (p ~0.41). A stability-passing partition can still
be a confound classifier; this gate is the check that catches it.

---

##### `genetic_anchoring_gate(...) -> ValidationResult`

**Gate 7 (feature-level).** Test whether the genes defining each subtype are
over-represented for disease genetic-risk genes. Runs from `run_all()` when
`genetic_anchoring` kwargs are supplied.

```python
result = validator.genetic_anchoring_gate(
    subtype_gene_sets: Dict[Any, Any],
    risk_genes: Any,
    background_universe: Any,
    reference_universe: Optional[Any] = None,
    alpha: float = 0.05,
    min_fold: float = 1.0
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `subtype_gene_sets` | `Dict[Any, Any]` | required | Subtype label → iterable of gene IDs defining that subtype. Identifiers must share a namespace with `risk_genes` and the universes |
| `risk_genes` | `Any` | required | Disease genetic-risk reference gene IDs |
| `background_universe` | `Any` | required | Background-matched universe of testable genes (the discriminating null, e.g. brain-expressed genes) |
| `reference_universe` | `Optional[Any]` | `None` | Optional genome-wide universe reported per subtype for contrast; never decides the gate |
| `alpha` | `float` | `0.05` | Significance level for the BH-adjusted hypergeometric p-values |
| `min_fold` | `float` | `1.0` | Minimum fold enrichment for a subtype to count as anchored |

**Logic:**
1. Hypergeometric over-representation test per subtype against `background_universe`
2. Benjamini–Hochberg-adjust across subtypes (a non-finite p-value is treated as 1.0)
3. A subtype is *anchored* if adjusted p < `alpha`, fold > `min_fold`, and it has at
   least one risk-gene hit
4. **PASS** if at least one subtype is anchored

Reported as `max_enrichment_fold` (comparison `>`); `details` carries
`anchored_subtypes`, `best_subtype`, universe sizes, the `per_subtype` breakdown, and
`gate_polarity: "positive_evidence"`.

**Polarity and caveats.** Unlike the negative controls, this is a positive-evidence
gate: germline variants are fixed at conception, upstream of any postmortem or technical
confound (PMI, RIN, dissection, batch, platform), so an enrichment cannot be
manufactured by one. The null must be background-matched, not genome-wide — brain-
expressed genes are already enriched for brain-disease risk, so a region-identity
subtype could show spurious genome-wide enrichment. It is a specific, low-power test: a
non-anchored result is weak evidence of absence, and a pass does not upgrade an unstable
or confounded partition.

---

##### `somatic_anchoring_gate(...) -> ValidationResult`

**Gate 7 (somatic mode).** The cancer counterpart to `genetic_anchoring_gate()`: test
whether the *tumors* in a subtype carry a somatic stratum — driver mutation, copy-number
alteration, or mutational-signature class — more than the other subtypes. Runs from
`run_all()` when `somatic_anchoring` kwargs are supplied.

```python
result = validator.somatic_anchoring_gate(
    cluster_labels: np.ndarray,
    somatic_strata: Dict[str, Any],
    cramers_v_min: float = 0.30,
    alpha: float = 0.05
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_labels` | `np.ndarray` | required | Cluster assignments (n_samples,) |
| `somatic_strata` | `Dict[str, Any]` | required | Stratum name → per-sample categorical value (e.g. driver carrier status, MSI-high/low). Missing values are dropped per stratum |
| `cramers_v_min` | `float` | `0.30` | Minimum Cramér's V for a stratum to count as an anchor |
| `alpha` | `float` | `0.05` | Significance level for the BH-adjusted p-values |

**Logic:**
1. Delegate to `genetics.somatic_anchoring.somatic_alignment()`
2. Chi-square + Bergsma-corrected Cramér's V per stratum, BH-adjusted — the confound
   gate's statistic with inverted polarity (a nuisance association fails Gate 6, a
   somatic-driver association passes Gate 7)
3. **PASS** if at least one stratum is both significant (adjusted p < `alpha`) and
   non-trivially associated (Cramér's V ≥ `cramers_v_min`)

Reported as `max_somatic_cramers_v` (comparison `>`); `details` carries
`anchored_strata`, `best_stratum`, the `per_stratum` breakdown, and
`gate_polarity: "positive_evidence"`.

**Confound caveat.** A somatic driver can co-vary with tissue of origin (pan-cancer) and
tumor purity (within a tumor type). Run the confound gate first so a "somatic anchor" is
not the tissue/purity axis re-labelled. As with the feature-level gate, a null is weak
evidence of absence.

---

## Utility Functions

### `cramers_v(contingency: np.ndarray, bias_correction: bool = True) -> float`

Cramér's V effect size for a contingency table — association strength between two
categorical variables (0 = independent, 1 = perfect association). This is the statistic
behind the confound association and somatic anchoring gates.

```python
import numpy as np
import pandas as pd
from pathway_subtyping.validation import cramers_v

table = pd.crosstab(cluster_labels, brain_region).values
v = cramers_v(table)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `contingency` | `np.ndarray` | required | r × c contingency table of counts |
| `bias_correction` | `bool` | `True` | Apply the Bergsma (2013) small-sample bias correction |

**Returns:** Cramér's V in [0, 1]; `0.0` for degenerate tables (a single row/column, or
no samples).

---

### `format_validation_report(result: ValidationGatesResult) -> str`

Format validation results as Markdown.

```python
from pathway_subtyping.validation import format_validation_report

markdown = format_validation_report(validation_result)
print(markdown)
```

**Output** (a run with no optional gate data, so only the three mandatory tests ran —
the count in the summary line is the number of gates actually run, not a fixed 3):
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

The report also emits an **Interpretation** section explaining the two negative controls,
the stability test, and cross-modal concordance.

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
