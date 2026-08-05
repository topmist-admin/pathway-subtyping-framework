# Discreteness API Reference

> **Subpackage**: `pathway_subtyping.discreteness`
> **Added in**: v0.8.0

Discreteness-aware subtype validation. `DiscretenessGateA` keeps the stability gate's
observed statistic (mean bootstrap-ARI at a fixed *k*) and replaces the *reference*
with a single-Gaussian (SigClust) null, so the gate tests **"are there discrete
clusters?"** rather than **"are the pathways mutually independent?"**.
`ReframedMembershipGate` reframes the conformal membership gate as assignment
sharpness *conditional on* a partition Gates A and B have already certified.

For the rationale — why the previous independence null was mis-specified, and what
that implies for prior "reproducible subtype" claims — see
[Discreteness-aware Gate A](../discreteness_gate.md). This page is the API reference.

---

## Decision Rule — read before citing this gate

**Only the single-Gaussian (SigClust) reference decides.** The rule is literally

```python
passed = bool(obs > sg_p95 and sg_p < self.alpha)
```

(`src/pathway_subtyping/discreteness/gate_a_discreteness_null.py:400`).

The **gap statistic and the Hartigan dip test are computed and reported but never
enter the decision.** Do not describe this gate as requiring a dataset to clear three
criteria — it requires one, and reports two more alongside it as diagnostics.

**Three outcomes, not two.** The gate returns certify, reject, or `not-testable` —
the last being an *abstention* triggered by the k-stability rule, not a rejection. On
the synthetic negative controls it abstained on **26 of 30 (87%)**, so any
false-positive rate must be quoted against the **testable denominator (4)**, never
against 30.

**It is conservative.** On the separation sweep it certifies nothing below δ = 2.5 SD:
good for refusing over-called subtypes, poor at detecting subtle real ones.

---

## Quick Example

```python
from pathway_subtyping.discreteness import DiscretenessGateA, ReframedMembershipGate

gate_a = DiscretenessGateA(seed=42, n_ref=200, n_bootstrap=50)
res = gate_a.run("my_cohort", pathway_scores_df, n_clusters=k)
print(res.verdict)   # "discrete structure" | "no discrete structure (continuum)"
                     # | "not-testable (no reproducible k)"
print(res.passed, res.sg_empirical_p, res.gap_optimal_k)

gate_c = ReframedMembershipGate(seed=42)
c = gate_c.run("my_cohort", pathway_scores_df, cluster_labels, k=k,
               gate_a_pass=res.passed, gate_b_pass=b_pass)
print(c.achieved_coverage, c.coverage_ci95, c.sharpness_confident_fraction)
```

---

## Classes

### `DiscretenessGateA`

```python
class DiscretenessGateA(ValidationGates):
    def __init__(
        self,
        seed: int = 42,
        n_ref: int = 200,
        n_bootstrap: int = 50,
        n_jobs: int = -1,
        alpha: float = 0.05,
        gap_n_ref: int = 50,
        k_range=(1, 2, 3, 4, 5, 6),
        **kw,
    )
```

Gate A v2 — discreteness null (single-Gaussian primary + gap + dip), with the
independence null demoted to a confound control. Subclasses `ValidationGates` and
reuses `ValidationGates.stability_test_bootstrap` verbatim for the observed,
reference, and confound-control statistics so all are strictly comparable.

**Constructor parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seed` | int | 42 | Random seed (forwarded to `ValidationGates`) |
| `n_ref` | int | 200 | Reference replicates for the single-Gaussian null (and for the confound control) |
| `n_bootstrap` | int | 50 | Bootstrap resamples per stability statistic |
| `n_jobs` | int | -1 | joblib parallelism (`prefer="threads"`) over reference draws |
| `alpha` | float | 0.05 | Significance level in the decision rule |
| `gap_n_ref` | int | 50 | Uniform reference draws per *k* in the gap statistic |
| `k_range` | tuple | `(1, 2, 3, 4, 5, 6)` | *k* values scanned by the gap statistic |
| `**kw` | — | — | Forwarded to `ValidationGates.__init__` (which is called with `show_progress=False`) |

#### `run()`

```python
def run(
    self,
    tumor: str,
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    gmm_seed: Optional[int] = None,
) -> GateAv2Result
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tumor` | str | required | Cohort label, echoed into the result |
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `n_clusters` | int | required | The fixed *k*, held identical for the observed statistic and every reference/gap replicate |
| `gmm_seed` | int or None | None | Mixture seed; defaults to `self.seed` |

**Returns:** `GateAv2Result`.

**What `run()` computes, in order:**

1. PCA-reduce the score matrix to `reduced_dim(n)` components (shared by observed data and every reference draw).
2. Fit `GaussianMixture(n_components=k, covariance_type="full", n_init=10, reg_covar=1e-6)` in the reduced space; take the observed mean bootstrap-ARI plus a 95% CI from `std_ari / sqrt(n_bootstrap)`.
3. **(A) Primary** — draw `n_ref` datasets from a single Gaussian `N(mean, diag(var))` with the retained per-PC variances, measure stability at the same fixed *k*, and form `sg_ref_p95` and the add-one empirical *p*.
4. **Confound control** — the demoted independence null: independently permute each pathway column, same statistic, reported only.
5. **(B)** Gap statistic over `k_range` with a uniform reference on the PCA-aligned bounding box; gap-optimal *k* by the first-standard-error rule.
6. **(C)** Hartigan dip test on PC1 and on the mixture discriminant axis (the direction between the two most separated component means).
7. k-stability routing, then the decision rule.

**Verdict routing:**

| Condition | `verdict` |
|-----------|-----------|
| `k_stability_frac < 0.5` | `"not-testable (no reproducible k)"` |
| testable and `passed` | `"discrete structure"` |
| testable and not `passed` | `"no discrete structure (continuum)"` |

Note the not-testable branch is checked **first**: an abstention overrides `passed`.

**Extra attributes attached to the returned result** (for plotting; deliberately *not*
in `to_dict()`): `sg_distribution`, `fp_distribution`, `gap_full`, `observed_labels`,
`Z`, `discriminant_axis_proj`.

---

### `ReframedMembershipGate`

```python
class ReframedMembershipGate:
    def __init__(
        self,
        seed: int = 42,
        n_splits: int = 80,
        primary_coverage: float = 0.90,
        tol: float = 0.05,
        min_n: int = 30,
    )
```

Gate C at a single 0.90 operating point, with a coverage CI, a min-n rule, and
explicit conditional-on-A/B framing. Wraps `ConformalMembershipGate`, evaluating only
at the primary target.

**Constructor parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seed` | int | 42 | Random seed (splits and the bootstrap CI) |
| `n_splits` | int | 80 | Repeated calibration/test splits |
| `primary_coverage` | float | `PRIMARY_COVERAGE` = 0.90 | The single, pre-agreed operating point |
| `tol` | float | `COVERAGE_TOL` = 0.05 | \|achieved − target\| ≤ tol ⇒ `coverage_achieved` |
| `min_n` | int | `MIN_N_GATE_C` = 30 | Below this, the 3-way split-conformal partition is reported as under-powered |

#### `run()`

```python
def run(
    self,
    tumor: str,
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    k: int,
    gate_a_pass: Optional[bool] = None,
    gate_b_pass: Optional[bool] = None,
) -> GateCReframedResult
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tumor` | str | required | Cohort label, echoed into the result |
| `pathway_scores` | pd.DataFrame | required | Pathway scores (samples x pathways) |
| `cluster_labels` | np.ndarray | required | The partition whose assignment sharpness is measured |
| `k` | int | required | Number of subtypes |
| `gate_a_pass` | bool or None | None | Gate A verdict, recorded in `conditional_on` |
| `gate_b_pass` | bool or None | None | Gate B verdict, recorded in `conditional_on` |

**Returns:** `GateCReframedResult`. When neither gate flag is supplied,
`conditional_on` is the string `"not conditioned"`.

**Interpretation caveat carried by the class itself:** conformal prediction guarantees
coverage of whatever label it is given; here that label is the mixture assignment, so
coverage is guaranteed by construction and cannot, on its own, distinguish a real
subtype from a sliced continuum. The result quantifies how sharply patients are
assigned **given** a partition — meaningful only when that partition is itself
certified.

---

### `ConformalMembershipGate`

```python
class ConformalMembershipGate:
    def __init__(
        self,
        seed: int = 42,
        n_splits: int = 50,
        target_coverages=(0.80, 0.90, 0.95),
        cal_fraction: float = 0.5,
        primary_coverage: float = 0.90,
        calib_tol: float = 0.05,
    )
```

The underlying per-patient conformal subtype-membership gate that
`ReframedMembershipGate` wraps. Its `membership_gate(tumor, pathway_scores,
cluster_labels, k) -> GateCResult` performs a proper 3-way split-conformal procedure
(train a `LogisticRegression` subtype classifier / calibrate the nonconformity
quantile / evaluate coverage), reporting empirical coverage, mean prediction-set size,
and the singleton ("confident") fraction at each target coverage, plus per-patient
confident flags at the primary coverage.

Prefer `ReframedMembershipGate` for reporting: it fixes the single 0.90 operating
point, avoiding the target-shopping the three-target output invites.

---

## Functions

### `reduced_dim()`

```python
def reduced_dim(n: int) -> int
```

`d = min(10, floor(n/10))`, floored at 2.

---

### `reduce_scores()`

```python
def reduce_scores(scores: np.ndarray, d: int, seed: int) -> tuple[np.ndarray, PCA]
```

PCA-reduce a samples x pathways matrix to `d` components (centred). `d` is
additionally capped at the number of available features and at `n_samples - 1`, so
narrow panels and very small *n* are handled without error.

**Returns:** `(Z, pca)` — the reduced score matrix and the fitted `PCA` object.

---

### `gap_statistic()`

```python
def gap_statistic(Z: np.ndarray, k_range, seed: int, n_ref: int = 50) -> Dict[str, Any]
```

Gap statistic over a uniform reference on the PCA-aligned bounding box (`Z` is already
in PCA space, so the principal axes are the coordinate axes and the box is
axis-aligned).

**Returns:** a dict with keys `k_range` (list), `gap` (list), `gap_se` (list), `logW`
(list), `gap_optimal_k` (int, by the first-standard-error rule: smallest *k* with
`gap(k) >= gap(k+1) - s(k+1)`).

> This is a **reported diagnostic**. It does not enter the Gate A decision rule.

---

## Dataclasses

### `GateAv2Result`

| Attribute | Type | Description |
|-----------|------|-------------|
| `tumor` | str | Cohort label |
| `n` | int | Number of samples |
| `n_clusters` | int | Fixed *k* used |
| `reduced_d` | int | Effective PCA dimension after capping |
| `observed_stability` | float | Observed mean bootstrap-ARI |
| `observed_ci95` | List[float] | 95% CI on the observed statistic |
| `sg_ref_mean` | float | Single-Gaussian reference mean |
| `sg_ref_p05` | float | Reference 5th percentile |
| `sg_ref_p50` | float | Reference median |
| `sg_ref_p95` | float | **Reference 95th percentile — half the decision rule** |
| `sg_empirical_p` | float | **Add-one empirical *p* — the other half** |
| `sg_effect_z` | float | (observed − reference mean) / reference SD; `nan` if SD = 0 |
| `passed` | bool | `obs > sg_ref_p95 and sg_empirical_p < alpha` |
| `gap_optimal_k` | int | Gap-optimal *k* (diagnostic) |
| `gap_supports_k` | bool | `gap_optimal_k >= 2 and gap_optimal_k == n_clusters` (diagnostic) |
| `gap_at_k` | float | Gap value at the fixed *k*; `nan` if *k* outside `k_range` (diagnostic) |
| `dip_pc1_p` | float | Dip-test *p* on PC1 (diagnostic; `nan` without `diptest`) |
| `dip_discriminant_p` | float | Dip-test *p* on the discriminant axis (diagnostic) |
| `unimodal_flag` | bool | Both dip *p* > 0.05 (diagnostic) |
| `k_stability_frac` | float | Fraction of bootstraps whose silhouette-optimal *k* equals the modal *k* |
| `modal_k` | int | Modal silhouette-optimal *k* across bootstraps |
| `modal_matches_fixed` | bool | `modal_k == n_clusters` |
| `testable` | bool | `k_stability_frac >= 0.5` |
| `verdict` | str | See the verdict routing table above |
| `fp_ref_mean` | float | Confound-control (independence null) mean |
| `fp_ref_p95` | float | Confound-control 95th percentile |
| `fp_empirical_p` | float | Confound-control empirical *p* |
| `fp_would_pass` | bool | Whether the demoted independence null would have passed |
| `n_ref` | int | Finite single-Gaussian reference draws actually used |
| `n_bootstrap` | int | Bootstrap resamples per statistic |
| `interpretation` | str | Human-readable summary of all four blocks |

**Methods:** `to_dict()`

### `GateCReframedResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `tumor` | str | Cohort label |
| `n` | int | Number of samples |
| `k` | int | Number of subtypes |
| `conditional_on` | str | Gate A / Gate B conditioning string, or `"not conditioned"` |
| `target_coverage` | float | The single operating point (0.90) |
| `achieved_coverage` | float | Mean empirical coverage across splits |
| `coverage_ci95` | List[float] | Bootstrap 95% CI on achieved coverage |
| `coverage_achieved` | bool | \|achieved − target\| ≤ `tol` |
| `sharpness_confident_fraction` | float | Singleton-set fraction at 0.90 |
| `sharpness_ci95` | List[float] | Bootstrap 95% CI on the sharpness |
| `mean_set_size` | float | Mean prediction-set size at 0.90 |
| `n_splits_used` | int | Splits that produced a usable coverage value |
| `under_powered` | bool | `n < min_n` |
| `interpretation` | str | Human-readable summary |

**Methods:** `to_dict()`

---

## Module Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `MIN_N_GATE_C` | 30 | Minimum *n* for a trustworthy 3-way split-conformal partition: train (≥ k·2, ~45%) + calibration (≥ 10) + test (≥ 5) |
| `PRIMARY_COVERAGE` | 0.90 | The single Gate C operating point |
| `COVERAGE_TOL` | 0.05 | Tolerance for "coverage achieved" |

Defined in `pathway_subtyping.discreteness.gate_c_reframed`.

---

## Small-*n* Hardening

Fixed in advance (not tuned) and applied **identically** to observed data and every
reference draw:

| Control | Rule |
|---------|------|
| Dimensionality | reduce to `d = min(10, floor(n/10))` PCA components, capped at n_pathways and n−1; the single Gaussian is simulated in this PCA space with a diagonal covariance equal to the retained per-PC variances |
| Fixed *k* | *k* is the inherited selected value, identical for the observed statistic and every reference/gap replicate; never re-selected per replicate |
| k-stability routing | per observed bootstrap, record the silhouette-optimal *k* (full-covariance BIC over-selects in the reduced space); if the modal *k* recurs in < 50% of resamples, route to `not-testable` |
| Reporting | full reference distribution plus a 95% CI on the observed statistic |

---

## Optional Dependency

The dip test requires the `diptest` extra:

```bash
pip install pathway-subtyping[discreteness]
```

`diptest` is GPLv2+ and is intentionally **not** a core dependency, so the default
(MIT) install stays MIT-clean. When it is absent the dip fields are returned as `nan`
and the gate proceeds unchanged — the decision rule never used them.

---

## References

- **SigClust (primary reference)**: Liu Y, Hayes DN, Nobel A, Marron JS. Statistical
  Significance of Clustering for High-Dimension, Low-Sample Size Data. *JASA*
  103(483):1281–1293, 2008. doi:10.1198/016214508000000454
- **Gap statistic (diagnostic)**: Tibshirani R, Walther G, Hastie T. Estimating the
  Number of Clusters in a Data Set Via the Gap Statistic. *JRSS B* 63(2):411–423,
  2001. doi:10.1111/1467-9868.00293
- **Dip test (diagnostic)**: Hartigan JA, Hartigan PM. The Dip Test of Unimodality.
  *Ann. Statist.* 13, 1985. doi:10.1214/aos/1176346577

---

## See Also

- [Discreteness-aware Gate A](../discreteness_gate.md) — rationale, small-*n* hardening table, and the implication for prior "reproducible subtype" claims
- [Validation Gates](validation.md) — the base `ValidationGates` class and the original stability gate
- [Uncertainty](uncertainty.md) — `ConformalPathwayPredictor`, the split-conformal machinery Gate C builds on
- [Deep-Learning Baselines](clustering_dl.md) — the gate is clusterer-agnostic; DL partitions need it too
