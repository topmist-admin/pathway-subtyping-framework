# `uncertainty` — Calibrated pathway-score uncertainty (F1)

Distribution-free prediction intervals, non-parametric bootstrap summaries, Bayesian posterior sampling, and calibration diagnostics for PSF pathway scores.

**Import path:** `pathway_subtyping.uncertainty`
**Source:** [`src/pathway_subtyping/uncertainty/`](../../src/pathway_subtyping/uncertainty/)
**Guide:** [uncertainty.md](../guides/uncertainty.md) · [Notebook 21](../../examples/notebooks/21_uncertainty.ipynb)
**Install:** core (no extra needed)

---

## Public surface

```python
from pathway_subtyping.uncertainty import (
    # Split-conformal prediction
    ConformalInterval, ConformalPathwayPredictor,
    # Bootstrap
    BootstrapMSV, BootstrapResult,
    # Bayesian
    BayesianPathwayGMM, PosteriorSample,
    # Calibration diagnostics
    CalibrationReport, brier_score, ece, reliability_curve,
)
```

---

## `ConformalPathwayPredictor`

Split-conformal wrapper around any score function. Returns finite-sample-oracle-valid prediction intervals for pathway-score regression.

```python
ConformalPathwayPredictor(
    score_fn: Callable[[np.ndarray], np.ndarray],
    coverage: float = 0.90,
)
```

- `score_fn` — any callable mapping `X (n, d)` to `y_hat (n,)` (a fitted linear model, sklearn estimator, `msv_head.transform`, etc.). Not differentiated or re-trained by the predictor; treat as a black box.
- `coverage` — nominal target in `(0, 1)`. `0.90` = 90% intervals.

**Methods:**

| Method | Returns | Notes |
|---|---|---|
| `.calibrate(X_cal, y_cal)` | `self` | Computes the conformity quantile on the held-out calibration split. Must be called before `.predict`. |
| `.predict(X)` | `List[ConformalInterval]` | One interval per row of `X`. Each has `.low`, `.high`, `.point`, `.width`, `.contains(y)`. |
| `.coverage_on(X_test, y_test)` | `float` | Empirical coverage on a held-out test set. Use to validate the nominal target. |
| `.quantile` | `float \| None` | The calibration quantile actually used. None before `.calibrate` is called. |
| `.n_calibration` | `int` | Size of the calibration split. |

**Finite-sample oracle coverage:** for `n_cal` calibration points, the *actual* achievable upper coverage is `ceil((n_cal + 1) * coverage) / (n_cal + 1)`. The roadmap gate compares against this oracle, not the nominal target, because discretisation dominates at small `n_cal`. See [validate_f1_real_data.py](../../scripts/validate_f1_real_data.py) for the exact calculation.

---

## `ConformalInterval`

Dataclass for a single prediction interval.

```python
interval = ConformalInterval(low=0.42, high=0.71, point=0.56)
interval.width         # 0.29
interval.contains(0.5) # True
interval.to_dict()     # {'low': 0.42, 'high': 0.71, 'point': 0.56, 'width': 0.29}
```

---

## `BootstrapMSV`

Non-parametric bootstrap over cells / samples. Produces per-element percentile intervals (pathway-wise or per-cell) and standard errors.

```python
BootstrapMSV(
    n_bootstrap: int = 200,
    coverage: float = 0.95,
    n_jobs: int = 1,
    seed: Optional[int] = None,
)
```

**Method:**

```python
result = BootstrapMSV(n_bootstrap=500, seed=42).run(
    X=pathway_scores,               # (n_cells, n_pathways)
    statistic=np.mean,              # or any function mapping (n, d) → (d,) or (n, d)
    mode="stat",                    # "stat" (per-pathway) or "per_cell" (per-cell × pathway)
)
```

**`BootstrapResult`:**

| Attribute | Shape | Description |
|---|---|---|
| `.point_estimate` | `(d,)` or `(n, d)` | The statistic applied to the full sample. |
| `.low` / `.high` | same as point | Bootstrap percentile bounds at `coverage`. |
| `.se` | same as point | Bootstrap standard error. |
| `.width` | same as point | `.high - .low`. |
| `.summary()` / `.to_dict()` | — | String / JSON-serialisable summary. |

Per-cell mode is memory-heavy on large cohorts; prefer `mode="stat"` when you only need pathway-level aggregate uncertainty.

---

## `BayesianPathwayGMM`

Posterior-sampling GMM for pathway-score subtyping. Returns draws over component assignments rather than a single MAP clustering.

```python
BayesianPathwayGMM(
    n_components: int,
    n_samples: int = 500,
    seed: Optional[int] = None,
)
```

**Methods:**

- `.fit(X)` — run variational inference + posterior sampling
- `.predict_proba(X)` — `(n_samples, n_components)` posterior probabilities
- `.posterior_samples` — `List[PosteriorSample]` with component means / weights per draw

Use the per-sample posterior probabilities as the `probs=` argument to [`ActiveSampleSelector`](active.md) for label-efficient subtype refinement.

---

## Calibration diagnostics

Functional API for probabilistic classifiers (any predictor returning `(n, n_classes)` probabilities).

```python
from pathway_subtyping.uncertainty import reliability_curve, ece, brier_score, CalibrationReport

curve = reliability_curve(y_true, y_prob, n_bins=10)   # (bin_centres, observed_freq)
err   = ece(y_true, y_prob, n_bins=10)                  # Expected Calibration Error
brier = brier_score(y_true, y_prob)                     # Brier score

report = CalibrationReport.from_probabilities(y_true, y_prob)
print(report)                                           # pretty-printed summary
```

`CalibrationReport` packages ECE + Brier + reliability curve + per-bin diagnostics so you can drop one call into a notebook and get the full calibration story.

---

## Acceptance gate (roadmap F1)

- Conformal oracle deviation within ±2% of nominal target on held-out TCGA-COAD + GSE28521 splits.
- Asserted in [`tests/test_uncertainty_real_data.py`](../../tests/test_uncertainty_real_data.py).
- Reproduce: `python scripts/validate_f1_real_data.py`.

---

## See also

- [Guide: uncertainty](../guides/uncertainty.md) — beginner-facing walkthrough
- [`active`](active.md) — consumes `BayesianPathwayGMM.predict_proba` for label selection
- [`perturb`](perturb.md) — wraps `ConformalPathwayPredictor` around the MSV head (see F5 conformal coverage)
