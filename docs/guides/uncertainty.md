# Uncertainty Quantification Guide

*PSF v0.6 Phase 1 — Rigor Layer*

Every MSV output in PSF (`maturation_score`, `stress_score`, `drift_score`,
etc.) is a single scalar. That leaves researchers unable to distinguish
"the model is confident this cell is in a stress state" from "the model is
uncertain; the true value could plausibly sit anywhere from 0.4 to 0.9."

The `pathway_subtyping.uncertainty` module closes this gap with three
complementary uncertainty channels plus a calibration assessment toolkit.

See [examples/notebooks/21_uncertainty.ipynb](../../examples/notebooks/21_uncertainty.ipynb)
for a runnable walkthrough of everything described here.

---

## When to reach for which tool

| Question you're trying to answer | Use |
|---|---|
| "What is the plausible range of this pathway score?" | `ConformalPathwayPredictor` |
| "How sensitive is the group-level score to the particular cells I sampled?" | `BootstrapMSV` |
| "How confident is this cell's subtype assignment?" | `BayesianPathwayGMM` |
| "Are the probabilities my classifier reports actually calibrated?" | `CalibrationReport` |

The four tools solve overlapping but distinct problems. Conformal gives
finite-sample coverage guarantees without any model assumption; bootstrap
gives sampling-variability intervals on any statistic; Bayesian GMM gives
subtype-assignment uncertainty under a generative mixture; calibration
metrics grade how truthful reported probabilities are.

---

## ConformalPathwayPredictor

Distribution-free prediction intervals that wrap any scoring function.
Requires a held-out **calibration set** — a second disjoint dataset
separate from the model's training data.

```python
from pathway_subtyping.uncertainty import ConformalPathwayPredictor

predictor = ConformalPathwayPredictor(
    score_fn=my_pathway_scorer,   # X -> np.ndarray
    coverage=0.9,                  # target marginal coverage
)
predictor.calibrate(X_calibration, y_calibration)

intervals = predictor.predict(X_test)  # list[ConformalInterval]
empirical = predictor.coverage_on(X_test, y_test)
```

### Guarantee

Under the exchangeability assumption, split-conformal gives
P(y_test in interval) ≥ 1 - α, where α = 1 - coverage. The guarantee is
**finite-sample** and **distribution-free** — no parametric assumptions
about the score function, no asymptotics.

### Modes

| `mode` | Score function returns | Interval format |
|---|---|---|
| `"regression"` (default) | `(n,)` float predictions | `(lower, upper)` bracket |
| `"classification"` | `(n, n_classes)` probabilities | `label_set` of plausible classes |

### Gotchas

- **Exchangeability matters.** Calibration and test data must be drawn
  exchangeably. If the test distribution shifts (new batch, new platform),
  conformal guarantees are void.
- **Over-coverage is expected on small calibration sets.** The
  finite-sample correction `⌈(n+1)(1-α)⌉ / n` is conservative.
- **A minimum of 10 calibration samples is required.** Below this
  `calibrate()` raises.

---

## BootstrapMSV

Non-parametric bootstrap over cells. Resamples the input with replacement,
recomputes any MSV score, reports percentile-based intervals.

```python
from pathway_subtyping.uncertainty import BootstrapMSV

bootstrap = BootstrapMSV(n_bootstrap=500, ci_level=0.95, seed=42, n_jobs=1)
result = bootstrap.run(expression_df, score_fn=compute_stress_score)

result.point   # point estimate on the full input
result.lower   # 2.5th percentile across replicates
result.upper   # 97.5th percentile
result.se()    # standard error
```

### Two modes

- **Default (`per_cell=False`)**: `score_fn` returns a scalar or vector
  (e.g., one score per pathway). Bootstrap intervals describe the
  sampling distribution of that aggregate statistic.
- **Per-cell (`per_cell=True`)**: `score_fn` returns one value per input
  cell. Intervals are reported cell-by-cell by aggregating all replicates
  in which that cell was included.

### Parallelism

`n_jobs > 1` dispatches resampling to a thread pool. Thread-based rather
than process-based so the parent score function's state (PyTorch, scanpy,
etc.) is available without pickling overhead. Switch to `joblib` for
process-based parallelism if your `score_fn` is CPU-bound pure-Python.

---

## BayesianPathwayGMM

Drop-in replacement for the point-estimate GMM used in
`pathway_subtyping.clustering` and `pathway_subtyping.pipeline`. Inference
is variational (sklearn's `BayesianGaussianMixture`). Same
`.fit() / .predict() / .predict_proba() / .score_samples()` API, plus
`.sample_assignments()` and `.sample_parameters()` for posterior draws.

```python
from pathway_subtyping.uncertainty import BayesianPathwayGMM

model = BayesianPathwayGMM(
    n_components=3,
    covariance_type="full",
    weight_concentration_prior_type="dirichlet_process",
    random_state=42,
).fit(X)

# Point-estimate-compatible interface:
labels = model.predict(X)           # mode-collapsed argmax
probs  = model.predict_proba(X)     # mean posterior per cell

# Posterior sampling:
draws  = model.sample_assignments(X, n_samples=100)   # (100, n_cells) ints
params = model.sample_parameters(n_samples=50)        # list[PosteriorSample]
```

### When to trust component count

The Dirichlet-process prior prunes unused components. Read
`model.n_effective_components` (components with weight > 1e-3) instead of
the `n_components` you passed in — the latter is an upper bound, not a
committed cardinality.

### Reproducibility with point-estimate GMM

On well-identified synthetic data, mode-collapsed `BayesianPathwayGMM`
reproduces sklearn `GaussianMixture` predictions within 1% (adjusted
Rand index > 0.99) — the acceptance target from the v0.6 roadmap.
Where the two diverge in practice: heavily overlapping components and
small sample sizes, where the Bayesian prior regularises the solution.

---

## CalibrationReport

Probability calibration diagnostics in one call:

```python
from pathway_subtyping.uncertainty import CalibrationReport

report = CalibrationReport.from_predictions(y_true, y_prob, n_bins=10)
report.ece          # expected calibration error
report.brier        # Brier score
report.reliability  # bin centers + observed frequencies
report.plot()       # reliability diagram (matplotlib)
```

### Interpreting the numbers

- **ECE**: weighted gap between mean predicted probability and observed
  frequency, binned over the probability axis. A perfectly calibrated
  classifier has ECE = 0. Typical "good" values are < 0.05.
- **Brier score**: mean squared error between probability and outcome.
  Combines discrimination and calibration; useful as a single summary.
  A random predictor gives 0.25 on balanced binary problems.
- **Reliability diagram**: a calibrated classifier's curve lies on the
  diagonal. Curves below the diagonal = over-confident; above =
  under-confident.

---

## Acceptance criteria (v0.6 roadmap)

The layer is considered complete when:

- [x] Conformal intervals achieve nominal coverage within ±2% on held-out
      synthetic data (mean across seeds).
- [x] Bootstrap intervals stable to < 5% relative width variance across
      3 independent resampling seeds.
- [x] `BayesianPathwayGMM` reproduces point-estimate GMM within 1%
      (ARI > 0.99) on well-identified synthetic data.
- [x] Conformal coverage within ±2% on TCGA-COAD (n=57) and GSE28521
      autism benchmark (n=79) against the finite-sample oracle.
- [x] No breaking changes to existing scoring APIs.

### Real-data results

`scripts/validate_f1_real_data.py` computes hallmark (TCGA-COAD) or curated
autism (GSE28521) pathway scores and, for each pathway, uses the other
pathway scores as features to predict it — then calibrates conformal
intervals and measures empirical coverage on a held-out split. Reported
across 10–20 seeds per pathway.

| Dataset | Samples | Pathways | Target | Oracle dev | Nominal dev | Runs |
|---|---|---|---|---|---|---|
| TCGA-COAD  | 57 | 50 | 0.80 | +0.0004 | +0.0152 | 1000 |
| TCGA-COAD  | 57 | 50 | 0.90 | −0.0031 | +0.0228 | 1000 |
| TCGA-COAD  | 57 | 50 | 0.95 | −0.0040 | +0.0089 | 1000 |
| GSE28521   | 79 | 15 | 0.80 | −0.0006 | +0.0102 |  300 |
| GSE28521   | 79 | 15 | 0.90 | +0.0009 | +0.0198 |  300 |
| GSE28521   | 79 | 15 | 0.95 | −0.0007 | +0.0223 |  300 |

**Oracle deviation** is the mean of (empirical coverage − finite-sample
oracle coverage), where the oracle is `⌈(n_cal+1)(1-α)⌉ / (n_cal+1)` — the
smallest realized coverage split-conformal can produce at a given `n_cal`.
The nominal-deviation excess of ~2% at 0.90/0.95 on these cohorts is
purely the discretization offset between the nominal and oracle levels at
small `n_cal` (~25–35). Both datasets pass the roadmap criterion against
the oracle.

Artifacts: `results/f1_validation/TCGA-COAD_conformal_coverage.json` and
`results/f1_validation/GSE28521_conformal_coverage.json`. Acceptance is
asserted by `tests/test_uncertainty_real_data.py`, which skips if the
artifacts are absent so CI without the underlying datasets still passes.

---

## See also

- [examples/notebooks/21_uncertainty.ipynb](../../examples/notebooks/21_uncertainty.ipynb)
  — runnable walkthrough
- [docs/roadmap-v06-codeberg.md](../roadmap-v06-codeberg.md) — v0.6 Rigor
  Layer roadmap
- Angelopoulos & Bates (2021). *A Gentle Introduction to Conformal
  Prediction and Distribution-Free Uncertainty Quantification.*
  arXiv:2107.07511.
