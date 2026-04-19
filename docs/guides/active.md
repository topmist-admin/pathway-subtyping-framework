# Active Sample Selection (F12)

> v0.6 Phase 3 F12 — given a pool of unlabelled samples and a fixed label budget, pick the samples whose labels would most reduce PSF subtype uncertainty.

**No extra required** — the `active` layer is part of the core install.

**When to use:** you have partially-labelled pathway-score data (most samples are unlabelled, a few have ground-truth labels) and a budget to send more samples to wet-lab / clinical review. You want to pick the next ones strategically rather than randomly.

Three strategies ship; pick based on what you care about:

| Strategy | Picks | Uses | Best when |
|---|---|---|---|
| `uncertainty` | Samples the current model is least sure about | `probs` (e.g. Bayesian-GMM posterior) or `uncertainty_scores` | You trust your current subtype probabilities and want to refine low-confidence regions |
| `diversity` | Samples farthest from the currently-labelled pool (greedy k-centre) | `features` only | You don't have a trustworthy uncertainty estimate; you want broad coverage |
| `hybrid` | Convex combination: α · uncertainty + (1 − α) · diversity | both | Usually the right default; `alpha=0.5` is a sensible start |

**API cheat sheet:**

| Class | Role |
|---|---|
| `ActiveSampleSelector` | User-facing entry point; `__init__(strategy, alpha, seed)` + `.select(features, budget, probs=None, uncertainty_scores=None)`. |
| `SelectionResult` | Return type. `.selected_indices`, `.scores` (per-sample), `.strategy`, `.budget`, `.to_dict()`. |
| `SelectionStrategy` | Enum: `UNCERTAINTY`, `DIVERSITY`, `HYBRID`. Strings accepted (`'uncertainty'`, etc). |

---

## 5-minute example

```python
import pandas as pd
from pathway_subtyping.active import ActiveSampleSelector
from pathway_subtyping.uncertainty import BayesianPathwayGMM

# Samples × pathways (unlabelled pool)
X = pd.read_csv("pool_pathway_scores.csv", index_col=0)

# Fit a posterior subtype model to derive uncertainty (one option)
bgmm = BayesianPathwayGMM(n_components=4, n_samples=500)
bgmm.fit(X.to_numpy())
posterior_probs = bgmm.predict_proba(X.to_numpy())  # (n, n_components)

# Pick 30 informative samples under a hybrid strategy
selector = ActiveSampleSelector(strategy="hybrid", alpha=0.5, seed=42)
result = selector.select(
    features=X,
    budget=30,
    probs=posterior_probs,
)

chosen = X.iloc[result.selected_indices]
chosen.to_csv("outputs/active_selection_round1.csv")
print(f"Selected {len(result.selected_indices)} samples using {result.strategy}")
```

---

## Bring your own uncertainty

If you already have a custom uncertainty signal (e.g., conformal prediction interval widths from F1, or a classifier's margin), skip `probs` and pass `uncertainty_scores` directly:

```python
from pathway_subtyping.uncertainty import ConformalPathwayPredictor

# Train a conformal predictor, get per-sample interval widths as uncertainty
cp = ConformalPathwayPredictor(score_fn=my_regressor, coverage=0.90)
cp.calibrate(X_cal.to_numpy(), y_cal.to_numpy())
intervals = cp.predict(X.to_numpy())
widths = [iv.width for iv in intervals]

selector = ActiveSampleSelector(strategy="uncertainty")
result = selector.select(features=X, budget=30, uncertainty_scores=widths)
```

Higher `uncertainty_scores` values = more uncertain. `ActiveSampleSelector` sorts descending, picks top-`budget`.

---

## Diversity-only (cold start)

When you have no posterior and no prior labels — e.g., your first batch selection — use diversity to cover the feature space evenly:

```python
selector = ActiveSampleSelector(strategy="diversity", seed=42)
result = selector.select(features=X, budget=30)
# Uses greedy k-centre on the features matrix; no probs required.
```

Greedy k-centre is deterministic given `seed` — the first sample is picked at index `seed % n`, then each subsequent pick is the sample farthest from the current selection.

---

## Hybrid weighting

The hybrid strategy combines uncertainty and diversity at each step:

```
score(sample) = alpha * normalised_uncertainty
              + (1 - alpha) * normalised_min_distance_to_selected
```

Tune `alpha`:
- `alpha=1.0` → equivalent to uncertainty-only
- `alpha=0.0` → equivalent to diversity-only but seeded from the most-uncertain sample
- `alpha=0.5` (default) → balanced

The min-distance term is re-normalised at every step so the two signals stay comparable even as the selected set grows.

---

## Roadmap acceptance (synthetic)

On the v0.6 synthetic cohort (`tests/test_active.py`), a **40% label budget** (select 40/100 samples) hits **≥ 90% accuracy** on the downstream GMM classifier — demonstrating the budget efficiency the layer claims. See [`tests/test_active.py`](../../tests/test_active.py) for the exact acceptance assertion.

---

## Known gotchas

- **Diversity requires a metric-sensible feature matrix.** Greedy k-centre uses Euclidean distance on the raw features. If your pathways are on wildly different scales, z-score first.
- **`uncertainty` strategy silently falls back to random if neither `probs` nor `uncertainty_scores` is supplied.** The fallback is documented in the docstring but easy to miss — if you expect "uncertainty" and don't pass a signal, you'll get uniform scores. Always pass one of the two.
- **`hybrid` picks are sensitive to feature scale through the distance term.** If some pathways dominate distance simply by having larger variance, the diversity component over-weights them. Standardise or PCA-project features first if this is a concern.
- **Budget can exceed pool size** — it's silently clipped to `n_samples`. No error.

---

## See also

- [API reference: `active`](../api/active.md) — full class + method signatures
- [Notebook 30 — active_learning](../../examples/notebooks/30_active_learning.ipynb) — budget-vs-accuracy curves
- [F1 uncertainty guide](uncertainty.md) — the source of calibrated `probs` / `uncertainty_scores`
- Settles, B. *Active Learning Literature Survey.* UW-Madison CS Tech Report 1648, 2009.
