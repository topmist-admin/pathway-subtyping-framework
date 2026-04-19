# `active` — Active-learning sample selection (F12)

Pick the next samples to label under a fixed budget. Three strategies: uncertainty (entropy on posterior), diversity (greedy k-centre), hybrid (convex combination).

**Import path:** `pathway_subtyping.active`
**Source:** [`src/pathway_subtyping/active/`](../../src/pathway_subtyping/active/)
**Guide:** [active.md](../guides/active.md) · [Notebook 30](../../examples/notebooks/30_active_learning.ipynb)
**Install:** core (no extra needed)

---

## Public surface

```python
from pathway_subtyping.active import (
    ActiveSampleSelector,
    SelectionStrategy,
    SelectionResult,
)
```

---

## `SelectionStrategy` (enum)

```python
class SelectionStrategy(str, Enum):
    UNCERTAINTY = "uncertainty"
    DIVERSITY = "diversity"
    HYBRID = "hybrid"
```

Strings accepted everywhere — `selector = ActiveSampleSelector(strategy="hybrid")` works.

---

## `ActiveSampleSelector`

```python
ActiveSampleSelector(
    strategy: SelectionStrategy | str = SelectionStrategy.UNCERTAINTY,
    alpha: float = 0.5,                  # hybrid weighting; alpha=1.0 == uncertainty-only
    seed: Optional[int] = None,          # determinism for diversity / hybrid
)
```

**Method:**

```python
result = selector.select(
    features: np.ndarray | pd.DataFrame,       # (n_samples, d)
    budget: int,                                # number of samples to return
    probs: Optional[np.ndarray] = None,         # (n_samples, n_components) posterior
    uncertainty_scores: Optional[np.ndarray] = None,  # (n_samples,) custom uncertainty (higher = more uncertain)
) -> SelectionResult
```

**Strategy contract:**

| Strategy | `features` required | `probs` or `uncertainty_scores` required |
|---|---|---|
| `uncertainty` | ✗ | ✓ (one of) |
| `diversity` | ✓ | ✗ |
| `hybrid` | ✓ | ✓ (one of) |

- `probs` is converted to per-sample entropy via `_categorical_entropy`. `uncertainty_scores` overrides `probs` if both supplied.
- Diversity uses greedy k-centre on raw `features` (Euclidean distance). Seed sample is `seed % n`.
- Hybrid renormalises the distance term at every step so uncertainty + distance stay on comparable scales.
- `budget > n_samples` is silently clipped.

---

## `SelectionResult`

| Attribute | Type | Description |
|---|---|---|
| `.selected_indices` | `np.ndarray` | Integer indices into the pool, in selection order |
| `.scores` | `pd.Series` | Per-sample score used for ranking (NaN for samples not scored by strategy) |
| `.strategy` | `SelectionStrategy` | Echo |
| `.budget` | `int` | Actual budget used (post-clipping) |
| `.to_dict()` | `dict` | JSON-serialisable |

---

## Roadmap acceptance (F12)

- On the v0.6 synthetic cohort with planted subtype structure, a **40% label budget** hits **≥ 90%** downstream GMM classifier accuracy. Test: [`tests/test_active.py`](../../tests/test_active.py).

---

## Gotchas (restated from guide)

- **Uncertainty strategy needs `probs` or `uncertainty_scores`** — falling through with neither gives uniform scores (silent random-like selection).
- **Diversity is Euclidean on raw features** — z-score first if pathway scales differ.
- **Hybrid sensitivity** to feature scale via the distance term — standardise if needed.
- **Deterministic given `seed`** for diversity + hybrid; uncertainty-only is deterministic given identical input.

---

## See also

- [Guide: active learning](../guides/active.md) — beginner-facing walkthrough with strategy decision table
- [`uncertainty`](uncertainty.md) — produces `probs` (via `BayesianPathwayGMM.predict_proba`) or `uncertainty_scores` (via conformal interval widths)
- Settles, B. (2009). *Active Learning Literature Survey.* UW-Madison CS Tech Report 1648.
