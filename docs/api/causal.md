# `causal` — Invariant Causal Prediction for pathway-level parent discovery (F11)

ICP-style identification of causal parents of a target pathway across multiple environments.

**Import path:** `pathway_subtyping.causal`
**Source:** [`src/pathway_subtyping/causal/`](../../src/pathway_subtyping/causal/)
**Guide:** [causal.md](../guides/causal.md) · [Notebook 29](../../examples/notebooks/29_causal_inference.ipynb)
**Install:** core (no extra needed)

---

## Public surface

```python
from pathway_subtyping.causal import (
    InvariantPathwayPredictor,
    CausalParentReport,
    invariance_pvalue,
)
```

---

## `InvariantPathwayPredictor`

```python
InvariantPathwayPredictor(
    alpha: float = 0.05,                    # invariance-test significance level
    max_subset_size: Optional[int] = 4,     # cap on subset enumeration
    include_empty_set: bool = True,         # test the null subset explicitly
)
```

**Method:**

```python
report = predictor.fit(
    X: pd.DataFrame,                  # samples × candidate parents (pathway scores)
    y: pd.Series,                     # target values, same index as X
    target_name: str,                 # for reporting; must NOT be a column of X
    environments: Sequence[Any],      # env label per sample (hashable)
) -> CausalParentReport
```

- **Requires ≥ 2 environments.** Single-environment fits raise.
- **`target_name` must not be a column of X.** The caller separates X and y explicitly.
- **Missing environment labels** (None / NaN) cause singleton-environment degeneracy — validate upstream.

**Subset enumeration:** with `p` candidates and `max_subset_size=k`, enumerates up to `sum(C(p, i) for i in 0..k)` subsets. Set `max_subset_size=None` for exhaustive search; see guide for practical sizing table.

---

## `CausalParentReport`

| Attribute / method | Returns | Description |
|---|---|---|
| `.target` | `str` | Echo of `target_name` |
| `.candidate_parents` | `List[str]` | All columns of X |
| `.invariant_subsets` | `List[FrozenSet[str]]` | Every subset whose residuals were invariant at level `alpha` |
| `.identifiable_parents` | `Set[str]` | Intersection of `.invariant_subsets` — the ICP identifiable set |
| `.pvalues` | `Dict[FrozenSet[str], float]` | Per-subset invariance p-values |
| `.n_invariant_subsets` | `int` | `len(.invariant_subsets)` |
| `.recall_against(ground_truth)` | `float` | Recall of `.identifiable_parents` vs a held-out gold set |
| `.precision_against(ground_truth)` | `float` | Precision of `.identifiable_parents` vs gold |
| `.summary()` | `str` | Human-readable multi-line summary |
| `.to_dict()` | `dict` | JSON-serialisable |

---

## `invariance_pvalue`

Low-level hook for testing a single subset directly.

```python
p = invariance_pvalue(
    X: np.ndarray,                    # (n, k) features for this subset (may be empty (n, 0))
    y: np.ndarray,                    # (n,)
    environments: np.ndarray,         # (n,) hashable labels
    mean_variance_split: bool = True, # combined mean + variance invariance test
) -> float
```

- Returns the combined p-value (mean + variance invariance). The small-subset default test is linear regression residuals; richer residual models live behind the `[causal]` extra (future work).
- Useful when you have a hypothesis about a specific subset and don't want the combinatorial search.

---

## Acceptance gate (roadmap F11)

- On a 2-environment synthetic benchmark with planted causal structure, recall = **1.0** on the identifiable parents; precision remains controlled at `alpha`.
- Test: [`tests/test_causal.py`](../../tests/test_causal.py).

---

## Known gotchas (restated from the guide for quick lookup)

- Need ≥ 2 environments. Require environment diversity in P(X) for power.
- `target_name not in X.columns` — caller must separate.
- Default residual model is linear; non-linear relationships may under-test.
- Empty `.identifiable_parents` is an honest "insufficient environment variation" answer, not a bug.

---

## See also

- [Guide: causal](../guides/causal.md) — beginner-facing walkthrough + sizing guidance
- [`genesets`](genesets.md) — complementary: similarity-based (correlational) expansion
- Peters, Bühlmann, Meinshausen (2016). *JRSS-B* — canonical ICP reference.
