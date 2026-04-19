# Causal Inference (F11)

> v0.6 Phase 3 F11 — Invariant Causal Prediction (ICP) for pathway-level parent-of-target discovery across multiple environments.

**No extra required** — the `causal` layer is part of the core install.

**When to use:** you have pathway scores across multiple *environments* (cell lines, batches, donors, stimulus conditions, tissues) and you want to know which pathways are **causal parents** of a target pathway / readout — not just correlated with it. If the conditional distribution P(Y | X) is stable when you change environments but P(X) shifts, X is a plausible causal parent; if P(Y | X) changes with environment, X is confounded.

This is the right tool when correlational cascade scoring (v0.5) is too noisy and you have the environment variation to exploit.

**API cheat sheet:**

| Class / function | Role |
|---|---|
| `InvariantPathwayPredictor` | Fit ICP; returns `CausalParentReport`. |
| `invariance_pvalue(...)` | Low-level: test a single subset for invariant residuals across environments. |
| `CausalParentReport` | Result. Has `.identifiable_parents` (intersection of invariant subsets), `.invariant_subsets` (full set), `.to_dict()`, `.recall_against(...)`, `.precision_against(...)`. |

---

## 5-minute example

```python
import pandas as pd
from pathway_subtyping.causal import InvariantPathwayPredictor

# samples × pathways DataFrame. Each row is one sample from one environment.
X = pd.read_csv("pathway_scores_multi_env.csv", index_col=0)
# Target is one pathway column (or a derived outcome)
y = X.pop("TARGET_PATHWAY")
# Environment label per sample (cell line, batch, donor, ...)
envs = pd.read_csv("sample_metadata.csv", index_col=0)["cell_line"]

predictor = InvariantPathwayPredictor(alpha=0.05, max_subset_size=3)
report = predictor.fit(
    X=X,
    y=y,
    target_name="TARGET_PATHWAY",
    environments=envs.loc[X.index],
)

print(report.summary())
print("Identifiable causal parents:", report.identifiable_parents)
# -> {'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_E2F_TARGETS'}
```

ICP reports the **intersection** of all statistically-invariant subsets — this is the only subset the data can *identify* with confidence under the ICP assumptions. If no subset is invariant, `identifiable_parents` is empty and the correct answer is "insufficient environment variation to identify the causal structure."

---

## Choosing `max_subset_size`

`max_subset_size` caps the combinatorial search. With `p` candidate pathways and `max_subset_size=k`, ICP enumerates `sum(C(p, i) for i in 0..k)` subsets. Rough guidance:

| Candidates (`p`) | Suggested `max_subset_size` | Subsets searched |
|---|---|---|
| ≤ 10 | `None` (exhaustive) | up to 2^10 = 1024 |
| 10–30 | 3 | ~1,000 (fast) |
| 30–100 | 2 | ~5,000 (seconds) |
| 100+ | 1 | linear — only tests individual parents |

If the "true" causal parent set is larger than `max_subset_size`, ICP will return an empty `identifiable_parents` (nothing passed the invariance test at that subset size). That's an explicit null result, not a false claim — prefer widening `max_subset_size` or reducing the candidate pool by a correlational pre-filter before assuming the answer is "no parents."

---

## Interpreting `identifiable_parents` vs `invariant_subsets`

- `report.invariant_subsets` — **every** subset whose residuals were invariant across environments at level `alpha`. Often includes subsets that mutually contradict each other.
- `report.identifiable_parents` — the **intersection** of all those subsets. A pathway appears here only if it's in *every* invariant subset. This is the subset ICP provably identifies.

A common workflow: if `identifiable_parents` is small, inspect `invariant_subsets` directly — it often reveals that two different sets of parents each explain the data (i.e., the environments don't disambiguate them). That's useful scientific information in its own right.

---

## Benchmarking against a ground truth

When you have a synthetic or known-causal ground truth (e.g., a published regulatory edge):

```python
gold = {"HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS"}
print(f"Recall:    {report.recall_against(gold):.2f}")
print(f"Precision: {report.precision_against(gold):.2f}")
```

Recall and precision are computed against `report.identifiable_parents`.

---

## Known gotchas

- **You need ≥ 2 environments.** One environment means no variation in P(X), no ICP signal. Two is minimum; more is better — ICP power scales with the diversity of environments.
- **Environments must vary P(X).** If every environment has the same marginal distribution over candidates, ICP has nothing to test. In practice this means: don't use a single disease's samples across technical replicates as "environments" — use different donors, batches, or stimulus conditions.
- **ICP assumes correct model specification.** The default residual model is linear regression with heteroskedasticity-adjusted invariance testing. If your true relationship is non-linear, ICP under-tests for invariance and may miss parents. The `[causal]` extra exposes richer residual models; see the API docs.
- **`target_name` must not be a column of `X`.** The caller is responsible for separating X from y (common bug: passing the full `pathway_scores` without dropping the target column → `ValueError`).
- **Missing environment labels are a silent failure mode.** If any sample's environment is `None` / `NaN`, ICP treats each as its own singleton environment and the test degenerates. Always validate that `environments` is dense and categorical before calling `.fit()`.

---

## See also

- [API reference: `causal`](../api/causal.md) — `InvariantPathwayPredictor` + `invariance_pvalue` signatures
- [Notebook 29 — causal_inference](../../examples/notebooks/29_causal_inference.ipynb) — synthetic 2-environment walkthrough
- [F7 gene-set expansion](genesets.md) — complementary: F7 ranks by similarity (correlational), F11 ranks by invariance (causal).
- Peters, Bühlmann, Meinshausen. *Causal inference using invariant prediction: identification and confidence intervals.* JRSS-B, 2016.
