# Time-sliced KG Sensitivity API (Gate K)

> **Module**: `pathway_subtyping.kg_sensitivity`

Tests whether a finding survives being recomputed against a **different knowledge-graph version**. Every other gate holds the KG fixed and varies the data; this one varies the KG. A subtype that appears under Reactome 2026 and vanishes under Reactome 2025 is a property of the knowledge base, not of the cohort — and it will pass bootstrap stability, the discreteness gate and confound association, because all three are computed against the single KG version that produced it.

For the rationale, the reading rule and the limits, see **[discreteness-style narrative doc: `../kg_sensitivity_gate.md`](../kg_sensitivity_gate.md)**. This page is the API reference.

---

## Quick Example

```python
from pathway_subtyping import kg_timeslice_sensitivity

def partition_fn(kg, cohort):
    """Must be deterministic — seed any internal clustering yourself."""
    scores = score_pathways_with_kg(kg, cohort)
    return cluster(scores, k=3, seed=42)

res = kg_timeslice_sensitivity(
    kg_reactome_2025,     # v1 — baseline
    kg_reactome_2026,     # v2 — comparison
    partition_fn,
    cohort,
    n_null=50,
    ari_min=0.80,
    seed=42,
    rewiring="degree",    # "uniform" | "degree" (default) | "module"
)

print(res.verdict)        # robust | kg-sensitive | generically-fragile | not-testable (...)
print(res.observed_ari, res.null_median, res.empirical_p)
print(res.rewiring)       # which null produced the verdict — always report it
print(res.diff.summary)   # what actually changed between the two graphs
```

---

## `kg_timeslice_sensitivity`

```python
kg_timeslice_sensitivity(
    v1, v2, partition_fn, cohort, *,
    n_null=50, ari_min=0.80, alpha=0.05, seed=42,
    rewiring="degree", score_fns=None, tolerance=0.05,
) -> KGSensitivityResult
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `v1` | `KnowledgeGraph` | — | Baseline graph (e.g. Reactome 2025) |
| `v2` | `KnowledgeGraph` | — | Comparison graph (e.g. Reactome 2026) |
| `partition_fn` | `PartitionFn` | — | `f(kg, cohort) -> labels`. Must be deterministic |
| `cohort` | `Any` | — | Opaque payload passed through to `partition_fn` |
| `n_null` | `int` | `50` | Size-matched random rewirings to draw |
| `ari_min` | `float` | `0.80` | Agreement at or above which the finding is `robust` |
| `alpha` | `float` | `0.05` | Significance level for the empirical p |
| `seed` | `int` | `42` | Seeds the rewiring generator; the run is reproducible |
| `rewiring` | `str` | `"degree"` | Null model — see [Null models](#null-models) |
| `score_fns` | `Iterable[ScoreFn]` \| `None` | `None` | Optional scalar scorers forwarded to `run_kg_regression` |
| `tolerance` | `float` | `0.05` | Relative-delta tolerance for that regression report |

### The decision rule, in full

```python
if observed_ari >= ari_min:    verdict = "robust"
elif empirical_p < alpha:      verdict = "kg-sensitive"
else:                          verdict = "generically-fragile"
```

**Both terms are live.** `ari_min` decides robustness; `empirical_p` decides how a non-robust result is *explained*. Nothing else enters the rule — the `KGDiff` and the optional `run_kg_regression` report are context, **not** criteria, and a regression test asserts that a flagged scalar score leaves the verdict unchanged.

`empirical_p` is an add-one empirical p in the **low** tail: how often a size-matched random rewiring disrupts the partition at least as much as the real KG change did.

### Abstentions

Returns a `not-testable (...)` verdict, with `testable=False`, when:

- the two graphs are identical (nothing to test);
- either partition has fewer than two clusters;
- the two partitions cover different sample counts;
- ARI is undefined on the partitions;
- the null distribution could not be constructed;
- `rewiring="module"` but the graph has no `GENE_IN_PATHWAY` edges.

---

## `KGSensitivityResult`

| Attribute | Type | Description |
|---|---|---|
| `verdict` | `str` | `robust`, `kg-sensitive`, `generically-fragile`, or `not-testable (...)` |
| `passed` | `bool` (property) | True only for `robust` |
| `testable` | `bool` | False on abstention — **use this as the denominator filter** |
| `observed_ari` | `float` | Agreement between the v1 and v2 partitions |
| `null_aris` | `List[float]` | Agreement between v1 and each rewired partition |
| `null_median` | `float` | Median of `null_aris` |
| `null_p05` | `float` | 5th percentile of `null_aris` |
| `empirical_p` | `float` | Add-one empirical p, low tail |
| `n_clusters_v1` / `n_clusters_v2` | `int` | Cluster counts under each graph |
| `ari_min` / `alpha` | `float` | The thresholds used |
| `rewiring` | `str` | Which null produced this verdict |
| `diff` | `KGDiff` \| `None` | Structural diff between the graphs |
| `regression` | `KGRegressionReport` \| `None` | Optional scalar-score report |

Methods: `summary() -> str`, `to_dict() -> Dict[str, Any]`.

> **`testable` is not cosmetic.** Gate A abstained on 28 of 30 synthetic negatives and its resulting "FPR 0.000" was quoted against n=30 rather than the testable n=2. Any rate computed over many Gate K runs must use the testable subset as its denominator.

---

## Null models

Set with `rewiring=`. What the null preserves sets how strict it is, and **the choice can change the verdict on identical data**.

| Mode | Preserves | Use |
|---|---|---|
| `"uniform"` | node set, edge-type counts | comparison baseline only — the weakest null |
| `"degree"` **(default)** | + exact in/out-degree sequence | general use |
| `"module"` | + the diff's within-/cross-module split | when a release concentrated edits in a few pathways |

`degree` uses Maslov–Sneppen double-edge swaps (`a→b`, `c→d` become `a→d`, `c→b`), so every node keeps its exact in- and out-degree. Each swap consumes one removal and one addition, covering the *balanced* part of the budget; any residual is degree-**weighted**, which preserves the shape of the degree distribution but not exact values — strict preservation is impossible there, since adding an edge necessarily raises somebody's degree.

`uniform` destroys the degree sequence and misreads changes at both ends of it: hub edges look specially targeted when losing one is ordinary, and peripheral edges look ordinary when losing one is specific. On identical data with an identical observed ARI, `uniform` returns `generically-fragile` (p = 0.066) where `degree` returns `kg-sensitive` (p = 0.016).

---

## `rewire_kg`

```python
rewire_kg(
    kg, n_remove_by_type, n_add_by_type, rng, *,
    rewiring="degree", module_map=None, within_module_frac=None,
) -> KnowledgeGraph
```

Builds one size-matched perturbation. The input graph is **not** modified.

| Parameter | Type | Description |
|---|---|---|
| `kg` | `KnowledgeGraph` | Baseline to perturb |
| `n_remove_by_type` | `Dict[str, int]` | `{edge_type_value: count}` edges to delete |
| `n_add_by_type` | `Dict[str, int]` | `{edge_type_value: count}` edges to introduce |
| `rng` | `np.random.Generator` | Seeded; the same seed reproduces the same rewiring |
| `rewiring` | `str` | `"uniform"`, `"degree"`, or `"module"` |
| `module_map` | `Dict[str, frozenset]` \| `None` | Required for `"module"` |
| `within_module_frac` | `Dict[str, float]` \| `None` | Per-edge-type target within-module fraction |

Raises `ValueError` on an unknown mode, or on `"module"` without a non-empty `module_map`.

Node set and edge-type composition are always preserved. Candidate additions are drawn from nodes already participating in that edge type, keeping every proposal schema-valid; proposals that already exist or that the schema rejects are skipped rather than retried indefinitely, and a shortfall is logged at DEBUG rather than raised.

---

## Module helpers

### `module_map_from_pathways(kg) -> Dict[str, frozenset]`

Derives modules from the KG's own `GENE_IN_PATHWAY` membership: two genes share a module when they share a pathway. Needs no external community detection. Genes in no pathway map to an empty `frozenset` and share a module with nobody.

### `within_module_fractions(diff, module_map) -> Dict[str, float]`

Per-edge-type fraction of an observed `KGDiff` that is within-module. This is what `rewiring="module"` matches, so a release that concentrated 80% of its edits inside shared pathways yields `{edge_type: 0.8}` and the null places 80% of its own changes within-module.

---

## Constants

```python
VERDICT_ROBUST        = "robust"
VERDICT_KG_SENSITIVE  = "kg-sensitive"
VERDICT_FRAGILE       = "generically-fragile"

REWIRING_UNIFORM = "uniform"
REWIRING_DEGREE  = "degree"
REWIRING_MODULE  = "module"
```

---

## Limits

Gate K tests the KG versions **you supply**. A `robust` verdict means the finding survived that specific swap, not that it is invariant to curation in general — two adjacent releases that barely differ return `robust` almost by construction, which is why `res.diff.summary` should be reported alongside every verdict.

The null holds the *number* of changed edges fixed and, in `degree`/`module` mode, the degree sequence and module split. It does not model which edges a curator would plausibly touch beyond that — releases also concentrate by evidence type, publication recency and organism.

Finally, the gate inherits `partition_fn`'s determinism. A partition function with unseeded internal randomness produces a null that measures its own noise floor, and the gate cannot detect that on the caller's behalf.

---

## See Also

- [`../kg_sensitivity_gate.md`](../kg_sensitivity_gate.md) — rationale and reading rule
- [`validation.md`](validation.md) — the gate battery that holds the KG fixed
- [`discreteness.md`](discreteness.md) — Gate A, whose null-design lesson this module follows
