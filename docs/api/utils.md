# Utilities API

> **Module**: `pathway_subtyping.utils`

Cross-cutting helpers: deterministic seeding, memory-bounded readers for large cohorts, progress reporting, and the **ARI validity guards** that prevent degenerate ground truth from producing a misleadingly perfect score.

---

## Metrics — the ARI guards

These exist because of a real defect. `sklearn.metrics.adjusted_rand_score` returns **1.0** when the ground truth has a single cluster, which is not a perfect recovery but a degenerate comparison. That behaviour produced 14 invalid rows in a published 47-dataset benchmark — 13 of them carrying a spurious ARI of 1.0 — which in turn drove a calibration model that was later retracted. See [`../../CORRECTION_2026-07/ERRATUM_2026-07-08.md`](../../CORRECTION_2026-07/ERRATUM_2026-07-08.md).

**Use these in place of `adjusted_rand_score` anywhere ground truth is not guaranteed well-formed.** They are already wired into `benchmark.py` and `pipeline.py`.

### `safe_adjusted_rand_score`

```python
safe_adjusted_rand_score(
    labels_true, labels_pred, *,
    min_samples=2, min_true_clusters=2, default=math.nan,
) -> float
```

Drop-in replacement returning `default` (NaN) on degenerate input instead of a misleading score. Raises `ValueError` on length mismatch, which is a real caller bug rather than a degenerate input.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `labels_true` | array-like | — | Ground-truth labels |
| `labels_pred` | array-like | — | Predicted labels |
| `min_samples` | `int` | `2` | Below this, return `default` |
| `min_true_clusters` | `int` | `2` | Below this many *distinct true* clusters, return `default` |
| `default` | `float` | `nan` | Value returned on degenerate input |

> The guard keys on ground-truth **structure** (`n_true_clusters < 2`), not on the output value. A value-based check (`== 1.0`) would miss the degenerate row that happens to return 0.0 — which is exactly what one row in the benchmark did.

### `ari_with_validity`

```python
ari_with_validity(
    labels_true, labels_pred, *, min_samples=2, min_true_clusters=2,
) -> Tuple[float, bool, Optional[str]]
```

Same computation, but returns `(value, valid, reason)` so a caller can record *why* a comparison was excluded rather than silently propagating NaN.

---

## Seeding

### `set_global_seed(seed: int) -> None`

Sets the process-wide seed across the framework's random sources. Call once at entry.

### `get_rng(seed: Optional[int], module_name: str = "") -> numpy.random.RandomState`

Returns a seeded generator. `module_name` is used for logging so a run's seeding is auditable.

> Determinism is a hard requirement for several gates. `kg_timeslice_sensitivity`, for example, builds a null by repeatedly re-partitioning; a partition function with unseeded internal randomness makes that null measure its own noise floor, and the gate cannot detect that on the caller's behalf.

---

## Performance

| Function | Purpose |
|---|---|
| `chunked_vcf_reader` | Stream a VCF in bounded-memory chunks |
| `compute_gene_burdens_chunked` | Gene-burden computation over chunked input |
| `parallel_pathway_scores` | Parallelised pathway scoring |
| `estimate_memory_usage` | Predict peak memory for a cohort before running |
| `downsample_cohort` | Reduce a cohort for smoke runs |
| `ProgressTracker` | tqdm-backed progress reporting |

---

## Exports

```python
__all__ = [
    "set_global_seed", "get_rng",
    "chunked_vcf_reader", "compute_gene_burdens_chunked",
    "parallel_pathway_scores", "estimate_memory_usage",
    "downsample_cohort", "ProgressTracker",
]
```

The metrics guards live at `pathway_subtyping.utils.metrics` and are imported from there directly:

```python
from pathway_subtyping.utils.metrics import safe_adjusted_rand_score, ari_with_validity
```

---

## See Also

- [`benchmark.md`](benchmark.md) and [`pipeline.md`](pipeline.md) — where the ARI guards are wired in
- [`../../KNOWN-ISSUES.md`](../../KNOWN-ISSUES.md) — the NB17 entry describing the original defect
