# Result 2 — empirical audit of the corrected 47-dataset benchmark

**Claim this package supports:** across independent public transcriptomic cohorts,
discrete-subtype partitions almost never clear a reproducibility bar — and the
silhouette-calibrated "adaptive threshold" that the withdrawn manuscript proposed
to set that bar does not hold at all.

Input is the **corrected** benchmark only: `corrected_benchmark_47datasets_v2.csv`,
released with the [2026-07-08 erratum](../../../CORRECTION_2026-07/ERRATUM_2026-07-08.md),
Zenodo **10.5281/zenodo.21262112** (v2.0). The pre-correction v1.0/v1.1 records are
superseded and must not be used.

## Reproduce

```bash
python scripts/audit_corrected_benchmark.py     # deterministic, no network
```

Reads the deposited CSV, writes `results/benchmark_audit.json`.

## 2a. The retracted threshold model does not hold under any row screen

The withdrawn manuscript reported an adaptive bootstrap threshold calibrated on
these 47 benchmarks at **R² = 0.889, slope +0.914, LOOCV RMSE 0.051**. Refitting
the same regression (bootstrap ARI ~ silhouette) on the released data:

| Row screen | n | R² | slope | p |
|---|---|---|---|---|
| all rows | 37 | 0.111 | −0.599 | 0.044 |
| valid rows | 22 | **0.015** | −0.037 | 0.589 |
| valid + stricter ground-truth screen (§2c) | 15 | **0.001** | +0.011 | 0.918 |

Silhouette does not predict reproducibility in these data. The slope even reverses
sign depending on the screen, which is what a null relationship looks like. This
reproduces the erratum's refit exactly (R² 0.0149, slope −0.0369 on valid rows).

## 2b. Reproducibility across cohorts is uniformly low — the Result 2 headline

Distribution of the benchmark's reproducibility statistic across valid cohorts:

| Screen | n | median | max | ≥0.8 | ≥0.5 | ≥0.2 |
|---|---|---|---|---|---|---|
| valid | 22 | −0.002 | 0.391 | **0/22** | **0/22** | 1/22 |
| valid + stricter (§2c) | 15 | 0.000 | 0.391 | **0/15** | **0/15** | 1/15 |

**No cohort clears 0.5.** A single cohort (GSE45291, immunology, n=805) exceeds
0.2, at 0.391. The median sits at zero. This holds under both screens.

The honest reading: the withdrawn manuscript's own benchmark, once its degenerate
rows are removed, does not contain a single public cohort on which a discrete
partition was reproducible by the framework's own criterion. That is a real
negative result about the field's central practice, and it is the empirical
backbone of the cautionary thesis — not a defect specific to this framework.

It also disposes of the adaptive-threshold idea from a second direction: there is
no spread of reproducibility values to calibrate a threshold *against*.

## 2c. The correction was itself incomplete (found by applying our own standard)

The erratum flagged "impossible labels" using the rule **`n_true_clusters >
n_samples`**. That catches GSE92332 (1,957 labels for 533 samples) and nothing
else. It does not catch rows where the label count *equals or approaches* the
sample count — one ground-truth cluster per sample, which is equally unusable as a
ground truth but fails a strict `>` test.

**Seven of the 22 rows still marked `valid=True` fail a stricter screen**
(`n_true_clusters / n_samples ≥ 0.5`):

| Dataset | n samples | n "true clusters" | ratio |
|---|---|---|---|
| GSE2109 | 2158 | 2158 | 1.00 |
| GSE5204 | 79 | 79 | 1.00 |
| GSE17537 | 55 | 54 | 0.98 |
| GSE44228 | 72 | 69 | 0.96 |
| GSE42127 | 176 | 145 | 0.82 |
| GSE66351 | 190 | 96 | 0.51 |
| GSE16759 | 16 | 8 | 0.50 |

GSE2109 and GSE5204 have *exactly one label per sample*. They sit one row away
from GSE92332 in kind, and were retained only because the erratum's rule used `>`
rather than `≥` or a ratio.

**This is disclosed, not buried.** A paper arguing that the field under-audits its
own subtype claims has no standing to skip auditing its own correction. The
material point is that **the conclusion is robust to the fix**: every headline
above is reported under both screens, and none of them move (0/22 → 0/15 at every
threshold ≥ 0.5; R² falls further, 0.015 → 0.001).

## ⚠️ Interpretive caveat — read before quoting a number

The benchmark's reproducibility column is `bootstrap_ari_5th_percentile`: the
**5th percentile** of a bootstrap ARI distribution, i.e. a deliberately
conservative lower bound, not a central estimate. Two consequences:

1. A near-zero 5th percentile does **not** license the claim "median
   reproducibility is zero." The deposited CSV carries only the 5th percentile, so
   the central tendency of each cohort's bootstrap distribution is not recoverable
   from this artifact. Do not write "subtypes are never reproducible."
2. The defensible claim is about the **gate criterion**: the framework compares
   this statistic to a threshold, and under that comparison no valid cohort passes.

There is also an unresolved question about the column's exact construction. The
erratum states that rows with `n_true_clusters = 1` make `adjusted_rand_score`
"return an arithmetically meaningless value rather than measured stability," which
implies the ground-truth labels entered the computation somewhere — yet a pure
bootstrap-stability statistic should not depend on ground truth at all. The
original benchmark code that produced the column was rewritten before the
correction, so this could not be settled from the repository. **Flag this as a
limitation in the manuscript rather than asserting a mechanism.** It does not
affect 2a (the refit is a property of the released numbers whatever they measure)
and it is the reason 2b is phrased as a statement about the gate criterion.

## What this package does NOT claim

- Not that these cohorts contain no biological structure. Reproducibility of a
  *partition* and existence of *structure* are different questions — that is the
  entire point of the discreteness gate (see `../gate_calibration/`, where real
  cohorts with certified discrete structure still fail the fixed stability bar).
- Not that the framework is uniquely bad. The benchmark measures a general
  practice — cluster, then check whether the partition survives resampling — and
  the finding is about that practice.
- Nothing about the pre-correction v1.0/v1.1 numbers, which are retracted.
