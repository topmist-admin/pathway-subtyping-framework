# Trajectory-validation Gate T — specification

> **Status: PROPOSED. Nothing in this document is implemented.** There is no
> `pathway_subtyping.trajectory` module, no `TrajectoryGateT` class, and no test
> file. This is a design spec written to the same standard as
> [`discreteness_gate.md`](discreteness_gate.md), which describes *shipped* code —
> do not cite this one as a capability. Target: v0.9 candidate, untriaged.
>
> **Prepared:** 2026-07-28 · **Prerequisite:** none of the code below exists yet;
> the reusable pieces it builds on (Gate 6, ICP, conformal) do.

## Why a reproducible trajectory is not a real one

This is the same mis-specification that produced Gate A, one dimension over.

A trajectory-inference method takes a sample-by-feature matrix and returns an
ordering — a pseudotime value per sample. The natural validation instinct is to
ask *"does the ordering recur under resampling?"*, and to compute something like
the mean rank correlation between orderings fitted on bootstrap resamples. That
observed statistic is fine. As with Gate A, the **reference distribution** decides
what "recurs more than chance" means, and the obvious references are the wrong
ones.

A single correlated Gaussian blob has **no trajectory**. It does, however, have a
dominant principal axis, and every trajectory method will happily project onto it
and return an ordering that is highly reproducible across resamples — because the
axis is a stable property of the covariance, not evidence of a process unfolding
in time. Rejecting a null of "no cross-feature dependence" tells you the data has
correlated structure. It does not tell you that structure is a trajectory.

Worse, and specific to this project's data: when the dominant axis is a
**confound**, the reproducible ordering is a confound gradient wearing a
trajectory's clothes. This is not hypothetical here. On GSE80655 a k=3 partition
passed every stability control at bootstrap ARI ≈ 0.92 while corresponding almost
exactly to **brain region** (Cramér's V ≈ 0.66, p ≈ 4e-26) and being independent
of **diagnosis** (p ≈ 0.41). A pseudotime axis fitted to that same matrix would
recover the region gradient, score extremely well on ordering persistence, and be
reported as a disease trajectory. The 141 samples also come from only **48
donors**, so any resampling that treats samples as independent inflates the
persistence statistic by counting pseudoreplicates.

> The ordering is reproducible. The question the gate must answer is what it is an
> ordering *of*.

## Observed statistic

**Trajectory persistence.** Fit the trajectory on *B* bootstrap resamples; for each
resample compute the Spearman rank correlation between its pseudotime and the
full-data pseudotime, restricted to shared samples; report the mean of the
**absolute** correlations.

The absolute value is not cosmetic. Pseudotime is identified only up to direction
unless a root is externally specified, so a resample that recovers the identical
ordering reversed is a success, not a failure. Averaging signed correlations would
score a perfectly stable trajectory near zero whenever the root flips.

This is the exact analog of Gate A's mean bootstrap ARI, and it inherits the same
property: **on its own it is uninformative**, because a Gaussian blob scores high.

## The trajectory null

> **What decides the verdict — and a deliberate departure from Gate A.**
>
> ```python
> passed = bool(
>     persistence > sg_p95        # T1 — exceeds the single-Gaussian reference
>     and sg_p < self.alpha       # T1 — add-one empirical p
>     and confound_max_v < self.v_max   # T2 — not a confound gradient
> )
> ```
>
> **Both T1 and T2 enter the rule.** This is stated explicitly because Gate A
> shipped with three named criteria of which only one decided — the gap statistic
> and dip test are computed, reported, and inert, and that gap between the
> documented framing and `passed = obs > sg_p95 and sg_p < alpha` was found by
> hostile review rather than by us. T3 and T4 below are diagnostics and are
> labelled as such; if a future revision wants them to decide, the rule above must
> change in the same commit as the prose.

### T1 — Single-Gaussian trajectory reference (primary, decides)

Fit one multivariate Gaussian to the dimension-reduced score matrix, draw *n_ref*
synthetic datasets from it, run the identical trajectory pipeline on each, and
compute each one's persistence. Pass iff observed persistence exceeds the
reference 95th percentile **and** the add-one empirical *p* < 0.05.

The logic is SigClust's, transposed: the null model of "no trajectory" is a single
Gaussian, and a genuine trajectory is a significant departure from what a
structureless blob's principal axis already gives you. Because a blob's ordering
*is* reproducible, this reference is typically much harder to beat than a
permutation null — which is the point.

*Liu Y, Hayes DN, Nobel A, Marron JS. Statistical Significance of Clustering for
High-Dimension, Low-Sample Size Data. JASA 103(483):1281–1293, 2008.
doi:10.1198/016214508000000454.*

### T2 — Confound attribution (decides, conjunctive)

Reuse the Gate 6 machinery in
[`validation.confound_association_gate`](../src/pathway_subtyping/validation.py)
against the inferred axis rather than against cluster labels:

- **categorical confounds** (brain region, batch, cohort, sex, platform) — bin
  pseudotime into deciles and compute bias-corrected Cramér's V against each
  confound, with BH-adjusted χ².
- **continuous confounds** (PMI, RIN, age at death, donor age) — Spearman ρ
  between the confound and pseudotime, BH-adjusted.

Fail if any confound is both statistically associated and non-trivially so
(default `v_max = 0.30`, matching Gate 6). As in Gate 6, **diagnosis is the
biology of interest, not a nuisance**: it is reported for context and never fails
the gate.

This is the term that would have caught the GSE80655 region gradient, and it is
the reason T2 decides rather than merely informing.

*Bergsma W. A bias-correction for Cramér's V and Tschuprow's T. J. Korean Stat.
Soc. 42(3):323–328, 2013. doi:10.1016/j.jkss.2012.10.002.*

### T3 — Invariance across environments (diagnostic, does not decide)

Where two or more environments exist (cohort, brain bank, platform, batch), test
whether the pseudotime→pathway-score relationship is invariant across them using
the existing
[`causal.InvariantPathwayPredictor`](../src/pathway_subtyping/causal/invariant.py)
and `invariance_pvalue`. A trajectory that reflects shared biology should hold its
shape across environments; one that reflects a site effect should not.

**Diagnostic only, for a reason worth stating.** Invariant causal prediction needs
genuinely multiple environments with overlap, and most postmortem cohorts supply
one. Making a criterion decide when it is usually inapplicable would convert
"could not test" into "failed," which is the abstention-counted-as-rejection error
this project already made once.

*Peters J, Bühlmann P, Meinshausen N. Causal inference by using invariant
prediction: identification and confidence intervals. JRSS B 78(5):947–1012, 2016.
doi:10.1111/rssb.12167.*

### T4 — Ordering-permutation null (diagnostic, does not decide)

Refit the trajectory with donor labels permuted across the ordering covariate
(age, illness duration) and recompute persistence. Reported as a sanity floor. It
does not decide because it is a weak null for the same reason the feature-permute
null was weak for Gate A: it destroys the covariate association while leaving the
covariance that generates the reproducible axis intact.

## Design hardening

Fixed in advance, applied **identically** to observed data and every reference draw:

| control | rule | rationale |
|---|---|---|
| **resampling unit** | bootstrap over **donors**, never samples | GSE80655 is 141 samples from 48 donors. Sample-level resampling treats pseudoreplicates as independent and inflates persistence. This is the single most important control in the table |
| minimum design | `n_donors >= 30` and `>= 3` distinct values on the ordering covariate, else **not-testable** | below this, persistence is dominated by resampling noise and no reference comparison is meaningful |
| dimensionality | PCA to *d* = min(10, ⌊n_donors/10⌋), capped at n_pathways and n_donors−1 | same HDLSS conditioning as Gate A; the Gaussian reference is simulated in this PCA space with diagonal covariance = retained per-PC variances |
| direction | mean **\|Spearman\|**; root, if any, chosen once on full data and held fixed across replicates | pseudotime is identified up to reversal; re-rooting per replicate destroys comparability |
| method fixed | the trajectory method and all hyperparameters are fixed across observed and every reference/bootstrap replicate | re-tuning per replicate lets the method find *some* axis in structureless data |
| direction stability routing | if no consistent orientation recurs in ≥ 50% of donor bootstraps → **not-testable** | analog of Gate A's k-stability routing; an unstable direction means there is nothing coherent to test |
| reporting | full reference distribution, 95% CI on observed persistence, per-confound V/ρ table, and the **donor count alongside every sample count** | the sample-vs-donor distinction must be impossible to lose downstream |

## Outcomes

Three, not two — matching Gate A, and for the same reason:

| verdict | meaning |
|---|---|
| `trajectory supported` | T1 and T2 both pass |
| `no trajectory support` | T1 fails — indistinguishable from a single-Gaussian axis |
| `confounded trajectory` | T1 passes but T2 fails — there is an axis, and it tracks a confound. Report **which** confound |
| `not-testable` | design minimum unmet, or direction unstable. **An abstention, not a rejection** |

Any false-positive rate computed from this gate must be quoted with its
**testable** denominator. Gate A abstained on 28/30 synthetic negatives and the
resulting "FPR 0.000" was misleading; Gate T will abstain at least as often,
because donor counts in postmortem psychiatric cohorts are small.

## What this gate does *not* do

It does not make a cross-sectional design longitudinal. If each donor contributes
one terminal snapshot, then recovering a within-individual trajectory from
between-individual variation requires an ergodicity-style assumption — that donors
at different ages trace the path one donor would follow — which is untestable in
that data and frequently false. Gate T can tell you an axis is not a
single-Gaussian artifact and not a measured confound. It **cannot** tell you the
axis is time. Passing it licenses "a reproducible axis exists that is not
attributable to the confounds we measured," and nothing stronger.

Designs with genuinely repeated measurement — longitudinal tumor biopsies
(pre-treatment / mid-therapy / relapse), organoid time courses — do not need this
caveat, because their time axis is observed rather than inferred. **Those are the
right first targets for both the gate and any trajectory work built on it.**

*Weinreb C, Wolock S, Tusi BK, Socolovsky M, Klein AM. Fundamental limits on
dynamic inference from single-cell snapshots. PNAS 115(10):E2467–E2476, 2018.
doi:10.1073/pnas.1714723115.*

## Synthetic validation (required before this ships)

`tests/test_trajectory_gate.py` must verify, on data with known ground truth, that
the gate:

1. **fails** a single correlated Gaussian (no trajectory, reproducible axis) —
   the case the naive persistence statistic passes;
2. **fails with `confounded trajectory`** on a matrix whose dominant axis is a
   synthetic batch/region label — the GSE80655 shape;
3. **passes** a genuine simulated trajectory with additive noise;
4. **abstains** below the donor minimum rather than returning a verdict;
5. returns **identical** results when samples are duplicated within donor —
   the pseudoreplication guard.

Test 5 is the one that most directly encodes the lesson from the flagship
reanalysis and should not be dropped for convenience.

## Proposed usage

```python
from pathway_subtyping.trajectory import TrajectoryGateT   # DOES NOT EXIST YET

gate_t = TrajectoryGateT(seed=42, n_ref=200, n_bootstrap=50, v_max=0.30)
res = gate_t.run(
    "gse80655",
    pathway_scores_df,
    pseudotime=pt,                  # from any trajectory method
    donor_ids=meta["donor_id"],     # REQUIRED — resampling unit
    confounds={"region": meta["region"], "batch": meta["batch"],
               "pmi": meta["pmi"], "rin": meta["rin"]},
    biology_of_interest={"diagnosis": meta["dx"]},
)
print(res.verdict)        # "trajectory supported" | "no trajectory support"
                          # | "confounded trajectory" | "not-testable (unstable direction)"
print(res.persistence, res.sg_empirical_p, res.attributed_confound)
print(res.n_donors, res.n_samples)   # both, always
```

## Implication for trajectory claims

Any pseudotime axis certified only by ordering reproducibility may be a
reproducibly-ordered blob, or a confound gradient, rather than a biological
trajectory. The failure transfers from Gate A intact: *stable is not discrete*
becomes *stable is not temporal*.

**Two limits to state whenever this gate is cited**, once it exists. It will
abstain frequently — donor counts, not sample counts, set the power, and 30 donors
is a real constraint on most postmortem cohorts. And like Gate A it is built to
refuse over-called structure, so it will be poor at detecting subtle real
trajectories. It is an instrument for saying no.

---

**References verified 2026-07-28** via the Crossref API (`api.crossref.org/works/<DOI>`),
per the practice used for the Monti 2003 reference in the cautionary manuscript.
All four resolve with matching author, journal, volume and year:

| DOI | Verified as |
|---|---|
| `10.1198/016214508000000454` | Liu et al., *JASA* 103, 2008 |
| `10.1016/j.jkss.2012.10.002` | Bergsma, *J. Korean Stat. Soc.* 42, 2013 |
| `10.1111/rssb.12167` | Peters et al., *JRSS B* 78, 2016 |
| `10.1073/pnas.1714723115` | Weinreb et al., *PNAS* 115, 2018 |
