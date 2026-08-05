# Discreteness-aware Gate A (v0.8.0)

## Why the previous stability gate was mis-specified

The stability gate compares an observed statistic — the mean bootstrap adjusted
Rand index (ARI) of a Gaussian-mixture partition at a fixed number of clusters *k*
— against a reference distribution. The observed statistic asks *"does the
partition recur under resampling?"* The **reference** decides what "recurs more
than chance" means, and the previous reference was the wrong one.

The old null permuted **each pathway-score column independently across samples**.
That destroys inter-pathway correlation while preserving every pathway's marginal
distribution, the sample size, and the dimensionality — so it is a null model of
*"the pathways are mutually independent."* Rejecting it answers *"is there
cross-pathway dependence?"*, which is **not** the same as *"are there discrete
clusters?"*

The two questions diverge exactly where it matters at small *n*. A single
correlated Gaussian blob, or a one-dimensional continuous gradient (tumor purity,
immune infiltration, proliferation), has **no discrete clusters**. Yet a mixture
model fit to such data reproducibly bisects it along the dominant axis on every
bootstrap resample — a high, stable observed ARI. The feature-permuted null,
having destroyed the correlation that produced that axis, collapses to a low ARI.
Observed ≫ null, and a smooth continuum is certified a "reproducible subtype." The
failure is structural, not stochastic: it recurs whenever the data is a correlated
blob, which is the typical shape of small-*n* expression scores.

> The observed statistic is fine; the reference distribution was wrong.

## The discreteness null

`DiscretenessGateA` keeps the identical observed statistic and changes what it is
compared against, adding two corroborating diagnostics.

> **What actually decides the verdict.** Only item 1 below. The decision rule is
> literally
> ```python
> passed = bool(obs > sg_p95 and sg_p < self.alpha)
> ```
> (`src/pathway_subtyping/discreteness/gate_a_discreteness_null.py:400`) — the
> single-Gaussian reference and nothing else. **The gap statistic (2) and the dip
> test (3) are computed and reported, but never enter the decision.** On the
> deposited ablation, thresholding the SigClust *p* alone reproduces the full
> gate's TPR and FPR exactly. Treat 2 and 3 as diagnostics you read alongside the
> verdict, not as criteria the data must also clear.

1. **Single-Gaussian reference (SigClust) — primary.** Fit **one** multivariate
   Gaussian to the (dimension-reduced) score matrix and ask whether the observed
   partition is more stable than partitions of data drawn from that single
   Gaussian. This is the SigClust principle: the null model of "no clusters" is a
   single Gaussian, and a genuine cluster is a statistically significant departure
   from it. Pass iff observed ARI exceeds the reference 95th percentile **and** the
   add-one empirical *p* < 0.05.
   *Liu Y, Hayes DN, Nobel A, Marron JS. Statistical Significance of Clustering for
   High-Dimension, Low-Sample Size Data. JASA 103(483):1281–1293, 2008.
   doi:10.1198/016214508000000454.*

2. **Gap statistic — complementary.** Compares log within-cluster dispersion to its
   expectation under a uniform reference over the PCA-aligned bounding box, and
   returns the gap-optimal *k* by the first-standard-error rule. A continuum yields
   no positive peaked gap at *k* ≥ 2.
   *Tibshirani R, Walther G, Hastie T. Estimating the Number of Clusters in a Data
   Set Via the Gap Statistic. JRSS B 63(2):411–423, 2001. doi:10.1111/1467-9868.00293.*

3. **Hartigan's dip test — corroborating, optional.** Tests unimodality of the
   score distribution on PC1 and on the mixture discriminant axis. A non-significant
   dip is evidence of a unimodal, non-discrete distribution. Requires the `diptest`
   extra (`pip install pathway-subtyping[discreteness]`); skipped gracefully if
   absent. *Hartigan JA, Hartigan PM. The Dip Test of Unimodality. Ann. Statist. 13,
   1985. doi:10.1214/aos/1176346577.*

The old independence (feature-permute) null is **retained and reported**, demoted
to a confound / marginal control: it now answers a secondary question — does the
structure depend on genuine cross-pathway dependence, or on a few marginal
distributions? — and no longer decides Gate A.

## Small-*n* hardening

Fixed in advance and applied **identically** to observed data and every reference
draw:

| control | rule | rationale |
|---|---|---|
| dimensionality | reduce to *d* = min(10, ⌊n/10⌋) PCA components (capped at n_pathways and n−1) | a 50-dim full-covariance mixture is singular at n in the tens; PCA to *d* ≪ n conditions it and exits the HDLSS regime. The single Gaussian is simulated in this PCA space with diagonal covariance = retained per-PC variances (PCs are uncorrelated by construction → faithful null) |
| fixed *k* | *k* = the selected value, identical for observed and every reference/gap replicate | *k* is never re-selected per replicate — structureless data would collapse to *k* = 1 and the ARI would degenerate |
| k-stability routing | per observed bootstrap, record the silhouette-optimal *k*; if no modal *k* recurs in ≥ 50% of resamples → **not-testable** | full-covariance BIC over-selects *k* in the reduced space; silhouette cleanly separates real clusters from continua |
| reporting | full reference distribution + a 95% CI on the observed statistic | transparency and old-vs-new comparison |

## Synthetic validation

`tests/test_discreteness_gate.py` verifies three properties on data with known
ground truth: the gate **fails** a single Gaussian and a 1-D continuum (the case
the old independence null wrongly passed) and **passes** two separated clusters.

## Gate C, reframed

Conformal prediction guarantees coverage of whatever label it is given; here that
label is the mixture assignment, so "coverage of the assignment" is guaranteed by
construction and cannot, on its own, distinguish a real subtype from a sliced
continuum. `ReframedMembershipGate` therefore reports **assignment sharpness
conditional on a partition already certified by Gates A and B** — at a single 0.90
operating point, with a bootstrap coverage CI, the singleton-set (sharpness)
fraction, and a minimum-n rule.

## Usage

```python
from pathway_subtyping.discreteness import DiscretenessGateA, ReframedMembershipGate

gate_a = DiscretenessGateA(seed=42, n_ref=200, n_bootstrap=50)
res = gate_a.run("my_cohort", pathway_scores_df, n_clusters=k)
print(res.verdict)          # "discrete structure" | "no discrete structure (continuum)" | "not-testable (no reproducible k)"
print(res.passed, res.sg_empirical_p, res.gap_optimal_k)

gate_c = ReframedMembershipGate(seed=42)
c = gate_c.run("my_cohort", pathway_scores_df, cluster_labels, k=k,
               gate_a_pass=res.passed, gate_b_pass=b_pass)
print(c.achieved_coverage, c.coverage_ci95, c.sharpness_confident_fraction)
```

## Implication for prior "reproducible subtype" claims

Any subtype previously certified only by the bootstrap-stability gate may be a
reproducibly-sliced continuum rather than a discrete subtype. A "reproducible
subtype" claim should now clear the discreteness gate — the **single-Gaussian
(SigClust) reference**, which is the whole decision rule — rather than the
independence null alone. Gap and dip are reported alongside as diagnostics; they
do not gate the verdict.

**Two limits to state whenever the gate is cited.** It frequently returns
`not-testable` rather than a ruling — on the synthetic negative controls it
abstained on 26/30 — so a low false-positive rate is largely abstention, not
rejection, and must be quoted with its testable denominator. And it is
conservative: on a separation sweep it certifies nothing below δ=2.5 SD, so it is
good for refusing over-called subtypes and poor at detecting subtle real ones.

## Related gates

This gate varies the **reference distribution** while holding the data and the
knowledge graph fixed. Two siblings vary something else:

- **[Gate K — time-sliced KG sensitivity](kg_sensitivity_gate.md)** *(implemented)*
  varies the **knowledge graph**: does the finding survive being recomputed
  against a different database release? Every gate here holds the KG fixed, so
  this is a class of artifact the rest of the battery cannot see. Its null design
  follows the lesson on this page — a raw agreement number between two KG
  versions is as uninterpretable as a raw bootstrap ARI, for the same reason.
- **[Gate T — trajectory validation](roadmap-trajectory-gate.md)** *(proposed, not
  implemented)* transfers the argument to pseudotime: *stable is not discrete*
  becomes *stable is not temporal*.
