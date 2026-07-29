# PSF Roadmap — F9 True Competition / Starvation Interaction Model

**Planned:** July 28, 2026
**Status:** **PROPOSED — not implemented.** No code in this document exists.
**Target release:** **v0.9.0** (next minor)
**Replaces:** the current `CrosstalkDetector._compute_competition`, soft-deprecated 2026-07-28
**Public/private posture:** Public Codeberg track, MIT, alongside the rest of the framework.

> ⚠️ Do not cite this as a capability. `pathway_subtyping.qc.crosstalk` currently
> ships the **broken** score described below, guarded by a `FutureWarning` and
> withheld from `qc.__all__`. This document specifies what should replace it.
>
> **Disambiguation:** "F9" is overloaded in this project. This concerns **F9 of the molecular QC layer** (`qc.crosstalk.CrosstalkDetector`). It is unrelated to **F9 of the v0.6 release** (`qc.offtarget_sequence`, Evo 2 CRISPR off-target scoring), which is unaffected.

---

## Why this note exists

F9's module docstring states the intent plainly: *"When multiple pathways share
intermediate transducers, one pathway can starve another."* The shipped score
does not test that claim, or anything near it.

`_compute_competition` subtracts the shared gene set from both pathways and
returns `mean(A-exclusive) − mean(B-exclusive)`. The shared genes never enter the
arithmetic — they act only as a presence gate. Demonstrated 2026-07-28:

- varying shared-node expression across four orders of magnitude (0.1 → 500)
  leaves the score **bit-identical** and `dominant` unchanged;
- holding shared nodes perfectly balanced and varying only the exclusive genes
  swings the score **+8.0 → −8.0** and flips `dominant`.

So the score is blind to the quantity it names and fully determined by genes that
are, by construction, not shared. A difference of means between two pathways'
private machinery is a statement about **relative pathway activity**, not about
competition for a shared resource. No amount of rescaling turns one into the other:
starvation is a claim about a *relationship between* A and B that *depends on* S,
and a difference of two means encodes no relationship and no dependence.

**One thing the current code gets right and the replacement must keep.** The
`- set(shared)` subtraction is correct and load-bearing. If shared genes were left
inside `E_A` and `E_B`, both would be mechanically correlated with `S`, and every
statistic below would be confounded by construction. The exclusion is not the bug.

---

## The model

For each observation *i* (cell or sample — see [Unit of analysis](#unit-of-analysis)):

| Symbol | Definition |
|---|---|
| `E_A(i)` | mean expression over genes **exclusive** to pathway A (A minus shared) |
| `E_B(i)` | mean expression over genes **exclusive** to pathway B (B minus shared) |
| `S(i)` | mean expression over the **shared** transducer genes |

### Stage 1 — Partial correlation (screen)

Competition implies A and B trade off once the common driver is removed:

```
ρ_AB·S  =  (ρ_AB − ρ_AS · ρ_BS) / sqrt( (1 − ρ_AS²)(1 − ρ_BS²) )
```

Screen criterion: **ρ_AB·S significantly negative** (BH-adjusted across pathway pairs).

Controlling for `S` is what makes this non-trivial. If the shared transducer drives
both pathways, the *raw* `ρ_AB` is positive and hides any competition underneath;
partialling `S` out removes that common-driver component. Use Spearman by default —
expression is heavy-tailed and the linearity assumption behind Pearson is not
safe here.

**Stage 1 alone is not sufficient**, and the spec should say so loudly, because it
is the tempting place to stop. A negative partial correlation is consistent with
competition but also with any unmeasured factor that pushes A up while pushing B
down. It does not test that `S` is the *mechanism*.

### Stage 2 — Moderation by shared-transducer level (confirmatory)

The starvation hypothesis is specifically that competition **appears when the
shared transducer is limiting and disappears when it is abundant**. That is an
interaction, and it is the actual test:

```
E_B  =  β₀ + β₁·Ẽ_A + β₂·S̃ + β₃·(Ẽ_A × S̃) + ε
```

with `Ẽ_A` and `S̃` mean-centred, so `β₁` is the A→B slope at *typical* `S` rather
than at the meaningless `S = 0`.

Starvation predicts **both**:

- **β₁ < 0** — at typical transducer levels, more A goes with less B;
- **β₃ > 0** — that negative slope *attenuates* as `S` rises, because an abundant
  transducer is no longer the binding constraint.

The conditional slope is `∂E_B/∂E_A = β₁ + β₃·S̃`. Report it evaluated at low,
median and high `S` (e.g. 10th / 50th / 90th percentile) rather than reporting
`β₃` alone — the sign flip across the range is the interpretable result and the
one a reviewer will want.

`β₃ > 0` with `β₁ ≥ 0` is **not** competition; it is some other moderation and must
not be reported as starvation. Requiring both signs is what stops the model from
finding "competition" in every dataset with a significant interaction term.

### Verdict

Emit competition only when Stage 1 **and** Stage 2 pass, and report a four-way
outcome consistent with Gates A, K and T:

| Verdict | Meaning |
|---|---|
| `competition` | ρ_AB·S < 0 significant, β₁ < 0, β₃ > 0 |
| `no-competition` | tested, criteria not met |
| `common-driver` | raw ρ_AB > 0 but ρ_AB·S ≈ 0 — S drives both; explicitly *not* competition |
| `not-testable` | too few observations, degenerate variance, `S` has no spread, or fewer than the minimum shared genes |

`not-testable` is an **abstention, not a rejection**. Carry a `testable` flag on the
result object, exactly as `KGSensitivityResult` does, so any rate computed over many
pathway pairs uses the testable subset as its denominator. This project has already
published one false-positive rate against the wrong denominator; the type should
make that mistake hard to repeat.

---

## Failure modes the implementation must handle

These are ranked by how likely they are to manufacture a false competition call.
The first two are the reason this is a rewrite rather than a patch.

### 1. Compositional normalisation manufactures anticorrelation ⚠️ highest risk

If expression is normalised to a fixed total per observation (CPM, TPM, relative
abundance), the features are compositionally constrained: one gene rising forces
others down **arithmetically**, with no biology involved. That artifact looks
exactly like starvation and will produce negative `ρ_AB·S` on data with no
competition whatsoever.

This is the same species of error as the artifacts behind this project's 2026-07
correction — a normalisation choice generating the very signal being claimed. It
must be handled explicitly, not assumed away:

- prefer a log-ratio / centred-log-ratio treatment, or an explicitly
  compositional association measure, over raw correlation on closed data;
- **require** a compositional negative control in the test suite (below);
- record the normalisation of the input matrix in the result object, and abstain
  (or at minimum warn) when it is closed and untransformed.

*Aitchison J. The Statistical Analysis of Compositional Data. JRSS B 44(2):139–177, 1982.
doi:10.1111/j.2517-6161.1982.tb01195.x. · Lovell D, Pawlowsky-Glahn V, Egozcue JJ, Marguerat S, Bähler J.
Proportionality: a valid alternative to correlation for relative data. PLoS Comput Biol
11(3):e1004075, 2015. doi:10.1371/journal.pcbi.1004075.*

### 2. Library size / sequencing depth induces positive correlation

Deeper-sequenced cells have higher counts for essentially every gene, inducing
positive correlation across all pairs and *masking* real competition. Normalise
per observation before computing any of the three summaries, and consider
including depth as an additional covariate in Stage 2.

### 3. `S` with no spread

If the shared transducer is uniformly expressed across observations, the
interaction term is unidentifiable — there is no variation in the moderator. Abstain
(`not-testable`) rather than fitting a model whose key coefficient has no support.

### 4. Averaging over heterogeneous shared genes

`S` as a plain mean assumes the shared genes behave as one transducer pool. If they
split into sub-modules with opposite behaviour, the mean cancels them. Worth
reporting the shared-gene set's internal coherence (e.g. mean pairwise correlation)
and abstaining or warning when it is low.

### 5. Multiple testing

All pathway pairs are tested. BH-adjust across pairs and report both raw and
adjusted p, as the existing gates do.

---

## Unit of analysis

Open decision, and it changes the statistics:

- **Per cell** (scRNA-seq) gives large *n* and real power for the interaction term,
  but dropout and technical zeros dominate the correlation structure, and cells
  from one donor are not independent.
- **Per sample / pseudobulk** is far cleaner but leaves *n* in the tens, which is
  thin for a four-parameter interaction model.

Whichever is chosen, **non-independence must be respected** — cells nested in
donors, samples nested in batches. Treating pseudoreplicates as independent is the
specific error that inflated significance in this project's flagship reanalysis
(141 samples from 48 donors). If the per-cell route is taken, the model needs a
donor-level random effect or donor-level bootstrap, not plain OLS.

**Recommendation:** implement per-sample first with an explicit minimum *n*, and
treat per-cell with mixed effects as a follow-on.

---

## Validation requirements

The current F9 tests could not have caught the existing defect: all three draw every
gene from the same distribution (`normal(5.0, 0.5)`), under which the score is ~0 for
*any* formula, and none asserts on `competition_score` at all. The replacement must
ship with tests that **discriminate**, on constructed ground truth:

| # | Fixture | Required outcome |
|---|---|---|
| 1 | Genuine competition: a fixed shared pool allocated between A and B | `competition` |
| 2 | A and B independent, `S` varying freely | `no-competition` |
| 3 | Common driver: `S` raises both A and B | `common-driver`, **never** `competition` |
| 4 | **Compositional artifact:** independent A and B, then CPM-normalised | **not** `competition` |
| 5 | Competition present only at low `S`, absent at high `S` | `competition`, and the reported conditional slope must flip sign across the `S` range |
| 6 | `S` constant across observations | `not-testable` |
| 7 | Below the minimum-*n* threshold | `not-testable`, with `testable=False` |

Test 4 is the one that matters most — it is the failure mode most likely to reach a
published claim, and the one no current test approaches.

Additionally, the inverse of the deprecation test: the score **must** respond to
shared-node expression. `TestCrosstalkSoftDeprecation::test_score_is_insensitive_to_shared_nodes`
currently pins the defect and is designed to fail once F9 is fixed — that failure is
the signal to proceed to the acceptance criteria below.

---

## Acceptance criteria — what un-deprecates F9

All must hold before `CrosstalkDetector`, `CrosstalkResult` and `InterferenceEdge`
are restored to `pathway_subtyping.qc.__all__` and the `FutureWarning` removed:

1. All seven validation fixtures pass.
2. The score demonstrably responds to shared-node expression (deprecation test
   inverted and deleted).
3. Four-way verdict implemented, with `testable` on the result object.
4. `dominant`, if retained at all, means something checkable — or is removed. Its
   current definition ("which pathway dominates shared nodes") has never been
   computed, and it currently names a winner off pure sampling noise because it
   tests `score > 0` with no dead zone.
5. `n_significant` / PASS–FAIL re-derived from the new verdict, not from
   `abs(score) > 0.3`.
6. Normalisation of the input matrix recorded in the result, with abstention or a
   warning on closed untransformed data.
7. Module docstring, `docs/api/qc.md` and the CHANGELOG updated in the **same
   commit** as the code — the current mismatch between comment and computation is
   exactly what this whole exercise is correcting.

---

## Scope note

Nothing here affects `KnowledgeGraph.get_pathway_crosstalk()` or
`get_shared_genes()`, which are unrelated topology helpers that happen to share the
word "crosstalk".

No published result depends on F9: `CrosstalkDetector` is invoked only in its own
unit tests, and no deposited artifact contains `competition_score`. This is a
forward-looking fix, not a correction to anything on the record.

---

**References verified 2026-07-28** via the Crossref API
(`api.crossref.org/works/<DOI>`), per the practice used for the Monti 2003
reference in the cautionary manuscript:

| DOI | Verified as |
|---|---|
| `10.1111/j.2517-6161.1982.tb01195.x` | Aitchison, *JRSS B* 44, 1982 — "The Statistical Analysis of Compositional Data" |
| `10.1371/journal.pcbi.1004075` | Lovell et al., *PLOS Comput Biol* 11, 2015 — "Proportionality: A Valid Alternative to Correlation for Relative Data" |

The DOI originally written for Aitchison from memory (`10.2307/2345821`) does
**not** resolve in Crossref; the JRSS-B DOI above is the one that does.
