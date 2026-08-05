# Gate ablation + clusterer-generality sweep (reviewer R3.10, R2.2)

Synthetic-ground-truth validation of the gate battery: what does each gate
actually *buy*, and what does the clusterer sweep actually establish? (Answer: only that the gate's interface accepts any clusterer — see the retraction below.) This is the
**Result 1 (methods validation)** evidence for the cautionary-framework paper.

## ⚠️ READ FIRST — the honest re-analysis supersedes the headline numbers below

A hostile-review pass showed the original ablation write-up was misleading in three
ways, all confirmed against the raw CSV. **The corrected numbers live in
`results/ablation_honest.json`** (from `scripts/reanalyze_ablation_honest.py`); use
those, not the "FPR 0.000 / rejected 100%" framing in the older sections:

1. **Three-way accounting.** The gate returns certify / reject / **not-testable**.
   On the 30 negative-control datasets it returned `not-testable` (an abstention) on
   **26 (87%)** — it did not *judge* them as continua, it declined. The old write-up
   scored abstentions as correct rejections, manufacturing FPR=0.000. Excluding
   abstentions the testable-negative denominator is **4** (Wilson 95% CI on that FPR:
   [0.00, 0.49] — nearly uninformative).
2. **The real, defensible result is the paired contrast.** The stability-only gate
   falsely certified **13/30** negatives (FPR 0.433, 95% CI [0.27, 0.61]); the
   recalibrated gate certified **0/30**. Exact McNemar on the negatives: b=13, c=0,
   **p=0.0002** — the reduction in false certification is significant. The mechanism is
   partly abstention, and we say so, but the improvement over stability-only is real.
3. **No measurable TPR cost on this benchmark — and that is not the same as none.**
   TPR is **1.000 for both** gate sets: the two classify all 30 genuine-structure
   datasets identically, so there are **0 discordant pairs** and exact McNemar gives
   **p=1.0**. With 30 positives the study still cannot exclude a small cost, so write
   "no cost was measurable here", not "the gate is free". ⚠️ Earlier revisions of this
   file reported a **1.000 → 0.967** cost from one lost `discrete_k3` replicate; that
   came from a run whose datasets were drawn under a per-process-salted `hash()` seed
   and could not be regenerated. It is retired, not merely updated.

**And the increment is a null recalibration, not a new instrument.** Head-to-head in
`ablation_honest.json`: the SigClust-style single-Gaussian p-value **alone**
(threshold α=0.05) reproduces the full gate's TPR and FPR exactly. The gap statistic
and Hartigan dip are computed but never enter the decision rule
(`passed = obs > sg_p95 and sg_p < alpha`); the k-stability route only adds
abstentions. Frame the contribution as replacing the feature-permutation null with a
single-Gaussian null (SigClust; Liu et al. 2008) on the existing bootstrap-ARI
statistic — not as a multi-test discreteness instrument.

**Resolution / operating characteristic** (`separation_sweep.py` /
`results/separation_sweep.json`, 20 reps/step, 7 δ points): sweeping component
separation δ from 0 to 3 SD (n=120, n_ref=100), the gate **does** resolve — certify
rate 0.00 through δ≤2.0, then **0.15 at δ=2.5** and **0.55 at δ=3.0**, while median
single-Gaussian p descends 0.78 → 0.23 → 0.04 — refuting "flat at the floor". But the
transition is **sharp and late** (confined to δ=2.5–3.0), so the recalibrated null is
**conservative** (low sensitivity to moderate structure), abstaining in the mid-range.
Good for rejecting over-called subtypes; not a sensitive detector of subtle real
structure. Reported as a limitation, not tuned away. (The script default is now
`--reps 20`, so a re-run reproduces the deposited table exactly.)

The scripts live in the framework repo root (they are framework tooling, not
package-local): `scripts/gate_ablation_study.py` and `scripts/plot_gate_ablation.py`.
The deposited outputs are here because the scripts' default `--out`
(`outputs/gate_ablation/`) is **gitignored** — these artifacts would otherwise
not be version-controlled at all.

## Reproduce

```bash
# from the framework repo root; deterministic (seed 42)
python scripts/gate_ablation_study.py --out consolidation-cautionary/cross-domain/gate_ablation/results
python scripts/plot_gate_ablation.py  --out consolidation-cautionary/cross-domain/gate_ablation/results
```

⚠️ **The full run took 4,623.9 s (~77 min)** on the deposited settings
(n=60, p=50, 15 reps/condition, n_ref=120, n_bootstrap=40). Use `--quick` for a
smoke run — `--quick` numbers will NOT match the deposited values.

Requires the **v0.8 line** of the framework (`pathway_subtyping.discreteness`,
`pathway_subtyping.clustering_dl`). v0.8.0 is released on PyPI:
`pip install pathway-subtyping==0.8.0`. Neither module exists in v0.7.0, so an
older pin cannot run this package. See [`../../RUNME.md`](../../RUNME.md).

## Result 1a — the ablation (R3.10)

R3.10 asked: *"What happens to the performance if some of the validation gates
are removed?"* Four synthetic conditions with known ground truth — two positive
(`discrete_k2`, `discrete_k3`: genuine well-separated components) and two
negative (`single_gaussian`: one correlated blob; `continuum_1d`: a 1-D gradient,
the tumor-purity / immune-infiltration analog).

| Gate subset | TPR | FPR † | `continuum_1d` cert rate |
|---|---|---|---|
| `stability_only` (pre-v0.8.0) | 1.000 | **0.433** | **0.867** |
| `discreteness_only` | 1.000 | 0.000 † | 0.000 |
| `stability+discreteness` (v0.8.0) | 1.000 | **0.000** † | 0.000 |

† **The 0.000 cells are computed against n=30 and are the over-generous
convention.** The gate abstains (`not-testable`) on 26 of those 30 negatives, so
the testable-negative denominator is **4** and the honest Wilson interval is
[0.00, 0.49] — nearly uninformative. Never quote "FPR 0.000" from this table
without that denominator. See the READ FIRST section above.

**Reading.** The pre-v0.8.0 stability-only gate false-certifies the 1-D continuum
in 87% of runs — this is the erratum finding reproduced under controlled
conditions: the old bootstrap null tested *pathway independence*, not
*discreteness*, so a reproducibly-bisected gradient passed. Adding Gate A removes
those false certifications — predominantly by declining to rule (abstaining) on
the negatives rather than by rejecting them.

**Honest note on TPR:** the ablation leaves true positives untouched on this
benchmark — TPR is **1.000 for both** gate sets, with **0 discordant pairs** and
exact McNemar **p=1.0**. That is an absence of *measurable* cost, not a proof of
none: 30 positives cannot exclude a small one. Write "no measurable cost", never
"the gate is free". (See item 3 in the READ FIRST section for the retired
1.000 → 0.967 figure and why it was withdrawn rather than updated.)

Source: `results/gate_ablation_results.json` (+ `_raw.csv` for per-replicate rows).
Figure: `results/gate_ablation_figure.{png,svg}`.

## Result 1b — clusterer generality (R2.2)

The same synthetic data clustered three ways — GMM, DEC (Xie 2016), VAE-GMM
(VaDE, Jiang 2017) — then passed to Gate A.

**Primary result:** Gate A **declines to certify** the 1-D continuum on **100%** of
runs, and does so identically regardless of which algorithm drew the partition.

⚠️ **"Declines to certify" is not "rejects."** The bulk of those non-certifications
are `not-testable` abstentions — the gate found no reproducible k and so declined to
rule, rather than ruling the data continuous. The earlier wording here ("rejects the
1-D continuum on 100% of runs") is the framing the READ FIRST section retracts; the
per-run breakdown is in `results/clusterer_sweep_results.md`, which labels each run
CERTIFY / REJECT / ABSTAIN.

⚠️ **The clusterer-agnostic part of the claim is ALSO retracted as an empirical finding.** The gate is called once per condition and its verdict copied into all three clusterer rows, so the arms are identical by construction, not by agreement. What survives is a statement about the **interface**: `DiscretenessGateA.run()` takes no partition argument and re-clusters internally, so it can be applied downstream of any clusterer — the framework is a validation layer rather than another clustering method to benchmark against DEC/VAE. That is definitional, not experimental, and must not be reported as a result.

**Secondary observation (reported honestly):** on the continuum, mean self-stability
ARI is GMM 0.88 vs DEC 0.34 vs VAE-GMM 0.22. GMM's near-0.8 self-stability is
exactly the classic continuum false positive Gate A was built to catch; the deep
methods are *additionally* unstable under resampling at this small n — itself a
caution against trusting algorithmic confidence (consistent with R3.5/R3.6/R3.9).
Competitive DL *recovery* of real subtypes is a separate claim belonging to the
large cohorts (see `../cancer_r38/`), not this small-n control.

Source: `results/clusterer_sweep_results.md` (+ `_raw.csv`).

## Scope limits (do not overstate)

- **Synthetic only.** These are controlled conditions with known ground truth.
  The real-data counterpart — does the gate certify genuine structure and reject a
  genuine continuum on public cohorts? — is `../gate_calibration/`, and that is
  the evidence that answers "is the gate just pessimistic?".
- **n=60 per dataset.** Deliberately small, to probe the small-n regime where the
  old null fails. Behaviour at n≈1000 is characterized on real data in
  `../cancer_r38/` (large cohorts, where Gate A is well-powered).
- DL baselines run at standard settings, not cohort-tuned.
