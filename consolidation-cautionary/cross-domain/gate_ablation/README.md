# Gate ablation + clusterer-generality sweep (reviewer R3.10, R2.2)

Synthetic-ground-truth validation of the gate battery: what does each gate
actually *buy*, and is the discreteness gate clusterer-agnostic? This is the
**Result 1 (methods validation)** evidence for the cautionary-framework paper.

## ⚠️ READ FIRST — the honest re-analysis supersedes the headline numbers below

A hostile-review pass showed the original ablation write-up was misleading in three
ways, all confirmed against the raw CSV. **The corrected numbers live in
`results/ablation_honest.json`** (from `scripts/reanalyze_ablation_honest.py`); use
those, not the "FPR 0.000 / rejected 100%" framing in the older sections:

1. **Three-way accounting.** The gate returns certify / reject / **not-testable**.
   On the 30 negative-control datasets it returned `not-testable` (an abstention) on
   **28 (93%)** — it did not *judge* them as continua, it declined. The old write-up
   scored abstentions as correct rejections, manufacturing FPR=0.000. Excluding
   abstentions the testable-negative denominator is **2** (Wilson 95% CI on that FPR:
   [0.00, 0.66] — nearly uninformative).
2. **The real, defensible result is the paired contrast.** The stability-only gate
   falsely certified **11/30** negatives (FPR 0.367, 95% CI [0.22, 0.54]); the
   recalibrated gate certified **0/30**. Exact McNemar on the negatives: b=11, c=0,
   **p=0.001** — the reduction in false certification is significant. The mechanism is
   partly abstention, and we say so, but the improvement over stability-only is real.
3. **No detectable TPR cost.** TPR 1.000 → 0.967 is a single discordant pair; exact
   McNemar **p=1.0**. Do not claim a "3% cost" — it is not distinguishable from zero
   (though the study has no power to exclude one).

**And the increment is a null recalibration, not a new instrument.** Head-to-head in
`ablation_honest.json`: the SigClust-style single-Gaussian p-value **alone**
(threshold α=0.05) reproduces the full gate's TPR and FPR exactly. The gap statistic
and Hartigan dip are computed but never enter the decision rule
(`passed = obs > sg_p95 and sg_p < alpha`); the k-stability route only adds
abstentions. Frame the contribution as replacing the feature-permutation null with a
single-Gaussian null (SigClust; Liu et al. 2008) on the existing bootstrap-ARI
statistic — not as a multi-test discreteness instrument.

**Resolution / operating characteristic** (`separation_sweep.py` /
`results/separation_sweep.json`): sweeping component separation δ from 0 to 3 SD
(n=120, n_ref=100), the gate **does** resolve — median p collapses from ≈0.8 to 0.03
and certify rate rises 0 → 0.75 between δ=2 and δ=3 — refuting "flat at the floor".
But the transition is **sharp and late**: it certifies nothing through δ≤2, so the
recalibrated null is **conservative** (low sensitivity to moderate structure), abstaining
in the mid-range. Good for rejecting over-called subtypes; not a sensitive detector of
subtle real structure. Reported as a limitation, not tuned away.

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
`pathway_subtyping.clustering_dl`) — see the v0.8.0 availability caveat in
[`../../RUNME.md`](../../RUNME.md).

## Result 1a — the ablation (R3.10)

R3.10 asked: *"What happens to the performance if some of the validation gates
are removed?"* Four synthetic conditions with known ground truth — two positive
(`discrete_k2`, `discrete_k3`: genuine well-separated components) and two
negative (`single_gaussian`: one correlated blob; `continuum_1d`: a 1-D gradient,
the tumor-purity / immune-infiltration analog).

| Gate subset | TPR | FPR | `continuum_1d` cert rate |
|---|---|---|---|
| `stability_only` (pre-v0.8.0) | 1.000 | **0.367** | **0.733** |
| `discreteness_only` | 0.967 | 0.000 | 0.000 |
| `stability+discreteness` (v0.8.0) | 0.967 | **0.000** | 0.000 |

**Reading.** The pre-v0.8.0 stability-only gate false-certifies the 1-D continuum
in 73% of runs — this is the erratum finding reproduced under controlled
conditions: the old bootstrap null tested *pathway independence*, not
*discreteness*, so a reproducibly-bisected gradient passed. Adding Gate A
collapses the false-positive rate to zero.

**Honest note on TPR:** the ablation does not leave true positives *untouched* —
TPR goes 1.000 → 0.967, i.e. one of the 30 genuine-structure datasets (a
`discrete_k3` replicate) is lost. The correct claim is that the discreteness gate
eliminates the continuum false positives **at a cost of ~3% of true positives**,
not that TPR is unchanged. Do not state "TPR held".

Source: `results/gate_ablation_results.json` (+ `_raw.csv` for per-replicate rows).
Figure: `results/gate_ablation_figure.{png,svg}`.

## Result 1b — clusterer generality (R2.2)

The same synthetic data clustered three ways — GMM, DEC (Xie 2016), VAE-GMM
(VaDE, Jiang 2017) — then passed to Gate A.

**Primary result:** Gate A rejects the 1-D continuum on **100%** of runs, and does
so identically regardless of which algorithm drew the partition. This is the
substantive R2.2 answer: the framework is **not another clustering method to be
benchmarked against DEC/VAE** — it is a validation layer that wraps any of them,
because Gate A tests the discreteness of the *data* at a given k, not the
confidence of the algorithm.

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
