# Discreteness-gate real-data calibration (within-study controls)

The reframed paper hinges on one claim a hostile reviewer will attack: *is the
gate calibrated, or just pessimistic?* This experiment shows on REAL public data
(cBioPortal, no auth) that the gate certifies genuine discrete structure and
rejects a genuine continuous gradient — with **both controls drawn from within a
single study each**, so study cannot stand in for the ground-truth label.

**Run:** `scripts/calibrate_within_study.py` → `results/gate_calibration_within_study.json`

## ⚠️ Supersedes the pooled-3-study version — read this

The earlier `calibrate_discreteness_gate.py` used a discrete positive control that
**pooled three cBioPortal studies** (COAD/GBM/LUAD), one tumour type per study. That
made study perfectly confounded with the label. The tell is in that script's own
run: within-study z-scores gave recovery ARI **0.05**; only across-study z-scores
(which reintroduce the batch axis) gave 0.92. A reviewer correctly flagged that the
"discrete" control was largely a **batch effect** — the very failure this paper
diagnoses in the flagship. **That result is withdrawn.** The pooled script and its
`gate_calibration.json` are retained in the repo only as the documented negative
example; do not cite their 0.921.

## Result — within-study controls, both called correctly

| Control | Ground truth | Recovery | Gate A verdict | Correct? |
|---|---|---|---|---|
| **DISCRETE** — IDH-mutant vs IDH-wildtype low-grade glioma (TCGA-LGG, n=507; 415 mut / 92 wt), **single study** | IDH status — the most discrete molecular dichotomy in adult glioma | ARI **0.418** (within-study z-scores) | **CERTIFIED** "discrete structure" (SigClust p ≤ 1/151) | ✓ |
| **CONTINUUM** — immune-infiltration gradient (TCGA-LUAD, n=510), **single study** | continuous cytolytic signature (unimodal, Hartigan dip p=0.98) | k=2 split just bisects the gradient (std-gap 1.09) | **REJECTED** "no discrete structure (continuum)" (SigClust p=0.17) | ✓ |

**Both correct, with no batch confound.** The gate says YES to a genuinely discrete
within-disease axis and NO to a within-disease continuum.

## Honest reading of the recovery number

IDH recovery is **ARI 0.418**, not the pooled version's inflated 0.921. That is the
*right* number and it is on-thesis: generic Hallmark-pathway subtyping only weakly
recovers even a genuinely discrete molecular axis (the same point Result 3 makes for
PAM50). The load-bearing fact is not the recovery magnitude — it is that the gate
**certifies** the discrete axis and **rejects** the continuum. Discreteness (does
structure exist?) and recovery (how well do generic pathways reconstruct known
labels?) are different questions; the gate answers the first correctly while
recovery of the second is honestly modest.

## Scope limit — this is two anchors, not an operating characteristic

These two controls sit near opposite ends of the difficulty range. They establish
that the gate is **not degenerate** (not always-yes, not always-no) on real data.
They do **not** establish that it resolves borderline within-disease structure. That
question — where does the verdict flip as separation decreases? — is answered by the
synthetic separation sweep in `../gate_ablation/` (`separation_sweep.py`), which
varies component separation from 0 to a clean split and reports the certify-rate and
p-value transition. Cite the sweep, not this pair, for resolution claims.

## p-value floor

The SigClust empirical p is `(#ref ≥ obs + 1)/(n_ref + 1)`, so with `n_ref=150` the
smallest reportable value is `1/151 ≈ 0.0066`. The discrete control sits at that
floor; report it as **p ≤ 1/(n_ref+1)**, not as a point estimate. The continuum
control's p=0.17 is above the floor, so it carries resolution.

## Provenance

cBioPortal studies `lgg_tcga_pan_can_atlas_2018` (IDH via patient-level `SUBTYPE`,
collapsed IDHmut/IDHwt) and `luad_tcga_pan_can_atlas_2018`; Hallmark 50-set panel;
within-study median z-scores for the discrete control (valid — one study);
deterministic (seed 42, n_ref=150). Network: cbioportal.org (public).
