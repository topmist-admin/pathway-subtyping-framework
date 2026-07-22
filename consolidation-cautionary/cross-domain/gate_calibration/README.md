# Discreteness-gate real-data calibration (the paper's defensive keystone)

The cautionary-framework paper hinges on one claim a hostile reviewer will
attack: *is the discreteness gate calibrated, or just pessimistic?* A synthetic
ablation is not enough. This experiment shows on REAL public data (cBioPortal, no
auth) that the gate certifies genuine discrete structure and rejects a genuine
continuous gradient. Run: `scripts/calibrate_discreteness_gate.py`.

## Result — the gate calls BOTH real controls correctly

| Control | Ground truth | Recovery | Gate A verdict | Correct? |
|---|---|---|---|---|
| **DISCRETE** — 3 pooled tumor types (COAD+GBM+LUAD, n=1262) | tumor-of-origin (unambiguous) | ARI **0.921** | **CERTIFIED** "discrete structure" (SigClust p=0.007) | ✓ |
| **CONTINUUM** — immune-infiltration gradient (LUAD, n=510) | continuous cytolytic/T-cell signature (unimodal, dip p=0.98) | k=2 split tracks the gradient (std-gap 1.09) | **REJECTED** "no discrete structure (continuum)" (SigClust p=0.17) | ✓ |

**Both correct.** The gate says YES to real discreteness and NO to a real
continuum — it is calibrated, not merely pessimistic. This is the direct answer
to the strongest attack on the reframe.

## Why each control is defensible
- **Discrete:** three biologically distinct cancers are unambiguously ≥3 groups.
  Recovery ARI 0.921 confirms the partition IS the tumor types; the gate then
  certifies that structure. (A gate blind to this would be broken.)
- **Continuum:** immune infiltration is a canonical CONTINUOUS axis of tumor
  variation (cytolytic signature, Rooney et al. 2015). The k=2 clustering merely
  bisects the gradient (immune-hi vs immune-lo, std-gap 1.09), and the gate
  correctly declines to call that a discrete subtype.

## Honesty notes / limitations
- **Purity was the first-choice continuum but unavailable** — the pan-can atlas
  cBioPortal studies do not expose tumor-purity/ESTIMATE as a clinical attribute
  (only aneuploidy). Immune infiltration is the substitute continuum: defensible,
  and the gate's PRIMARY single-Gaussian (SigClust) test does the rejecting, not
  a purity assumption.
- **Unimodality of the immune axis** confirmed: Hartigan dip p=**0.98** (cannot
  reject unimodality → a genuine gradient, not two groups). Requires the `diptest`
  extra (installed).
- **Informative nuance — why SigClust, not the dip test, is the primary gate.**
  On the DISCRETE control, the dip test on PC1 is *unimodal* (p=0.81) even though
  the data is genuinely 3-group — the tumor types don't separate along PC1 alone.
  The dip-on-one-axis test would MISS this structure; Gate A's primary
  single-Gaussian (SigClust) test, which uses the full reduced space, correctly
  certifies it (p=0.007, recovery 0.92). This real-data case validates Gate A's
  design choice: SigClust as the decision test, dip as corroborating only. On the
  CONTINUUM both agree (dip unimodal AND SigClust not-discrete).
- **A bug in the first version of THIS experiment** pooled within-study z-scores,
  which erased the between-tumor-type signal (recovery collapsed to 0.05).
  Corrected by pooling continuous expression and z-scoring across the combined
  set → recovery 0.921. (Noted because a calibration experiment that was itself
  miscalibrated is exactly the kind of error the paper's thesis warns about.)

## Role in the paper
This is **Result 1 (methods validation), real-data half** — pairs with the
synthetic ablation figure. Together: the discreteness gate is validated on
synthetic ground truth AND on real cohorts with known discrete/continuous
structure. Only after this is established do the audit (Result 2) and the
cautionary worked examples (cancer, psychiatry) carry weight.

## Provenance
cBioPortal continuous RNA-seq (`*_rna_seq_v2_mrna`), Hallmark 50-set panel;
immune signature CD8A/GZMA/GZMB/GZMK/PRF1/IFNG/NKG7/CCL5/CXCL9-11/IDO1/STAT1/
GBP1/HLA-DRA. Deterministic (seed 42), Gate A n_ref=150. `results/gate_calibration.json`.
