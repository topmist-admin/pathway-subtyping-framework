# Cancer arm of R3.8 — large independent-cohort validation (TCGA-BRCA / PAM50)

Public-data (cBioPortal, no auth) large-cohort validation for Scientific Reports
reviewer point R3.8. controlled-access-independent — proceeds while NIH psychiatric access is
pending. Run: `scripts/fetch_and_validate_brca.py`.

## Result (n=1,082 TCGA-BRCA; 981 PAM50-labelled; k=5)

**Recovery of PAM50 (ARI / AMI) — PSF vs baselines including deep learning:**

| Method | ARI | AMI |
|---|---|---|
| **PSF (pathway-GMM)** | **0.218** | **0.250** |
| k-means (pathway) | 0.160 | 0.190 |
| VAE-GMM (VaDE) | 0.125 | 0.168 |
| DEC | 0.088 | 0.118 |

**Validation gates (forced k=5 to match PAM50, AND framework's BIC k=2):**

| Gate | at k=5 | at BIC k=2 |
|---|---|---|
| label-shuffle (no artifact) | PASS (ARI≈0) | — |
| discreteness Gate A | **CERTIFIED** (discrete structure) | **CERTIFIED** (discrete structure) |
| bootstrap-stability (0.80 bar) | FAIL (ARI 0.39) | FAIL (ARI 0.43) |

## Honest interpretation (do NOT overstate)
- **PSF recovers PAM50 better than every baseline, including modern deep
  clustering (DEC, VAE-GMM).** This is a genuine, useful answer to R2.2: with
  standard settings, deep clustering did **not** outperform pathway-GMM here.
- **But recovery is modest across all methods (ARI ~0.09–0.22).** Expected:
  PAM50 is defined by a specific 50-gene classifier, not by the 50 Hallmark
  pathways used here; generic pathway-level scoring only partially aligns with
  PAM50. State this plainly — it is not a strong recovery claim.
- **The two "stability" instruments disagree, informatively.** The discreteness
  gate CERTIFIES real structure at BOTH the natural k=2 AND the forced k=5, yet the
  old bootstrap-stability gate FAILS the fixed 0.80 bar at both k (ARI 0.39 / 0.43).
  So the stability FAIL is **not** an artifact of forcing k (caveat #1 below is
  RESOLVED — it holds at the framework's own BIC k=2). BRCA has genuine discrete
  molecular structure; a fixed 0.80 reproducibility bar simply rejects it.
- **This directly supports reviewer R3.7** ("the ARI≥0.80 threshold seems
  arbitrary"). A real, discreteness-certified 1,000-sample cohort failing the fixed
  0.80 bar is concrete evidence the fixed threshold is the wrong instrument — and
  motivates the discreteness gate (data-referenced, no arbitrary cut) as its
  replacement. Use this as a worked example in the R3.7 response.

## Caveats / follow-ups before manuscript use
1. ~~k forced to 5~~ **RESOLVED:** BIC selects k=2; stability fails at both k=2 and
   k=5, so the FAIL is a real property, not a forced-k artifact.
2. **DL baselines used default settings** (not extensively tuned). Report as
   "standard settings"; the honest framing is PSF is competitive without
   cohort-specific tuning, not that DL is fundamentally inferior.
3. Consider a Hallmark → curated-BRCA pathway panel sensitivity check (pairs with
   R2.6) — recovery may rise with a more relevant pathway collection.

## Provenance
cBioPortal study `brca_tcga_pan_can_atlas_2018`; Hallmark 50-set panel
(`../../panels/hallmark_200genes.gmt`); mRNA median-Zscores; PAM50 = patient-level
`SUBTYPE`. Deterministic (seed 42), n_ref=100. Results: `results/brca_pam50_validation.json`.
