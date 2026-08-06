# Cancer arm of R3.8 — large independent-cohort validation (TCGA-BRCA / PAM50)

Public-data (cBioPortal, no auth) large-cohort validation for Scientific Reports
reviewer point R3.8. Uses public data only; large-N psychiatric validation on
controlled-access cohorts is future work. Run: `scripts/fetch_and_validate_brca.py`.

> ⏱️ **This is the slowest package in the bundle: expect well over 50 minutes**, and it
> prints nothing between the fetch line and the final result. A long silence is expected,
> not a hang. Do **not** run it at the same time as another cBioPortal package — two
> concurrent `molecular-data/fetch` POSTs throttle each other to a standstill (measured:
> two jobs that each complete alone made no progress for 45 minutes when run together).
>
> The DL baselines additionally require `torch`, pinned to **2.11.0** in
> `requirements.txt`. They are **not deterministic across torch releases**, so the
> acceptance bar for this package is the *conclusion* — PSF competitive with (not better
> than) the DL baselines, k=5 stability failing the 0.80 bar, recovery metric-dependent —
> not the digits. The torch version used is recorded under
> `provenance.environment.torch`.

## Result (n=1,082 TCGA-BRCA; 981 PAM50-labelled; k=5)

**Recovery of PAM50 — TWO metrics (the metric choice is itself a finding):**

| Method | k-way ARI (strict) | best-subtype enrichment (CMS4-style, generous) |
|---|---|---|
| **PSF (pathway-GMM)** | **0.218** | LumA **87.6%** (OR 11.6, p 8e-47) |
| VAE-GMM (VaDE) | 0.125 | LumA **89.1%** (OR 12.6) |
| k-means (pathway) | 0.160 | LumA 76.6% (OR 4.4) |
| DEC | 0.088 | LumA 75.7% (OR 4.2) |

**The metric choice is a cautionary finding in itself.** The SAME PSF partition
scores ARI 0.218 (looks weak) and LumA-enrichment 87.6% (looks strong). The
manuscript's headline "75.9% CMS4" was single-subtype enrichment — the generous
metric; measured that way PSF's BRCA recovery (87.6%) is comparable-to-stronger.
Single-subtype enrichment FLATTERS subtyping; the k-way ARI is the honest test.
This is a worked example of a meta-point for the paper: reported recovery
numbers depend enormously on the metric, and the commonly-reported one overstates.

**Method comparison is NOT metric-robust:** PSF leads on ARI; VAE-GMM edges it on
enrichment (89.1% vs 87.6%). Do NOT claim "PSF beats deep learning" — the honest
statement is that pathway-GMM is *competitive* with DL, with the ranking
depending on the metric.

**Validation gates (forced k=5 to match PAM50, AND framework's BIC k=2):**

| Gate | at k=5 | at BIC k=2 |
|---|---|---|
| label-shuffle (no artifact) | PASS (ARI≈0) | — |
| discreteness Gate A | **CERTIFIED** (discrete structure) | **CERTIFIED** (discrete structure) |
| bootstrap-stability (0.80 bar) | FAIL (ARI 0.39) | FAIL (ARI 0.43) |

## Honest interpretation (do NOT overstate)
- **Pathway-GMM is COMPETITIVE with deep clustering, ranking depends on metric**
  (PSF leads on ARI; VAE-GMM edges it on enrichment). The honest R2.2 answer is
  that standard-setting DL does not clearly beat pathway-GMM — not that PSF wins.
  (Superseded the earlier "PSF beats every baseline" claim, which was ARI-only.)
- **Recovery is metric-dependent, not simply "modest."** By strict k-way ARI it
  is weak (~0.09–0.22, because PAM50 is a specific gene classifier not Hallmark
  pathways); by single-subtype enrichment it is strong (LumA 87.6%). Report BOTH;
  the gap between them is the cautionary point.
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

---

## CPTAC-BRCA multi-omic (R1.6 — protein modality; n=122)

Matched mass-spec **protein** + mRNA for the same tumors (CPTAC breast, Krug et
al. Cell 2020), public via cBioPortal. Run: `scripts/fetch_and_validate_cptac_brca.py`.

**Recovery of PAM50 (k=5) — both metrics:**
| Modality | k-way ARI (strict) | best-subtype enrichment (generous) |
|---|---|---|
| mRNA (pathway) | 0.189 | Basal **88.2%** (OR 48.8, p 1e-9) |
| **protein (pathway)** | **0.172** | LumA **87.0%** (OR 11.2, p 1e-5) |

The metric-dependence seen on TCGA-BRCA recurs here across BOTH modalities: weak
by ARI (~0.18), strong by enrichment (~87–88%). Protein enrichment (87%) matches
mRNA — the R1.6 protein-modality point holds on the generous metric too.

**Gates (both modalities):** discreteness Gate A **CERTIFIED** (discrete
structure); bootstrap-stability **FAIL** at 0.80 bar (mRNA 0.47, protein 0.36).
**Expression↔protein subtype concordance:** ARI **0.166** (n=122 shared).

**Honest interpretation:**
- **Directly answers R1.6.** PSF's pathway subtyping operates on the PROTEIN
  modality and recovers PAM50 *comparably to expression* (protein AMI 0.245 is
  actually higher than mRNA's 0.199). The framework is not expression-specific —
  the validation the reviewer said was "mainly gene expression" now extends to
  mass-spec proteomics.
- **Gates are modality-general.** The discreteness gate certifies real structure
  in both modalities; the same fixed-0.80 stability behaviour recurs — reinforcing
  the R3.7 "fixed threshold is the wrong instrument" point across modalities.
- **Modest expression↔protein concordance (ARI 0.17) is expected, not a defect.**
  mRNA–protein correspondence is known to be imperfect; the two modalities agree
  with each other about as much as each agrees with PAM50. State plainly; it
  motivates careful (not naive-concatenation) multi-omic fusion.
- Caveat: n=122 is small — CPTAC's value here is the *modality*, not size (TCGA-BRCA
  supplies the large-cohort size). Recovery modest for the same Hallmark≠PAM50
  reason as TCGA-BRCA.

## Provenance
cBioPortal study `brca_tcga_pan_can_atlas_2018`; Hallmark 50-set panel
(`../../panels/hallmark_200genes.gmt`); mRNA median-Zscores; PAM50 = patient-level
`SUBTYPE`. n_ref=200 (the BRCA fetcher's default; the CPTAC fetcher uses 100). ⚠️ The gate is deterministic (seed 42), but the deposited **bootstrap-stability** figures were produced before `gmm_seed` was passed to `stability_test_bootstrap`, so their per-bootstrap GMM initialisations were unseeded — which is why two deposited runs give k=5 ARI 0.4076 and 0.3918. The call sites are now seeded; re-running will produce one stable value that will differ slightly from both deposited figures. Results: `results/brca_pam50_validation.json`.
