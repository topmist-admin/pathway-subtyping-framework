# Expanding PSF beyond autism: other psychiatric disorders + cancer

**Type:** Research-direction roadmap (design + demonstrated confound remap). Speculative items ⚠️.
**Date:** 2026-07-16 · **Standalone bundle** — separate from the autism reproduction done earlier.
**Builds on:** the PSF v0.7.0 autism/SCZ reproduction (Gate 6 caught a brain-region artifact; feature-level
Gate 7 anchoring validated). This document asks: does that machinery generalize to (a) other psychiatric
disorders and (b) cancer — and what changes?

> **Research use only.** Named drugs/trials are cited as real development history for methodological
> illustration, not treatment recommendations. Clinical decisions belong to qualified clinicians.

---

## 0. The one-sentence thesis

**PSF is really a *confounded-subtype detector*; the disease only sets what the confound is.** The autism
work showed Gate 6 flagging a partition that was secretly a **brain-region** classifier. Expanding to other
diseases is mostly a matter of **re-specifying the confound registry** — not re-architecting the framework.
Cancer is the more valuable addition of the two, for reasons that do *not* apply to psychiatry (§2).

A runnable demonstration ships with this bundle: `scripts/confound_remap_demo.py` plants a nuisance-aligned
partition in three domains and shows the *same* v0.7.0 gate flag the domain-appropriate axis each time
(results in `out/confound_remap/`):

| domain | dominant nuisance | measured Cramér's V (planted alignment) | gate verdict |
|---|---|---|---|
| psychiatry (multi-region brain) | brain region | 0.93 (align 0.92) | FAIL (worst = brain region) |
| cancer (pan-cancer pooled) | tissue of origin | 0.88 (align 0.90) | FAIL (worst = tissue of origin) |
| cancer (single tumor type) | tumor purity | 0.48 (align 0.55) | FAIL (worst = tumor purity) |

Identical gate, identical logic, remapped registry.

---

## 1. Other psychiatric disorders — the work already half-covers this

The autism reproduction was **not** SCZ-only. GSE80655 is a four-diagnosis cortex cohort, and the
cross-diagnosis region gate (Table 2 of the reproduction) already ran the full panel:

| contrast | region V (Bergsma) | diagnosis p |
|---|---|---|
| SCZ + Control | 0.671 | 0.203 (n.s.) |
| BD + Control | 0.659 | 0.466 (n.s.) |
| MDD + Control | 0.702 | 0.940 (n.s.) |
| All-4 pooled | 0.716 | 0.001 |

So the cautionary result **already replicates across schizophrenia, bipolar, and major depression**: region
dominates (V ≈ 0.66–0.72), diagnosis is n.s. in every pairwise contrast. Psychiatric expansion is therefore
about **breadth and the harder affirmative test**, not a new capability:

1. **More disorders / cohorts** — PTSD, OCD, ASD-in-other-banks; any multi-region cohort with region metadata.
   Effort: low. Reuse the autism scripts unchanged.
2. **The affirmative counterpart — within-region, single-diagnosis subtyping.** The autism audit's F8 notes
   the one design that *worked*: Bowen 2019 found two SCZ subtypes from a **single region**, replicated in an
   independent cohort. Running Gate 6 → Gate 7 on a within-region design and asking whether a genuine
   within-disorder subtype survives is where psychiatric expansion earns its keep. This is the positive case
   to complement the negative one.
3. **Disorder-specific confound remaps.** Region is not the only nuisance. For mood disorders, medication
   (lithium, antipsychotics) and post-mortem interval carry large transcriptomic signal; for any disorder,
   the §5 registry must be re-derived (a confound is a claim about causal structure, not a variable name).

**Ceiling:** psychiatry inherits the autism study's core limitation — germline genetic risk is heterogeneous
and low per-subject, so genetic anchoring (Gate 7) is a *low-power* test. Feature-level anchoring works
(demonstrated for autism); subject-level stays access- and power-limited.

---

## 2. Cancer — a different, and in three ways stronger, fit

Cancer is worth adding for reasons that psychiatry cannot supply:

**2.1 The framework already ingests cancer data (architecturally).** The plumbing that reads a psychiatric
expression matrix reads a tumor expression matrix — no re-architecting; the v0.7.0 47-dataset benchmark
included TCGA oncology rows, so tumor data flows through the pipeline unchanged. **What this bundle does NOT
claim is a validated cancer-subtype recovery.** The one TCGA row in that benchmark — **TCGA-COAD** — was
flagged **`valid=False`** in the 2026-07 correction (`corrected_benchmark_47datasets_v2.csv`:
`degenerate_ground_truth(n_true_clusters=1)→empty_input_ARI_artifact; NOT_independent:TCGA-COAD_is_manuscript_primary_cohort`),
i.e. its apparent "perfect" recovery was the empty-input ARI artifact, not a real result — and TCGA-COAD is
non-independent of the manuscript anyway. So no recovery statistic from that benchmark can be cited as evidence
here. Demonstrating that the battery **recovers a known cancer subtype (CMS / PAM50)** on a *valid, independent*
cohort is exactly **Stage A** below — a test to be run, not one already passed. (Earlier drafts cited a
"TCGA CMS4 recovery at OR ≈ 16.5"; that number is not in the corrected benchmark and pointed at the invalidated
TCGA-COAD row — removed.)

**2.2 Cancer is the positive-control domain the psychiatric story lacks.** Colorectal **CMS** (consensus
molecular subtypes), breast **PAM50**, and TCGA pan-cancer classes are *real, replicated, clinically-actionable*
molecular subtypes. Running the battery on them should return **subtype = REAL** — precisely the ground-truth
calibration the audit says the instrument needs (§0.5). A detector that only ever returns "artifact" is as
uninformative as one that only returns "subtype." Cancer proves the gates can *pass* a true positive, which
psychiatry (where the honest answer keeps being "region artifact") never demonstrates.

**2.3 Genetic anchoring (Gate 7) is far more powerful in cancer.** Autism's germline risk is diffuse; cancer
carries **somatic drivers, copy-number alterations, and mutational signatures** — high per-tumor genetic
signal, causally upstream of the transcriptome, and directly druggable. A transcriptomic subtype that aligns
with a driver-mutation stratum is anchored with orders of magnitude more power than any autism analog. In
cancer, Gate 7's anchor changes from *germline GWAS enrichment* to *somatic-driver / CNA / signature
alignment* — a stronger and easier test.

### 2.4 The catch: the confound structure inverts

In psychiatry the dominant nuisance was **brain region**. In cancer it becomes **tissue of origin** (pan-cancer)
and **tumor purity / stromal-immune infiltration** (within a tumor type). A pan-cancer partition will recover
*organ*, exactly as the psychiatric one recovered *anatomy* (demonstrated in the §0 table). So Gate 6 does not
transfer unchanged — its **nuisance registry must be re-specified**:

| | psychiatry (autism work) | cancer |
|---|---|---|
| dominant nuisance | brain region | tissue of origin; **tumor purity** (ESTIMATE / ABSOLUTE score) |
| secondary nuisances | PMI, RIN, pH, batch, medication | stromal & immune fraction, batch, FFPE vs frozen, sequencing centre |
| biology-of-interest (gate-exempt) | diagnosis | cancer type, or known subtype (CMS / PAM50) |
| Gate-7 anchor | germline GWAS enrichment (low power) | **somatic driver / CNA / mutational signature (high power)** |
| replication axis | subject-disjoint cohort | independent cohort (TCGA → ICGC / CPTAC) |

This remap is a small, well-defined code change (add a purity/infiltration confound module feeding Gate 6's
nuisance set), not a redesign — see `confound_remap_demo.py`.

---

## 3. Target nomination in cancer — grounded in fetched evidence

The colorectal-cancer CMS framework is the cleanest worked example: each molecular stratum maps to an
approved, subtype-stratified therapy. **All identifiers below were fetched live this session** (ChEMBL for
drug/MoA, ClinicalTrials.gov for Phase-3 stratified trials; see `evidence/`):

| CRC stratum (≈CMS) | druggable node | approved agent | ChEMBL id / MoA | Phase-3 subtype-stratified trials |
|---|---|---|---|---|
| BRAF V600E | BRAF kinase | encorafenib | CHEMBL3301612 — "Serine/threonine-protein kinase B-raf inhibitor" | 4 (e.g. **NCT02928224**, the BEACON-CRC design, completed) |
| MSI-H / dMMR (≈CMS1 immune) | PD-1 | pembrolizumab | CHEMBL3137343 — "Programmed cell death protein 1 inhibitor" | 8 (e.g. NCT03755739, NCT04776148) |
| RAS-wt (≈CMS2 canonical) | EGFR | cetuximab | EGFR-targeting mAb (biologic) | 72 (e.g. NCT00208546, NCT03391934) |

The existence of Phase-3 **subtype-stratified** trials for each stratum is the external validation the autism
work could only aspire to: in cancer, the molecular-subtype → drug-target → stratified-trial chain is already
built and clinically proven. This is why cancer is the domain where PSF's therapeutic-ranking / knowledge-graph
modules (`autism/therapeutic/`, `knowledge_graph/`) pay off most directly — the oncology target space is dense
and validated.

**The five-filter nomination pipeline** (from the autism roadmap §3) carries over, with the anchoring filter
strengthened: subtype-specificity → **somatic genetic support** (driver/CNA alignment — far stronger than
germline PRS) → directionality (in-silico perturbation) → druggability (ChEMBL/DGIdb) → cell-type/purity
grounding (single-cell or purity-adjusted, the cancer analog of the interneuron-class check).

---

## 4. Staged plan (each stage has a publishable negative)

- **Stage A — Cancer positive control.** Cluster TCGA-CRC in a validated feature space; confirm the battery
  **recovers CMS** (subtype = REAL). This calibrates the instrument on a known truth — the check the
  psychiatric side never got to pass.
- **Stage B — Confound remap (Gate 6).** Add tumor-purity (ESTIMATE/ABSOLUTE) + tissue-of-origin to the
  nuisance registry; verify a pan-cancer pooled partition FAILS on tissue (as the §0 demo shows synthetically),
  and a within-CRC partition is not merely a purity axis.
- **Stage C — Somatic Gate 7.** Test whether a transcriptomic subtype aligns with somatic driver strata
  (BRAF/KRAS/MSI), CNAs, or mutational signatures — the high-power cancer anchor.
- **Stage D — Subject-disjoint replication.** TCGA → ICGC/CPTAC (independent-cohort, and unlike the legacy
  psychiatric GEO sets, genuinely disjoint donors).
- **Stage E — Target nomination.** Five-filter pipeline; the CRC table in §3 is the template.
- **Psychiatry track (parallel, low-effort):** add PTSD/OCD cohorts; run the Bowen-style within-region
  affirmative test.

---

## 5. What this bundle contains

```
docs/01_cross_disorder_cancer_roadmap.md   <- this document
scripts/confound_remap_demo.py             <- runnable: same Gate 6, three domains (psychiatry/cancer x2)
out/confound_remap/                        <- its output (region V=0.93, tissue V=0.88, purity V=0.48)
evidence/
  crc_chembl_targets_drugs.json            <- fetched: encorafenib/pembrolizumab ChEMBL ids + MoA
  crc_subtype_stratified_trials.json       <- fetched: Phase-3 stratified CRC trial counts + NCT ids
logs/                                      <- run logs
```

Requires `pathway-subtyping==0.7.0` + numpy for the demo (same env as the autism bundle; see that bundle's
HANDOFF_README §1 for pins). Evidence JSONs are static snapshots — re-fetch via the ChEMBL / ClinicalTrials
connectors.

---

## 6. Bottom line

- **Psychiatry:** expansion = breadth + the affirmative within-region test. The negative (region artifact)
  already replicates across SCZ/BD/MDD. Low effort, inherits the low-power anchoring ceiling.
- **Cancer:** the higher-value expansion. It gives PSF (1) a domain where the gates should *pass* a true
  positive (CMS/PAM50), (2) a much stronger, somatic genetic anchor, and (3) a dense, clinically-validated
  target space where the therapeutic modules pay off — at the cost of one well-defined confound remap
  (region → purity/tissue), which this bundle demonstrates is a registry change, not a redesign.
