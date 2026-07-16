# Gate 7 (somatic) — real TCGA-CRC positive control

Wires `ValidationGates.somatic_anchoring_gate` against **real** colorectal-cancer
data (TCGA PanCancer Atlas, public cBioPortal API) with actual **BRAF-V600E /
KRAS / MSI** strata — the roadmap's Stage C, on real driver calls rather than
synthetic ones, and the somatic analog of the Voineagu positive control that
validated the feature-level anchoring gate.

Two transcriptomic partitions are tested against the **same** real strata, and
the contrast is the point.

## Reproduce

```bash
# needs the framework (>=0.8.0 line) + numpy pandas scikit-learn requests
python scripts/fetch_and_run_tcga_crc.py       # fetches live from cBioPortal
```

Deposited snapshot: `results/tcga_crc_somatic_result.json`. Real strata recovered:
**BRAF-V600E = 48, KRAS-mut = 218, MSI-H = 85** (MANTIS > 0.4).

## A. Canonical CMS partition — the strong positive control (PASS)

Partition = the **published Consensus Molecular Subtypes** (Guinney et al. 2015,
CRC Subtyping Consortium public-tier labels, vendored in `inputs/`). CMS is
expression-defined and independent of the somatic strata, so this is a genuine,
non-circular test — and it reproduces textbook CRC biology.

n = 430 tumors (CMS-labelled ∩ sequenced): CMS1 57 · CMS2 188 · CMS3 62 · CMS4 123.

| gate | result |
|---|---|
| Gate 6 (confound: tissue) | **pass**, V = 0.175 — mildly organ-associated, not a tissue classifier |
| Gate 7 (somatic) | **PASS**, anchored = **BRAF-V600E + MSI**, max **V = 0.656** (q ≈ 1e-40) |

| stratum | Cramér's V | BH q | anchored |
|---|---|---|---|
| BRAF-V600E | **0.656** | 6e-40 | **yes** |
| MSI | **0.656** | 1e-39 | **yes** |
| KRAS | 0.280 | 6e-8 | no (significant, just below the 0.30 bar) |

Row-fraction crosstabs (deposited) — the biology the gate is anchoring to:

| CMS | % MSI-H | % BRAF-V600E | % KRAS-mut |
|---|---|---|---|
| **CMS1** (MSI-immune) | **80%** | **60%** | 25% |
| CMS2 (canonical) | 2% | 0% | 31% |
| **CMS3** (metabolic) | 23% | 3% | **71%** |
| CMS4 (mesenchymal) | 11% | 5% | 37% |

The gate PASSES, anchoring the transcriptomic subtype to BRAF-V600E and MSI at a
large effect size — the somatic counterpart of the feature-level gate reproducing
Voineagu's 16.5× enrichment. (KRAS is real CMS3 biology — 71% mutant — but diluted
across the 4-way partition to V = 0.28, so it lands just under the anchor bar;
significant, correctly not over-called.)

## B. Unbiased PSF clustering — the honest contrast (not anchored)

Partition = built by the framework from the disease-agnostic Hallmark panel
(pathway means → Gaussian mixture, **k = 2 by BIC**, not tuned to the outcome).

| gate | result |
|---|---|
| Gate 6 (confound: tissue) | pass, V = 0.000 — not a tissue classifier |
| Gate 7 (somatic) | **not anchored**, max V = 0.13 — KRAS/MSI significant (q < 0.05) but effect-size-modest (V < 0.30) |

A generic whole-transcriptome partition couples to the drivers only weakly.
Recovering the driver-aligned structure de novo needs the **dedicated CMS
classifier**, not generic clustering — and the clustering was **not** tuned to
raise the association (that would be the outcome-fitting the framework's 2026-07
correction warns against).

## Together

Same gate, same real strata, two partitions:

- a **validated transcriptomic subtype (CMS)** anchors strongly → gate **PASS**
  (V = 0.66);
- a **generic unbiased partition** couples weakly → gate **withholds** the anchor
  (V < 0.30), correctly refusing to over-call a significant-but-weak association.

That is exactly the behaviour a positive-evidence gate needs: it fires when the
partition captures real driver-aligned biology, and stays silent when it does not.

## Provenance / caveats

- **Somatic + expression:** cBioPortal `coadread_tcga_pan_can_atlas_2018` (public,
  no auth). Expression = `rna_seq_v2_mrna_median_all_sample_Zscores`; mutations =
  `..._mutations`; MSI = `MSI_SCORE_MANTIS` (> 0.4 → MSI-H); tissue =
  `CANCER_TYPE_DETAILED`. Wild-type = in the `_sequenced` list but not in the
  gene's mutation records (a real negative, not "untested").
- **CMS labels (`inputs/cms_labels_public_all.txt`):** CRC Subtyping Consortium
  **public-tier** labels from Guinney JJ et al., *Nat Med* 2015 (21:1350–1356,
  doi:10.1038/nm.3967). Column `CMS_final_network_plus_RFclassifier_in_nonconsensus_samples`;
  `NOLBL` (unclassified) dropped; patient barcodes mapped to the `-01` primary-tumor
  sample. Obtained from the public `biobakery/crc-subtyping-paper` mirror.
- **Panel (result B):** `../../panels/hallmark_200genes.gmt` (the repo's Hallmark
  panel, ~4,384 genes / 50 sets; disease-agnostic).
- Deterministic (seed 42). cBioPortal is live — re-fetched counts can drift by a
  few tumors over time.

**Research use only.** Not clinical guidance.
