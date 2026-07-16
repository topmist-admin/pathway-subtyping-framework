# Gate 7 (somatic) — real TCGA-CRC positive control

Wires `ValidationGates.somatic_anchoring_gate` against **real** colorectal-cancer
data (TCGA PanCancer Atlas, via the public cBioPortal API) with actual
**BRAF-V600E / KRAS / MSI** strata. This is the somatic analog of the Voineagu
positive control that validated the feature-level anchoring gate — the roadmap's
Stage C, on real driver calls rather than synthetic ones.

No CMS labels are borrowed. The transcriptomic partition is built by the
framework from a **disease-agnostic** Hallmark gene panel, so its alignment with
the somatic strata is a genuine, non-circular test.

## Reproduce

```bash
# needs the framework (>=0.8.0 line) + numpy pandas scikit-learn requests
python scripts/fetch_and_run_tcga_crc.py           # fetches live from cBioPortal
```

Deposited snapshot: `results/tcga_crc_somatic_result.json` (cBioPortal is a live
service — re-fetched counts can drift slightly over time).

## Two results, both on the same real cohort (n = 532 sequenced tumors)

Real strata recovered: **BRAF-V600E = 48**, **KRAS-mut = 218**, **MSI-H = 85**
(MANTIS > 0.4) — all realistic CRC frequencies.

### A. The wiring — unbiased transcriptomic partition vs drivers

Score the Hallmark panel into 50 pathway means, cluster with a Gaussian mixture
(**k = 2 selected by BIC**, not tuned to the somatic outcome), then run the gates:

| gate | result | reading |
|---|---|---|
| Gate 6 (confound: tissue COAD/READ) | **pass**, V = 0.000 | the partition is **not** a tissue classifier |
| Gate 7 (somatic: BRAF/KRAS/MSI) | **not anchored**, max V = 0.13 | KRAS (q≈4e-3) and MSI (q≈1e-2) are *significant* but *effect-size-modest* (V < 0.30) |

**This is the honest whole-transcriptome result.** An unbiased Hallmark partition
of CRC couples to somatic drivers only modestly — recovering the driver-aligned
subtypes at a medium effect size (V ≥ 0.30) needs the dedicated **CMS classifier**,
not generic pathway clustering. The clustering was **not** tuned to raise the
association (that would be exactly the kind of outcome-fitting the framework's
2026-07 correction warns against). The value here is that Gate 7 shows the right
**effect-size discipline**: it does not over-call a significant-but-weak
association as an anchor — a naive "q < 0.05 → anchored" rule would have wrongly
flagged KRAS and MSI.

### B. Gate-PASS positive control — real strong association

BRAF-V600E and MSI-high strongly co-occur in CRC (serrated / MLH1-methylation
pathway) — a textbook strong genomic association. Using MSI status as the
partition and BRAF/KRAS as strata:

| stratum | Cramér's V | BH q | anchored |
|---|---|---|---|
| **BRAF-V600E** | **0.549** | **2.5e-36** | **yes** |
| KRAS | 0.000 | 0.34 | no (correctly — KRAS is depleted in MSI-H) |

**Gate 7 passes**, anchored on BRAF-V600E, and correctly rejects KRAS. This
confirms the gate's PASS branch fires on real data when the alignment is genuinely
strong — the somatic counterpart to the feature-level gate reproducing Voineagu's
16.5× enrichment.

## Together

The two results validate the somatic gate end-to-end on real TCGA-CRC: it **fires
PASS** on a real strong driver association (B, BRAF↔MSI, V = 0.55), and it
**withholds** an anchor call when a real transcriptomic partition couples to
drivers only weakly (A, V < 0.30) — the exact behaviour a positive-evidence gate
needs to be trustworthy.

## Provenance / caveats

- **Source:** cBioPortal `coadread_tcga_pan_can_atlas_2018`, public API, no auth.
  Expression = `rna_seq_v2_mrna_median_all_sample_Zscores`; mutations =
  `..._mutations`; MSI = `MSI_SCORE_MANTIS` (> 0.4 → MSI-H); tissue =
  `CANCER_TYPE_DETAILED`.
- **Panel:** `../../panels/hallmark_200genes.gmt` (the repo's Hallmark panel —
  ~4,384 genes across 50 sets; disease-agnostic, chosen so the partition is not
  defined by any CRC-subtype marker).
- Wild-type = present in the `_sequenced` sample list but not in the gene's
  mutation records (so "wt" is a real negative, not "untested").
- Deterministic (seed 42, BIC k-selection). Live re-fetch may shift counts by a
  few tumors as cBioPortal updates.

**Research use only.** Not clinical guidance.
