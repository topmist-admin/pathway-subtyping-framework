# Public psychiatric meta-cohort enumeration (Track-A scoping)

Reproducible enumeration tooling that asks one question: **can a large-N
psychiatric postmortem-brain cohort be assembled from public, uniformly
reprocessed data alone** (no controlled access), to answer the reviewers'
small-N criticism (R1.5 / R3.5 / R3.6 / R3.9)?

**The answer is no.** This package exists to document that negative result
reproducibly, and to leave reusable tooling behind. No analysis was run on the
enumerated cohort — the enumeration itself was the deliverable, and it redirected
the strategy.

## Reproduce

```bash
# 1. Systematic GEO/SRA enumeration (NCBI E-utilities; set NCBI_API_KEY to raise the rate ceiling)
python scripts/enumerate_candidates.py      # -> results/candidate_studies.tsv
python scripts/resolve_srp.py               # fills missing GSE->SRP in place

# 2. Flag which SRA studies are actually in recount3 (needs R + recount3)
Rscript scripts/check_recount3.R            # -> results/track_a_recount3.tsv
```

Network-dependent and against live NCBI/recount3 indices, so counts can drift as
GEO grows. The deposited TSVs are the 2026-07-22 snapshot the finding is based on.

## Method

Systematic NCBI E-utilities search of GEO (`db=gds`) — not hand-picking:

```
(schizophrenia|bipolar|"major depress"|psychiatric)[Title] AND brain
  AND "Expression profiling by high throughput sequencing"[DataSet Type] AND human
```

Curation keeps **bulk postmortem RNA-seq only** — the DataSet-Type filter still
admits methylation, ChIP/ATAC, iPSC-derived, and single-nucleus studies, which
are dropped by title/summary keyword (`iPSC`, `hiPSC`, `organoid`, `scRNA`,
`snRNA`, `single-cell`, `single-nucleus`).

**57 candidate series → 30 bulk-RNA-seq keepers after curation.**

## FINDING — Track A is not large-N

`recount3` is a **fixed ~2021 release**, so uniform reprocessing (the thing that
would kill cross-study batch) is available only for pre-~2021 studies. Of the
curated keepers, only **6 studies / 243 samples** are in recount3:

| GSE | SRP | n (recount3) |
|---|---|---|
| GSE53239 | SRP033725 | 62 |
| GSE107638 | SRP126056 | 51 |
| GSE119290 | SRP159190 | 44 |
| GSE80336 | SRP073382 | 36 |
| GSE112523 | SRP136794 | 34 |
| GSE81396 | SRP074904 | 16 |
| | **subtotal** | **243** |

Adding the anchor cohort **GSE80655 / SRP073813 (n=352)** — which the
title-restricted search *missed*, and which is already the primary cohort of the
"stable but confounded" analysis — brings Track A to **~595 samples across ~7
studies**. Eight further curated keepers (e.g. GSE159487 n=150, GSE311978 n=97)
are **not** in recount3's 2021 release; see `results/track_a_recount3.tsv`.

**Conclusion:** Track A reaches roughly 2× the largest existing single cohort
(n=281), not the thousands the small-N criticism calls for — and it does so by
pooling seven small studies that differ in brain region, disorder, and design.
For a paper whose thesis is that subtype claims are frequently confound
artifacts, manufacturing a maximally heterogeneous pooled cohort is the wrong
instrument: cross-study batch would plausibly dominate any structure found.

**Strategy shift:** a single well-characterized large cohort beats aggregating
seven tiny heterogeneous ones. **CommonMind (n=986)** via controlled access is
dramatically better and became the priority. Track B (direct GEO download of
post-2021 studies) adds N but with mixed pipelines and the same heterogeneity —
low value here.

Track A is retained only as a possible modest *supplementary public replication*;
the enumeration tooling is reusable for any future public-cohort scoping.

## Honesty note

This is a negative result reported as such. The `~595` figure is the honest
ceiling (243 enumerated + 352 anchor), not an achieved cohort — nothing was
pulled, harmonized, or analyzed. Phenotype harmonization across these studies
(free-text `colData` → common disease/region/age/sex/PMI/RIN schema) would be the
substantive work if Track A were ever pursued, and it was not done.
