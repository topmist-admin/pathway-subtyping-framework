# GTEx brain at scale (n=2,931) — large-N region structure + gate behavior

Public (recount3, no controlled access needed) large-N analysis: 2,931 GTEx brain
samples across 13 subregions, Hallmark pathway scores. Fetch:
`scripts/fetch_gtex_brain_recount3.R`;
analysis: `scripts/analyze_gtex_brain.py`.

**Motivation:** proposed as the large-N rebuttal of the reviewers' small-N
criticism (region-confound at scale). **The result is more nuanced than that
and is reported honestly — it does NOT cleanly deliver "region dominates."**

## Result

| Measure | Value |
|---|---|
| 13-way region recovery (ARI) | **0.151** (AMI 0.27) — WEAK |
| best single-region enrichment (k=13) | Cerebellum 58.9% (OR 22, p 9e-77) |
| region recovery @ BIC k=3 (ARI) | 0.074 |
| best enrichment @ k=3 | Spinal cord 24.7% (OR 147, p 9e-106) |
| discreteness gate @ k=13 | **CERTIFIED** discrete structure |
| discreteness gate @ BIC k=3 | **REJECTED** — "no discrete structure (continuum)" |

## Honest interpretation (do NOT overclaim)
- **Brain region is NOT a clean 13-way discrete axis in Hallmark-pathway space.**
  Generic pathways weakly recover the 13 regions (ARI 0.15). This is partly a
  PANEL limitation: regional identity is driven by specific neuronal/glial genes,
  not cancer-flavored Hallmark sets. A region-specific panel would likely recover
  region far better; do not read ARI 0.15 as "region has no structure."
- **The dominant discrete structure is cerebellum/hindbrain vs forebrain.**
  Cerebellum (OR 22) and spinal cord (OR 147) are molecularly distinct; the
  cortical/subcortical regions psychiatric studies actually use form more of a
  CONTINUUM. The gate captured exactly this: fine (k=13) structure is discrete,
  the coarse (k=3) split is a continuum.
- **This does NOT cleanly rebut the small-N criticism.** The projected
  "region-confound at scale" claim came out weaker and more caveated. Report it
  as nuance, not a headline. The strong region signal (cerebellum) is also
  typically EXCLUDED from psychiatric postmortem cohorts.

## What genuinely holds
1. **Large-N calibration point.** At n=2,931 the discreteness gate still
   discriminates — certifies genuine fine structure, rejects a coarse continuum.
   Extends the `gate_calibration/` result to a much larger cohort.
2. **Metric-dependence recurs.** Weak ARI (0.15) yet strong single-region
   enrichment (OR 22–147) — the same pattern as the cancer cohorts (see
   `cancer_r38/`). Reinforces the meta-finding that enrichment flatters while ARI
   is strict.

## Implication for the rebuild
GTEx brain is a supporting calibration/metric-dependence result, NOT the large-N
psychiatric headline. The reviewers' small-N criticism is better answered by a
large-N DISEASE cohort (CommonMind n=986 via controlled access, or an aggregated
psychiatric SRA meta-cohort via recount3). A region-specific (not Hallmark) panel on
GTEx would strengthen the region-axis point if pursued.

## Provenance
recount3 GTEx BRAIN (uniform reprocessing), Hallmark 50-set panel, gene-wise
z-score → pathway means. Deterministic (seed 42), Gate A n_ref=60.
`results/gtex_brain_region_confound.json`.
