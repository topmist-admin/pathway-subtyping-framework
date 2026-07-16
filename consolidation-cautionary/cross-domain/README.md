# Cross-domain expansion — Gate 6 confound remap (psychiatry + cancer)

Reproduction package for the **cross-disorder / cancer** expansion of the
"Stable but confounded" work. It demonstrates one claim: **the confound gate
(Gate 6) is domain-agnostic — only its *confound registry* changes, not the
code.** The autism reproduction (sibling `../` and `../genetic-anchoring/`) showed
Gate 6 catching a partition that was really a *brain-region* classifier; this
shows the same gate flag the domain-appropriate nuisance axis in cancer too
(tissue of origin; tumor purity).

Unlike the genetic-anchoring package, this one adds **no new framework code** on
the confound side — `confound_remap_demo.py` calls the shipped
`ValidationGates.confound_association_gate` unchanged. The genuine new code the
roadmap calls for (a somatic anchor for Gate 7) **is** now in the framework:
`ValidationGates.somatic_anchoring_gate` (see `../../src/pathway_subtyping/genetics/`).

---

## Layout

```
cross-domain/
  README.md                       <- you are here
  scripts/confound_remap_demo.py  <- runnable: the SAME Gate 6, three domains
  results/confound_remap_results.json  <- deposited reference output
  evidence/                       <- fetched provenance for the CRC target table (see caveats)
    crc_chembl_targets_drugs.json
    crc_subtype_stratified_trials.json
  docs/cross_disorder_cancer_roadmap.md  <- the design roadmap (psychiatry breadth + cancer)
```

## Reproduce

```bash
# needs the framework (any version with Gate 6) + numpy; deterministic (seeded)
python scripts/confound_remap_demo.py --out results
```

## Results → claims

| domain (scenario) | dominant nuisance | Cramér's V (planted align) | gate verdict | claim it supports |
|---|---|---|---|---|
| psychiatry, multi-region brain | brain region | **0.9329** (0.92) | FAIL (worst = brain region) | reproduces the autism/SCZ region-artifact case |
| cancer, pan-cancer pooled | tissue of origin | **0.8797** (0.90) | FAIL (worst = tissue of origin) | a pan-cancer partition recovers *organ*, as psychiatry recovered *anatomy* |
| cancer, single tumor type | tumor purity | **0.4841** (0.55) | FAIL (worst = tumor purity) | within a tumor type the nuisance is purity / stromal-immune fraction |

Identical gate, identical logic, remapped registry — diagnosis is exempt
(biology-of-interest) in every scenario. Numbers are deposited in
`results/confound_remap_results.json` and reproduce exactly from a fresh run.

The **confound remap** the roadmap prescribes (psychiatry → cancer):

| | psychiatry | cancer |
|---|---|---|
| dominant nuisance | brain region | tissue of origin; **tumor purity** |
| Gate-7 anchor | germline GWAS enrichment (low power) — `genetic_anchoring_gate` | **somatic driver / CNA / signature** (high power) — `somatic_anchoring_gate` |
| gate-exempt biology | diagnosis | cancer type or known subtype (CMS / PAM50) |

## Honesty notes (read before reusing the evidence)

- **No validated cancer-subtype recovery is claimed here.** Whether the battery
  *recovers* CMS/PAM50 on a valid, independent cohort is **Stage A** of the
  roadmap (a test to be run), not something this bundle demonstrates. The one TCGA
  row in the v0.7.0 47-dataset benchmark (TCGA-COAD) was flagged **`valid=False`**
  (degenerate ARI artifact + non-independent) in the 2026-07 correction, so no
  recovery statistic from it can be cited. An earlier roadmap draft cited a "TCGA
  CMS4 recovery at OR ≈ 16.5" — removed, because that number is not in the
  corrected benchmark and pointed at the invalidated row. See roadmap §2.1.
- **The `evidence/` fetch is partial.** `crc_chembl_targets_drugs.json` records
  live-fetch **errors** for the target lookups (`'str' object has no attribute
  'get'`) and a cetuximab timeout — only encorafenib (CHEMBL3301612) and
  pembrolizumab (CHEMBL3137343) mechanisms came through cleanly. The CRC
  target-nomination table is therefore illustrative, not a clean end-to-end fetch;
  re-run the ChEMBL/ClinicalTrials connectors before relying on it.
- **Trial counts are search-approximate.** The Phase-3 "subtype-stratified" counts
  (BRAF 4 / MSI-H 8 / RAS-wt 72) come from keyword search, not a curated
  stratified-trial set (e.g. the MSI-H list includes a DEBIRI liver-mets trial that
  is not obviously MSI-stratified). Treat the counts as order-of-magnitude.

## Scope

Feature-level (germline) and somatic subject-level anchoring are both shipped
(`pathway_subtyping.genetics`). **Germline** subject-level anchoring (rare-variant
burden / PRS on same-donor genotypes) remains access-gated (dbGaP/SFARI/Synapse)
and is not part of this package.

**Research use only.** Named drugs/trials are cited as development history for
methodological illustration, not treatment recommendations.
