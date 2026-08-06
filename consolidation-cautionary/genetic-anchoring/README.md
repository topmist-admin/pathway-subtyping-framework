# Genetic anchoring extension — Gate 6 + Gate 7

Reproduction package for the **genetic-anchoring** results that extend the paper
*"Stable but confounded: brain region, not diagnosis, drives a validation-passing
molecular subtype in postmortem psychiatric brain"* (Chauhan 2026, Research Square
[`10.21203/rs.3.rs-8913089`](https://doi.org/10.21203/rs.3.rs-8913089), under
revision). It is the sequel to the parent harness one level up
([`../`](..)): where the parent shows a stable partition can be a **region
artifact**, this package shows (a) the artifact is now **caught** by the confound
gate, and (b) the framework can **recover a known genetic signal** (the Gate 7
calibration) — the external, causally-prior check the original battery lacked.

Both gates are now **first-class framework code**, not bolt-on scripts:
`ValidationGates.confound_association_gate` (Gate 6) and
`ValidationGates.genetic_anchoring_gate` (Gate 7, backed by
`pathway_subtyping.genetics`). The drivers here only do I/O and call the package —
they exercise exactly what a third party installs.

---

## Layout

```
genetic-anchoring/
  README.md            <- you are here (results -> figures/tables map)
  provenance.json      <- input SHA-256, versions, seeds, result-file map, headline numbers
  scripts/
    run_gate6_confound.py    <- driver -> ValidationGates.confound_association_gate
    run_gate7_anchoring.py   <- driver -> pathway_subtyping.genetics + genetic_anchoring_gate
    fetch_asd_risk_genes.py  <- (re)build the risk set from Open Targets  [network]
    fetch_universes.py       <- (re)build the genome-wide null + symbol map from BioMart [network]
  inputs/              <- deposited fetched inputs, so the gates run OFFLINE
    asd_risk_ensembl.json         (1,020 ASD risk genes; Open Targets MONDO_0005258)
    opentargets_asd_targets.json  (full Open Targets pull, for audit)
    protein_coding_ensembl.json   (genome-wide reference null, N=23,272)
    testset_symbol2ens.json       (neuronal+glial symbol -> Ensembl map)
  results/             <- deposited reference outputs (a correct run must match these)
    gate6_confound_association.json
    gate7_genetic_anchoring_poscontrol.json
```

Two inputs are **not** vendored here (they live with the parent harness / are
large):
- `../panels/schizophrenia_pathways.gmt` — the cell-type/pathway panels (parent dir).
- `GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz` — the brain-expressed
  universe is derived from it at run time; fetch from GEO
  (`GSE80655`, standard supplementary download) and pass `--gse80655-expr`.

---

## Reproduce

Environment: the published `pathway-subtyping` release (**v0.8.0 line**, for the
`genetics` subpackage) + `numpy pandas scipy scikit-learn`. See the parent
harness `../requirements.txt` for the reproduction-sensitive pins (scikit-learn
in particular).

```bash
# Gate 6 — confound gate on the deposited partition (needs the parent partition CSV)
python scripts/run_gate6_confound.py \
  --labels ../data/partition/sample_metadata_with_subtypes.csv

# Gate 7 — genetic anchoring positive control (offline; uses ./inputs)
python scripts/run_gate7_anchoring.py \
  --gse80655-expr /path/to/GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz
```

To rebuild the deposited inputs from source instead of using `./inputs`:

```bash
python scripts/fetch_asd_risk_genes.py --out inputs           # Open Targets
python scripts/fetch_universes.py --gmt ../panels/schizophrenia_pathways.gmt --out inputs  # BioMart
```

---

## Results → figures / tables

Headline numbers are also in [`provenance.json`](provenance.json). "Paper element"
points at where each result lands in the **next version** of the manuscript; the
Layer A/B figures it builds on are reproduced by the parent harness (`../scripts`,
`../figures`).

| Result file | Metric (deposited) | What it establishes | Paper element |
|---|---|---|---|
| [`results/gate6_confound_association.json`](results/gate6_confound_association.json) | `passed=false`; max **nuisance Cramér's V = 0.66** on `brain region` (χ²=125.1, p≈4e-26); **diagnosis exempt** (V=0.0, p=0.408) | The same partition that passed the original stability battery (bootstrap ARI ≈0.914) is **caught as a region classifier** by the confound gate; diagnosis is treated as biology-of-interest and never fails it. | Methods/Results — the new **Gate 6**; ties directly to the Fig 1 region crosstab (parent Layer A: χ²=125.12, region V=0.666, diagnosis p=0.408). |
| [`results/gate7_genetic_anchoring_poscontrol.json`](results/gate7_genetic_anchoring_poscontrol.json) | **NEURONAL/SYNAPTIC** fold **16.51×, p=1.04e-23**; **GLIAL/IMMUNE** fold 2.65×, p=0.17 (n.s.) — under the brain-expressed-matched null | **Discrimination confirmed**: the pipeline recovers Voineagu (2011) F5 from public data — a neuronal axis carries ASD genetic risk, a glial axis does not. This is the §0.5 ground-truth calibration the original battery lacked; passing it is the precondition for trusting Gate 7 on any *novel* subtype. | New figure/table — the **Gate 7** genetic-anchoring calibration. |
| (same file, both nulls) | neuronal fold **10.14× genome-wide → 16.51× brain-matched** | ⚠️ **CORRECTED — the null matters, but NOT in the direction stated here.** The brain-expressed universe is not intersected with the protein-coding set, so it is *larger* (36,971 vs 23,272) and risk-gene density *falls* (4.38% → 2.69%); the 10.14×→16.51× rise is a denominator effect, not concentration. Do not cite 16.51× as a matched-background enrichment. The direction of the finding (neuronal significant, glial not) is unaffected. ~~Risk genes concentrate in the brain-expressed background, so a genome-wide null *understates* true enrichment for a real neuronal axis — and would *overstate* it for a region-identity artifact. This is why the matched null is the default and decides the gate. | Methods — null-choice justification for Gate 7. |

---

## Scope — feature-level only (read before extending the claim)

- **Feature-level anchoring (what is here): executable now, zero data-use
  agreement.** Tests whether a subtype's *defining genes* are enriched for ASD
  genetic risk. Confound-immune (germline variants are upstream of every
  postmortem/technical confound) but **necessary, not sufficient** — the claim
  ceiling is *"genetically-implicated subtype axis,"* not *"aetiologically
  distinct donors."* It is a **specific, low-power** test: a null is weak evidence
  of absence, never proof.
- **Subject-level anchoring (not here): access-gated.** Tests whether subtype
  *donors* carry differential rare-variant burden / polygenic score. Needs
  transcriptome + genotype on the **same donors** (PsychENCODE / CommonMind /
  SPARK) behind dbGaP / SFARI Base / Synapse data-use agreements. The barrier is
  structural, not temporal. This is a later framework addition, not part of this
  package. See `../../docs/` roadmap material.

## Reproducibility notes

- The **deposited reference outputs** were produced with `pathway-subtyping 0.7.0`
  plus the Gate-7 prototype now folded into the package; the integrated framework
  gate reproduces them to the decimal (neuronal 16.51× / p=1.04e-23).
- Open Targets and Ensembl are **live services** — a re-fetch months later can
  shift the risk-gene count slightly as evidence accrues; the deposited
  `inputs/` freeze the exact set used, and `provenance.json` carries their
  SHA-256. Prefer the offline path for exact reproduction.
- Gate 6 needs the parent harness's deposited partition
  (`../data/partition/sample_metadata_with_subtypes.csv`). Use the 3-level
  `brain region` column (AnCg/DLPFC/nAcc), not a 2-level coding — see the parent
  README's region-coding trap.

## Citation / integrity

- Framework: `pathway-subtyping`, RRID:SCR_028051 (concept DOI
  10.5281/zenodo.18638048; v0.7.0 DOI 10.5281/zenodo.21279842).
- Risk genes: Open Targets Platform (autism spectrum disorder, MONDO_0005258,
  `genetic_association` evidence).
- Manuscript under extension: Chauhan 2026, Research Square
  `10.21203/rs.3.rs-8913089` (under revision).

**Research use only.** Nothing here is clinical guidance.
