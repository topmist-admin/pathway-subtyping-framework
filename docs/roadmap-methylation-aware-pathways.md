# PSF Roadmap — Methylation-Aware Pathway Feature Module

**Planned:** April 29, 2026
**Status:** **Prioritized backlog** — pre-implementation, not yet scheduled for a specific release.
**Target release window:** v0.7.x or v0.8.x
**Public/private posture:** Public Codeberg track. MIT-licensed alongside the rest of the framework.
**Author:** Rohit Chauhan

---

## Why This Note Exists

PSF v0.6 ingests gene expression (RNA-seq, microarray) and scores 50–100 biological pathways via ssGSEA, then discovers patient subtypes via Gaussian mixture modeling with five mandatory validation gates. **This pipeline is transcriptomics-only.** That restriction has been acceptable through v0.5–v0.6 because every validation case (sepsis GSE65682, autism, schizophrenia, colorectal cancer cross-validation) is RNA-driven.

That single-modality restriction is a real limitation:

- Integrated methylome-plus-variant assays (e.g. Illumina's 5-Base DNA Prep) produce methylation and high-accuracy variant calls from one library, and PSF cannot ingest the methylation channel today.
- Cancer epigenetics in particular is an entire analysis surface (methylation atlases, HNSCC and other methylation-rich cancers, TCGA reanalysis tracks) where a transcriptomics-only tool either competes weakly or cannot participate.
- Methylation is a *second genomics modality* that does not commit PSF to any single disease — it broadens the tool's applicability rather than narrowing it.

This roadmap formalizes the methylation-aware pathway feature module as a **prioritized backlog item** with clear architecture, dependencies, and validation strategy.

---

## Strategic Rationale

PSF's pathway-subtyping value proposition rests on three claims:

1. Pathway-level features are more biologically interpretable than individual genes.
2. Cross-platform reproducibility is achievable when scoring at the pathway level (already demonstrated microarray ↔ RNA-seq).
3. Statistical validation gates ensure subtype calls survive replication (the five mandatory gates).

**Methylation is the natural next modality** because:

- Promoter methylation is mechanistically linked to gene expression repression — there is a biological path from methylation pattern to pathway activity that pathway-level scoring is uniquely positioned to capture.
- The same ssGSEA + GMM + validation-gates machinery that powers the transcriptomic path can be reused with a methylation-scored feature vector. The architecture extends, not replaces.
- Cross-platform reproducibility translates: 450k arrays, EPIC, EPIC v2, and bisulfite/5-base sequencing all produce CpG-level methylation calls. PSF's existing reproducibility story extends to methylation data.

---

## Architectural Approach

**Design principle:** Extend, don't replace. Existing transcriptomic pipeline stays intact; methylation is a parallel data path that produces an ssGSEA-compatible feature vector at the same shape and column convention as the existing pathway-score matrix.

```
                                      ┌──────────────────────────────────┐
   RNA-seq / microarray ──ssGSEA──►   │  Pathway × Sample score matrix   │
                                      │  (existing v0.6 path)            │
                                      └────────────────┬─────────────────┘
                                                       │
                                                       ▼
                                         Existing GMM + 5-gate validation
                                                       ▲
                                                       │
                                      ┌────────────────┴─────────────────┐
   450k / EPIC / 5-Base ─►            │  Pathway × Sample methylation    │
   CpG-level methylation ─►           │  score matrix (NEW path)         │
   (β-values or M-values)             └──────────────────────────────────┘
              │
              │
              ▼
    CpG → gene aggregation
    (promoter / gene body / TFBS-weighted)
              │
              ▼
    Gene-level methylation matrix
              │
              ▼
    Methylation-aware pathway scoring
    (mGSEA — methylation analog of ssGSEA;
     directionality flipped: hypermethylation
     of promoters → reduced pathway activity)
```

Two operating modes both supported:

- **Parallel modality:** methylation pathway scores join transcriptomic pathway scores at the GMM step → joint multi-omic subtyping.
- **Single modality:** methylation pathway scores feed GMM directly → methylation-only subtyping (use case: cohorts with methylation but no matched RNA).

---

## Features

### F1: CpG-to-Gene Aggregation

**Module:** `pathway_subtyping.methylation.aggregation`
**Priority:** HIGH — foundational for everything downstream.
**Depends on:** New module + reference annotation files (Illumina 450k / EPIC / EPIC v2 manifests, GENCODE gene models).

**Problem:** Methylation arrays and bisulfite/5-Base sequencing produce CpG-level (probe-level or position-level) β-values or methylation calls. PSF's pathway scoring operates at the gene level. We need a principled, configurable, well-tested aggregation step.

**Design:** Three aggregation strategies, user-selectable:

1. **Promoter mean:** average β across CpGs in TSS ± 2kb. Default. Strongest gene-expression linkage.
2. **Gene body mean:** average β across CpGs in gene body. Useful for active-gene methylation patterns (paradoxically positive correlation with expression for some genes).
3. **TFBS-weighted:** weight CpGs by transcription-factor binding site overlap. Higher resolution; depends on KG annotation.

**Deliverables:**
- `aggregate_cpg_to_gene(meth_matrix, manifest, mode="promoter_mean")` API
- Manifest loaders for Illumina 450k, EPIC, EPIC v2; bisulfite-call BED ingest
- Unit tests across all three modes + edge cases (CpG islands, shores, shelves)
- Smoke benchmark: TCGA-HNSC methylation → gene matrix → known-gene methylation profile reproduction

### F2: Methylation-Aware Pathway Scoring (mGSEA)

**Module:** `pathway_subtyping.methylation.pathway_scoring`
**Priority:** HIGH — this is the actual pathway scoring step.
**Depends on:** F1 + existing pathway gene set definitions (MSigDB, Reactome, KEGG already loaded by v0.6).

**Problem:** ssGSEA assumes higher input value = higher pathway activity. Methylation flips that for promoter-aggregated input — *higher* β = repression. The scoring step must invert directionality cleanly without breaking the pathway-set machinery.

**Design:**

```
mGSEA(meth_gene_matrix, pathway_gene_sets, direction="repressive"):
    - For each pathway, score the gene set against the methylation rank distribution
    - direction="repressive": invert β before scoring (high β → low score)
    - direction="activating": pass β through (rare; for gene-body methylation use cases)
    - Output: pathway × sample score matrix, same shape/convention as existing ssGSEA output
```

**Deliverables:**
- `mgsea(meth_gene_matrix, gene_sets, direction)` API
- Joint API `score_pathways(modality="rna"|"methylation"|"joint", ...)` that routes to ssGSEA / mGSEA / both
- Validation: TCGA-HNSC methylation-derived MAPK / WNT / immune pathway scores must correlate with expression-derived scores at biologically expected sign (negative for promoter mode)

### F3: 5-Base / Bisulfite Sequencing Ingest

**Module:** `pathway_subtyping.methylation.ingest`
**Priority:** MEDIUM — required for the sequencing-based methylation use case but not for any array-based backlog work.
**Depends on:** F1 (aggregation re-uses BED-style position input).

**Problem:** Sequencing-derived methylation calls come in BED-like format (chrom, start, end, β, coverage) rather than array probe IDs. PSF needs a clean ingest path that produces the same intermediate as array data.

**Design:**
- BED ingest with coverage filter (default ≥10x)
- Position → manifest probe-ID matching (for sites overlapping array positions)
- Position → genomic feature mapping (for non-array CpGs) to reuse F1's annotation logic

**Deliverables:**
- `ingest_5base(bed_file, coverage_min=10)` returning the same gene-level matrix shape as F1's output
- Round-trip test: TCGA-HNSC methylation array → simulated 5-Base BED → ingest → identical pathway scores within tolerance

### F4: Multi-Omic GMM (joint clustering)

**Module:** `pathway_subtyping.clustering.multiomic`
**Priority:** MEDIUM — strong scientific selling point but not blocking for single-modality subtyping.
**Depends on:** F2 + existing GMM machinery in v0.6.

**Problem:** When transcriptomic and methylation data are both available for a cohort, joint multi-omic subtyping is the highest-impact use case. Naïve concatenation of feature vectors over-weights whichever modality has more pathways. Need feature-balanced or modality-weighted clustering.

**Design:**
- Feature-balanced concatenation (z-score within modality before joining)
- Optional MOFA-style latent-factor pre-step for very high-dimensional cohorts
- Existing 5-gate validation runs unchanged on the joint feature space

**Deliverables:**
- `joint_gmm(rna_scores, meth_scores, weighting="balanced"|"modality_weighted")` API
- Cross-validation test: TCGA-HNSC joint subtypes must recover known basal/mesenchymal/atypical/classical with ARI ≥ 0.7

### F5: Validation Suite

**Module:** `tests/integration/test_methylation_validation.py` + benchmark notebook
**Priority:** HIGH — modules that ship without validation are a liability.

**Validation cohorts (public, no new data needed):**
- TCGA-HNSC (methylation 450k + RNA-seq matched, n ≈ 530)
- TCGA-COAD (matched cohort, n ≈ 460)
- TCGA-BRCA (matched cohort, n ≈ 870; PAM50 reproduction is gold-standard test)

**Validation criteria:**
1. F1 + F2 single-modality methylation subtypes recover known molecular subtypes at ARI ≥ 0.5 on at least 2 of the 3 TCGA cohorts
2. F4 joint subtypes outperform single-modality on at least 1 cohort (ARI improvement ≥ 0.1)
3. Cross-platform: 450k → EPIC reproduction at within-cohort cluster stability ≥ 0.85
4. Cross-tech: array vs simulated-5-Base reproduction at within-cohort cluster stability ≥ 0.80

**Deliverables:**
- `BENCHMARK_METHYLATION_TCGA.ipynb` — paper-ready validation
- Continuous-integration test suite covering F1–F4

---

## Dependencies and Out-of-Scope

### New runtime dependencies

- `methylprep` or equivalent (Illumina array manifest parsing)
- `pyBigWig` or `pybedtools` (BED ingest for F3)
- Existing PSF stack handles everything else (numpy, pandas, scikit-learn, scipy)

### Reference annotations to vendor or download

- Illumina 450k manifest, EPIC, EPIC v2 manifests (vendor-provided, redistributable)
- GENCODE gene model (for promoter / gene body coordinates)

### Out of scope (do NOT take on in this module)

- Raw bisulfite / 5-Base sequencing pipeline (alignment, methylation calling, DMR detection) — PSF accepts methylation-call BED as input. Upstream pipelines (Bismark, methylKit, ChAMP) are the user's responsibility.
- Allele-specific methylation
- Histone-mark integration (H3K27me3, H3K4me3, etc.) — separate roadmap if ever pursued
- Single-cell methylation — keep focused on bulk first

---

## Prioritization

This is a **prioritized backlog item**, not active work. The commitments it represents:

1. **The module is in the public roadmap** — its architecture and scope are decided in the open.
2. **Architecture is decided** — F1–F5 are the work, in this order, with these deliverables. No re-architecting from scratch under deadline pressure.
3. **Validation cohorts are pre-selected** — TCGA-HNSC / COAD / BRCA. No scrambling for benchmark data later.
4. **Backstop scheduled** — this lands in a v0.7.x or v0.8.x release, and F1–F2–F3 (the ingest + scoring core) can be pulled forward as a scoped sprint if a matched methylation cohort becomes available sooner.

---

## Related Files

- [roadmap-molecular-qc-v05.md](roadmap-molecular-qc-v05.md) — v0.5 QC roadmap (retrospective; shipped) — pattern this file follows
- [roadmap-v06-codeberg.md](roadmap-v06-codeberg.md) — v0.6 public roadmap (retrospective)

---

**Last updated:** April 29, 2026
