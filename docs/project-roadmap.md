# Pathway Subtyping Framework: v0.1 → v0.4 Roadmap

This document outlines the project milestones for the Pathway Subtyping Framework,
modeled after the successful Autism Pathway Framework 90-day plan.

## Phase 1: v0.1 Completion (Weeks 1-2)

### [Week 1] Zenodo DOI + GitHub Release

**Theme:** Finalize v0.1 with citable DOI

**Deliverables:**
- [ ] Link GitHub repo to Zenodo
- [ ] Configure Zenodo metadata: title, description, keywords, communities
- [ ] Create GitHub release `v0.1.0` with release notes
- [ ] Zenodo DOI minted and added to README

**Acceptance Criteria:**
- GitHub release page exists with clear notes
- Zenodo DOI resolves correctly
- README badges include DOI

---

### [Week 2] Collaborator Outreach Package

**Theme:** External collaborator adoption materials

**Deliverables:**
- [ ] One-pager PDF/MD for genomics labs (`docs/one-pager.md` → PDF)
- [ ] Email template for lab outreach
- [ ] "How to run on your cohort" appendix
- [ ] Target list: 5-10 neurogenetics labs/consortia
- [ ] Send to at least 1 external collaborator

**Acceptance Criteria:**
- Outreach bundle complete and peer-reviewed
- At least 1 introductory email sent
- Follow-up scheduled for +7 days

---

## Phase 2: v0.2 Development (Weeks 3-10)

### [Week 3] User Feedback Integration

**Theme:** Incorporate early adopter feedback

**Deliverables:**
- [ ] Collect feedback from Week 2 outreach
- [ ] Prioritize top 3 user-requested features
- [ ] Create v0.2 scope document
- [ ] GitHub Issues created for v0.2 features

**Acceptance Criteria:**
- At least 1 user feedback session completed
- v0.2 scope approved
- Backlog prioritized

---

### [Week 4] Real-World Data Support

**Theme:** Handle messy real-world VCFs

**Deliverables:**
- [ ] Support for multi-allelic variants
- [ ] Graceful handling of missing annotations
- [ ] VEP/ANNOVAR annotation helper validated on real data
- [ ] Error messages that guide users to fixes

**Acceptance Criteria:**
- Pipeline runs on user-provided VCF without crashing
- Warnings logged for data quality issues
- Documentation updated with troubleshooting

---

### [Week 5] Additional Disease Pathways

**Theme:** Expand pathway library

**Deliverables:**
- [ ] Validate schizophrenia pathways with literature
- [ ] Validate epilepsy pathways with literature
- [ ] Validate intellectual disability pathways with literature
- [ ] Add Parkinson's disease pathways
- [ ] Add bipolar disorder pathways

**Acceptance Criteria:**
- Each pathway file has literature citations
- GMT files validated with parser
- Example configs for each disease

---

### [Week 6] Subtype Characterization Tools

**Theme:** Help users interpret subtypes

**Deliverables:**
- [ ] Subtype pathway enrichment analysis
- [ ] Gene-level contribution scores
- [ ] Subtype comparison heatmaps
- [ ] Export to standard formats (CSV, Excel)

**Acceptance Criteria:**
- Users can identify which pathways drive each subtype
- Visualizations are publication-ready
- Outputs documented in dictionary

---

### [Week 7] Cross-Cohort Validation

**Theme:** Compare subtypes across datasets

**Deliverables:**
- [x] Document cross-cohort validation workflow
- [x] Add example with synthetic cohorts
- [x] Validation report template
- [x] Guidance on interpreting transfer ARI

**Acceptance Criteria:**
- Cross-cohort tutorial in notebook
- Clear criteria for "replicable subtypes"
- Example report generated

**Status: COMPLETED** (Feb 14, 2026) — Added `to_dict()`, `format_report()`, `get_citations()` to dataclasses, `generate_synthetic_cohort_pair()`, CLI example script, user guide, API reference. 12 new tests.

---

### [Week 8] Performance Optimization

**Theme:** Scale to 10,000+ samples

**Deliverables:**
- [x] Benchmark on large synthetic dataset
- [x] Optimize memory usage for chunked processing
- [x] Add progress bars for long operations
- [x] Document recommended hardware specs

**Acceptance Criteria:**
- 10,000 sample cohort processes in <30 minutes
- Memory stays under 8GB
- Clear guidance on resource requirements

**Status: COMPLETED** (Feb 14, 2026) — Added tqdm progress bars to 4 modules, chunked VCF processing in pipeline, benchmark script, hardware guide. 9 new tests. 620 total tests passing.

---

### [Week 9] Advanced Visualization

**Theme:** Interactive exploration

**Deliverables:**
- [x] Interactive HTML report with Plotly
- [x] UMAP/t-SNE visualization option
- [x] Subtype trajectory visualization
- [x] Export figures in multiple formats

**Acceptance Criteria:**
- HTML report opens in browser
- Figures are customizable
- Publication-quality exports

**Status: COMPLETED** (Feb 14, 2026) — New `visualization.py` module with Plotly interactive HTML reports, PCA/t-SNE/UMAP scatter plots, pathway heatmaps, radar charts, and multi-format export (PNG, SVG, PDF, HTML). Optional `[viz]` dependency group. Static matplotlib fallbacks. Pipeline integration via `generate_interactive_report` config. 53 new tests. 716 total tests passing.

---

### [Week 10] v0.2 Release

**Theme:** Ship v0.2

**Deliverables:**
- [ ] Update CHANGELOG.md for v0.2
- [ ] Create GitHub release `v0.2.0`
- [ ] Update PyPI package
- [ ] Update Zenodo DOI
- [ ] Announce on relevant channels

**Acceptance Criteria:**
- All v0.2 features merged
- Tests passing on CI
- Documentation current
- Announcement posted

---

## Phase 3: Community Growth (Weeks 11-13)

### [Week 11] Tutorial Series

**Theme:** Educational content

**Deliverables:**
- [ ] Video walkthrough: Getting Started (5-10 min)
- [ ] Tutorial: Adapting for Your Disease
- [ ] Tutorial: Interpreting Validation Gates
- [ ] Blog post: "Pathway-Based Subtyping Explained"

**Acceptance Criteria:**
- Videos hosted on YouTube/platform
- Tutorials linked from README
- Blog post published

---

### [Week 12] Community Contributions

**Theme:** Enable external contributors

**Deliverables:**
- [ ] Create "Good First Issues" labels
- [ ] Document contribution workflow
- [ ] Set up GitHub Discussions
- [ ] Respond to community questions within 48h

**Acceptance Criteria:**
- At least 3 "Good First Issues" open
- First external contribution merged
- Discussion forum active

---

### [Week 13] Conference/Publication Prep

**Theme:** Academic dissemination

**Deliverables:**
- [ ] Draft abstract for conference submission
- [ ] Identify target conferences/journals
- [ ] Collect preliminary results from collaborators
- [ ] Prepare poster/presentation materials

**Acceptance Criteria:**
- Abstract submitted to at least 1 venue
- Collaborator data collected (with consent)
- Materials ready for presentation

---

## Phase 4: v0.3 Multi-Omic Integration (Weeks 14-21)

**Community mandate:** 67% poll vote for multi-omic support. This phase transforms the framework
from a genomic-variant-only tool into a multi-modal transcriptomics and genomics platform.

### [Week 14] Bulk RNA-seq Pathway Scoring

**Theme:** Extend pathway aggregation engine to expression data

**Deliverables:**
- [ ] Accept gene expression matrix (counts/TPM) as input alongside VCF
- [ ] Pathway scoring from expression data (GSVA, ssGSEA, or mean-Z methods)
- [ ] Unified pathway score matrix regardless of input modality
- [ ] Tests with synthetic expression data

**Acceptance Criteria:**
- `pipeline.run()` accepts `--input-type expression` flag
- Pathway scores from expression data pass existing validation gates
- Documentation updated with expression workflow

**GitHub Issue:** #29

---

### [Week 15] Single-Cell Pathway Scoring

**Theme:** Per-cell and per-cell-type pathway score computation

**Deliverables:**
- [ ] Accept AnnData (h5ad) or cell-by-gene matrix as input
- [ ] Pathway scoring per cell and aggregated per cell-type
- [ ] Cell-type-level pathway score matrix compatible with existing GMM clustering
- [ ] Handle sparsity and dropout common in scRNA-seq

**Acceptance Criteria:**
- Framework clusters cell types into pathway-defined subtypes
- Works on standard 10X Genomics-format data
- Memory-efficient for datasets up to 50K cells

**GitHub Issue:** #30

---

### [Week 16] Bulk Deconvolution Integration

**Theme:** Estimate cell-type proportions from bulk RNA-seq using single-cell references

**Deliverables:**
- [ ] Implement or wrap deconvolution method (MuSiC or NNLS-based)
- [ ] Add estimated cell-type proportions as features alongside pathway scores
- [ ] Combined feature matrix (pathway scores + cell-type fractions) → GMM clustering
- [ ] Subtypes become both pathway-defined AND cell-type-aware

**Acceptance Criteria:**
- Deconvolution runs with synthetic bulk + reference single-cell data
- Combined clustering produces more biologically interpretable subtypes
- Validation gate results reported for combined vs pathway-only features

**GitHub Issue:** #31

---

### [Week 17] Signaling Pathway Databases

**Theme:** Add intercellular communication pathway sources

**Deliverables:**
- [ ] Import ligand-receptor interaction pathways from CellChatDB or CellPhoneDB
- [ ] Score signaling pathways using existing aggregation engine
- [ ] New pathway category: `signaling` alongside `metabolic`, `regulatory`, etc.
- [ ] Example: signaling pathway scores reveal immune-neuronal crosstalk subtypes

**Acceptance Criteria:**
- At least 200 signaling pathway gene sets loaded
- Signaling scores pass validation gates
- Documentation with biological interpretation examples

**GitHub Issue:** #32

---

### [Week 18] Cross-Modal Validation Gate

**Theme:** Validation gate #5 — do subtypes replicate across data modalities?

**Deliverables:**
- [ ] Cross-modal concordance test: subtypes from bulk expression vs subtypes from variants
- [ ] Single-cell validation: do bulk-defined subtypes correspond to distinct cell-type compositions?
- [ ] Concordance metrics (ARI, NMI) across modalities
- [ ] Gate pass/fail threshold calibrated on synthetic multi-modal data

**Acceptance Criteria:**
- New validation gate integrated into pipeline report
- Gate detects when subtypes are modality-specific artifacts vs biologically real
- Methods documentation updated

**GitHub Issue:** #33

---

### [Week 19] Spatial Transcriptomics Support (Experimental)

**Theme:** Pathway scores per spatial spot for tissue-level subtype maps

**Deliverables:**
- [ ] Accept Visium/spatial expression matrix with coordinates
- [ ] Pathway scoring per spatial spot
- [ ] Spatial subtype visualization (colored tissue maps)
- [ ] Identify spatially coherent subtype regions (spatial autocorrelation)

**Acceptance Criteria:**
- Works on standard 10X Visium format
- Spatial subtype maps are visually interpretable
- Marked as experimental/beta in documentation

**GitHub Issue:** #34

---

### [Week 20] Multi-Omic Pipeline Integration

**Theme:** End-to-end pipeline combining all modalities

**Deliverables:**
- [ ] Unified config accepting multiple input types (VCF, expression, scRNA-seq, spatial)
- [ ] Pipeline auto-detects available modalities and builds combined feature matrix
- [ ] Nextflow/Snakemake templates updated for multi-omic workflows
- [ ] Colab notebook: multi-omic quick demo

**Acceptance Criteria:**
- Single `pathway-subtyping run --config multi_omic.yaml` processes all available data
- Report includes per-modality and combined results
- Colab demo runs in <5 minutes

**GitHub Issue:** #35

---

### [Week 21] v0.3 Release

**Theme:** Ship v0.3

**Deliverables:**
- [ ] Update CHANGELOG.md for v0.3
- [ ] Create GitHub release `v0.3.0`
- [ ] Update PyPI package
- [ ] Update Zenodo DOI
- [ ] Announce on LinkedIn, Substack, Twitter/X

**Acceptance Criteria:**
- All v0.3 features merged and tested
- Multi-omic documentation complete
- Announcement posted with Colab demo link

**GitHub Issue:** #36

---

## Phase 5: v0.4 Deep Multi-Omics & Metabolic Disease (Weeks 22-29)

**Rationale:** Phase 4 adds transcriptomics and spatial data. Phase 5 completes the multi-omics vision by adding proteomics, metabolomics, and advanced integration methods — enabling the framework to handle all five major data modalities (DNA, RNA, protein, metabolite, functional screens). This phase also adds domain-specific pathway databases for metabolic and mitochondrial diseases.

### [Week 22] MitoCarta 3.0 Mitochondrial Pathway Database

**Theme:** First-class support for mitochondrial disease research via the gold-standard MitoCarta 3.0 database

**Deliverables:**
- [ ] Parser for MitoCarta 3.0 MitoPathways (149 hierarchical pathways, 1,136 genes)
- [ ] GMT export from 3-tier hierarchy (top-level → sub-pathway → individual complexes)
- [ ] Pre-built `mitochondria_pathways.gmt` (OXPHOS I-V, TCA, ISC biogenesis, CoQ, mt-translation, dynamics)
- [ ] Sub-organelle localization annotations (matrix, inner membrane, intermembrane space, outer membrane)
- [ ] Example config and documentation

**Acceptance Criteria:**
- MitoCarta GMT files load in existing pipeline
- At least 100 pathway gene sets from MitoCarta hierarchy
- Documentation includes mitochondrial disease use case

**GitHub Issue:** #37

---

### [Week 23] Proteomics Pathway Scoring Module

**Theme:** Protein abundance as a pathway scoring input modality

**Deliverables:**
- [ ] Accept protein abundance matrix (label-free, TMT, DIA) as input type
- [ ] Pathway scoring via existing aggregation engine (mean-Z, rank-based)
- [ ] Missing value handling: MinProb, KNN, zero-fill with flagging
- [ ] Normalization: median centering, quantile normalization, variance stabilization
- [ ] Unified pathway score matrix across all modalities

**Acceptance Criteria:**
- `pipeline.run()` accepts `--input-type proteomics`
- Handles 30-50% missingness without crashing
- Pathway scores pass existing validation gates
- Documentation with proteomics workflow

**GitHub Issue:** #38

---

### [Week 24] Metabolomics Feature Integration

**Theme:** Metabolite abundance as a pathway scoring input modality

**Deliverables:**
- [ ] Accept metabolite abundance matrix (LC-MS/MS, NMR) as input type
- [ ] Map metabolites to pathways via HMDB or KEGG compound IDs
- [ ] Metabolic pathway activity scoring (aggregate metabolite levels per pathway)
- [ ] Support for targeted panels: TCA intermediates, redox cofactors, acylcarnitines, amino acids
- [ ] Combined feature matrix: gene + protein + metabolite pathway scores

**Acceptance Criteria:**
- Pipeline accepts `--input-type metabolomics`
- Metabolite-to-pathway mapping works for HMDB and KEGG identifiers
- Combined feature matrix produces valid clustering
- Documentation with metabolomics workflow example

**GitHub Issue:** #39

---

### [Week 25] Multi-Omics Factor Analysis (MOFA) Integration

**Theme:** Joint latent factor discovery across data modalities

**Deliverables:**
- [ ] Simplified latent factor model (or MOFA+ wrapper via mofapy2)
- [ ] Accept 2+ modalities → extract shared and modality-specific factors
- [ ] Use latent factors as features for GMM clustering
- [ ] Variance decomposition: which modality drives which subtype
- [ ] Compare factor-based vs pathway-score-based subtypes (ARI concordance)

**Acceptance Criteria:**
- Factor-based clustering passes bootstrap validation gate
- Variance decomposition report shows per-modality contribution
- Falls back gracefully with single modality
- Documentation with multi-omics integration tutorial

**GitHub Issue:** #40

---

### [Week 26] CRISPR Screen Functional Overlay

**Theme:** Bridge computational subtypes with experimental gene essentiality

**Deliverables:**
- [ ] Import CRISPR screen results (MAGeCK, BAGEL, custom formats)
- [ ] Annotate subtype-driving pathways with essentiality scores
- [ ] Enrichment test: subtype pathways vs CRISPR hits (Fisher's exact, FDR)
- [ ] Visualization: pathway enrichment vs essentiality scatterplot

**Acceptance Criteria:**
- Accepts standard MAGeCK gene_summary.txt format
- Produces interpretable, FDR-corrected enrichment statistics
- Documentation with example workflow

**GitHub Issue:** #41

---

### [Week 27] Experimental Design Module

**Theme:** Study planning tools for multi-omics subtyping experiments

**Deliverables:**
- [ ] Power calculator: samples needed per modality to detect k subtypes
- [ ] Batch design optimizer: minimize confounding in sample-to-batch allocation
- [ ] Modality cost-benefit analysis: incremental ARI gain per additional modality
- [ ] Exportable study design report

**Acceptance Criteria:**
- Power calculator handles 1-5 input modalities
- Integrates with existing simulation and batch correction modules
- Documentation with study planning tutorial

**GitHub Issue:** #42

---

### [Week 28] Multi-Omics QC Dashboard

**Theme:** Unified quality control across all data modalities

**Deliverables:**
- [ ] RNA-seq QC: library size, gene detection rate, mitochondrial read fraction
- [ ] Proteomics QC: missing value patterns, batch effects, normalization
- [ ] Metabolomics QC: CV per metabolite, batch drift, left-censoring
- [ ] Cross-modality: sample concordance, multi-modality outlier flagging
- [ ] Self-contained HTML QC report

**Acceptance Criteria:**
- Runs on any combination of available modalities
- Clear pass/warn/fail per sample per modality
- Conservative outlier flagging (warns, does not auto-exclude)
- Integration with existing variant QC module

**GitHub Issue:** #43

---

### [Week 29] v0.4 Release

**Theme:** Ship v0.4

**Deliverables:**
- [ ] Update CHANGELOG.md for v0.4
- [ ] Create GitHub release `v0.4.0`
- [ ] Update PyPI package
- [ ] Update Zenodo DOI
- [ ] Colab notebook: full multi-omics demo (5+ modalities)
- [ ] Announce on LinkedIn, Substack, Twitter/X

**Acceptance Criteria:**
- All Phase 5 features merged and tested
- 6 input modalities supported (VCF, bulk RNA, scRNA, spatial, proteomics, metabolomics)
- Factor analysis, CRISPR overlay, experimental design operational
- MitoCarta mitochondrial pathways included
- Full documentation and examples

**GitHub Issue:** #44

---

## Community Feedback: LinkedIn Poll (2026-02-02)

**Question:** Preferred feature for pathway-subtyping framework?
**Total votes:** 33

| Feature | Votes | % |
|---------|-------|---|
| **Multi-omic support (RNAseq + WES)** | 22 | 67% |
| Integrate with other pipelines | 8 | 24% |
| More pre-built disease pathways | 2 | 6% |
| Cloud ready deployment | 1 | 3% |

### Key Insights

1. **Multi-omic is the clear priority** — users want to combine genomics with transcriptomics (RNAseq + WES integration)
2. **Pipeline integration is secondary** — interoperability with existing bioinformatics workflows matters
3. **Pre-built pathways and cloud deployment are low priority** — users are willing to curate their own pathways and run locally

### Roadmap Implications

| Poll Result | Roadmap Impact |
|-------------|----------------|
| Multi-omic (67%) | **New priority**: Add RNAseq integration module in v0.3 |
| Pipeline integration (24%) | Strengthen Nextflow/Snakemake examples in Week 6 |
| Pre-built pathways (6%) | Week 5 scope is sufficient |
| Cloud deployment (3%) | Defer to post-v0.2 |

### Validated Use Cases

- **Naresh Doni Jayavelu (Benaroya Research)**: Asthma research — heterogeneous disease, interested in multi-omic approach
---

## Success Metrics

| Metric | v0.1 Target | v0.2 Actual | v0.3 Target | v0.4 Target |
|--------|-------------|-------------|-------------|-------------|
| GitHub Stars | 10 | — | 100 | 250 |
| PyPI Downloads | 100 | — | 2,000 | 10,000 |
| External Collaborators | 1 | 4 responding | 10 | 20 |
| Disease Pathways | 4 | 6 | 6+ | 7+ (+ MitoCarta) |
| Issues Closed | 80% | ~70% | 80% | 85% |
| Test Coverage | 64 tests | 716 tests | 750 tests | 900 tests |
| Input Modalities | 1 (VCF) | 1 (VCF) | 4 (VCF, bulk RNA, scRNA, spatial) | 6 (+ proteomics, metabolomics) |
| Validation Gates | 3 | 4 | 5 (+ cross-modal) | 5 (+ multi-omics QC) |
| Integration Methods | 1 (GMM) | 1 (GMM) | 1 (GMM) | 2 (GMM + MOFA) |

---

## Project Board Setup

To create this project in GitHub:

```bash
# Create project
gh project create "Pathway Subtyping Framework: v0.2 Roadmap" \
  --owner topmist-admin

# Add issues for each milestone
# (See individual issue templates below)
```

---

*Last updated: 2026-02-15*
