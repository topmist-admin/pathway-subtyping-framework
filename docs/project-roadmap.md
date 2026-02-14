# Pathway Subtyping Framework: v0.1 → v0.3 Roadmap

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
- [ ] Document cross-cohort validation workflow
- [ ] Add example with synthetic cohorts
- [ ] Validation report template
- [ ] Guidance on interpreting transfer ARI

**Acceptance Criteria:**
- Cross-cohort tutorial in notebook
- Clear criteria for "replicable subtypes"
- Example report generated

---

### [Week 8] Performance Optimization

**Theme:** Scale to 10,000+ samples

**Deliverables:**
- [ ] Benchmark on large synthetic dataset
- [ ] Optimize memory usage for chunked processing
- [ ] Add progress bars for long operations
- [ ] Document recommended hardware specs

**Acceptance Criteria:**
- 10,000 sample cohort processes in <30 minutes
- Memory stays under 8GB
- Clear guidance on resource requirements

---

### [Week 9] Advanced Visualization

**Theme:** Interactive exploration

**Deliverables:**
- [ ] Interactive HTML report with Plotly
- [ ] UMAP/t-SNE visualization option
- [ ] Subtype trajectory visualization
- [ ] Export figures in multiple formats

**Acceptance Criteria:**
- HTML report opens in browser
- Figures are customizable
- Publication-quality exports

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

| Metric | v0.1 Target | v0.2 Actual | v0.3 Target |
|--------|-------------|-------------|-------------|
| GitHub Stars | 10 | — | 100 |
| PyPI Downloads | 100 | — | 2,000 |
| External Collaborators | 1 | 4 responding | 10 |
| Disease Pathways | 4 | 6 | 6+ |
| Issues Closed | 80% | ~70% | 80% |
| Test Coverage | 64 tests | 599 tests | 700 tests |
| Input Modalities | 1 (VCF) | 1 (VCF) | 4 (VCF, bulk RNA, scRNA, spatial) |
| Validation Gates | 3 | 4 | 5 (+ cross-modal) |

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

*Last updated: 2026-02-13*
