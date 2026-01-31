# Pathway Subtyping Framework: v0.1 → v0.2 Roadmap

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

## Success Metrics

| Metric | v0.1 Target | v0.2 Target |
|--------|-------------|-------------|
| GitHub Stars | 10 | 50 |
| PyPI Downloads | 100 | 500 |
| External Collaborators | 1 | 5 |
| Contributed Pathways | 4 | 8 |
| Issues Closed | 80% | 80% |
| Test Coverage | 64 tests | 100 tests |

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

*Last updated: 2026-01-29*
