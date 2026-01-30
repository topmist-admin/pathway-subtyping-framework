# Project Proposal Template

Use this template to propose a new disease study using the Pathway Subtyping Framework.

---

## Project Proposal: [Disease Name] Molecular Subtyping

**Date**: [Date]  
**Proposer**: [Your Name]  
**Status**: Draft / Under Review / Approved

---

### 1. Executive Summary

*Provide a 2-3 sentence summary of the proposed project.*

```
[Example]
We propose to identify molecular subtypes of [Disease] by analyzing 
rare variant burden in biological pathways using data from [Repository]. 
This study will be the [first/second] disease application of the 
Pathway Subtyping Framework beyond autism.
```

---

### 2. Background and Rationale

#### 2.1 Disease Overview

| Attribute | Description |
|-----------|-------------|
| Disease | [Full name] |
| Prevalence | [X per Y people] |
| Genetic architecture | [Describe: monogenic, oligogenic, complex] |
| Known genes | [Number and examples] |
| Current challenges | [What problems does subtyping address?] |

#### 2.2 Scientific Rationale

*Why is pathway-based subtyping appropriate for this disease?*

- [ ] Multiple genes implicated (genetic heterogeneity)
- [ ] Genes converge on identifiable biological pathways
- [ ] Current treatments are not targeted to mechanisms
- [ ] Subtyping could inform precision medicine approaches

#### 2.3 Prior Work

*Summarize relevant prior studies:*

| Study | Year | Finding |
|-------|------|---------|
| [Reference] | [Year] | [Key finding] |
| [Reference] | [Year] | [Key finding] |

---

### 3. Specific Aims

**Aim 1**: [Pathway burden analysis]
- Compute pathway-level rare variant burden scores
- Pathways: [List proposed pathways]
- Data: [Source and sample size]

**Aim 2**: [Subtype discovery]
- Apply GMM clustering to identify subtypes
- Validation: label shuffle, random genes, bootstrap
- Expected subtypes: [Number, if hypothesized]

**Aim 3**: [Phenotypic characterization]
- Correlate subtypes with clinical features
- Features: [List key phenotypes to analyze]

**Aim 4**: [Optional - validation/replication]
- Replicate in independent cohort
- Compare to existing subtyping approaches

---

### 4. Data Sources

#### 4.1 Primary Dataset

| Attribute | Value |
|-----------|-------|
| Repository | [SFARI, UK Biobank, dbGaP, etc.] |
| Accession | [If known] |
| Sample size | [N cases, N controls if applicable] |
| Data types | [WES, WGS, phenotypes] |
| Access status | [Open, Applied, To apply] |
| Estimated access date | [Date] |

#### 4.2 Secondary/Validation Dataset (if applicable)

| Attribute | Value |
|-----------|-------|
| Repository | |
| Sample size | |
| Purpose | [Replication, validation] |

---

### 5. Pathway Selection

*List the biological pathways to be analyzed:*

| Pathway | Genes (n) | Rationale |
|---------|-----------|-----------|
| [Pathway 1] | [N] | [Why relevant to disease] |
| [Pathway 2] | [N] | [Why relevant to disease] |
| [Pathway 3] | [N] | [Why relevant to disease] |
| [Pathway 4] | [N] | [Why relevant to disease] |

**Source of pathway definitions**: [Literature, GO, KEGG, custom curation]

---

### 6. Methods

#### 6.1 Variant Filtering

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| MAF threshold | [0.001] | [Rare variants] |
| Consequences | [LoF, missense] | [Functional impact] |
| CADD threshold | [25 for missense] | [Deleteriousness] |

#### 6.2 Clustering Approach

| Parameter | Value |
|-----------|-------|
| Method | Gaussian Mixture Model |
| K range | [2-8] |
| Selection criterion | BIC |
| Covariance type | Full |

#### 6.3 Validation Strategy

| Test | Threshold | Purpose |
|------|-----------|---------|
| Label shuffle | ARI < 0.15 | Null hypothesis |
| Random genes | ARI < 0.15 | Pathway biology |
| Bootstrap | ARI > 0.80 | Stability |

---

### 7. Expected Outcomes

*What do you expect to find and deliver?*

- [ ] Identification of [N] molecular subtypes
- [ ] Characterization of dominant pathway per subtype
- [ ] Phenotypic correlations for each subtype
- [ ] GMT file with disease-specific pathways
- [ ] Config file for framework
- [ ] Peer-reviewed publication
- [ ] Open-source analysis code

---

### 8. Timeline

| Phase | Duration | Activities |
|-------|----------|------------|
| Planning | [X weeks] | Literature review, pathway curation |
| Data access | [X weeks] | Application, IRB, approval |
| Analysis | [X weeks] | Pipeline execution, validation |
| Interpretation | [X weeks] | Phenotype analysis, manuscript |
| Publication | [X weeks] | Submission, revision |

**Total estimated duration**: [X months]

---

### 9. Resources Needed

#### 9.1 Compute

| Resource | Specification | Cost |
|----------|---------------|------|
| Cloud compute | [AWS/GCP instance type] | [$X/month] |
| Storage | [X TB] | [$X/month] |

#### 9.2 Data Access

| Repository | Access fee | Status |
|------------|------------|--------|
| [Repository] | [$X] | [Applied/Pending] |

#### 9.3 Personnel

| Role | Effort | Person |
|------|--------|--------|
| Analysis lead | [X%] | [Name] |
| Domain expert | [X%] | [Name/TBD] |
| Reviewer | [X%] | [Name] |

---

### 10. Risks and Mitigation

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Data access delay | Medium | High | Apply early, have backup plan |
| Insufficient sample size | Low | High | Target large repositories |
| No clear subtypes | Medium | Medium | Document null result, publish |
| Pathway definitions incomplete | Medium | Medium | Iterative curation |

---

### 11. Success Criteria

The project will be considered successful if:

- [ ] All validation gates pass (label shuffle, random genes, bootstrap)
- [ ] Subtypes show interpretable pathway profiles
- [ ] At least one phenotypic association identified
- [ ] Results published or preprinted
- [ ] Code and configs contributed to framework

---

### 12. Approvals

| Role | Name | Date | Signature |
|------|------|------|-----------|
| Proposer | | | |
| Project Lead | | | |
| Reviewer | | | |

---

### Appendix A: References

*List key references supporting the proposal.*

1. [Reference 1]
2. [Reference 2]

### Appendix B: Proposed GMT File

*Attach or link draft pathway definitions if available.*

---

*Submit completed proposal to project lead for review.*
