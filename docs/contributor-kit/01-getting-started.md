# Getting Started Guide

Welcome to the genetic research team! This guide will help you get up and running with the Pathway Subtyping Framework.

## Prerequisites

### Required Background
- [ ] Basic understanding of genetics (genes, variants, inheritance patterns)
- [ ] Familiarity with Python programming
- [ ] Understanding of statistical concepts (clustering, validation)

### Recommended Training
- [Coursera: Introduction to Genomics](https://www.coursera.org/learn/introduction-genomics)
- [EdX: Data Analysis for Life Sciences](https://www.edx.org/professional-certificate/harvardx-data-analysis-for-life-sciences)
- Internal training: "Genomics 101 for Engineers" (contact project lead)

## Onboarding Checklist

### Week 1: Setup and Access

- [ ] **GitHub Access**
  - Request access to `topmist-admin/pathway-subtyping-framework`
  - Clone the repository locally
  - Review the README and documentation

- [ ] **Development Environment**
  - Install Python 3.9+ 
  - Set up virtual environment
  - Install framework dependencies
  - Run the demo notebook successfully

- [ ] **Team Communication**
  - Join the #genetic-research Slack/Teams channel
  - Schedule intro meeting with project lead
  - Add project meetings to calendar

### Week 2: Learning the Framework

- [ ] **Read Core Documentation**
  - Framework architecture overview
  - Validation gates explanation
  - Pathway curation guide

- [ ] **Complete Tutorials**
  - Run `01_getting_started.ipynb` notebook
  - Understand input data formats (VCF, phenotypes, GMT)
  - Review example configurations

- [ ] **Shadow a Research Session**
  - Attend a working session with experienced team member
  - Ask questions about workflow and decisions

### Week 3: Data Access Preparation

- [ ] **Identify Your Disease Focus**
  - Review available disease areas
  - Discuss with project lead
  - Document rationale for selection

- [ ] **Begin Data Access Application**
  - Select appropriate data repository
  - Gather required documentation
  - Draft application using templates

- [ ] **Ethics Training**
  - Complete CITI Human Subjects Research training
  - Review data use agreement requirements
  - Understand PHI handling procedures

## Key Concepts to Understand

### 1. Genetic Heterogeneity
Many diseases (autism, schizophrenia, epilepsy) are caused by mutations in hundreds of different genes. Our framework groups patients by which *biological pathways* are disrupted, not individual genes.

### 2. Pathway-Based Analysis
Instead of analyzing individual genes, we aggregate variants into biological pathways (e.g., "synaptic signaling", "chromatin remodeling"). This reduces the dimensionality problem and provides biologically interpretable subtypes.

### 3. Validation Gates
We use three statistical tests to ensure our subtypes are real, not artifacts:
- **Label Shuffle**: Are clusters better than random?
- **Random Gene Sets**: Does pathway biology matter?
- **Bootstrap Stability**: Are clusters reproducible?

### 4. GMT Files
Gene Matrix Transposed format defines pathways:
```
PATHWAY_NAME<tab>DESCRIPTION<tab>GENE1<tab>GENE2<tab>GENE3...
```

## Getting Help

| Issue | Contact |
|-------|---------|
| Technical setup | #genetic-research channel |
| Data access questions | Project lead |
| Framework bugs | GitHub Issues |
| Research direction | Weekly team meeting |

## Next Steps

After completing this onboarding:

1. Read [02-data-access-guide.md](02-data-access-guide.md) to understand data sources
2. Choose a disease area and review relevant literature
3. Begin your data access application process
4. Set up your analysis environment

---

*Questions? Reach out to the project lead or post in the team channel.*
