# UK Biobank Data Access Application Template

**UK Biobank** is a large-scale biomedical database containing genetic and health information from 500,000 UK participants.

## Portal Information

- **URL**: https://www.ukbiobank.ac.uk/
- **Data Available**: WES, WGS, genotyping, phenotypes, health records
- **Cohort Size**: ~500,000 participants
- **Typical Review Time**: 8-12 weeks
- **Access Fee**: Yes (varies by data requested)

## Pre-Application Checklist

- [ ] Register for UK Biobank Access Management System (AMS)
- [ ] Confirm institutional eligibility
- [ ] Identify collaborating UK-based researcher (if required)
- [ ] Prepare research proposal
- [ ] Obtain institutional sign-off
- [ ] Budget for access fees

## Application Template

### Section 1: Project Information

**Project Title:**
```
Pathway-Based Molecular Subtyping of Neuropsychiatric Conditions 
Using Rare Variant Burden Analysis in UK Biobank
```

**Lay Summary** (for public posting):
```
Many brain conditions like schizophrenia, epilepsy, and depression 
are caused by changes in many different genes. We want to find 
groups of patients who share similar biological problems, even if 
they have changes in different genes. This could help doctors 
understand these conditions better and find more targeted treatments.

We will analyze genetic data from UK Biobank participants to 
identify these patient groups based on which biological pathways 
are affected by their genetic changes.
```

### Section 2: Scientific Rationale

**Background:**
```
Neuropsychiatric conditions including schizophrenia, bipolar disorder, 
major depression, and epilepsy show substantial genetic heterogeneity. 
Rare coding variants across hundreds of genes contribute to disease 
risk, making it challenging to identify coherent biological mechanisms.

Pathway-based approaches aggregate variant effects across functionally 
related genes, potentially revealing shared mechanisms among patients 
with variants in different genes. UK Biobank's large sample size and 
rich phenotyping enable well-powered subtype discovery.
```

**Research Questions:**
```
1. Can neuropsychiatric conditions be stratified into molecular 
   subtypes based on pathway-level rare variant burden?

2. Do identified subtypes correlate with clinical features, 
   treatment response, or disease trajectory?

3. Are subtypes reproducible across independent cohorts?
```

### Section 3: Specific Aims

```
Aim 1: Pathway Burden Analysis
- Compute rare variant burden scores aggregated to curated 
  neurobiological pathways
- Pathways: synaptic transmission, ion channels, chromatin 
  regulation, immune signaling, mTOR/autophagy

Aim 2: Subtype Discovery
- Apply unsupervised clustering (GMM) to pathway scores
- Validate using negative controls and bootstrap resampling
- Target conditions: schizophrenia, bipolar, MDD, epilepsy

Aim 3: Phenotypic Characterization
- Correlate subtypes with clinical features
- Assess medication response patterns
- Examine comorbidity profiles

Aim 4: Cross-Disease Analysis
- Compare subtype structure across conditions
- Identify shared vs. condition-specific subtypes
```

### Section 4: Data Requested

**Genetic Data:**
- [x] Whole-exome sequencing (450K participants)
- [ ] Whole-genome sequencing (200K participants)
- [ ] Genotyping array data
- [x] Variant annotations (VEP, CADD scores)

**Phenotype Data:**
- [x] ICD-10 diagnosis codes (mental health)
- [x] Self-reported conditions
- [x] Medication history
- [x] Cognitive assessments
- [x] Mental health questionnaires
- [ ] Imaging data

**Specific Fields** (UK Biobank field IDs):
```
Diagnoses: 41202, 41204 (ICD-10 codes)
Medications: 20003 (treatment codes)
Mental health: 20126, 20127, 20421 (PHQ, GAD, etc.)
Cognitive: 20016, 20023 (fluid intelligence, reaction time)
Demographics: 31, 34, 52 (sex, year of birth, ethnicity)
```

**Sample Size:**
```
Estimated participants with neuropsychiatric diagnoses:
- Schizophrenia: ~2,000
- Bipolar disorder: ~4,000
- Major depression: ~40,000
- Epilepsy: ~8,000
```

### Section 5: Methods

**Analysis Pipeline:**
```
1. Variant Filtering
   - Rare variants: MAF < 0.1% in gnomAD
   - Functional: LoF + missense (CADD > 25)
   - Quality: PASS filter, DP > 10, GQ > 20

2. Gene Burden Calculation
   - Count qualifying variants per gene per individual
   - Weight by predicted impact (LoF = 2x, missense = 1x)

3. Pathway Aggregation
   - Sum gene burdens within curated pathways
   - Z-normalize across cohort
   - Pathways from literature + GO/KEGG

4. Clustering
   - Gaussian mixture models (k = 2-10)
   - Model selection via BIC
   - Stability assessment via bootstrap

5. Validation
   - Label permutation test (null ARI < 0.15)
   - Random gene set test
   - Holdout replication

Software: Pathway Subtyping Framework (open-source)
```

### Section 6: Investigator Information

**Principal Investigator:**
```
Name: [Your Name]
Institution: OpenText Corporation
Role: [Your Role]
Email: [your.email]@opentext.com
ORCID: [Your ORCID]
```

**Co-Investigators:**
```
[List any collaborators, especially UK-based if required]
```

### Section 7: Data Security

**Approved Research Environment:**
```
[x] UK Biobank Research Analysis Platform (RAP)
[ ] Institutional secure environment (describe)
[ ] Other approved cloud (describe)
```

**Security Measures:**
```
- Analysis performed within UK Biobank RAP
- No data download to local systems
- Results reviewed before export
- Access limited to approved personnel
- Compliance with UK Biobank security requirements
```

### Section 8: Publication and Data Return

```
Publication Plans:
- Submit findings to peer-reviewed journal
- Return summary statistics to UK Biobank
- Acknowledge UK Biobank per requirements

Data Sharing:
- Summary statistics will be shared
- Individual-level data will not be redistributed
- Analysis code will be open-sourced
```

### Section 9: Ethics and Approvals

```
UK Biobank operates under ethics approval from the 
North West Multi-centre Research Ethics Committee.

This project falls under UK Biobank's generic ethics 
approval for health-related research.

Institutional approval: [Describe any additional approvals needed]
```

## Fee Structure (Approximate)

| Data Type | Cost |
|-----------|------|
| Application fee | £3,000 |
| Exome data access | £1,500/year |
| Phenotype data | Included |
| RAP compute | Variable (usage-based) |

*Fees subject to change - verify current rates on UK Biobank website*

## Submission Checklist

- [ ] Complete AMS registration
- [ ] Submit research application
- [ ] Upload research proposal
- [ ] List all team members
- [ ] Agree to Material Transfer Agreement
- [ ] Pay application fee
- [ ] Await review decision

---

*Template based on UK Biobank requirements as of January 2026*
