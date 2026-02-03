# SFARI Data Access Application Template

**Simons Foundation Autism Research Initiative (SFARI)**

SFARI provides access to genetic and phenotypic data from autism research cohorts including the Simons Simplex Collection (SSC) and SPARK.

## Portal Information

- **URL**: https://base.sfari.org/
- **Data Available**: WES, WGS, phenotypes, medical history
- **Cohort Size**: ~30,000 families
- **Typical Review Time**: 4-6 weeks

## Pre-Application Checklist

- [ ] Create SFARI Base account
- [ ] Verify institutional affiliation
- [ ] Prepare research proposal
- [ ] Prepare investigator biosketch
- [ ] Confirm IRB status (exemption likely sufficient)

## Application Template

### Section 1: Project Information

**Project Title:**
```
Pathway-Based Molecular Subtyping of Autism Spectrum Disorder Using 
Rare Variant Burden Analysis
```

**Keywords:**
```
autism, molecular subtyping, rare variants, pathway analysis, 
precision medicine, genetic heterogeneity
```

### Section 2: Research Proposal

**Background and Rationale** (adapt as needed):
```
Autism spectrum disorder (ASD) is genetically heterogeneous, with 
rare variants in over 100 genes implicated. This heterogeneity 
presents challenges for understanding disease mechanisms and 
developing targeted treatments.

We propose to identify molecular subtypes of ASD by aggregating 
rare variant burden into biologically coherent pathways (synaptic, 
chromatin remodeling, etc.) and applying unsupervised clustering. 
This approach has potential to:

1. Identify patient subgroups with shared molecular mechanisms
2. Improve genotype-phenotype correlations
3. Guide development of pathway-targeted therapies
```

**Specific Aims:**
```
Aim 1: Compute pathway-level rare variant burden scores for ASD 
probands using curated gene sets (synaptic transmission, chromatin 
regulation, Wnt signaling, mTOR pathway).

Aim 2: Apply Gaussian mixture model clustering to identify molecular 
subtypes based on pathway burden profiles.

Aim 3: Validate subtypes using: (a) negative control tests to rule 
out spurious clustering, (b) bootstrap resampling to assess stability, 
(c) phenotypic correlation analysis.

Aim 4: Characterize subtypes by clinical features, co-occurring 
conditions, and developmental trajectories.
```

**Methods:**
```
Data: Whole-exome sequencing data from SSC and/or SPARK cohorts.

Variant Filtering: Rare variants (MAF < 0.1%) with predicted 
functional impact (LoF, missense with CADD > 25).

Pathway Scoring: Gene-level burden aggregated to curated pathways 
(4-8 pathways based on ASD literature). Scores z-normalized across 
cohort.

Clustering: Gaussian mixture models with BIC-based model selection 
(k = 2-8 clusters tested).

Validation: 
- Label shuffle test (ARI < 0.15 expected under null)
- Random gene set test (pathway biology vs. random)
- Bootstrap stability test (ARI > 0.8 required)

Analysis Software: Pathway Subtyping Framework (open-source, 
available at github.com/topmist-admin/pathway-subtyping-framework)
```

**Expected Outcomes:**
```
1. Identification of 3-5 molecular subtypes of ASD
2. Characterization of each subtype's dominant pathway disruption
3. Correlation of subtypes with phenotypic features
4. Open-source analysis pipeline for replication
5. Peer-reviewed publication of findings
```

### Section 3: Data Requested

**Datasets:**
- [ ] Simons Simplex Collection (SSC)
- [ ] SPARK
- [ ] Other: _______________

**Data Types:**
- [x] Whole-exome sequencing (VCF)
- [ ] Whole-genome sequencing (VCF)
- [x] Phenotype data (ADI-R, ADOS, IQ, medical history)
- [ ] Other: _______________

**Estimated Sample Size Needed:**
```
Minimum: 1,000 probands with WES + phenotypes
Ideal: 5,000+ probands for robust subtype discovery
```

### Section 4: Investigator Information

**Principal Investigator:**
```
Name: [Your Name]
Title: [Your Title]
Institution: [Your Institution]
Email: [your.email]
```

**Biosketch/CV:**
- Attach 2-page NIH-style biosketch or CV
- Highlight relevant computational biology or genetics experience

### Section 5: Institutional Information

**Institution:**
```
[Your Institution]
[Address]
```

**IRB Status:**
```
[ ] IRB approval obtained (attach documentation)
[x] IRB exemption (de-identified data)
[ ] IRB review in progress
```

### Section 6: Data Security

**Computing Environment:**
```
Analysis will be performed on:
[ ] Institutional secure server
[x] Cloud environment (AWS/GCP with encryption)
[ ] Other: _______________

Security measures:
- Data encrypted at rest (AES-256)
- Data encrypted in transit (TLS 1.3)
- Access restricted to approved personnel
- Audit logging enabled
- Data deleted after project completion per DUA
```

### Section 7: Publication Plans

```
We plan to publish findings in a peer-reviewed journal 
(e.g., Nature Genetics, American Journal of Human Genetics).

SFARI will be acknowledged per data use agreement requirements.
Analysis code will be made publicly available.
```

## Submission Checklist

- [ ] Complete all sections of online application
- [ ] Upload research proposal (PDF)
- [ ] Upload investigator biosketch (PDF)
- [ ] Upload IRB documentation
- [ ] Agree to data use terms
- [ ] Submit application

## Post-Approval Steps

1. Sign Data Use Agreement
2. Receive SFARI Base download credentials
3. Download approved datasets
4. Transfer to secure analysis environment
5. Begin analysis per approved protocol

---

*Template based on SFARI application requirements as of January 2026*
