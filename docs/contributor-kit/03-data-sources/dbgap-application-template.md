# dbGaP Data Access Application Template

**Database of Genotypes and Phenotypes (dbGaP)** is the NIH repository for controlled-access genetic studies.

## Portal Information

- **URL**: https://dbgap.ncbi.nlm.nih.gov/
- **Data Available**: GWAS, WES, WGS, phenotypes from NIH-funded studies
- **Studies**: 1,000+ studies across many diseases
- **Typical Review Time**: 6-8 weeks

## Pre-Application Checklist

- [ ] Obtain eRA Commons account (required for PI)
- [ ] Confirm institutional signing official (SO)
- [ ] Identify specific dbGaP studies to request
- [ ] Prepare research proposal
- [ ] Obtain IRB approval or exemption
- [ ] Complete required NIH training

## Key dbGaP Studies for This Framework

### Autism/Neurodevelopment
| Study | Accession | Samples | Data Types |
|-------|-----------|---------|------------|
| Autism Sequencing Consortium | phs000298 | ~35,000 | WES |
| SPARK Autism | phs001908 | ~50,000 | WES, WGS |
| Simons Simplex Collection | phs000294 | ~2,500 | WES |

### Schizophrenia
| Study | Accession | Samples | Data Types |
|-------|-----------|---------|------------|
| PGC Schizophrenia | phs000021 | ~150,000 | GWAS, WES |
| COGS | phs000979 | ~2,000 | WES |

### Epilepsy
| Study | Accession | Samples | Data Types |
|-------|-----------|---------|------------|
| Epi25 | phs001489 | ~25,000 | WES |
| EPGP | phs000280 | ~4,000 | WES |

### Intellectual Disability
| Study | Accession | Samples | Data Types |
|-------|-----------|---------|------------|
| DDD Study | phs001232 | ~14,000 | WES |
| PCGC (CHD) | phs001194 | ~10,000 | WES |

## Application Template

### Section 1: Research Use Statement

**Project Title:**
```
Pathway-Based Molecular Subtyping of [DISEASE] Using Rare 
Variant Burden Analysis
```

**Research Use Statement** (2,500 character limit):
```
BACKGROUND: [Disease] is a genetically heterogeneous condition with 
rare variants in hundreds of genes contributing to disease risk. 
This heterogeneity complicates efforts to understand disease 
mechanisms and develop targeted treatments. Pathway-based approaches 
that aggregate variant effects across functionally related genes 
may reveal coherent molecular subtypes.

OBJECTIVES: We propose to identify molecular subtypes of [Disease] 
by:
1) Computing pathway-level rare variant burden scores using curated 
   biological pathways
2) Applying unsupervised clustering to identify patient subgroups
3) Validating subtypes using rigorous statistical controls
4) Characterizing subtypes by phenotypic features

METHODS: We will analyze whole-exome sequencing data from the 
[Study Name] (phs#####). Rare variants (MAF < 0.1%) with predicted 
functional impact will be identified and aggregated into biological 
pathways. Gaussian mixture model clustering will identify subtypes, 
validated by label permutation, random gene set controls, and 
bootstrap resampling. 

EXPECTED OUTCOMES: Identification of biologically meaningful 
molecular subtypes that may inform patient stratification and 
treatment development. Analysis code will be open-sourced. 
Findings will be published in peer-reviewed journals with 
appropriate acknowledgment of dbGaP data.

This research represents a non-commercial academic use of the data.
```

### Section 2: Non-Technical Summary

```
[Disease] affects many people but can be caused by changes in many 
different genes. We want to find groups of patients who have similar 
biological problems even though they have changes in different genes. 
This could help researchers understand the disease better and 
potentially help develop more targeted treatments in the future.
```

### Section 3: Requested Studies

```
Study 1: [Study Name]
Accession: phs######.v#.p#
Requested data: [x] Genotypes  [x] Phenotypes  [ ] Expression
Consent groups: [Specify or "All available"]

Study 2: [If requesting multiple studies]
Accession: phs######.v#.p#
...
```

### Section 4: Investigator Information

**Principal Investigator:**
```
Name: [Your Name]
eRA Commons ID: [Required]
Institution: [Your Institution]
Title: [Your Title]
Email: [your.email]
```

**Signing Official:**
```
Name: [Institutional SO]
Email: [SO email]
Institution: [Your Institution]
```

**Institutional Contact:**
```
Name: [IT Security contact if required]
Email: [contact email]
Role: [e.g., Information Security Officer]
```

### Section 5: IRB Information

```
IRB Name: [Your IRB or "Exempt"]
IRB Number: [Protocol number or N/A]
IRB Status: 
[ ] Approved (attach documentation)
[x] Exempt - using de-identified data
[ ] Pending

Expiration Date: [If applicable]
```

### Section 6: IT Security

**dbGaP IT Security Plan:**
```
1. Data Storage
   - Location: [Institutional server / Cloud / Both]
   - Encryption: AES-256 at rest
   - Access controls: Role-based, minimum necessary

2. Data Transmission
   - Download via dbGaP secure download manager
   - Encrypted transfer (TLS 1.3)
   - No email transmission of data

3. Computing Environment
   - Description: [e.g., AWS GovCloud, institutional HPC]
   - Security certifications: [e.g., SOC2, HIPAA]
   - Firewall protection: Yes

4. User Management
   - Access limited to approved project personnel
   - Unique user accounts with strong passwords
   - Multi-factor authentication enabled

5. Incident Response
   - Security incidents reported within 24 hours
   - Contact: [IT Security contact]

6. Data Destruction
   - Data deleted within 30 days of project end
   - Secure deletion methods used
   - Destruction certified to dbGaP
```

### Section 7: Data Use Certification

The following certifications must be agreed to:

- [ ] I will not attempt to identify individual participants
- [ ] I will not share data with unauthorized individuals
- [ ] I will not use data for purposes outside approved research
- [ ] I will acknowledge dbGaP and contributing studies
- [ ] I will follow all data use limitations specified by studies
- [ ] I will report any data security incidents immediately
- [ ] I will destroy data at end of approved research period

### Section 8: Acknowledgment Language

```
Data for this study were obtained from dbGaP (accession phs######). 
We thank the [Study Name] investigators and participants for 
contributing data. The [Study Name] was supported by [funding source].
```

## Submission Process

1. **PI creates request** in dbGaP portal using eRA Commons login
2. **SO approves request** electronically
3. **dbGaP Data Access Committee (DAC)** reviews (6-8 weeks)
4. **Approval notification** sent via email
5. **Download data** using dbGaP download tools

## Data Access Request (DAR) Checklist

- [ ] Log into dbGaP with eRA Commons
- [ ] Create new DAR
- [ ] Select requested studies
- [ ] Enter Research Use Statement
- [ ] Enter Non-Technical Summary
- [ ] Enter PI and SO information
- [ ] Complete IT security attestation
- [ ] Upload IRB documentation
- [ ] Submit to Signing Official
- [ ] SO submits to dbGaP

## Post-Approval Requirements

1. **Annual renewal** - Requests must be renewed annually
2. **Project closeout** - Certify data destruction at project end
3. **Publications** - Submit to dbGaP bibliography
4. **Amendments** - Request approval for scope changes

---

*Template based on dbGaP requirements as of January 2026. 
Requirements may vary by study - check study-specific documentation.*
