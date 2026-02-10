# Data Access Guide

This guide explains how to obtain access to genetic data for research using the Pathway Subtyping Framework.

## Overview of Data Access Tiers

| Tier | Access Type | Examples | Timeline |
|------|-------------|----------|----------|
| **Open** | No application needed | gnomAD, ClinVar, OMIM | Immediate |
| **Registered** | Simple registration | Gene Expression Omnibus | 1-2 days |
| **Controlled** | Full application + IRB | SFARI, UK Biobank, dbGaP | 2-6 months |

## Recommended Data Sources by Disease

### Autism Spectrum Disorder
| Repository | Data Type | Size | Access |
|------------|-----------|------|--------|
| **SFARI** | WES/WGS + phenotypes | ~30,000 families | Controlled |
| **MSSNG** | WGS + phenotypes | ~12,000 individuals | Controlled |
| **Autism Sequencing Consortium** | WES | ~35,000 cases | Controlled |

### Schizophrenia
| Repository | Data Type | Size | Access |
|------------|-----------|------|--------|
| **UK Biobank** | WES/WGS + phenotypes | ~500,000 (subset) | Controlled |
| **dbGaP - PGC** | GWAS + WES | ~150,000 | Controlled |
| **CLOZUK** | WES + clinical | ~12,000 | Collaboration |

### Epilepsy
| Repository | Data Type | Size | Access |
|------------|-----------|------|--------|
| **Epi25** | WES + phenotypes | ~25,000 cases | Controlled |
| **UK Biobank** | WES/WGS + phenotypes | ~500,000 (subset) | Controlled |
| **dbGaP - EPGP** | WES + clinical | ~4,000 families | Controlled |

### Intellectual Disability
| Repository | Data Type | Size | Access |
|------------|-----------|------|--------|
| **DDD Study** | WES + phenotypes | ~14,000 trios | Controlled |
| **UK Biobank** | WES/WGS + cognitive | ~500,000 | Controlled |
| **dbGaP** | Various | Variable | Controlled |

### Parkinson's Disease
| Repository | Data Type | Size | Access |
|------------|-----------|------|--------|
| **AMP-PD** | WGS + longitudinal clinical | ~10,000 cases | Controlled |
| **PPMI** | WGS + biomarkers + imaging | ~1,500 PD + prodromal | Controlled |
| **GP2** | WGS + phenotypes | ~150,000 (growing) | Controlled |
| **UK Biobank** | WES/WGS + phenotypes | ~500,000 (subset) | Controlled |
| **dbGaP - IPDGC** | GWAS + WES | ~40,000 | Controlled |

**Notes:**
- **AMP-PD** (Accelerating Medicines Partnership - Parkinson's Disease): https://amp-pd.org/ - Best resource for longitudinal WGS data with deep phenotyping
- **PPMI** (Parkinson's Progression Markers Initiative): https://www.ppmi-info.org/ - Includes prodromal/at-risk cohorts
- **GP2** (Global Parkinson's Genetics Program): https://gp2.org/ - Largest global PD genetics initiative

### Bipolar Disorder
| Repository | Data Type | Size | Access |
|------------|-----------|------|--------|
| **UK Biobank** | WES/WGS + phenotypes | ~500,000 (subset) | Controlled |
| **dbGaP - PGC Bipolar** | GWAS + WES | ~40,000 cases | Controlled |
| **STEP-BD** | WES + clinical + treatment | ~4,000 | Controlled |
| **Stanley Brain Collection** | WES + postmortem | ~600 brains | Controlled |
| **CommonMind Consortium** | WES + brain expression | ~1,000 | Controlled |

**Notes:**
- **PGC Bipolar**: https://pgc.unc.edu/for-researchers/resources/bipolar/ - Largest BD GWAS consortium
- **CommonMind**: https://www.synapse.org/cmc - Combines genomics with brain transcriptomics

## Application Process Overview

### Step 1: Prepare Documentation (2-4 weeks)

**Required for most controlled-access repositories:**

1. **Research Proposal** (1-2 pages)
   - Scientific rationale
   - Specific aims
   - Methodology overview
   - Expected outcomes

2. **IRB Approval or Exemption**
   - May need institutional IRB review
   - Some repositories accept exemption for de-identified data
   - Allow 4-8 weeks for IRB review if needed

3. **Data Use Agreement (DUA)**
   - Legal agreement between your institution and data provider
   - Typically handled by legal/compliance team
   - Allow 2-4 weeks for review

4. **PI/Investigator Credentials**
   - Biosketch or CV
   - List of relevant publications
   - Institutional affiliation letter

### Step 2: Submit Application (1-2 weeks)

1. Create account on repository portal
2. Complete online application forms
3. Upload required documents
4. Submit for review

### Step 3: Review Period (4-12 weeks)

- **SFARI**: 4-6 weeks typical
- **UK Biobank**: 8-12 weeks typical
- **dbGaP**: 6-8 weeks typical
- **AMP-PD**: 2-4 weeks typical (streamlined process)
- **PPMI**: 2-4 weeks typical
- **GP2**: 4-6 weeks typical

### Step 4: Data Access (1-2 weeks)

Once approved:
1. Sign final data use agreement
2. Receive access credentials
3. Download data to approved computing environment
4. Begin analysis

## Open-Access Resources

These resources don't require applications and can be used immediately:

### gnomAD (Genome Aggregation Database)
- **URL**: https://gnomad.broadinstitute.org/
- **Content**: Variant frequencies from ~140,000 exomes and ~76,000 genomes
- **Use Case**: Filter variants by population frequency, assess pathogenicity

### ClinVar
- **URL**: https://www.ncbi.nlm.nih.gov/clinvar/
- **Content**: Variant-disease associations with clinical interpretations
- **Use Case**: Identify known pathogenic variants

### OMIM (Online Mendelian Inheritance in Man)
- **URL**: https://omim.org/
- **Content**: Gene-disease relationships, phenotype descriptions
- **Use Case**: Curate disease-relevant pathways

### Gene Ontology (GO)
- **URL**: http://geneontology.org/
- **Content**: Standardized gene function annotations
- **Use Case**: Build pathway definitions, functional enrichment

### KEGG Pathways
- **URL**: https://www.genome.jp/kegg/pathway.html
- **Content**: Curated biological pathways
- **Use Case**: Reference pathways for GMT files

## Data Security Requirements

All controlled-access genetic data must be:

1. **Stored securely** - Encrypted at rest and in transit
2. **Accessed from approved environments** - Cloud or institutional compute
3. **Not shared** - Cannot redistribute raw data
4. **Deleted after project** - Per DUA terms
5. **Audited** - Maintain access logs

See [04-research-compliance.md](04-research-compliance.md) for detailed requirements.

## Application Templates

Ready-to-use application templates are available in:
- [03-data-sources/sfari-application-template.md](03-data-sources/sfari-application-template.md)
- [03-data-sources/ukbiobank-application-template.md](03-data-sources/ukbiobank-application-template.md)
- [03-data-sources/dbgap-application-template.md](03-data-sources/dbgap-application-template.md)

## Timeline Planning

| Activity | Duration |
|----------|----------|
| Prepare documentation | 2-4 weeks |
| IRB review (if needed) | 4-8 weeks |
| Application submission | 1-2 weeks |
| Repository review | 4-12 weeks |
| Data access setup | 1-2 weeks |
| **Total** | **3-6 months** |

**Recommendation**: Start data access applications early while learning the framework with synthetic data.

---

## Important: Data Provenance of This Framework

**This framework does not ship with, embed, or include any data from the repositories listed above.** The data sources in this guide are provided as a reference for users who wish to apply the framework to real cohorts for their own research.

### What This Repository Contains

- **Synthetic test data only** — generated by `SyntheticDataGenerator` with random numbers and fixed seeds
- **Pathway gene lists** — standard HGNC gene symbols curated from public scientific literature
- **No controlled-access data** — no dbGaP, SFARI, UK Biobank, or other restricted data was used in development

### Independence Statement

This framework was developed as an independent, open-source personal research project. No data from any employer, client, institution, or commercial entity was used at any stage — not in development, testing, validation, or documentation. Users must obtain their own data through the proper channels described above.

For full provenance details, see [DISCLAIMER.md](../../DISCLAIMER.md).

---

*Questions about data access? Contact the project lead.*
