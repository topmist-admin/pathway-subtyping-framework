# Research Compliance Guide

This guide covers ethical, legal, and regulatory requirements for genetic research using the Pathway Subtyping Framework.

## Overview

Working with human genetic data requires compliance with:
- Institutional Review Board (IRB) requirements
- Data Use Agreements (DUAs)
- Data security standards
- Privacy regulations (HIPAA, GDPR where applicable)

## IRB Requirements

### What is an IRB?

An Institutional Review Board reviews research involving human subjects to ensure ethical conduct and participant protection.

### When is IRB Approval Needed?

| Scenario | IRB Required? |
|----------|---------------|
| De-identified public data (gnomAD, ClinVar) | No |
| De-identified controlled data (dbGaP, SFARI) | Usually exempt* |
| Data with indirect identifiers | Review needed |
| Direct participant contact | Full review |

*Many institutions require IRB exemption documentation even for de-identified data

### IRB Exemption Categories

For most controlled-access genetic data (dbGaP, SFARI, UK Biobank):

**Category 4 Exemption** typically applies:
> Research involving secondary use of identifiable private information 
> or identifiable biospecimens, if the information is publicly available 
> OR recorded by the investigator in such a manner that the identity of 
> the human subjects cannot be readily ascertained.

### Sample IRB Exemption Request

```
Protocol Title: Pathway-Based Molecular Subtyping of [Disease]

Principal Investigator: [Name]

Description: This study involves secondary analysis of existing 
de-identified genetic and phenotypic data from [Repository]. 
No direct contact with participants will occur. Data are coded 
with no direct identifiers, and the research team will not have 
access to the code key.

Exemption Category: 4 (Secondary research, de-identified)

Data Source: [Repository name and accession]

Data Elements: Genetic variants (VCF), phenotype codes (no PHI)

Consent: Original participants consented to research use of data.
```

## Data Use Agreements

### What is a DUA?

A legal contract between your institution and the data provider specifying:
- Permitted uses of data
- Security requirements
- Publication and sharing rules
- Duration of access
- Destruction requirements

### Common DUA Terms

| Requirement | Typical Terms |
|-------------|---------------|
| **Use limitations** | Research purposes only, no commercial use |
| **Sharing** | Cannot redistribute individual-level data |
| **Security** | Encrypted storage, access controls |
| **Publication** | Aggregate results only, acknowledge source |
| **Duration** | 1-3 years, renewable |
| **Destruction** | Delete data within 30 days of project end |

### DUA Process

1. **Review terms** - Read DUA carefully before requesting data
2. **Institutional sign-off** - Your legal/compliance team reviews
3. **Sign agreement** - Authorized institutional official signs
4. **Receive data** - Access granted after DUA execution
5. **Maintain compliance** - Follow terms throughout project

## Data Security Requirements

### Minimum Security Standards

All controlled-access genetic data must meet these requirements:

#### Storage Security
- [ ] Encrypted at rest (AES-256 or equivalent)
- [ ] Access-controlled storage location
- [ ] No storage on personal devices
- [ ] No cloud storage without approval

#### Access Controls
- [ ] Unique user accounts for all personnel
- [ ] Strong passwords (12+ characters)
- [ ] Multi-factor authentication where available
- [ ] Role-based access (minimum necessary)

#### Network Security
- [ ] Encrypted data transmission (TLS 1.2+)
- [ ] Firewall protection
- [ ] No email transmission of data
- [ ] VPN for remote access (if permitted)

#### Audit and Monitoring
- [ ] Access logs maintained
- [ ] Regular access reviews
- [ ] Incident detection capabilities
- [ ] Log retention per policy

### Approved Computing Environments

| Environment | Acceptable? | Notes |
|-------------|-------------|-------|
| Institutional HPC | Yes | With security approval |
| AWS GovCloud | Yes | For government data |
| AWS/GCP/Azure | Usually | Check DUA terms |
| Personal laptop | No | Never for individual data |
| Shared drives | No | Unless specifically secured |

### Cloud Security Checklist

If using cloud computing:
- [ ] Region restrictions honored (e.g., UK data stays in UK)
- [ ] Encryption enabled for storage and transit
- [ ] IAM policies configured (least privilege)
- [ ] VPC/network isolation configured
- [ ] Audit logging enabled
- [ ] Approved by data provider

## Privacy Considerations

### Re-identification Risk

Even "de-identified" genetic data carries re-identification risk:
- Genetic variants are individually unique
- Combined with demographics, identification possible
- Family data can reveal relationships

### Mitigation Strategies

1. **Aggregate reporting** - Report group statistics, not individuals
2. **Minimum cell sizes** - No groups smaller than 5 individuals
3. **No individual genotypes** - Don't publish individual variant lists
4. **Geographic masking** - Avoid precise location data

### GDPR Considerations

For data involving EU participants:
- Genetic data is "special category" under GDPR
- Additional consent requirements may apply
- Data subject rights must be respected
- Cross-border transfer restrictions may apply

## Incident Response

### What is a Security Incident?

- Unauthorized access to data
- Data loss or theft
- Accidental disclosure
- System compromise
- Policy violation

### Reporting Requirements

| Repository | Reporting Timeframe |
|------------|---------------------|
| dbGaP | Within 24 hours |
| UK Biobank | Within 24 hours |
| SFARI | As soon as possible |
| General | Per institutional policy |

### Incident Response Steps

1. **Contain** - Stop ongoing breach, isolate affected systems
2. **Report** - Notify data provider and institutional security
3. **Investigate** - Determine scope and cause
4. **Remediate** - Fix vulnerabilities
5. **Document** - Record incident and response

### Incident Report Template

```
Date/Time of Discovery: 
Description of Incident:
Data Potentially Affected:
Systems Involved:
Immediate Actions Taken:
Root Cause (if known):
Remediation Plan:
Reporter Contact:
```

## Training Requirements

### Required Training

| Training | Provider | Frequency |
|----------|----------|-----------|
| Human Subjects Research | CITI Program | Initial + every 3 years |
| Data Security Awareness | Institutional | Annual |
| Repository-specific | Data provider | As required |

### CITI Training

Register at https://www.citiprogram.org/

Required modules (typical):
- Belmont Report and Basic Principles
- History and Ethics
- IRB Regulations
- Informed Consent
- Privacy and Confidentiality

## Compliance Checklist

Before starting analysis:

- [ ] IRB approval or exemption obtained
- [ ] Data Use Agreement signed
- [ ] Required training completed
- [ ] Computing environment approved
- [ ] Security controls implemented
- [ ] Access limited to approved personnel
- [ ] Audit logging enabled

During analysis:

- [ ] Access only approved data
- [ ] Follow DUA restrictions
- [ ] Maintain security controls
- [ ] Report any incidents immediately

At project end:

- [ ] Certify data destruction
- [ ] Submit required reports
- [ ] Archive approved results only
- [ ] Close out IRB protocol

## Data Provenance of This Repository

**Important**: The compliance requirements above apply to **users' own data** when they apply this framework to real cohorts. The framework itself **does not contain any controlled-access, proprietary, or patient data**.

### What This Repository Contains

| Content | Source | Contains Patient Data? |
|---------|--------|----------------------|
| Python source code | Original implementation | No |
| Synthetic VCF files | `SyntheticDataGenerator` (random numbers, fixed seeds) | No |
| Synthetic phenotype CSVs | `SyntheticDataGenerator` (random numbers, fixed seeds) | No |
| Pathway GMT files | Curated from public literature (SFARI Gene, KEGG, Reactome, MSigDB, GO) | No |
| Test fixtures | Programmatically generated in test setup | No |
| Configuration examples | Manually authored YAML | No |

### What This Repository Does NOT Contain

- **No real patient or clinical data** from any source
- **No data from any employer, client, or commercial entity** — this project was developed as independent personal research
- **No data obtained under any Data Use Agreement (DUA)** — no dbGaP, SFARI, UK Biobank, or other controlled-access data was used
- **No data from any Life Sciences, pharmaceutical, or healthcare company**
- **No proprietary databases** — only open-access, publicly available references are cited
- **No data requiring IRB approval** — all data is computationally generated or from public literature

### Independence Statement

This framework was developed as an independent, open-source research project. It is not associated with, derived from, or based on any work product from any employer or commercial entity. The framework is designed so that **users supply their own data**; it does not ship with, embed, or depend on any private, restricted, or commercially licensed datasets.

### Verification

Anyone can verify these claims by:
1. Inspecting `data/sample/` — all VCF/CSV files are generated by `SyntheticDataGenerator` with documented random seeds
2. Inspecting `data/pathways/*.gmt` — gene symbols are standard HGNC identifiers from public scientific literature
3. Running `git log --all --stat` — the full commit history shows no external data imports
4. Searching the repository for any database connection strings, API keys, or data download scripts — none exist (the optional `validation_datasets.py` downloads only from public FTP/HTTP endpoints: ClinVar and Reactome)

For more details, see [DISCLAIMER.md](../../DISCLAIMER.md).

## Resources

- **CITI Program**: https://www.citiprogram.org/
- **OHRP (HHS)**: https://www.hhs.gov/ohrp/
- **NIH Data Sharing**: https://sharing.nih.gov/
- **GDPR**: https://gdpr.eu/

---

*Questions about compliance? Contact your institutional IRB or research compliance office.*
