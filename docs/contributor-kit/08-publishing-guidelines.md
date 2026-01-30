# Publishing Guidelines

This guide covers publication requirements, co-authorship policies, and acknowledgment standards for research using the Pathway Subtyping Framework.

## Publication Requirements

### Data Use Agreement Compliance

Before publishing, ensure compliance with all DUAs:

- [ ] Only aggregate results shared (no individual-level data)
- [ ] Minimum cell sizes respected (typically N ≥ 5)
- [ ] Data source properly acknowledged
- [ ] Publication reviewed per DUA terms (if required)

### Required Acknowledgments

All publications must acknowledge:

1. **Data sources** (as specified by each DUA)
2. **Framework** (Pathway Subtyping Framework)
3. **Funding** (if applicable)

**Example acknowledgment:**
```
Data for this study were obtained from [Repository Name] (accession: XXX). 
We thank the study participants and investigators for contributing data.
Analysis was performed using the Pathway Subtyping Framework 
(github.com/topmist-admin/pathway-subtyping-framework). 
[Additional funding acknowledgments as required]
```

### Repository-Specific Requirements

| Repository | Requirements |
|------------|--------------|
| **SFARI** | Acknowledge SFARI, cite data descriptor paper |
| **UK Biobank** | Include UK Biobank application number, standard acknowledgment |
| **dbGaP** | Acknowledge funding source, cite contributing studies |
| **gnomAD** | Cite gnomAD paper |
| **ClinVar** | Cite ClinVar paper |

## Co-Authorship Guidelines

### ICMJE Criteria

Authorship requires meeting ALL of these criteria:
1. Substantial contributions to conception/design OR data acquisition/analysis
2. Drafting the work OR revising it critically
3. Final approval of the published version
4. Agreement to be accountable for all aspects

### Typical Author Roles

| Role | Contribution | Author Position |
|------|--------------|-----------------|
| Analysis lead | Primary analysis, writing | First author |
| Project lead | Conception, supervision | Senior/corresponding |
| Domain expert | Clinical interpretation | Middle author |
| Data generation | Original data collection | Middle author |
| Technical support | Pipeline development | Middle author or acknowledgment |

### Framework Contribution

Contributors to the framework itself may be acknowledged:
- **Acknowledgment**: For general framework use
- **Authorship**: If substantial custom development was needed

### Consortium Guidelines

If using consortium data (e.g., PGC, DDD):
- Follow consortium authorship policies
- May require consortium author list
- Check with data access committee

## Journal Selection

### Recommended Journals by Impact

**High-impact (IF > 15):**
- Nature Genetics
- Cell
- American Journal of Human Genetics
- Nature Medicine

**Mid-impact (IF 5-15):**
- Genome Medicine
- Molecular Psychiatry
- Biological Psychiatry
- Human Molecular Genetics

**Specialized:**
- Autism Research
- Schizophrenia Research
- Epilepsia
- Journal of Child Psychology and Psychiatry

### Open Access Considerations

- Many funders require open access publication
- Check DUA for open access requirements
- Budget for article processing charges (APCs)

## Data Availability Statement

### Template Statement

```
Data Availability:
The genetic data used in this study are available through 
[Repository] under accession [XXX] to qualified researchers 
who agree to the data use terms. Summary statistics and 
analysis code are available at [GitHub URL]. The Pathway 
Subtyping Framework is available at 
github.com/topmist-admin/pathway-subtyping-framework.
```

### What CAN Be Shared

- Aggregate statistics (allele frequencies, cluster sizes)
- Summary scores (mean pathway scores per subtype)
- Analysis code and configurations
- Pathway definitions (GMT files)

### What CANNOT Be Shared

- Individual-level genetic data
- Sample-level assignments (without permission)
- Re-identification risk data

## Code and Reproducibility

### Required Sharing

For reproducibility, share:
1. Analysis code (GitHub repository)
2. Configuration files (YAML)
3. Pathway definitions (GMT files)
4. Software versions (requirements.txt)

### Code Repository Setup

```bash
# Create publication repository
mkdir disease-subtyping-paper
cd disease-subtyping-paper

# Include these files
├── README.md              # Overview and instructions
├── configs/
│   └── analysis_config.yaml
├── pathways/
│   └── disease_pathways.gmt
├── scripts/
│   ├── run_analysis.py
│   └── generate_figures.py
├── requirements.txt
└── LICENSE                # MIT or Apache 2.0
```

### README Template

```markdown
# [Disease] Molecular Subtyping Analysis

Code and data for "[Paper Title]" (Authors, Journal, Year)

## Requirements
- Python 3.9+
- Pathway Subtyping Framework v0.1.0+

## Reproduce Analysis
1. Obtain data access from [Repository]
2. Install dependencies: `pip install -r requirements.txt`
3. Run analysis: `python scripts/run_analysis.py`

## Citation
[Full citation here]
```

## Preprints

### Preprint Servers

- **bioRxiv**: https://www.biorxiv.org/ (most common for genomics)
- **medRxiv**: https://www.medrxiv.org/ (clinical focus)

### Benefits of Preprints

- Establish priority
- Get early feedback
- Accelerate dissemination
- Link to final publication

### Preprint Policy Check

Before posting:
- [ ] Confirm DUA allows preprints
- [ ] Check target journal's preprint policy
- [ ] Ensure all co-authors approve

## Presentations

### Conference Presentations

Allowed under most DUAs:
- Oral presentations
- Poster presentations
- Abstract submissions

**Check DUA for:**
- Pre-approval requirements
- Acknowledgment requirements
- Embargo periods

### Recommended Conferences

| Disease | Conference |
|---------|------------|
| General genetics | ASHG (American Society of Human Genetics) |
| Autism | INSAR (International Society for Autism Research) |
| Psychiatry | SOBP (Society of Biological Psychiatry) |
| Neurology | AAN (American Academy of Neurology) |

## Citation Guidelines

### Citing the Framework

```bibtex
@software{pathway_subtyping_framework,
  author = {Chauhan, Rohit},
  title = {Pathway Subtyping Framework},
  year = {2026},
  url = {https://github.com/topmist-admin/pathway-subtyping-framework}
}
```

### Citing Data Sources

Follow each repository's citation guidelines. Example:

**SFARI:**
```bibtex
@article{sfari_ssc,
  title={The Simons Simplex Collection...},
  author={Fischbach, Gerald D and Lord, Catherine},
  journal={Neuron},
  year={2010}
}
```

## Review and Approval Process

### Internal Review

Before submission:
1. All co-authors review and approve manuscript
2. Data use compliance verified
3. Figures checked for individual data exposure
4. Acknowledgments verified

### External Review (if required)

Some DUAs require:
- Data access committee review
- Funder review
- Embargo compliance

**Allow 2-4 weeks for external review.**

## Checklist Before Submission

- [ ] All DUA requirements met
- [ ] Co-author approval obtained
- [ ] Acknowledgments complete and accurate
- [ ] Data availability statement included
- [ ] Code repository prepared and linked
- [ ] Preprint policy checked (if applicable)
- [ ] Figures contain only aggregate data
- [ ] Supplementary materials prepared
- [ ] Target journal selected
- [ ] Cover letter drafted

---

*Questions about publishing? Contact the project lead.*
