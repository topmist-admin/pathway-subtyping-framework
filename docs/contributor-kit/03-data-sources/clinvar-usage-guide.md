# ClinVar Usage Guide

**ClinVar** is NCBI's public archive of variant-disease relationships with clinical significance interpretations. This is an **open-access** resource.

## Overview

- **URL**: https://www.ncbi.nlm.nih.gov/clinvar/
- **Access**: Open - no application required
- **Content**: >2 million variant submissions
- **Use Case**: Identify known pathogenic variants, validate findings

## How We Use ClinVar

### 1. Identify Known Pathogenic Variants

Flag variants previously reported as pathogenic:

```python
# Clinical significance categories
PATHOGENIC = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
BENIGN = ['Benign', 'Likely benign', 'Benign/Likely benign']
VUS = ['Uncertain significance']

def get_clinical_significance(variant):
    """Check ClinVar classification."""
    clnsig = variant.get('CLNSIG', '')
    if any(p in clnsig for p in PATHOGENIC):
        return 'pathogenic'
    elif any(b in clnsig for b in BENIGN):
        return 'benign'
    else:
        return 'vus_or_unknown'
```

### 2. Gene-Disease Associations

ClinVar provides evidence linking genes to diseases:
- Gene symbols with disease variants
- Condition names (MedGen, OMIM)
- Inheritance patterns

### 3. Variant Evidence

For each variant, ClinVar reports:
- Review status (★★★★ = expert panel)
- Number of submitters
- Date of last evaluation
- Functional evidence

## Accessing ClinVar Data

### Web Interface
Search variants at https://www.ncbi.nlm.nih.gov/clinvar/

### Download Files

**Full ClinVar VCF:**
```bash
# Download latest ClinVar VCF
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

**Gene-specific summary:**
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt
```

**Variant summary:**
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
```

### Programmatic Access (E-utilities)

```python
from Bio import Entrez

Entrez.email = "your.email@example.com"

def search_clinvar_gene(gene_symbol):
    """Search ClinVar for variants in a gene."""
    handle = Entrez.esearch(
        db="clinvar",
        term=f"{gene_symbol}[gene] AND pathogenic[clinsig]",
        retmax=1000
    )
    record = Entrez.read(handle)
    return record['IdList']
```

## Integration with Framework

### Pre-Annotation
Annotate your VCF with ClinVar before running the pipeline:

```bash
# Using bcftools
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN \
    input.vcf.gz \
    -o annotated.vcf.gz
```

### In Analysis

```python
# Weight variants by ClinVar status
def get_variant_weight(variant, clinvar_boost=1.5):
    """Increase weight for ClinVar pathogenic variants."""
    base_weight = get_consequence_weight(variant)
    
    if is_clinvar_pathogenic(variant):
        return base_weight * clinvar_boost
    return base_weight
```

### For Pathway Curation

Use ClinVar to identify disease-relevant genes:

```python
def get_disease_genes(disease_term, min_pathogenic=3):
    """Get genes with multiple pathogenic variants for a disease."""
    # Query ClinVar for disease
    # Filter genes with sufficient pathogenic evidence
    # Return gene list for pathway curation
    pass
```

## Clinical Significance Levels

| Classification | Stars | Meaning |
|---------------|-------|---------|
| Practice guideline | ★★★★ | Expert panel review |
| Reviewed by expert panel | ★★★★ | Expert consensus |
| Criteria provided, multiple submitters | ★★★ | Consistent evidence |
| Criteria provided, single submitter | ★★ | One lab with criteria |
| Criteria provided, conflicting | ★ | Disagreement |
| No assertion criteria | - | No criteria provided |

**Best Practice**: Prioritize variants with ≥★★★ review status.

## Citation

When using ClinVar, cite:
```
Landrum MJ, et al. (2024). "ClinVar: improvements to accessing data." 
Nucleic Acids Research 52(D1):D733-D741.
```

## Resources

- **FTP**: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
- **Documentation**: https://www.ncbi.nlm.nih.gov/clinvar/docs/
- **Submission**: https://www.ncbi.nlm.nih.gov/clinvar/docs/submit/
- **API**: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/

---

*No application needed - ClinVar is freely available!*
