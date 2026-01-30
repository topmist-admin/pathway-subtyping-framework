# gnomAD Usage Guide

**Genome Aggregation Database (gnomAD)** provides variant frequency data from ~140,000 exomes and ~76,000 genomes. This is an **open-access** resource.

## Overview

- **URL**: https://gnomad.broadinstitute.org/
- **Access**: Open - no application required
- **Current Version**: gnomAD v4.x
- **Use Case**: Filter rare variants, assess pathogenicity

## How We Use gnomAD

### 1. Variant Frequency Filtering

Filter cohort variants to rare variants (MAF < 0.1%):

```python
# Example: Filter variants by gnomAD frequency
def is_rare_variant(variant, max_af=0.001):
    """Check if variant is rare based on gnomAD frequency."""
    gnomad_af = variant.get('gnomad_af', 0)
    return gnomad_af < max_af or gnomad_af is None
```

### 2. Population-Specific Frequencies

gnomAD provides frequencies for multiple populations:
- `AF_afr` - African/African American
- `AF_amr` - Latino/Admixed American
- `AF_asj` - Ashkenazi Jewish
- `AF_eas` - East Asian
- `AF_fin` - Finnish
- `AF_nfe` - Non-Finnish European
- `AF_sas` - South Asian

**Best Practice**: Use population-matched frequencies when available.

### 3. Constraint Metrics

gnomAD provides gene-level constraint metrics:

| Metric | Description | Use |
|--------|-------------|-----|
| pLI | Probability of LoF intolerance | Higher = more constrained |
| LOEUF | LoF observed/expected upper bound | Lower = more constrained |
| mis_z | Missense Z-score | Higher = fewer missense than expected |

**Use in framework**: Weight genes by constraint when computing burden.

## Accessing gnomAD Data

### Web Interface
Browse variants at https://gnomad.broadinstitute.org/

### Programmatic Access

**Download VCFs:**
```bash
# gnomAD v4 exomes (example - check current URLs)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz
```

**Query API:**
```python
import requests

def get_variant_info(chrom, pos, ref, alt):
    """Query gnomAD API for variant information."""
    url = f"https://gnomad.broadinstitute.org/api"
    query = """
    query {
      variant(variantId: "%s-%s-%s-%s", dataset: gnomad_r4) {
        exome {
          ac
          an
          af
        }
        genome {
          ac
          an
          af
        }
      }
    }
    """ % (chrom, pos, ref, alt)
    
    response = requests.post(url, json={'query': query})
    return response.json()
```

### Pre-Annotated Files

For pipeline use, download pre-computed annotation files:
- Sites VCF (all variants with frequencies)
- Constraint metrics by gene
- Coverage files

## Integration with Framework

### In Config File
```yaml
variant_filter:
  max_maf: 0.001  # 0.1% in gnomAD
  gnomad_version: "4.0"
  population: "nfe"  # or "all" for global AF
```

### In Pipeline
The framework automatically filters variants by gnomAD frequency if VCF is annotated with gnomAD fields.

## Citation

When using gnomAD, cite:
```
Chen, S., et al. (2024). "A genomic mutational constraint map using variation 
in 76,156 human genomes." Nature 625, 92-100.
```

## Resources

- **Documentation**: https://gnomad.broadinstitute.org/help
- **Downloads**: https://gnomad.broadinstitute.org/downloads
- **GitHub**: https://github.com/broadinstitute/gnomad_methods
- **Blog**: https://gnomad.broadinstitute.org/news

---

*No application needed - start using gnomAD today!*
