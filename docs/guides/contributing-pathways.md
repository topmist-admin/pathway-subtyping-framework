# Contributing Pathway Definitions

This guide explains how to contribute new disease-specific pathway definitions to the Pathway Subtyping Framework.

## Overview

Pathway definitions are the foundation of this framework. High-quality, well-curated pathways enable meaningful subtype discovery across genetically heterogeneous diseases.

## GMT File Format

Pathway definitions use the Gene Matrix Transposed (GMT) format, a tab-delimited text format:

```
PATHWAY_NAME<tab>SOURCE_URL<tab>GENE1<tab>GENE2<tab>GENE3...
```

### Example

```
SYNAPTIC_TRANSMISSION	https://example.org/synapse	SHANK2	SHANK3	SYNGAP1	GRIN2A	GRIN2B
CHROMATIN_REMODELING	https://example.org/chromatin	CHD8	ARID1B	KMT2A	SETD5	ASH1L
```

### File Header

Include a header with metadata:

```
# Disease Name Pathway Definitions
# Sources: [List primary sources]
# Version: X.Y (YYYY-MM-DD)
# Genes curated from: [Key publications with citations]
```

## Quality Standards

### Gene Selection Criteria

1. **Evidence-based**: Include genes with strong genetic evidence
   - Genome-wide significant associations (p < 5×10⁻⁸)
   - Exome-wide significant burden (p < 2.5×10⁻⁶)
   - Replicated findings across independent cohorts

2. **Functional relevance**: Genes should have biological connection to the pathway
   - Literature-supported pathway membership
   - Protein-protein interaction evidence
   - Expression in relevant tissues

3. **Avoid overlap bias**: Minimize excessive gene overlap between pathways
   - Some overlap is expected (e.g., GRIN2A in both glutamate and synaptic)
   - Flag pathways with >50% overlap for review

### Pathway Size Guidelines

| Category | Gene Count | Notes |
|----------|------------|-------|
| Minimum | 10 genes | Smaller pathways lack statistical power |
| Recommended | 15-30 genes | Optimal balance of specificity and power |
| Maximum | 50 genes | Larger pathways may be too broad |

### Naming Conventions

- Use `SCREAMING_SNAKE_CASE` for pathway names
- Be specific but concise: `VOLTAGE_GATED_SODIUM_CHANNELS` not `SODIUM`
- Avoid abbreviations unless universally understood (e.g., `GABA`, `MTOR`)

## Primary Sources

### Recommended Databases

| Database | URL | Best For |
|----------|-----|----------|
| KEGG | https://www.genome.jp/kegg/ | Metabolic & signaling pathways |
| Reactome | https://reactome.org/ | Detailed pathway hierarchies |
| Gene Ontology | http://geneontology.org/ | Functional annotations |
| MSigDB | https://www.gsea-msigdb.org/ | Curated gene sets |
| DisGeNET | https://www.disgenet.org/ | Disease-gene associations |
| SysID | https://sysid.cmbi.umcn.nl/ | Intellectual disability |
| SFARI Gene | https://gene.sfari.org/ | Autism spectrum disorder |

### Disease-Specific Resources

**Neuropsychiatric Disorders**
- PGC GWAS: https://pgc.unc.edu/
- SCHEMA: https://schema.broadinstitute.org/
- Epi25: https://epi25.org/

**Neurodevelopmental Disorders**
- DDD Study: https://www.ddduk.org/
- DECIPHER: https://www.deciphergenomics.org/

## Contribution Workflow

### 1. Fork and Clone

```bash
git clone https://github.com/YOUR_USERNAME/pathway-subtyping-framework
cd pathway-subtyping-framework
git checkout -b add-DISEASE-pathways
```

### 2. Create Pathway File

```bash
# Create your pathway file
touch data/pathways/your_disease_pathways.gmt

# Use an existing file as template
cp data/pathways/autism_pathways.gmt data/pathways/your_disease_pathways.gmt
```

### 3. Populate Pathways

Follow this process for each pathway:

1. **Identify core genes**: Start with GWAS/exome hits for your disease
2. **Expand with biology**: Add functionally related genes from pathway databases
3. **Validate**: Check that genes are expressed in relevant tissue
4. **Document**: Record the source for each gene addition

### 4. Validate Your Pathways

Run the validation script:

```bash
# Check GMT format
python -c "
from pathway_subtyping.utils.gmt_parser import load_gmt
pathways = load_gmt('data/pathways/your_disease_pathways.gmt')
print(f'Loaded {len(pathways)} pathways')
for name, genes in pathways.items():
    print(f'  {name}: {len(genes)} genes')
"
```

### 5. Test with Synthetic Data

Create synthetic test data for your pathways:

```bash
# Generate synthetic cohort
python scripts/generate_synthetic_cohort.py \
    --pathways data/pathways/your_disease_pathways.gmt \
    --output data/sample/your_disease_synthetic.vcf
```

### 6. Submit Pull Request

```bash
git add data/pathways/your_disease_pathways.gmt
git commit -m "feat(pathways): add YOUR_DISEASE pathway definitions"
git push origin add-DISEASE-pathways
```

Include in your PR:
- [ ] Pathway GMT file following format standards
- [ ] Header with sources and version
- [ ] At least 8 pathways with 10+ genes each
- [ ] No syntax errors (validated with parser)
- [ ] Brief description of curation methodology

## Review Checklist

Maintainers will review contributions against:

- [ ] GMT format validity
- [ ] Gene symbol correctness (HGNC approved)
- [ ] Source documentation
- [ ] Pathway size within guidelines
- [ ] Biological coherence of gene groupings
- [ ] Minimal excessive overlap between pathways
- [ ] Successful test run with synthetic data

## Example: Adding Bipolar Disorder Pathways

```gmt
# Bipolar Disorder Pathway Definitions
# Sources: PGC-BD3, BipEx consortium
# Version: 1.0 (2026-01-29)
# Genes curated from: Mullins et al. 2021 (Nature Genetics)

CALCIUM_SIGNALING	https://pgc.unc.edu/bipolar	CACNA1C	CACNA1D	CACNB2	ANK3	CAMK2A...
CIRCADIAN_RHYTHM	https://pgc.unc.edu/bipolar	CLOCK	ARNTL	PER1	PER2	CRY1	CRY2...
LITHIUM_RESPONSE	https://pgc.unc.edu/bipolar	GSK3B	CREB1	BDNF	NTRK2	IMPA1...
```

## Getting Help

- **Questions**: Open a GitHub Discussion
- **Issues**: File a bug report for pathway errors
- **Slack**: Join #pathway-curation channel

## Acknowledgments

Contributors who add high-quality pathway definitions will be:
- Listed in CONTRIBUTORS.md
- Acknowledged in publications using the framework
- Invited to collaborate on disease-specific analyses

---

Thank you for contributing to open science and helping advance precision medicine research!
