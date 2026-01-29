# Pathway Curation Guide

How to create high-quality pathway definitions for your disease.

---

## Why Pathway Quality Matters

The framework's ability to find meaningful subtypes depends entirely on pathway quality:
- **Good pathways** → Capture convergent biology → Meaningful subtypes
- **Poor pathways** → Random noise → Overfitting or no signal

This guide helps you build robust pathway definitions.

---

## Pathway Sources

### Tier 1: Disease-Specific Databases (Preferred)

| Disease Area | Database | URL |
|--------------|----------|-----|
| Autism | SFARI Gene | https://gene.sfari.org |
| Epilepsy | Epi25 | https://epi25.org |
| Schizophrenia | PGC | https://pgc.unc.edu |
| Intellectual Disability | SysID | https://sysid.cmbi.umcn.nl |
| Parkinson's | PDGene | http://www.pdgene.org |
| Cancer | COSMIC | https://cancer.sanger.ac.uk/cosmic |

### Tier 2: General Pathway Databases

| Database | Best For | GMT Available |
|----------|----------|---------------|
| KEGG | Metabolic, signaling | Yes (via MSigDB) |
| Reactome | Signaling, cellular processes | Yes |
| Gene Ontology | Biological process, molecular function | Yes (via MSigDB) |
| WikiPathways | Community-curated | Yes |

### Tier 3: Literature Mining

For diseases without curated databases:
1. Search PubMed for "[disease] gene" or "[disease] pathway"
2. Extract genes from supplementary tables of GWAS/WES papers
3. Group by functional annotation

---

## Pathway Design Principles

### 1. Optimal Pathway Size

| Size | Pros | Cons |
|------|------|------|
| <5 genes | High specificity | Too sparse, noisy |
| 5-30 genes | Good balance | — |
| 30-100 genes | Captures broad biology | May dilute signal |
| >100 genes | — | Too generic, low specificity |

**Recommendation:** Aim for 10-30 genes per pathway.

### 2. Pathway Independence

Pathways should capture distinct biological processes:

**Bad:**
```
SYNAPSE_1: SHANK3, NRXN1, DLG4, GRIN2B
SYNAPSE_2: SHANK3, NRXN1, NLGN1, SYNGAP1  # 50% overlap!
```

**Good:**
```
SYNAPTIC_SCAFFOLDING: SHANK3, DLG4, HOMER1, GRIP1
SYNAPTIC_ADHESION: NRXN1, NLGN1, CNTNAP2, LRRC4
```

### 3. Biological Coherence

Each pathway should represent a unified biological process:

**Bad (mixing unrelated processes):**
```
MISC_GENES: CHD8, SCN2A, PTEN, BRCA1, TP53
```

**Good (coherent process):**
```
CHROMATIN_REMODELING: CHD8, ARID1B, ASH1L, KMT2A, SETD5
```

### 4. Evidence Requirements

For each gene in a pathway:

| Evidence Level | Description | Weight |
|----------------|-------------|--------|
| Strong | Multiple GWAS hits, functional validation | Include |
| Moderate | Single significant study, biological plausibility | Include |
| Weak | Nominal association only | Consider excluding |
| Candidate | Hypothesized, no genetic evidence | Exclude |

---

## Building a GMT File

### Step-by-Step Process

1. **Create a spreadsheet** with columns: Gene, Pathway, Evidence, Source

2. **Populate genes** from your sources

3. **Assign pathways** based on function

4. **Review for quality:**
   - Remove genes with weak evidence
   - Merge overlapping pathways
   - Split pathways >50 genes

5. **Export to GMT format**

### Python Script for GMT Creation

```python
import pandas as pd

# Load your curated gene list
df = pd.read_csv("curated_genes.csv")
# Expected columns: gene, pathway, evidence

# Filter for sufficient evidence
df = df[df['evidence'].isin(['strong', 'moderate'])]

# Group by pathway and export
with open("my_pathways.gmt", "w") as f:
    for pathway, genes in df.groupby('pathway'):
        gene_list = "\t".join(genes['gene'].tolist())
        f.write(f"{pathway}\thttp://example.org\t{gene_list}\n")
```

---

## Validation Checklist

Before using your pathways:

- [ ] Each pathway has 5-50 genes
- [ ] Total pathways: 5-20
- [ ] Minimal overlap between pathways (<20% shared genes)
- [ ] Each gene has documented disease relevance
- [ ] Gene symbols are standard (HGNC)
- [ ] No duplicate genes within a pathway
- [ ] At least one "housekeeping" pathway as negative control

---

## Testing Pathway Quality

Run the framework with validation gates enabled:

```yaml
validation:
  run_gates: true
  random_genes_iterations: 100
```

**Interpretation:**

| Random Genes ARI | Meaning |
|------------------|---------|
| << 0.15 | Your pathways capture disease biology well |
| ≈ 0.15 | Borderline — pathways may need refinement |
| >> 0.15 | Random genes work as well — pathways not informative |

If random genes perform well, your pathways may be:
- Too generic (e.g., "all brain-expressed genes")
- Poorly curated (noise genes included)
- Not relevant to the disease

---

## Common Mistakes

| Mistake | Why It's Bad | Fix |
|---------|--------------|-----|
| Using all GWAS hits | GWAS hits ≠ pathways | Group by function |
| Too many pathways | Overfitting risk | Merge similar ones |
| Single-gene "pathways" | No aggregation benefit | Minimum 5 genes |
| Including non-coding genes | VCF typically codes protein-coding | Use protein-coding genes |
| Outdated gene symbols | Won't match VCF annotations | Use HGNC symbols |

---

## Resources

### Gene Symbol Mapping
- **HGNC**: https://www.genenames.org/ — Official gene nomenclature
- **Ensembl BioMart**: https://www.ensembl.org/biomart — Bulk symbol conversion

### Pathway Downloads (GMT format)
- **MSigDB**: https://www.gsea-msigdb.org/gsea/msigdb/
- **Enrichr**: https://maayanlab.cloud/Enrichr/#libraries

### Literature Search
- **PubMed**: https://pubmed.ncbi.nlm.nih.gov/
- **GWAS Catalog**: https://www.ebi.ac.uk/gwas/

---

## Example: Curating Epilepsy Pathways

### Sources Used
1. Epi25 collaborative gene list
2. EuroEPINOMICS publications
3. ClinVar epilepsy-associated genes
4. ILAE genetic testing guidelines

### Resulting Pathways
| Pathway | Gene Count | Evidence Basis |
|---------|------------|----------------|
| Sodium channels | 11 | Epi25, ClinVar |
| Potassium channels | 14 | Epi25, functional |
| GABA signaling | 14 | Epi25, ClinVar |
| mTOR pathway | 13 | TSC literature |
| Synaptic vesicle | 12 | Epi25 |

### Quality Metrics
- Average pathway size: 12.8 genes
- Max overlap: 8% (sodium/GABA share SCN1A context)
- Random gene ARI: 0.04 (PASS)

---

> **Remember:** Pathway curation is iterative. Start with literature, test with data, refine based on validation results.
