# Pathway Definition Files

This directory contains GMT (Gene Matrix Transposed) files that define biological pathways for each disease.

## GMT Format

Each line represents one pathway:

```
PATHWAY_NAME<tab>DESCRIPTION_OR_URL<tab>GENE1<tab>GENE2<tab>GENE3<tab>...
```

**Example:**
```
SYNAPTIC_TRANSMISSION	http://example.org/synaptic	SHANK3	NRXN1	NLGN1	SYNGAP1	GRIN2B
CHROMATIN_REMODELING	http://example.org/chromatin	CHD8	ARID1B	ASH1L	KMT2A	SETD5
```

## Available Pathway Files

| File | Disease | Status | Gene Count |
|------|---------|--------|------------|
| `autism_pathways.gmt` | Autism Spectrum Disorder | Validated | ~200 genes |
| `schizophrenia_pathways.gmt` | Schizophrenia | Template | ~250 genes |
| `epilepsy_pathways.gmt` | Epilepsy | Template | ~200 genes |
| `intellectual_disability_pathways.gmt` | Intellectual Disability | Template | ~350 genes |
| `parkinsons_pathways.gmt` | Parkinson's Disease | Template | ~280 genes |
| `bipolar_pathways.gmt` | Bipolar Disorder | Template | ~290 genes |

## Creating Your Own Pathway File

### Step 1: Identify Disease-Relevant Pathways

Sources for pathway curation:
- **KEGG**: https://www.kegg.jp/kegg/pathway.html
- **Reactome**: https://reactome.org/
- **MSigDB**: https://www.gsea-msigdb.org/gsea/msigdb/
- **Gene Ontology**: http://geneontology.org/
- **Disease-specific databases** (e.g., SFARI Gene for autism, Epi25 for epilepsy)

### Step 2: Select Genes for Each Pathway

Criteria for gene inclusion:
1. **Literature support**: Published association with the disease
2. **Functional relevance**: Gene function relates to pathway biology
3. **Expression data**: Expressed in relevant tissues (e.g., brain for neuropsychiatric)

### Step 3: Format as GMT

```bash
# Example: creating a simple GMT file
cat > my_disease_pathways.gmt << 'EOF'
PATHWAY_A	description	GENE1	GENE2	GENE3
PATHWAY_B	description	GENE4	GENE5	GENE6
EOF
```

### Step 4: Validate Your Pathways

Run the framework with validation gates enabled:
- If **random gene sets** perform as well as your pathways → your pathways may not be capturing real biology
- If **bootstrap stability** is low → consider merging similar pathways or removing noisy genes

## Pathway Curation Tips

1. **Start broad, refine later** — Begin with established pathway databases, then customize
2. **Balance pathway sizes** — Very large pathways (>100 genes) may dilute signal; very small (<5 genes) may be noisy
3. **Consider overlap** — Some overlap between pathways is expected, but excessive overlap can inflate correlations
4. **Document your sources** — Keep track of where each gene came from for reproducibility

## Converting from Other Formats

### From Gene List (one gene per line)

```python
# Convert gene list to single-pathway GMT
genes = open("gene_list.txt").read().strip().split("\n")
print(f"MY_PATHWAY\thttp://example.org\t" + "\t".join(genes))
```

### From Reactome/KEGG Export

Most pathway databases offer GMT export. Look for "Download" → "GMT format" options.

### From MSigDB

1. Go to https://www.gsea-msigdb.org/gsea/msigdb/
2. Select relevant gene sets (e.g., C2: curated, C5: ontology)
3. Download as GMT

## Recommended Pathways by Disease

### Autism / Intellectual Disability
- Synaptic transmission
- Chromatin remodeling
- Transcriptional regulation
- Wnt signaling
- mTOR signaling
- Neuronal migration

### Schizophrenia
- Synaptic transmission
- Dopamine signaling
- Glutamate signaling
- Immune/complement
- Calcium signaling
- Voltage-gated channels

### Epilepsy
- Ion channels (Na+, K+, Ca2+)
- GABA signaling
- Glutamate signaling
- mTOR pathway
- Synaptic vesicle cycle

### Parkinson's Disease
- Alpha-synuclein aggregation
- Mitochondrial function
- Autophagy-lysosomal pathway
- Dopamine metabolism
- Endolysosomal trafficking
- Immune/inflammation
- Oxidative stress

### Bipolar Disorder
- Calcium signaling
- Circadian rhythm
- Glutamate/GABA signaling
- WNT/GSK3 signaling (lithium target)
- Inositol phosphate pathway
- Synaptic transmission
- HPA stress response

---

See also: [Pathway Curation Guide](../docs/guides/pathway-curation-guide.md)
