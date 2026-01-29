# Adapting the Framework for Your Disease

This guide walks you through adapting the Pathway Subtyping Framework for a disease other than autism.

---

## Overview

The framework is disease-agnostic at its core. To adapt it for your disease, you need:

1. **Pathway definitions** — A GMT file with disease-relevant gene sets
2. **Configuration** — A YAML file pointing to your data and pathways
3. **Cohort data** — VCF + phenotype files for your patient population

Everything else (variant processing, clustering, validation) works the same.

---

## Step 1: Curate Disease-Relevant Pathways

### 1.1 Identify Key Biological Pathways

Start by reviewing the literature for your disease:

| Source Type | Examples | How to Use |
|-------------|----------|------------|
| GWAS catalogs | PGC, NHGRI-EBI | Find implicated loci → map to genes → group by function |
| Gene databases | OMIM, ClinVar, disease-specific DBs | Curated gene lists with phenotype associations |
| Pathway databases | KEGG, Reactome, MSigDB | Pre-defined pathways, filter for tissue relevance |
| Review articles | Nature Reviews, Cell | Expert synthesis of pathway biology |

### 1.2 Create Your GMT File

Format:
```
PATHWAY_NAME<tab>URL_OR_DESCRIPTION<tab>GENE1<tab>GENE2<tab>GENE3...
```

Example for a cardiac disease:
```
SARCOMERE	https://example.org	MYH7	MYBPC3	TNNT2	TNNI3	TPM1	ACTC1	MYL2	MYL3
ION_CHANNELS	https://example.org	SCN5A	KCNQ1	KCNH2	KCNE1	KCNE2	KCNJ2	HCN4
CALCIUM_HANDLING	https://example.org	RYR2	CASQ2	PLN	ATP2A2	CALM1	CALM2	CALM3
```

### 1.3 Pathway Design Tips

| Guideline | Reason |
|-----------|--------|
| 5-50 genes per pathway | Too small = noisy; too large = diluted signal |
| 5-20 pathways total | Enough for subtype separation without overfitting |
| Minimize overlap | Highly overlapping pathways inflate correlations |
| Include controls | Add 1-2 "housekeeping" pathways as negative controls |

---

## Step 2: Prepare Your Data

### 2.1 VCF File Requirements

Your VCF must include annotations:

```vcf
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CONSEQUENCE,Number=.,Type=String,Description="Variant consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD phred score">
```

**If your VCF lacks annotations**, run VEP:
```bash
vep -i your_cohort.vcf -o annotated.vcf \
    --symbol --canonical --pick \
    --plugin CADD,/path/to/cadd_scores.tsv.gz
```

### 2.2 Phenotype CSV Requirements

Minimum:
```csv
sample_id,diagnosis
SAMPLE001,AFFECTED
SAMPLE002,AFFECTED
SAMPLE003,CONTROL
```

Recommended (for subtype correlation):
```csv
sample_id,diagnosis,age_onset,severity,subtype_clinical,outcome
SAMPLE001,AFFECTED,25,mild,typeA,stable
SAMPLE002,AFFECTED,18,severe,typeB,progressive
```

---

## Step 3: Create Your Configuration

Copy an example config and modify:

```bash
cp configs/example_autism.yaml configs/my_disease.yaml
```

Key fields to change:

```yaml
pipeline:
  name: "my_disease_subtyping"  # Change this
  output_dir: "outputs/my_disease"  # Change this

data:
  vcf_path: "path/to/your/cohort.vcf"  # Your data
  phenotype_path: "path/to/your/phenotypes.csv"  # Your data
  pathway_db: "data/pathways/my_disease_pathways.gmt"  # Your pathways

clustering:
  n_clusters_range: [2, 6]  # Adjust based on expected subtypes
```

---

## Step 4: Run the Pipeline

```bash
python -m pathway_subtyping --config configs/my_disease.yaml
```

### Expected Output

```
outputs/my_disease/
├── pathway_scores.csv       # Pathway disruption scores per sample
├── subtype_assignments.csv  # Cluster assignments
├── validation_results.json  # Gate pass/fail status
├── report.md                # Human-readable summary
└── figures/
    └── summary.png
```

---

## Step 5: Interpret Validation Gates

### Gate 1: Label Shuffle
**What it tests:** Do clusters arise from real structure or noise?

| Result | Interpretation |
|--------|----------------|
| ARI < 0.15 | PASS — Clusters are not random |
| ARI ≥ 0.15 | FAIL — May be overfitting |

### Gate 2: Random Genes
**What it tests:** Do your pathways capture biology better than random gene sets?

| Result | Interpretation |
|--------|----------------|
| ARI < 0.15 | PASS — Pathways add signal |
| ARI ≥ 0.15 | FAIL — Pathways may not be informative |

### Gate 3: Bootstrap Stability
**What it tests:** Are clusters stable under resampling?

| Result | Interpretation |
|--------|----------------|
| ARI > 0.8 | PASS — Robust subtypes |
| ARI ≤ 0.8 | FAIL — Consider merging clusters or increasing N |

---

## Step 6: Iterate

If validation fails:

| Issue | Possible Fix |
|-------|--------------|
| Random genes perform well | Refine pathway definitions, add more specific genes |
| Bootstrap unstable | Reduce n_clusters, increase cohort size |
| No clear clusters | Disease may not have pathway-based subtypes (valid finding!) |

---

## Step 7: Validate Externally

Before publishing:

1. **Split your cohort** — Train on 70%, validate on 30%
2. **Cross-cohort replication** — Test on independent cohorts if available
3. **Clinical correlation** — Do subtypes associate with phenotypes?

---

## Example: Adapting for Parkinson's Disease

### Pathways (parkinson_pathways.gmt)
```
LYSOSOMAL	https://pdgenetics.org	GBA	SMPD1	CTSD	ATP13A2	TMEM175	GALC	ARSA
MITOCHONDRIAL	https://pdgenetics.org	PINK1	PRKN	DJ1	CHCHD2	VPS13C	POLG	HTRA2
VESICLE_TRAFFICKING	https://pdgenetics.org	LRRK2	VPS35	RAB29	RAB39B	SYNJ1	DNAJC6	AUXILIN
ALPHA_SYNUCLEIN	https://pdgenetics.org	SNCA	SNCAIP	UCHL1	FBXO7	GIGYF2
IMMUNE_INFLAMMATION	https://pdgenetics.org	HLA-DRB1	BST1	GPNMB	TMEM163	LAMP3
```

### Config (configs/example_parkinson.yaml)
```yaml
pipeline:
  name: "parkinson_subtyping"
  output_dir: "outputs/parkinson"
  seed: 42

data:
  vcf_path: "path/to/pd_cohort.vcf"
  phenotype_path: "path/to/pd_phenotypes.csv"
  pathway_db: "data/pathways/parkinson_pathways.gmt"

clustering:
  n_clusters_range: [2, 5]  # PD may have fewer genetic subtypes
```

### Phenotypes to correlate
- Age of onset
- Motor subtype (tremor-dominant vs PIGD)
- Progression rate
- Cognitive decline
- Drug response (L-DOPA)

---

## Getting Help

- **Issues**: https://github.com/topmist-admin/pathway-subtyping-framework/issues
- **Discussions**: https://github.com/topmist-admin/pathway-subtyping-framework/discussions
- **Email**: info@topmist.com

---

> **RESEARCH USE ONLY** — This framework is for hypothesis generation. Subtypes require clinical validation before any diagnostic or therapeutic application.
