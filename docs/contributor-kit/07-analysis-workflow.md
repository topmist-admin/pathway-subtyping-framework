# Analysis Workflow Guide

This guide provides a step-by-step workflow for running molecular subtyping analysis using the Pathway Subtyping Framework.

## Workflow Overview

```
┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│  1. Data    │───▶│ 2. Pathway  │───▶│ 3. Config   │───▶│ 4. Execute  │
│  Preparation│    │  Curation   │    │  Setup      │    │  Pipeline   │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
                                                                │
                                                                ▼
┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│  8. Report  │◀───│ 7. Phenotype│◀───│ 6. Interpret│◀───│ 5. Validate │
│  & Publish  │    │  Analysis   │    │  Results    │    │  Subtypes   │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
```

---

## Phase 1: Data Preparation

### 1.1 Obtain Data Access

Complete data access application per [02-data-access-guide.md](02-data-access-guide.md).

### 1.2 Download and Organize Data

```bash
# Create project directory
mkdir -p ~/projects/my_disease_study/{data,outputs,configs}
cd ~/projects/my_disease_study

# Download data (method varies by repository)
# Example for dbGaP:
prefetch SRR12345678
fasterq-dump SRR12345678
```

### 1.3 Prepare VCF File

**Required VCF annotations:**
- `GENE` - Gene symbol (HGNC)
- `CONSEQUENCE` - Variant effect (VEP terms)
- `CADD` - Deleteriousness score

**Annotate with VEP:**
```bash
vep --input_file raw_variants.vcf \
    --output_file annotated.vcf \
    --format vcf --vcf \
    --symbol --canonical \
    --plugin CADD,cadd_snvs.tsv.gz \
    --fields "SYMBOL,Consequence,CADD_PHRED"
```

**Convert to required format:**
```bash
# Extract needed fields to INFO
bcftools annotate \
    -c INFO/GENE:=INFO/SYMBOL \
    -c INFO/CONSEQUENCE:=INFO/CSQ \
    annotated.vcf.gz -Oz -o final.vcf.gz
```

### 1.4 Prepare Phenotype File

Create CSV with required columns:

```csv
sample_id,age,sex,diagnosis,severity_score
SAMPLE_001,8,M,ASD,moderate
SAMPLE_002,12,F,ASD,severe
SAMPLE_003,6,M,ASD,mild
```

**Required columns:**
- `sample_id` - Must match VCF sample names

**Optional columns:**
- Demographics (age, sex)
- Clinical scores
- Subtype labels (if known, for validation)

---

## Phase 2: Pathway Curation

### 2.1 Review Existing Pathways

Check if pathways exist for your disease:

```bash
ls data/pathways/
# autism_pathways.gmt
# schizophrenia_pathways.gmt
# ...
```

### 2.2 Create Disease-Specific GMT

If new pathways needed, follow the [pathway-curation-guide.md](../guides/pathway-curation-guide.md).

**GMT format:**
```
PATHWAY_NAME    Description    GENE1    GENE2    GENE3    ...
SYNAPTIC    Synaptic transmission    SHANK3    NRXN1    SYNGAP1    ...
CHROMATIN    Chromatin remodeling    CHD8    ARID1B    KMT2A    ...
```

### 2.3 Validate Pathways

```python
# Check pathway file
from pathway_subtyping.utils import load_gmt

pathways = load_gmt("data/pathways/my_disease_pathways.gmt")

for name, genes in pathways.items():
    print(f"{name}: {len(genes)} genes")
    
# Verify genes are in VCF
# Check for overlap between pathways
```

---

## Phase 3: Configuration Setup

### 3.1 Create Config File

Copy and modify a template:

```bash
cp configs/example_autism.yaml configs/my_disease_study.yaml
```

### 3.2 Edit Configuration

```yaml
# configs/my_disease_study.yaml

pipeline:
  name: "my_disease_subtyping"
  output_dir: "outputs/my_disease_study"
  seed: 42
  verbose: true

data:
  vcf_path: "data/my_cohort.vcf.gz"
  phenotype_path: "data/my_phenotypes.csv"
  pathway_db: "data/pathways/my_disease_pathways.gmt"

variant_filter:
  min_depth: 10
  min_gq: 20
  max_maf: 0.001
  consequences:
    - "frameshift_variant"
    - "stop_gained"
    - "stop_lost"
    - "splice_acceptor_variant"
    - "splice_donor_variant"
    - "missense_variant"
  min_cadd: 25

burden:
  lof_weight: 2.0
  missense_weight: 1.0
  normalize: true

clustering:
  method: "gmm"
  n_clusters: null  # Auto-select via BIC
  n_clusters_range: [2, 8]
  covariance_type: "full"

validation:
  run_gates: true
  label_shuffle_iterations: 100
  random_genes_iterations: 100
  bootstrap_iterations: 100
  stability_threshold: 0.8
  negative_control_threshold: 0.15

output:
  save_pathway_scores: true
  save_cluster_assignments: true
  generate_report: true
```

### 3.3 Validate Configuration

```bash
# Dry run to check config
psf --config configs/my_disease_study.yaml --dry-run
```

---

## Phase 4: Execute Pipeline

### 4.1 Run Full Pipeline

```bash
# Activate environment
source venv/bin/activate

# Run pipeline
psf --config configs/my_disease_study.yaml --verbose

# Or with Python API
python -c "
from pathway_subtyping import DemoPipeline, PipelineConfig
config = PipelineConfig.from_yaml('configs/my_disease_study.yaml')
pipeline = DemoPipeline(config)
pipeline.run()
"
```

### 4.2 Monitor Progress

The pipeline outputs progress:
```
[INFO] Loading input data...
[INFO] Loaded: 5000 variants, 500 samples, 6 pathways
[INFO] Computing gene burdens...
[INFO] Computing pathway scores...
[INFO] Running GMM clustering (k=2-8)...
[INFO]   k=2: BIC=-1234.5
[INFO]   k=3: BIC=-1567.8
[INFO]   k=4: BIC=-1823.2 <- optimal
[INFO] Running validation gates...
[INFO] Generating outputs...
[INFO] Pipeline completed in 45.2 seconds
```

### 4.3 Check Outputs

```bash
ls -la outputs/my_disease_study/

# Expected files:
# pathway_scores.csv      - Pathway scores per sample
# subtype_assignments.csv - Cluster labels
# report.json             - Machine-readable results
# report.md               - Human-readable report
# figures/summary.png     - Visualization
```

---

## Phase 5: Validate Subtypes

### 5.1 Review Validation Gates

Open `outputs/my_disease_study/report.md`:

```markdown
## Validation Gates

| Test | Status | Value | Threshold |
|------|--------|-------|-----------|
| Label Shuffle | PASS | 0.03 | < 0.15 |
| Random Genes | PASS | 0.05 | < 0.15 |
| Bootstrap | PASS | 0.92 | >= 0.80 |
```

**All three gates must pass for valid subtypes.**

### 5.2 Troubleshooting Failed Gates

**Label Shuffle Failed (ARI > 0.15):**
- Clusters may be driven by confounders (batch, ancestry)
- Check for population stratification
- Review variant QC

**Random Genes Failed (ARI > 0.15):**
- Subtypes not specific to pathway biology
- May need different pathways
- Check pathway gene overlap

**Bootstrap Failed (ARI < 0.80):**
- Subtypes are unstable
- May need larger sample size
- Try different k range

---

## Phase 6: Interpret Results

### 6.1 Examine Cluster Assignments

```python
import pandas as pd

assignments = pd.read_csv("outputs/my_disease_study/subtype_assignments.csv")

# Cluster distribution
print(assignments['cluster_label'].value_counts())

# Confidence scores
print(assignments['confidence'].describe())
```

### 6.2 Analyze Pathway Profiles

```python
import seaborn as sns
import matplotlib.pyplot as plt

scores = pd.read_csv("outputs/my_disease_study/pathway_scores.csv", index_col=0)
assignments = pd.read_csv("outputs/my_disease_study/subtype_assignments.csv")

# Merge and compute means
scores['cluster'] = assignments.set_index('sample_id')['cluster_label']
cluster_means = scores.groupby('cluster').mean()

# Heatmap
sns.heatmap(cluster_means, cmap='RdBu_r', center=0, annot=True)
plt.title("Mean Pathway Scores by Subtype")
plt.show()
```

### 6.3 Identify Dominant Pathways

For each subtype, identify the pathway with highest burden:

```python
for cluster in cluster_means.index:
    dominant = cluster_means.loc[cluster].idxmax()
    score = cluster_means.loc[cluster, dominant]
    print(f"Subtype {cluster}: Dominant pathway = {dominant} (z={score:.2f})")
```

---

## Phase 7: Phenotype Analysis

### 7.1 Correlate with Clinical Features

```python
# Merge assignments with phenotypes
phenotypes = pd.read_csv("data/my_phenotypes.csv")
merged = assignments.merge(phenotypes, on='sample_id')

# Compare features by subtype
for feature in ['age', 'severity_score']:
    print(f"\n{feature} by subtype:")
    print(merged.groupby('cluster_label')[feature].describe())
```

### 7.2 Statistical Testing

```python
from scipy import stats

# ANOVA for continuous variables
for feature in ['age', 'severity_score']:
    groups = [merged[merged['cluster_label'] == c][feature] 
              for c in merged['cluster_label'].unique()]
    f_stat, p_val = stats.f_oneway(*groups)
    print(f"{feature}: F={f_stat:.2f}, p={p_val:.4f}")

# Chi-square for categorical variables
for feature in ['sex', 'diagnosis']:
    contingency = pd.crosstab(merged['cluster_label'], merged[feature])
    chi2, p_val, dof, expected = stats.chi2_contingency(contingency)
    print(f"{feature}: χ²={chi2:.2f}, p={p_val:.4f}")
```

### 7.3 Visualize Phenotype Differences

```python
# Box plots for continuous features
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sns.boxplot(data=merged, x='cluster_label', y='severity_score', ax=axes[0])
axes[0].set_title("Severity Score by Subtype")

sns.boxplot(data=merged, x='cluster_label', y='age', ax=axes[1])
axes[1].set_title("Age by Subtype")

plt.tight_layout()
plt.savefig("outputs/my_disease_study/phenotype_comparison.png")
```

---

## Phase 8: Report and Publish

### 8.1 Generate Final Report

The framework generates `report.md` automatically. Supplement with:
- Phenotype analysis results
- Clinical interpretation
- Comparison to prior work

### 8.2 Prepare Figures

Standard figures for publication:
1. PCA/UMAP colored by subtype
2. Pathway score heatmap by subtype
3. Phenotype comparisons
4. Validation gate results

### 8.3 Code and Data Sharing

```bash
# Create reproducibility package
mkdir -p publication/
cp configs/my_disease_study.yaml publication/
cp data/pathways/my_disease_pathways.gmt publication/
cp -r outputs/my_disease_study/figures publication/

# Document software versions
pip freeze > publication/requirements.txt
```

### 8.4 Submit for Publication

See [08-publishing-guidelines.md](08-publishing-guidelines.md) for:
- Journal recommendations
- Co-authorship guidelines
- Acknowledgment requirements
- Data availability statements

---

## Quick Reference

### CLI Commands

```bash
# Run pipeline
psf --config config.yaml

# Run with verbose output
psf --config config.yaml --verbose

# Dry run (validate config only)
psf --config config.yaml --dry-run

# Specify output directory
psf --config config.yaml --output outputs/custom/

# Set random seed
psf --config config.yaml --seed 123
```

### Common Issues

| Issue | Solution |
|-------|----------|
| VCF parsing error | Check annotation format, ensure GENE/CONSEQUENCE/CADD in INFO |
| Memory error | Use smaller chunks, increase instance size |
| No clusters found | Check variant counts, may need lower MAF threshold |
| Validation fails | Review data quality, sample size, pathway definitions |

---

*Questions about the workflow? Contact the project lead or post in the team channel.*
