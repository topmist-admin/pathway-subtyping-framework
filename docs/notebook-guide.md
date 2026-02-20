# Notebook Guide

> Complete guide to the 16 analysis notebooks, their execution order, and how outputs flow between them.

---

## Notebook Inventory

### Tier 1: Tutorials (Synthetic Data)

These notebooks use synthetic data bundled with the framework. They run independently in any order and require no external downloads.

| # | Notebook | Topic | Runtime |
|---|----------|-------|---------|
| 00 | [quick_demo](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/00_quick_demo.ipynb) | 60-second end-to-end demo | <1 min |
| 01 | [getting_started](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/01_getting_started.ipynb) | Installation, basic pipeline, validation gates | 5 min |
| 02 | [expression_scoring](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/02_expression_scoring.ipynb) | ssGSEA, GSVA, mean-Z scoring comparison | 5 min |
| 03 | [multi_omic_fusion](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/03_multi_omic_fusion.ipynb) | Multi-modal data fusion (VCF + expression + scRNA-seq) | 5 min |
| 04 | [deconvolution](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/04_deconvolution.ipynb) | Cell-type deconvolution from bulk RNA-seq | 5 min |
| 05 | [visualization](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/05_visualization.ipynb) | PCA, t-SNE, UMAP, interactive Plotly reports | 5 min |
| 06 | [ancestry_batch_correction](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/06_ancestry_batch_correction.ipynb) | Population stratification and batch effects | 5 min |
| 07 | [sensitivity_analysis](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/07_sensitivity_analysis.ipynb) | Parameter robustness testing | 5 min |
| 08 | [characterization](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/08_characterization.ipynb) | Subtype profiling, heatmaps, gene contributions | 5 min |
| 09 | [signaling_databases](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/09_signaling_databases.ipynb) | CellPhoneDB and CellChatDB integration | 5 min |

### Tier 2: Real Data Validation (GEO Datasets)

These notebooks download real transcriptomics data from NCBI GEO and reproduce the manuscript analyses. They must be run in a specific order because later notebooks consume outputs from earlier ones.

| # | Notebook | Dataset | Tissue | N | Platform | Depends On |
|---|----------|---------|--------|---|----------|------------|
| 10 | [geo_autism_bulk](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/10_geo_autism_bulk.ipynb) | GSE28521 | Brain (frontal + temporal cortex) | 79 | Affymetrix microarray | None |
| 11 | [geo_autism_validation](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/11_geo_autism_validation.ipynb) | GSE64018 | Brain (temporal cortex) | 24 | Illumina RNA-seq | **10** |
| 12 | [geo_schizophrenia_bulk](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/12_geo_schizophrenia_bulk.ipynb) | GSE80655 | Brain (3 regions) | 281 | Illumina RNA-seq | **10** (optional) |
| 12b | [null_ari_permutation](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/12b_null_ari_permutation.ipynb) | GSE80655 | Brain (3 regions) | 141 | (reuses 12 outputs) | **12** |
| 13 | [geo_blood_ados](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/13_geo_blood_ados.ipynb) | GSE111175 | Blood (leukocytes) | 141 | Illumina BeadChip | **10** (optional) |
| 14 | [geo_blood_large_cohort](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/14_geo_blood_large_cohort.ipynb) | GSE18123 | Blood | 285 | Affymetrix (2 platforms) | **13**, **10** (optional) |

---

## Dependency Diagram

```
Tier 1 (Tutorials 00-09)          Tier 2 (Real Data 10-14)
  All standalone                   Must follow arrows
  Any order                        ─────────────────────

                                   10: GSE28521 Autism Brain
                                       (DISCOVERY - run first)
                                        │
                          ┌─────────────┼─────────────┬──────────────┐
                          ▼             ▼             ▼              ▼
                    11: GSE64018  12: GSE80655  13: GSE111175  (optional
                    Cross-cohort  Schizophrenia Blood + ADOS    inputs)
                    validation    Brain              │
                                       │             ▼
                                       ▼        14: GSE18123
                                  12b: Null ARI Blood Large Cohort
                                  Permutation   (needs 13 + 10)
```

**Arrow meaning:** The target notebook loads output files produced by the source notebook. If the source has not been run, the target will skip that analysis section (with a warning) or attempt to download pre-computed results from GitHub.

---

## Recommended Execution Sequence

To reproduce all manuscript analyses from scratch, run the Tier 2 notebooks in this order:

### Step 1: Notebook 10 -- GSE28521 Autism Brain (Foundation)

This is the discovery analysis. All subsequent notebooks reference its outputs.

```
Run first. No prerequisites.
```

**Produces:**
```
outputs/gse28521/
├── pathway_scores.csv                    # 79 samples x 15 pathways
├── gene_expression_processed.csv         # Full processed expression matrix
└── frontal_cortex/
    ├── fc_pathway_scores.csv             # 32 frontal cortex samples x 15 pathways
    ├── fc_sample_metadata_with_subtypes.csv  # Subtype assignments (k=2)
    └── fc_results_summary.json           # Validation gate results
```

**Key result:** GABA-Collapsed subtype (56% of ASD, 0% controls), Cohen's d = 3.21

### Step 2: Notebook 11 -- GSE64018 Cross-Cohort Validation

Cross-platform (RNA-seq vs microarray) and cross-region (temporal vs frontal cortex) replication.

**Requires from Notebook 10:**
```
outputs/gse28521/frontal_cortex/fc_pathway_scores.csv
outputs/gse28521/frontal_cortex/fc_sample_metadata_with_subtypes.csv
```

**Produces:**
```
outputs/gse64018/frontal_cortex/
├── gse64018_fc_pathway_scores.csv
├── gse64018_fc_sample_metadata_with_subtypes.csv
└── gse64018_fc_results_summary.json
```

**Key result:** Disease-enriched subtype (83% ASD) with same top 3 disrupted pathways

### Step 3: Notebook 12 -- GSE80655 Schizophrenia Cross-Disease

The primary validated result. Largest brain dataset, all 3 validation gates pass.

**Requires from Notebook 10 (optional, for cross-disease comparison):**
```
outputs/gse28521/frontal_cortex/fc_pathway_scores.csv
```

**Produces:**
```
outputs/gse80655/
├── pathway_scores_scz.csv               # 141 SCZ+CTL samples x 14 SCZ pathways
├── pathway_scores_asd.csv               # 141 SCZ+CTL samples x 15 ASD pathways
├── gene_expression_processed.csv         # Full processed expression matrix
├── sample_metadata_with_subtypes.csv
└── results_summary.json
```

**Key result:** ARI = 0.870 cross-disease convergence, bootstrap ARI = 0.923, all 3 gates PASS

### Step 4: Notebook 12b -- Null ARI Permutation Test

Statistical significance test for the cross-disease convergence finding.

**Requires from Notebook 12:**
```
outputs/gse80655/pathway_scores_scz.csv       # or research-results/GSE80655/
outputs/gse80655/pathway_scores_asd.csv
outputs/gse80655/gene_expression_processed.csv
```

**Produces:**
```
outputs/null_ari_permutation/
├── null_ari_results.json                 # p-value, percentile, null distribution stats
├── null_ari_histogram.png                # Figure for manuscript
└── manuscript_paragraph.txt              # Draft text for Section 6.5
```

**Key result:** Observed ARI exceeds 100% of null distribution, p < 0.001

### Step 5: Notebook 13 -- GSE111175 Blood + ADOS Clinical Correlation

First blood-based analysis. Tests whether pathway subtypes are detectable in peripheral tissue and correlate with clinical symptom severity.

**Requires from Notebook 10 (optional, for cross-tissue comparison):**
```
outputs/gse28521/frontal_cortex/fc_pathway_scores.csv
outputs/gse28521/frontal_cortex/fc_sample_metadata_with_subtypes.csv
```

If Notebook 10 outputs are not found, the cross-tissue comparison section is skipped. The notebook also checks `research-results/GSE28521/frontal-cortex/` and a GitHub URL as fallbacks.

**Produces:**
```
outputs/gse111175/
├── results_summary.json
├── pathway_scores_matrix.csv             # 141 samples x 15 pathways (full)
├── pathway_scores_asd.csv               # ASD-only subset
├── sample_metadata_with_subtypes.csv
├── pathway_enrichment.csv
├── gene_contributions.csv
├── subtype_summary.csv
├── pca_scatter.png
├── pathway_heatmap.png
├── gene_heatmap.png
├── benchmark_comparison.png
├── ados_by_subtype_boxplot.png
├── ados_correlation_pathways.png
└── blood_vs_brain_comparison.png
```

**Key result:** Synaptic transmission correlates with ADOS Social Affect (rho = -0.52, FDR p = 0.032)

### Step 6: Notebook 14 -- GSE18123 Large Blood Cohort

Largest cohort (n=285). Cross-cohort projection tests whether blood subtypes from Notebook 13 replicate.

**Requires from Notebook 13:**
```
outputs/gse111175/pathway_scores_asd.csv
outputs/gse111175/sample_metadata_with_subtypes.csv
```

**Requires from Notebook 10 (optional, for cross-tissue comparison):**
```
outputs/gse28521/frontal_cortex/fc_pathway_scores.csv
outputs/gse28521/frontal_cortex/fc_sample_metadata_with_subtypes.csv
```

If prior outputs are not found, cross-cohort/cross-tissue sections are skipped. The notebook checks `research-results/`, `../../research-results/`, `outputs/`, and GitHub URLs as fallbacks.

**Produces:**
```
outputs/gse18123/
├── results_summary.json
├── pathway_scores_matrix.csv             # 285 samples x 15 pathways (full)
├── pathway_scores_asd.csv               # ASD-only subset (72 samples)
├── sample_metadata_with_subtypes.csv
├── pathway_enrichment.csv
├── gene_contributions.csv
├── subtype_summary.csv
├── cross_cohort_projection_ari.json
├── pca_scatter.png
├── pathway_heatmap.png
├── gene_heatmap.png
├── benchmark_comparison.png
├── cross_cohort_comparison.png
└── blood_vs_brain_comparison.png
```

**Key result:** Cross-cohort projection ARI = 0.374 (exceeds 0.3 replication threshold)

---

## How to Run

### Option A: Google Colab (Recommended for First Run)

Each notebook has a "Open in Colab" badge at the top. Click it to run in your browser with no local setup.

**Requirements:** Google account, Colab Pro recommended for Tier 2 notebooks (15+ min runtime, 8+ GB RAM for GEO downloads).

**Output handling on Colab:** Notebooks save outputs to `./outputs/<dataset>/`. To preserve outputs between sessions, mount Google Drive:
```python
from google.colab import drive
drive.mount('/content/drive')
# Then copy outputs to Drive after execution
```

### Option B: Local Execution

```bash
# Clone and set up
git clone https://github.com/topmist-admin/pathway-subtyping-framework
cd pathway-subtyping-framework
python3 -m venv pathwayenv && source pathwayenv/bin/activate
pip install -e ".[dev,viz]"
pip install GEOparse mygene jupyter

# Run a notebook (headless)
jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.timeout=1800 \
  examples/notebooks/10_geo_autism_bulk.ipynb

# Or open in Jupyter
jupyter notebook examples/notebooks/
```

**Working directory matters:** Notebooks assume the working directory is `examples/notebooks/`. When running with `nbconvert`, outputs go to `examples/notebooks/outputs/<dataset>/`.

### Option C: VS Code

Open any `.ipynb` file in VS Code. Select the `pathwayenv` kernel. Run cells interactively.

---

## Output Management

### During Execution

Each notebook writes outputs to a local `outputs/<dataset>/` directory relative to the notebook's working directory:

```
examples/notebooks/
├── outputs/
│   ├── gse28521/           # Notebook 10
│   │   └── frontal_cortex/
│   ├── gse64018/           # Notebook 11
│   │   └── frontal_cortex/
│   ├── gse80655/           # Notebook 12
│   ├── null_ari_permutation/ # Notebook 12b
│   ├── gse111175/          # Notebook 13
│   └── gse18123/           # Notebook 14
└── data/                   # Downloaded GEO SOFT files (cached)
```

### For Archival

After execution, copy outputs to `research-results/` at the repository root for permanent storage:

```bash
# From the repository root:
cp -r examples/notebooks/outputs/gse28521/frontal_cortex/ research-results/GSE28521/frontal-cortex/
cp -r examples/notebooks/outputs/gse64018/frontal_cortex/ research-results/GSE64018/
cp -r examples/notebooks/outputs/gse80655/ research-results/GSE80655/
cp -r examples/notebooks/outputs/null_ari_permutation/ research-results/GSE80655/null_ari_permutation/
cp -r examples/notebooks/outputs/gse111175/ research-results/GSE111175/
cp -r examples/notebooks/outputs/gse18123/ research-results/GSE18123/
```

Later notebooks check both `outputs/` and `research-results/` when loading prior results.

### What Is Gitignored

The `outputs/` and `data/` directories under `examples/notebooks/` are gitignored. The `research-results/` directory is also gitignored. Only the clean (unexecuted) notebook `.ipynb` files are tracked in git.

---

## Manuscript Mapping

Each Tier 2 notebook corresponds to a section in the manuscript:

| Notebook | Manuscript Section | Key Claim |
|----------|-------------------|-----------|
| 10 | Section 5 (ASD Application) | GABA-Collapsed subtype, 56/44 split |
| 11 | Section 5.9 (Cross-Cohort Validation) | Cross-platform replication |
| 12 | Section 6 (Cross-Disease Validation) | ARI = 0.870, all 3 gates pass |
| 12b | Section 6.5 (Supplementary) | Permutation p < 0.001 |
| 13 | Section 6.10 (Blood Validation) | ADOS-pathway correlation |
| 14 | Section 6.10 (Blood Validation) | Cross-cohort replication ARI = 0.374 |

---

## Troubleshooting

### GEO Download Failures

NCBI's FTP server is unreliable. Notebooks 13 and 14 set `GEOPARSE_USE_HTTP_FOR_FTP=yes` and include retry logic with exponential backoff. If downloads still fail:

1. Download the SOFT file manually from the GEO page (e.g., https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28521)
2. Place the `.soft.gz` file in `examples/notebooks/data/`
3. Re-run the notebook -- it will detect the cached file and skip the download

### Multi-Platform Datasets (Notebook 14)

GSE18123 uses two Affymetrix platforms (GPL570 + GPL6244). The notebook:
1. Extracts expression data per-platform
2. Maps probes to gene symbols using platform-specific annotation columns
3. Merges on common gene symbols (8,937 genes)

If you see warnings about gene symbol mapping, they are expected -- not all probes map cleanly.

### Memory Issues

| Notebook | Peak RAM |
|----------|----------|
| 10 | ~4 GB |
| 11 | ~2 GB |
| 12 | ~6 GB |
| 12b | ~4 GB |
| 13 | ~4 GB |
| 14 | ~6 GB |

Use Colab Pro (25 GB RAM) or a machine with 16+ GB for notebooks 12 and 14.

### Missing Cross-Reference Data

If a notebook cannot find outputs from a prior notebook, it will:
1. Check multiple fallback paths (`outputs/`, `research-results/`, `../../research-results/`)
2. Attempt to download from the GitHub repository URL
3. If all fail, skip the dependent analysis section with a message

This means every notebook can run independently -- cross-references are optional enhancements, not hard requirements.

---

*Last updated: February 2026*
