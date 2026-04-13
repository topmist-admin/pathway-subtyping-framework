# Notebook Guide

> Complete guide to the 21 analysis notebooks (00-19 + 12b), their execution order, and how outputs flow between them.

---

## Notebook Inventory

### Tier 1: Tutorials (Synthetic Data)

These notebooks use synthetic data bundled with the framework. They run independently in any order and require no external downloads.

| # | Notebook | Topic | Runtime |
|---|----------|-------|---------|
| 00 | [quick_demo](examples/notebooks/00_quick_demo.ipynb) | 60-second end-to-end demo | <1 min |
| 01 | [getting_started](examples/notebooks/01_getting_started.ipynb) | Installation, basic pipeline, validation gates | 5 min |
| 02 | [expression_scoring](examples/notebooks/02_expression_scoring.ipynb) | ssGSEA, GSVA, mean-Z scoring comparison | 5 min |
| 03 | [multi_omic_fusion](examples/notebooks/03_multi_omic_fusion.ipynb) | Multi-modal data fusion (VCF + expression + scRNA-seq) | 5 min |
| 04 | [deconvolution](examples/notebooks/04_deconvolution.ipynb) | Cell-type deconvolution from bulk RNA-seq | 5 min |
| 05 | [visualization](examples/notebooks/05_visualization.ipynb) | PCA, t-SNE, UMAP, interactive Plotly reports | 5 min |
| 06 | [ancestry_batch_correction](examples/notebooks/06_ancestry_batch_correction.ipynb) | Population stratification and batch effects | 5 min |
| 07 | [sensitivity_analysis](examples/notebooks/07_sensitivity_analysis.ipynb) | Parameter robustness testing | 5 min |
| 08 | [characterization](examples/notebooks/08_characterization.ipynb) | Subtype profiling, heatmaps, gene contributions | 5 min |
| 09 | [signaling_databases](examples/notebooks/09_signaling_databases.ipynb) | CellPhoneDB and CellChatDB integration | 5 min |

### Tier 2: Real Data Validation (GEO Datasets)

These notebooks download real transcriptomics data from NCBI GEO and reproduce the manuscript analyses. They must be run in a specific order because later notebooks consume outputs from earlier ones.

| # | Notebook | Dataset | Tissue | N | Platform | Depends On |
|---|----------|---------|--------|---|----------|------------|
| 10 | [geo_autism_bulk](examples/notebooks/10_geo_autism_bulk.ipynb) | GSE28521 | Brain (frontal + temporal cortex) | 79 | Affymetrix microarray | None |
| 11 | [geo_autism_validation](examples/notebooks/11_geo_autism_validation.ipynb) | GSE64018 | Brain (temporal cortex) | 24 | Illumina RNA-seq | **10** |
| 12 | [geo_schizophrenia_bulk](examples/notebooks/12_geo_schizophrenia_bulk.ipynb) | GSE80655 | Brain (3 regions) | 281 | Illumina RNA-seq | **10** (optional) |
| 12b | [null_ari_permutation](examples/notebooks/12b_null_ari_permutation.ipynb) | GSE80655 | Brain (3 regions) | 141 | (reuses 12 outputs) | **12** |
| 13 | [geo_blood_ados](examples/notebooks/13_geo_blood_ados.ipynb) | GSE111175 | Blood (leukocytes) | 141 | Illumina BeadChip | **10** (optional) |
| 14 | [geo_blood_large_cohort](examples/notebooks/14_geo_blood_large_cohort.ipynb) | GSE18123 | Blood | 285 | Affymetrix (2 platforms) | **13**, **10** (optional) |
| 15 | [geo_scz_replication](examples/notebooks/15_geo_scz_replication.ipynb) | GSE53987 | Brain (3 regions) | 205 | Affymetrix microarray | **12** |
| 16 | [knowledge_graph_analysis](examples/notebooks/16_knowledge_graph_analysis.ipynb) | Multi-dataset | KG (STRING + DGIdb) | 1,075 | Network analysis | **10-15** |
| 17 | [tcga_cancer_validation](examples/notebooks/17_tcga_cancer_validation.ipynb) | TCGA-COAD | Colon adenocarcinoma | 452 | RNA-seq (NCI GDC API) | None (standalone) |
| 18 | [geo_clinical_phenotype](examples/notebooks/18_geo_clinical_phenotype.ipynb) | GSE15402 | LCL (lymphoblastoid) | 116 | TIGR 40K two-channel | None (standalone) |
| 19 | [scz_blood_multi_cohort](examples/notebooks/19_scz_blood_multi_cohort.ipynb) | Multi-GEO (5) | Blood (SCZ) | 407 | Multi-platform | None (standalone) |

---

## Dependency Diagram

```
Tier 1 (Tutorials 00-09)          Tier 2 (Real Data 10-17)
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
                                    │    │           ▼
                                    │    ▼      14: GSE18123
                                    │  12b:    Blood Large Cohort
                                    │  Null ARI (needs 13 + 10)
                                    │  Permutation
                                    ▼
                               15: GSE53987
                               SCZ Replication
                               (Affymetrix, needs 12)


                    16: Knowledge Graph Analysis
                    (uses results from 10-15;
                     STRING PPI + DGIdb drug targets)

                    17: TCGA-COAD Cancer Validation
                    (standalone; NCI GDC API download;
                     MSigDB Hallmark pathways;
                     CMS validation + survival analysis)

                    18: GSE15402 Clinical Phenotype
                    (standalone; ADI-R severity subgroups;
                     TIGR 40K EST microarray)

                    19: SCZ Blood Multi-Cohort
                    (standalone; 5 GEO datasets;
                     Hertzberg replication)
```

**Arrow meaning:** The target notebook loads output files produced by the source notebook. If the source has not been run, the target will skip that analysis section (with a warning) and proceed with available data.

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

If Notebook 10 outputs are not found, the cross-tissue comparison section is skipped. The notebook also checks `research-results/GSE28521/frontal-cortex/` as a fallback.

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

If prior outputs are not found, cross-cohort/cross-tissue sections are skipped. The notebook checks `research-results/`, `../../research-results/`, and `outputs/` as fallbacks.

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

### Step 7: Notebook 15 -- GSE53987 SCZ Replication (Affymetrix)

Independent replication of SCZ pathway subtypes on a different platform (Affymetrix vs RNA-seq) and different brain regions. Tests cross-cohort projection from GSE80655 and cross-disease ARI with ASD pathways.

**Requires from Notebook 12:**
```
outputs/gse80655/pathway_scores_scz.csv
outputs/gse80655/sample_metadata_with_subtypes.csv
```

If Notebook 12 outputs are not found, checks `research-results/GSE80655/` as fallback.

**Produces:**
```
outputs/gse53987/
├── results_summary.json
├── pathway_scores_scz_ctl.csv          # 103 SCZ+CTL samples x 14 SCZ pathways
├── pathway_scores_asd.csv              # 103 samples x 15 ASD pathways
├── sample_metadata_with_subtypes.csv
├── pathway_enrichment.csv
├── gene_contributions.csv
├── subtype_summary.csv
├── benchmark_results.csv
├── subtype_heatmap.png
├── gene_heatmap.png
├── pca_scatter_subtypes.png
├── benchmark_comparison.png
├── algorithm_comparison.png
├── model_selection.png
├── regions/                            # Per-region subtyping
│   ├── {region}_results.json
│   ├── {region}_pathway_heatmap.png
│   ├── cross_region_ari.json
│   └── cross_region_ari.png
├── cross_cohort/                       # GSE80655 → GSE53987 projection
│   ├── cross_cohort_projection_results.json
│   ├── centroid_correlation.csv
│   ├── centroid_correlation.png
│   └── cross_platform_comparison.csv
└── cross_disease/                      # ASD pathway comparison
    ├── multi_diagnosis_heatmap.png
    ├── multi_diagnosis_results.json
    └── shared_pathway_correlation.csv
```

**Key results:**
- Cross-cohort projection ARI = 0.319 (PASS, threshold 0.3)
- Cross-disease ARI = 0.792 (strong SCZ/ASD pathway convergence)
- Multi-diagnosis pooled silhouette = 0.450 (k=5, diagnosis-independent)

---

### Step 8: Notebook 16 -- Knowledge Graph Analysis

Cross-disease knowledge graph integrating gene contributions from all 6 datasets with STRING PPI network and DGIdb drug targets. Produces hub gene ranking, drug repurposing table, and manuscript-ready network figures.

**Requires from Notebooks 10-15:**
```
research-results/GSE28521/gene_contributions.csv
research-results/GSE64018/frontal-cortex/gene_contributions.csv
research-results/GSE80655/gene_contributions.csv
research-results/GSE53987/gene_contributions.csv
research-results/GSE111175/gene_contributions.csv
research-results/GSE18123/gene_contributions.csv
data/pathways/autism_pathways.gmt
data/pathways/schizophrenia_pathways.gmt
```

**Requires network access:** STRING API (string-db.org), DGIdb API (dgidb.org). Results are cached locally after first run.

**Produces:**
```
outputs/knowledge_graph/
├── kg_results_summary.json             # Complete metrics summary
├── hub_genes_ranked.csv                # Top 20 hub genes by betweenness
├── drug_repurposing_table.csv          # Ranked drug candidates
├── community_analysis.csv             # Louvain communities with disease enrichment
├── pathway_crosstalk_matrix.csv       # 21×21 PPI-based pathway crosstalk
├── string_ppi_edges.csv               # Cached STRING PPI data
├── dgidb_interactions.csv             # Cached DGIdb drug interactions
├── cross_disease_network.graphml      # Full KG (importable to Cytoscape)
├── convergence_subnetwork.graphml     # ASD↔SCZ bridge genes subnetwork
├── network_figure.png / .pdf          # Manuscript-ready network (300 DPI)
├── drug_network_figure.png            # Drug targets overlay
├── pathway_crosstalk_heatmap.png      # Pathway crosstalk heatmap
└── ppi_distributions.png              # PPI score and degree distributions
```

**Key results:**
- 336 genes, 4,378 PPI edges, 97.9% coverage
- 20 hub genes (11 cross-disease bridges: AKT1, CTNNB1, GSK3B, PTEN, etc.)
- 1,546 unique drug candidates across 44 target genes
- 6 Louvain communities (all cross-disease)

---

### Step 9: Notebook 17 -- TCGA-COAD Cancer Validation

**Standalone notebook** — demonstrates disease-agnostic applicability by applying the framework to colorectal adenocarcinoma from TCGA. Downloads RNA-seq data from NCI GDC REST API, scores 50 MSigDB Hallmark cancer pathways, discovers molecular subtypes, and validates against CMS1-4 using a pure Python NTP classifier.

**No prior notebooks required.** Can be run independently.

**Requires network access:** NCI GDC REST API (expression + clinical data), MSigDB Broad Institute server (Hallmark GMT). Results are cached locally after first run.

**Produces:**
```
research-results/tcga/
├── results_summary.json               # Complete metrics (k, silhouette, CMS ARI, etc.)
├── sample_metadata_with_subtypes.csv  # Subtype assignments + all clinical variables
├── pathway_scores.csv                 # 50 Hallmark pathway scores per sample
├── gene_expression_processed.csv      # Cleaned RNA-seq matrix (log2 TPM+1)
├── pathway_enrichment.csv             # Enriched pathways per subtype (FDR < 0.05)
├── gene_contributions.csv             # Top genes per subtype (Cohen's d)
├── subtype_summary.csv                # Subtype profile summary
├── cms_ntp_results.csv                # CMS predictions (NTP, per-sample FDR)
├── hallmark_score_distributions.png   # ssGSEA score distributions (top 12 pathways)
├── model_selection.png                # BIC + silhouette curves for k selection
├── pca_scatter.png                    # 3-panel PCA: our subtypes / CMS / MSI
├── pathway_heatmap.png                # Subtype × pathway heatmap (top 25 pathways)
├── subtype_pathway_heatmap.png        # Characterization pathway heatmap
├── subtype_gene_heatmap.png           # Top genes per subtype
├── benchmark_comparison.png           # Pathway-GMM vs NMF/PCA-Kmeans/Gene-Kmeans
├── cms_msi_comparison.png             # Comparison to CMS1-4 and MSI status
└── cms_validation_heatmap.png         # CMS NTP validation contingency heatmap
```

**Key results:**
- 452 TCGA-COAD primary tumors (6 TCGA-A6 batch outliers removed)
- 50 MSigDB Hallmark gene sets, k=3 (forced after degenerate BIC k=2)
- Subtype 0 (19%): Stromal/EMT → CMS4 at 76% (Fisher OR=16.7, p=1.4e-25)
- Subtype 1 (29%): Immune-cold → CMS2 at 60%
- Subtype 2 (52%): Proliferative/Metabolic → spread across CMS1-3
- CMS classification rate: 436/452 (96.5%, FDR ≤ 0.05)
- Benchmark: pathway_gmm wins; 2/3 validation gates pass
- Survival analysis: Kaplan-Meier curves, log-rank test, Cox PH regression

---

### Step 10: Notebook 18 -- GSE15402 Clinical Phenotype Validation

**Standalone notebook** — validates whether pathway-derived subtypes correlate with ADI-R clinical severity subgroups in lymphoblastoid cell lines.

**No prior notebooks required.** Can be run independently.

**Requires network access:** NCBI GEO (GSE15402 download), MyGene.info (probe-to-gene mapping).

**Produces:**
```
outputs/gse15402/
├── results_summary.json
├── pathway_scores_matrix.csv
├── sample_metadata_with_subtypes.csv
├── subtypes_vs_adi_r.png
└── ... (additional figures and CSVs)
```

**Key results:**
- 116 LCL samples (87 ASD, 29 controls), TIGR 40K platform (532 genes mapped)
- k=6 subtypes, silhouette=0.057, 2/3 validation gates pass
- ADI-R subgroup association: chi-square p=0.001
- S5→L subgroup OR=13.2 (raw p=0.008)

---

### Step 11: Notebook 19 -- SCZ Blood Multi-Cohort (Hertzberg Replication)

**Standalone notebook** — applies PSF to 5 schizophrenia blood GEO datasets for multi-cohort validation and Hertzberg concordance mapping.

**No prior notebooks required.** Can be run independently.

**Requires network access:** NCBI GEO (5 datasets: GSE38484, GSE27383, GSE38481, GSE18312, GSE48072).

**Produces:**
```
outputs/scz_blood_multi_cohort/
├── results_summary.json
├── merged_pathway_scores.csv
├── merged_metadata_with_subtypes.csv
├── cross_cohort_projection_results.json
├── hertzberg_concordance.json
├── per-dataset figures and CSVs
└── ... (20 output files total)
```

**Key results:**
- 407 total samples (177 SCZ) across 5 GEO datasets, 5 microarray platforms
- Per-dataset subtyping: k=2-4, silhouette 0.10-0.25
- Merged analysis (177 SCZ): k=7, silhouette=0.088, 2/3 gates
- Cross-cohort projection: mean ARI=0.205 (GSE38484→GSE18312 ARI=0.469 PASS)
- Hertzberg concordance: 4 Immune-like (A) + 3 Neuro-like (B) subtypes

---

## How to Run

### Option A: Binder (Recommended for First Run)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F00_quick_demo.ipynb)

Click the badge above (or any Binder link in the notebook tables) to launch a notebook in your browser with no local setup. Binder provides a free JupyterLab environment with all dependencies pre-installed.

**Requirements:** None — just a web browser. Binder sessions are temporary; download any outputs before the session expires.

**Binder links for individual notebooks:**

| # | Binder Link |
|---|-------------|
| 00 | [quick_demo](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F00_quick_demo.ipynb) |
| 01 | [getting_started](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F01_getting_started.ipynb) |
| 02 | [expression_scoring](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F02_expression_scoring.ipynb) |
| 03 | [multi_omic_fusion](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F03_multi_omic_fusion.ipynb) |
| 04 | [deconvolution](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F04_deconvolution.ipynb) |
| 05 | [visualization](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F05_visualization.ipynb) |
| 06 | [ancestry_batch_correction](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F06_ancestry_batch_correction.ipynb) |
| 07 | [sensitivity_analysis](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F07_sensitivity_analysis.ipynb) |
| 08 | [characterization](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F08_characterization.ipynb) |
| 09 | [signaling_databases](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F09_signaling_databases.ipynb) |

> **Note:** Tier 2 notebooks (10-19) require GEO/GDC data downloads and 4-6 GB RAM, which may exceed Binder's free resource limits. Run those locally (Option B) or via Docker (Option D).

### Option B: Local Execution

```bash
# Clone and set up
git clone https://codeberg.org/pathways/pathway-subtyping-framework
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

### Option D: Docker

```bash
# Pull the Jupyter image with all dependencies
docker pull rohitdataops/pathway-subtyping:0.5.0-jupyter

# Launch Jupyter server
docker run -p 8888:8888 \
  -v $(pwd)/examples/notebooks:/home/psf/examples/notebooks \
  rohitdataops/pathway-subtyping:0.5.0-jupyter

# Or use docker compose
docker compose up jupyter
# Open http://localhost:8888
```

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
│   ├── gse18123/           # Notebook 14
│   └── gse53987/           # Notebook 15
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
cp -r examples/notebooks/outputs/gse53987/ research-results/GSE53987/
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
| 15 | Section 6.11 (SCZ Replication) | Cross-platform ARI = 0.319, cross-disease ARI = 0.792 |
| 16 | Section 7 (Knowledge Graph) | 20 hub genes, 1,546 drug candidates, 6 cross-disease communities |
| 17 | Section 7 (Cancer Validation) | CMS4 recovery at 76% (Fisher p=1.4e-25), survival analysis |
| 18 | Section 6.10 (Clinical Validation) | ADI-R subgroup association chi-sq p=0.001 |
| 19 | Section 6.11 (SCZ Replication) | Multi-cohort SCZ blood, Hertzberg concordance |

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
| 15 | ~4 GB |

Use a machine with 16+ GB RAM for notebooks 12 and 14.

### Missing Cross-Reference Data

If a notebook cannot find outputs from a prior notebook, it will:
1. Check multiple fallback paths (`outputs/`, `research-results/`, `../../research-results/`)
2. If all local paths fail, skip the dependent analysis section with a message

This means every notebook can run independently -- cross-references are optional enhancements, not hard requirements.

---

*Last updated: March 2026*
