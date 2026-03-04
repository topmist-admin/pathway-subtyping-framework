# Pathway Subtyping Framework

**A Disease-Agnostic Tool for Pathway-Based Molecular Subtype Discovery**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18442426.svg)](https://doi.org/10.5281/zenodo.18442426)
[![PyPI version](https://badge.fury.io/py/pathway-subtyping.svg)](https://pypi.org/project/pathway-subtyping/)
[![CI](https://codeberg.org/pathways/pathway-subtyping-framework/badges/workflows/ci.yml/badge.svg)](https://codeberg.org/pathways/pathway-subtyping-framework)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![RRID:SCR_028051](https://img.shields.io/badge/RRID-SCR__028051-blue.svg)](https://scicrunch.org/resolver/RRID:SCR_028051)

---

## Overview

The **Pathway Subtyping Framework** is an open-source computational tool for identifying molecular subtypes in genetically heterogeneous diseases. Instead of analyzing individual genes, it aggregates rare variant burden at the **biological pathway level**, enabling:

- Better signal detection across genetically diverse cohorts
- Identification of biologically coherent patient subgroups
- Cross-cohort validation of discovered subtypes

Originally developed for [autism research](https://codeberg.org/pathways/autism-pathway-framework), this generalized version can be adapted for any disease with:
- Genetic heterogeneity (many implicated genes)
- Convergent pathway biology
- Available exome/genome sequencing data

## Supported Disease Areas

| Disease | Status | Pathway File |
|---------|--------|--------------|
| Autism Spectrum Disorder | Validated | `autism_pathways.gmt` |
| Schizophrenia | Template | `schizophrenia_pathways.gmt` |
| Epilepsy | Template | `epilepsy_pathways.gmt` |
| Intellectual Disability | Template | `intellectual_disability_pathways.gmt` |
| Parkinson's Disease | Template | `parkinsons_pathways.gmt` |
| Bipolar Disorder | Template | `bipolar_pathways.gmt` |
| *Your disease* | [Adapt it →](docs/guides/adapting-for-your-disease.md) | `your_pathways.gmt` |

## Key Features

| Feature | Description |
|---------|-------------|
| **Pathway Scoring** | Aggregate gene burdens across biological pathways |
| **Expression Scoring** | Bulk RNA-seq pathway scoring via ssGSEA, GSVA, or mean-Z methods |
| **Single-Cell Scoring** | Per-cell and pseudobulk pathway scoring from scRNA-seq (h5ad/CSV) |
| **Multi-Omic Fusion** | Fuse VCF + expression + single-cell scores (concatenate, weighted, intersection) |
| **Bulk Deconvolution** | Estimate cell-type proportions from bulk RNA-seq via NNLS; cell-type-aware subtypes |
| **Multiple Clustering** | GMM, K-means, Hierarchical, Spectral with cross-validation |
| **Ancestry Correction** | PCA-based population stratification correction with independence testing |
| **Batch Correction** | ComBat-style batch effect detection and correction |
| **Sensitivity Analysis** | Parameter robustness testing across algorithms, features, normalization |
| **Threshold Calibration** | Data-driven validation thresholds that adjust for sample size and cluster count |
| **Variant QC** | QUAL, call rate, HWE, MAF filters before burden computation |
| **Validation Gates** | 5 gates: negative controls, bootstrap stability, ancestry independence, cross-modal concordance |
| **Statistical Rigor** | FDR correction, effect sizes, confidence intervals |
| **Power Analysis** | Sample size recommendations, Type I error estimation |
| **Simulation** | Synthetic data generation with ground truth for validation |
| **Cross-Cohort Validation** | Transfer learning and projection-based replication testing |
| **Visualization** | Interactive Plotly HTML reports, UMAP/t-SNE scatter plots, radar charts, multi-format export |
| **Performance** | tqdm progress bars, chunked VCF processing, 10K+ sample support |
| **Reproducibility** | Deterministic execution, pinned dependencies, Docker |
| **Config-Driven** | YAML configuration for all parameters |

## Quick Start

### Installation

```bash
pip install pathway-subtyping
```

Optional extras:

```bash
pip install pathway-subtyping[vcf]   # VCF file processing (pysam)
pip install pathway-subtyping[viz]   # Interactive visualizations (Plotly, UMAP)
pip install pathway-subtyping[sc]    # Single-cell support (AnnData)
pip install pathway-subtyping[graph] # Network analysis (NetworkX, py4cytoscape)
```

#### Network Visualization (Cytoscape)

The `[graph]` extra enables publication-ready network figure generation via
[Cytoscape](https://cytoscape.org), an open-source desktop application for
network visualization. The `py4cytoscape` library communicates with Cytoscape's
CyREST API on `localhost:1234`.

**Setup:**

1. Download and install [Cytoscape desktop](https://cytoscape.org/download.html) (v3.10+)
2. Launch Cytoscape and wait for it to fully load
3. Install the Python extra: `pip install pathway-subtyping[graph]`
4. Verify the connection:
   ```bash
   python -c "import py4cytoscape as p4c; print(p4c.cytoscape_ping())"
   ```
5. Run the figure generator:
   ```bash
   python scripts/generate_cytoscape_figures.py
   ```

### Try in Browser (No Installation)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F00_quick_demo.ipynb)

**60-second demo** — generates a synthetic cohort, discovers subtypes, validates them, and visualizes results. No data needed. Click the badge above to launch in Binder.

**Full tutorial**: [01_getting_started.ipynb](https://mybinder.org/v2/git/https%3A%2F%2Fcodeberg.org%2Fpathways%2Fpathway-subtyping-framework.git/main?labpath=examples%2Fnotebooks%2F01_getting_started.ipynb) (Binder) | [Local](examples/notebooks/01_getting_started.ipynb)

### Notebooks

21 Jupyter notebooks covering tutorials through full manuscript reproduction. See [docs/notebook-guide.md](docs/notebook-guide.md) for execution order and dependencies.

**Tutorials (00-09)** -- Synthetic data, standalone, any order:

| # | Topic |
|---|-------|
| 00 | Quick demo (60 seconds) |
| 01 | Getting started (installation, pipeline, validation) |
| 02 | Expression scoring (ssGSEA, GSVA, mean-Z) |
| 03 | Multi-omic fusion |
| 04 | Cell-type deconvolution |
| 05 | Visualization (PCA, t-SNE, UMAP, Plotly) |
| 06 | Ancestry and batch correction |
| 07 | Sensitivity analysis |
| 08 | Subtype characterization |
| 09 | Signaling database integration |

**Real Data Validation (10-19)** -- GEO/TCGA datasets, run in order (later notebooks use earlier outputs):

| # | Dataset | Tissue | N | Manuscript Section |
|---|---------|--------|---|-------------------|
| 10 | GSE28521 | Brain (frontal + temporal cortex) | 79 | Section 5 |
| 11 | GSE64018 | Brain (temporal cortex, RNA-seq) | 24 | Section 5.9 |
| 12 | GSE80655 | Brain (3 regions, multi-diagnosis) | 281 | Section 6 |
| 12b | GSE80655 | Null ARI permutation test | 141 | Section 6.5 |
| 13 | GSE111175 | Blood + ADOS clinical scores | 141 | Section 6.10 |
| 14 | GSE18123 | Blood (largest cohort, 2 platforms) | 285 | Section 6.10 |
| 15 | GSE53987 | Brain (3 regions, Affymetrix, 4 diagnoses) | 205 | Section 6.11 |
| 16 | Multi-dataset | Knowledge graph (STRING PPI + DGIdb) | 1,075 | Section 7 |
| 17 | TCGA-COAD | Colon adenocarcinoma (452 tumors, CMS + survival) | 452 | Section 7 |
| 18 | GSE15402 | LCL (ADI-R clinical phenotype) | 116 | Section 6.10 |
| 19 | Multi-GEO | SCZ blood multi-cohort (Hertzberg replication) | 407 | Section 6.11 |

### Run with Sample Data

```bash
# Clone for sample data and configs
git clone https://codeberg.org/pathways/pathway-subtyping-framework.git
cd pathway-subtyping-framework

# Run the pipeline with synthetic test data
psf --config configs/test_synthetic.yaml

# View results
cat outputs/synthetic_test/report.md
```

### Run with Your Data

```bash
# Copy and customize a config
cp configs/example_autism.yaml configs/my_analysis.yaml

# Edit paths in my_analysis.yaml, then run
psf --config configs/my_analysis.yaml
```

### Docker

```bash
# Pull from Docker Hub
docker pull rohitdataops/pathway-subtyping:latest        # CLI runtime
docker pull rohitdataops/pathway-subtyping:0.4.0-jupyter  # Jupyter + notebooks

# Run pipeline
docker compose run pipeline

# Run tests (1003 tests)
docker compose run test

# Start Jupyter notebook
docker compose up jupyter
# Open http://localhost:8888
```

## Adapting for Your Disease

1. **Create a pathway GMT file** with disease-relevant gene sets
2. **Copy an example config** and point to your data
3. **Run the pipeline** — validation gates will tell you if subtypes are meaningful

See the full guide: [Adapting for Your Disease](docs/guides/adapting-for-your-disease.md)

## How It Works

```
VCF / Expression / scRNA-seq → Pathway Scoring → [Ancestry Correction] → [Batch Correction] → GMM Clustering → [Sensitivity Analysis] → Validation → Report
```

### 1. Pathway Scoring
Multiple input modalities are supported, each producing the same Z-normalized pathway score matrix:
- **VCF**: Rare damaging variants aggregated with LoF/CADD weights
- **Expression**: Bulk RNA-seq scored via ssGSEA, GSVA, or mean-Z
- **Single-cell**: Pseudobulk or per-cell scoring from scRNA-seq
- **Multi-omic**: Fuse scores from multiple modalities for unified subtype discovery

### 2. Subtype Discovery
Multiple clustering algorithms identify patient subgroups:
- **GMM** (default): Soft assignments, automatic selection via BIC
- **K-means**: Fast, spherical clusters
- **Hierarchical**: Dendogram-based, no K required
- **Spectral**: Nonlinear boundaries
- Cross-validation for stability assessment
- Algorithm comparison with pairwise ARI

### 3. Validation Gates
Built-in tests prevent overfitting:
- **Label shuffle**: Randomized labels should NOT cluster (ARI < 0.15)
- **Random genes**: Fake pathways should NOT work (ARI < 0.15)
- **Bootstrap**: Clusters should be stable under resampling (ARI > 0.8)
- **Ancestry independence**: Clusters should not correlate with ancestry PCs (when provided)
- **Cross-modal concordance**: Subtypes should replicate across data modalities (when multi-omic)

### 4. Statistical Rigor
Publication-quality statistics:
- **FDR correction**: Benjamini-Hochberg for multiple testing
- **Effect sizes**: Cohen's d with 95% bootstrap confidence intervals
- **Power analysis**: Sample size recommendations for target effect sizes
- **Type I error**: Estimation via null simulations

See [docs/METHODS.md](docs/METHODS.md) for full statistical methodology.

## Data Requirements

| Input | Format | Notes |
|-------|--------|-------|
| Variants | VCF | Annotated with gene symbols, consequences |
| Bulk Expression | CSV/TSV | Gene expression matrix (samples x genes) |
| Single-Cell | h5ad/CSV | AnnData or cell-by-gene matrix with cell type annotations |
| Phenotypes | CSV | Sample IDs + clinical features |
| Pathways | GMT | Gene sets for your disease |

**Your data stays on your infrastructure.** The framework runs locally or in your cloud environment.

## Data Provenance and Integrity

**This project contains zero proprietary, commercial, or third-party customer data.**

Every data file in this repository was either:

1. **Computationally generated** — The synthetic VCF and phenotype files in `data/sample/` were created by our `SyntheticDataGenerator` using random number generators with fixed seeds. They contain no real patient or clinical data whatsoever.
2. **Curated from public scientific literature** — The pathway GMT files in `data/pathways/` contain gene symbol lists assembled exclusively from publicly available, peer-reviewed sources: [SFARI Gene](https://gene.sfari.org/), [KEGG](https://www.kegg.jp/), [Reactome](https://reactome.org/), [MSigDB](https://www.gsea-msigdb.org/), and [Gene Ontology](http://geneontology.org/). Gene symbols (e.g., SHANK3, CHD8) are standard scientific identifiers published in thousands of research papers.
3. **Open-source code only** — All algorithms are original implementations or standard open-source libraries (scikit-learn, scipy, numpy, pandas). No proprietary software, commercial code, or licensed algorithms were used.

**No data from any employer, client, institution, or commercial entity was used at any stage of this project** — not in development, testing, validation, or documentation. The framework is designed so that users supply their own data; it does not ship with, embed, or depend on any private or restricted datasets.

For full details, see [DISCLAIMER.md](DISCLAIMER.md) and [docs/contributor-kit/04-research-compliance.md](docs/contributor-kit/04-research-compliance.md).

## Project Structure

```
pathway-subtyping-framework/
├── src/pathway_subtyping/       # Core Python package
│   ├── pipeline.py              # Main pipeline orchestrator
│   ├── clustering.py            # GMM, K-means, Hierarchical, Spectral
│   ├── validation.py            # Validation gates (5 gates)
│   ├── statistical_rigor.py     # FDR, effect sizes, burden weights
│   ├── simulation.py            # Synthetic data & power analysis
│   ├── expression.py            # Bulk RNA-seq pathway scoring (ssGSEA, GSVA, mean-Z)
│   ├── single_cell.py           # Single-cell scRNA-seq scoring (pseudobulk + per-cell)
│   ├── multi_omic.py            # Multi-omic pathway score fusion
│   ├── deconvolution.py         # Bulk deconvolution (NNLS cell-type proportions)
│   ├── cross_modal_validation.py # Cross-modal concordance gate (Gate 5)
│   ├── visualization.py         # Interactive Plotly reports, UMAP/t-SNE, export
│   ├── characterization.py      # Subtype profiling, heatmaps, gene contributions
│   ├── ancestry.py              # Population stratification correction
│   ├── batch_correction.py      # Batch effect detection & correction
│   ├── sensitivity.py           # Parameter sensitivity analysis
│   ├── benchmark.py             # Method comparison benchmarks
│   ├── cross_cohort.py          # Cross-cohort validation
│   ├── threshold_calibration.py # Data-driven threshold calibration
│   ├── variant_qc.py            # Variant QC (QUAL, HWE, MAF, call rate)
│   ├── validation_datasets.py   # ClinVar/Reactome integration
│   ├── data_quality.py          # VCF quality checks
│   └── utils/                   # Performance, seeding, progress tracking
├── configs/                     # Example YAML configurations
├── data/
│   ├── pathways/                # Pathway GMT files (6 diseases)
│   └── sample/                  # Synthetic test data
├── docs/
│   ├── METHODS.md               # Statistical methods documentation
│   ├── api/                     # API reference (13 modules)
│   └── guides/                  # User guides
├── scripts/                     # Utility scripts
│   ├── generate_cytoscape_figures.py  # Publication-ready network figures (requires [graph] + Cytoscape desktop)
│   ├── validate_with_public_data.py   # ClinVar/Reactome validation
│   └── benchmark_performance.py       # Performance benchmarks
├── examples/notebooks/          # Jupyter tutorials
├── tests/                       # Test suite (1003 tests)
├── Dockerfile                   # Container support
└── docker-compose.yml           # Easy orchestration
```

## Development

```bash
# Install with dev dependencies (from cloned repo)
pip install -e ".[dev,vcf,viz,sc]"

# Run tests
pytest tests/ -v

# Run linting
black src/ tests/
isort src/ tests/
flake8 src/ tests/

# Set up pre-commit hooks
pre-commit install
```

## Related Projects

- **[Autism Pathway Framework](https://codeberg.org/pathways/autism-pathway-framework)** — The original autism-focused implementation with SFARI cohort validation

## Contributing

Contributions welcome! Areas where help is needed:
- Additional disease pathway definitions
- Multi-omic integration (spatial transcriptomics, proteomics)
- Documentation and tutorials

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Citation

If you use this framework in your research, please cite the preprint:

```
Chauhan R. Pathway-Based Molecular Subtyping Reveals a GABA-Collapsed Autism Subgroup
and Cross-Disease Convergence with Schizophrenia in Human Cerebral Cortex.
Research Square. 2026. DOI: 10.21203/rs.3.rs-8913089/v1
https://www.researchsquare.com/article/rs-8913089/v1
```

To cite the software specifically:
```
Chauhan R. Pathway Subtyping Framework. Zenodo. 2026.
DOI: 10.5281/zenodo.18442426
https://codeberg.org/pathways/pathway-subtyping-framework
```

For autism-specific work, also cite:
```
Chauhan R. Autism Pathway Framework. Zenodo. 2026.
DOI: 10.5281/zenodo.18403844
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

**Rohit Chauhan**
- Email: info@topmist.com
- Codeberg: [@pathways](https://codeberg.org/pathways)
- ORCID: [0009-0003-9895-4629](https://orcid.org/0009-0003-9895-4629)

---

> **RESEARCH USE ONLY** — This framework is for hypothesis generation. Not for clinical diagnosis or treatment decisions.
