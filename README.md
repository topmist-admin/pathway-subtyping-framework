# Pathway Subtyping Framework

> ## ⚠️ CORRECTION NOTICE (2026-07-08)
> The **adaptive bootstrap-threshold model** previously reported at R²=0.889 does **not reproduce** and is **retracted**; the **47-dataset benchmark** contained an empty-input ARI artifact (14 invalid rows), an incorrect independence claim (TCGA-COAD was present despite being described as excluded), and count discrepancies. See **[`CORRECTION_2026-07/ERRATUM_2026-07-08.md`](CORRECTION_2026-07/ERRATUM_2026-07-08.md)** for the full notice and corrected artifacts. A **corrected Zenodo dataset v2.0 (`10.5281/zenodo.21262112`)** supersedes the prior versions. The associated preprint (Research Square rs-9284565) has been **withdrawn** from journal review; rs-8913089 is being revised. Do not use the retracted calibration model or the uncorrected benchmark.

**A Disease-Agnostic Tool for Pathway-Based Molecular Subtype Discovery**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18638048.svg)](https://doi.org/10.5281/zenodo.18638048)
[![PyPI version](https://badge.fury.io/py/pathway-subtyping.svg)](https://pypi.org/project/pathway-subtyping/)
[![CI](https://codeberg.org/pathways/pathway-subtyping-framework/badges/workflows/ci.yml/badge.svg)](https://codeberg.org/pathways/pathway-subtyping-framework)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
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
| **Validation Gates** | Discreteness (Gate A, SigClust), negative controls, bootstrap stability (demoted v0.8), ancestry independence, cross-modal concordance, confound association (Gate 6), genetic/somatic anchoring (Gate 7), KG-version sensitivity (Gate K) |
| **Statistical Rigor** | FDR correction, effect sizes, confidence intervals |
| **Power Analysis** | Sample size recommendations, Type I error estimation |
| **Simulation** | Synthetic data generation with ground truth for validation |
| **Cross-Cohort Validation** | Transfer learning and projection-based replication testing |
| **Visualization** | Interactive Plotly HTML reports, UMAP/t-SNE scatter plots, radar charts, multi-format export |
| **Molecular QC** | 12-feature manufacturing QC: cascade detection, dosage gating, off-target activation, drift detection, stress fingerprinting |
| **Knowledge Graph** | Topology-aware pathway scoring, hierarchical queries, cross-omics entity resolution, pathway crosstalk quantification |
| **Uncertainty Quantification** *(v0.6)* | Split-conformal prediction intervals, bootstrap, Bayesian pathway-GMM drop-in, calibration diagnostics |
| **Cross-Platform Harmonization** *(v0.6)* | UCE-anchored aligner for 10x / Smart-seq2 / bulk / spatial, with per-cell confidence scoring |
| **KG Refresh Infrastructure** *(v0.6)* | Pinned source manifest (OmniPath / SIGNOR / Reactome), diff + regression tools, reproducibility hashing |
| **AlphaMissense-modulated Cascade** *(v0.6)* | Variant-carrier down-weighting for pathway-cascade scoring |
| **In-silico Perturbation** *(v0.6)* | Geneformer wrapper + MSV-from-embedding head + batch screens with F1 uncertainty intervals |
| **scGPT / Nicheformer Embeddings** *(v0.6)* | Backend-pluggable cell embedders with content-hashed cache + joint dissociated/spatial analysis |
| **Regulatory Gene-Set Expansion** *(v0.6)* | Borzoi or co-expression-backed suggestion tool for custom seed sets |
| **CRISPR Sequence Off-Target** *(v0.6)* | Evo 2 sequence-level scoring complementing pathway-level off-target detection |
| **Multi-Omics Fusion** *(v0.6)* | RNA + ATAC + proteomics weighted fusion with RNA/protein discordance flagger |
| **Causal Inference** *(v0.6)* | Invariant causal prediction identifies pathway parents across environments |
| **Active Learning** *(v0.6)* | Uncertainty / diversity / hybrid sample selection under a fixed label budget |
| **Discreteness Gate** *(v0.8)* | `DiscretenessGateA` — asks whether structure is *discrete* or a reproducibly-sliced continuum, using a single-Gaussian (SigClust) reference in place of the independence null. Reports gap statistic + Hartigan dip as diagnostics. Three outcomes: certify / reject / not-testable |
| **Deep-Learning Baselines** *(v0.8)* | `clustering_dl` — DEC (Xie 2016) and VAE-GMM (VaDE, Jiang 2017) baselines for benchmarking pathway-GMM against deep clustering |
| **Genetic / Somatic Anchoring** *(v0.8)* | Gate 7 tests whether a partition tracks real genetic strata (validated against TCGA-CRC BRAF-V600E / KRAS / MSI) |
| **KG Sensitivity (Gate K)** | Tests whether a finding survives a knowledge-graph version swap, against a size-matched degree-preserving null. Catches subtypes that are artifacts of the database release they were derived under |
| **GNN Embeddings** | TransE/RotatE KG embeddings, OntologyAwareGNN with biological attention, gene risk classification *(experimental)* |
| **Autism Interpretation** | Neuro-symbolic rules (R1-R7), therapeutic hypothesis ranking with safety flags *(autism-only)* |
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
pip install pathway-subtyping[vcf]           # VCF file processing (pysam)
pip install pathway-subtyping[viz]           # Interactive visualizations (Plotly, UMAP)
pip install pathway-subtyping[sc]            # Single-cell support (AnnData)
pip install pathway-subtyping[graph]         # Network analysis (NetworkX, py4cytoscape)
pip install pathway-subtyping[qc]            # Molecular QC layer for manufacturing
pip install pathway-subtyping[harmonize]     # v0.6 F2 — UCE cross-platform harmonization
pip install pathway-subtyping[perturb]       # v0.6 F5 — Geneformer in-silico perturbation
pip install pathway-subtyping[embed]         # v0.6 F6/F8 — scGPT + Nicheformer embedders
pip install pathway-subtyping[genesets]      # v0.6 F7 — Borzoi regulatory gene-set expansion
pip install pathway-subtyping[qc-sequence]   # v0.6 F9 — Evo 2 CRISPR off-target sequence scoring
pip install pathway-subtyping[gnn]           # Graph neural networks (PyTorch) — experimental
pip install pathway-subtyping[autism]        # Autism-specific interpretation (pure Python)
pip install pathway-subtyping[discreteness]  # v0.8 — Hartigan dip test for Gate A (diptest, GPLv2+)
pip install pathway-subtyping[all]           # Everything
```

Also available: `[docs]` and `[notebooks]` (documentation and notebook tooling), and
`[dev]` (test/lint stack — see Development below).

The v0.6 foundation-model extras (`harmonize`, `perturb`, `embed`,
`genesets`, `qc-sequence`) each install the PyTorch substrate that the
corresponding `Official*Backend` needs (`torch>=2.0`; `perturb` additionally
pulls `transformers>=4.35`). The model-specific upstream package
(`geneformer`, `scgpt`, `nicheformer`, `borzoi`, `evo2`, `uce`) and its
checkpoint must be installed separately because not all of them ship
stable PyPI wheels yet. Every wrapper also ships a deterministic PCA-based
fallback that works with *none* of these extras, so tests and local smoke
runs don't depend on heavyweight checkpoints.

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

31 Jupyter notebooks covering tutorials, full manuscript reproduction, and v0.6 feature walkthroughs (21-30). See [docs/notebook-guide.md](docs/notebook-guide.md) for execution order and dependencies.

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

**v0.6 Feature Walkthroughs (21-30)** -- Standalone, any order:

| # | Feature | Topic | Extra |
|---|---------|-------|-------|
| 21 | F1 uncertainty | Conformal intervals + bootstrap MSV + calibration | (core) |
| 22 | F2 cross-platform | UCE harmonization across 10x / Smart-seq2 / bulk | `[harmonize]` |
| 23 | F5 perturbation | Geneformer in-silico knockouts + MSV head | `[perturb]` |
| 24 | F6 embeddings | scGPT + content-hashed cache | `[embed]` |
| 25 | F7 gene sets | Regulatory expansion (Borzoi or co-expression) | `[genesets]` |
| 26 | F8 spatial | Nicheformer joint dissociated + spatial embed | `[embed]` |
| 27 | F9 Evo 2 | Sequence-level off-target scoring | `[qc-sequence]` |
| 28 | F10 multi-omics | RNA + ATAC + protein fusion + discordance | (core) |
| 29 | F11 causal | Invariant causal prediction | (core) |
| 30 | F12 active learning | Sample-selection strategies | (core) |

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
docker pull rohitdataops/pathway-subtyping:latest         # CLI runtime
docker pull rohitdataops/pathway-subtyping:0.6.3-jupyter  # Jupyter + notebooks

# Run pipeline
docker compose run pipeline

# Run tests (1,675 tests on the public edition)
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

**5-Layer Architecture:**

| Layer | Extra | What It Adds |
|-------|-------|-------------|
| **Core** | *(none)* | GMT scoring, GMM clustering, validation gates, simulation |
| **Graph** | `[graph]` | KG builder, network propagation, topology-weighted scoring |
| **QC** | `[qc]` | 12-feature manufacturing QC (cascade, dosage, drift, off-target, stress) |
| **GNN** | `[gnn]` | TransE/RotatE embeddings, OntologyAwareGNN, biological attention *(experimental)* |
| **Autism** | `[autism]` | Neuro-symbolic rules (R1-R7), therapeutic hypothesis ranking *(autism-only)* |

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
- **Discreteness (Gate A, v0.8)**: is the structure *discrete*, or a reproducibly-sliced
  continuum? Compares the observed bootstrap ARI against a **single-Gaussian (SigClust)
  reference** — `DiscretenessGateA` in `pathway_subtyping.discreteness`. This is the gate
  that decides discreteness.
- **Bootstrap stability** (ARI > 0.8): retained and reported, but **demoted in v0.8** to a
  marginal/confound control. It tested *pathway independence*, not discreteness, and could
  certify a continuous gradient as a subtype — the finding behind the 2026-07 correction.
  It no longer decides discreteness on its own.
- **Ancestry independence**: Clusters should not correlate with ancestry PCs (when provided)
- **Cross-modal concordance**: Subtypes should replicate across data modalities (when multi-omic)
- **Confound association (Gate 6)** and **genetic / somatic anchoring (Gate 7)**: does the
  partition track a known confound, or a real genetic stratum?
- **Time-sliced KG sensitivity (Gate K)**: does the finding survive being recomputed
  against a *different knowledge-graph version*? Every gate above holds the KG fixed and
  varies the data; this one varies the KG, and catches subtypes that are artifacts of the
  database release they were derived under —
  `pathway_subtyping.kg_sensitivity.kg_timeslice_sensitivity`.

⚠️ **Two properties of Gate A to know before quoting it.** It has three outcomes —
certify, reject, and `not-testable` (an abstention) — and on synthetic negative controls
it abstained on 28/30, so any false-positive rate must be quoted with its *testable*
denominator. It is also conservative: on a separation sweep it certifies nothing below
δ=2.5 SD. It is built to refuse over-called subtypes, not to detect subtle real ones.
See [`docs/discreteness_gate.md`](docs/discreteness_gate.md).

### 4. Statistical Rigor
Publication-quality statistics:
- **FDR correction**: Benjamini-Hochberg for multiple testing
- **Effect sizes**: Cohen's d with 95% bootstrap confidence intervals
- **Power analysis**: Sample size recommendations for target effect sizes
- **Type I error**: Estimation via null simulations

See [docs/how-it-works.md](docs/how-it-works.md) for a plain-language conceptual guide, or [docs/METHODS.md](docs/METHODS.md) for full statistical methodology.

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
│   ├── network_propagation.py   # PPI network propagation (RWR, PageRank)
│   ├── variant_qc.py            # Variant QC (QUAL, HWE, MAF, call rate)
│   ├── validation_datasets.py   # ClinVar/Reactome integration
│   ├── data_quality.py          # VCF quality checks
│   ├── utils/                   # Performance, seeding, progress tracking
│   ├── knowledge_graph/         # [graph] KG builder, schema, exporters
│   ├── qc/                      # [qc] 12-feature molecular QC layer
│   ├── gnn/                     # [gnn] GNN embeddings, attention, model (experimental)
│   └── autism/                  # [autism] Neuro-symbolic rules, therapeutic ranking
├── configs/                     # Example YAML configurations
├── data/
│   ├── pathways/                # Pathway GMT files (6 diseases)
│   └── sample/                  # Synthetic test data
├── docs/
│   ├── METHODS.md               # Statistical methods documentation
│   ├── api/                     # API reference (22+ modules)
│   └── guides/                  # User guides
├── scripts/                     # Utility scripts
│   ├── generate_cytoscape_figures.py  # Publication-ready network figures (requires [graph] + Cytoscape desktop)
│   ├── validate_with_public_data.py   # ClinVar/Reactome validation
│   └── benchmark_performance.py       # Performance benchmarks
├── examples/notebooks/          # Jupyter tutorials
├── tests/                       # Test suite (1,675 tests)
├── Dockerfile                   # Container support
└── docker-compose.yml           # Easy orchestration
```

## Development

```bash
# Install with dev dependencies (from cloned repo)
pip install -e ".[dev,vcf,viz,sc,graph,qc,autism]"

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
DOI: 10.5281/zenodo.18638048
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
