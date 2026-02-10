# Pathway Subtyping Framework

**A Disease-Agnostic Tool for Pathway-Based Molecular Subtype Discovery**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18442427.svg)](https://doi.org/10.5281/zenodo.18442427)
[![PyPI version](https://badge.fury.io/py/pathway-subtyping.svg)](https://pypi.org/project/pathway-subtyping/)
[![CI](https://github.com/topmist-admin/pathway-subtyping-framework/actions/workflows/ci.yml/badge.svg)](https://github.com/topmist-admin/pathway-subtyping-framework/actions/workflows/ci.yml)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## Overview

The **Pathway Subtyping Framework** is an open-source computational tool for identifying molecular subtypes in genetically heterogeneous diseases. Instead of analyzing individual genes, it aggregates rare variant burden at the **biological pathway level**, enabling:

- Better signal detection across genetically diverse cohorts
- Identification of biologically coherent patient subgroups
- Cross-cohort validation of discovered subtypes

Originally developed for [autism research](https://github.com/topmist-admin/autism-pathway-framework), this generalized version can be adapted for any disease with:
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
| **Multiple Clustering** | GMM, K-means, Hierarchical, Spectral with cross-validation |
| **Ancestry Correction** | PCA-based population stratification correction with independence testing |
| **Batch Correction** | ComBat-style batch effect detection and correction |
| **Sensitivity Analysis** | Parameter robustness testing across algorithms, features, normalization |
| **Validation Gates** | Negative controls + bootstrap stability + ancestry independence testing |
| **Statistical Rigor** | FDR correction, effect sizes, confidence intervals |
| **Power Analysis** | Sample size recommendations, Type I error estimation |
| **Simulation** | Synthetic data generation with ground truth for validation |
| **Reproducibility** | Deterministic execution, pinned dependencies, Docker |
| **Config-Driven** | YAML configuration for all parameters |

## Quick Start

### Installation

```bash
pip install pathway-subtyping
```

For VCF file processing, install with the `vcf` extra:

```bash
pip install pathway-subtyping[vcf]
```

### Try in Browser (No Installation)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/00_quick_demo.ipynb)

**60-second demo** — generates a synthetic cohort, discovers subtypes, validates them, and visualizes results. No data needed.

**Full tutorial**: [01_getting_started.ipynb](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/01_getting_started.ipynb)

### Run with Sample Data

```bash
# Clone for sample data and configs
git clone https://github.com/topmist-admin/pathway-subtyping-framework
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
# Run pipeline
docker-compose run pipeline

# Run tests
docker-compose run test

# Start Jupyter notebook
docker-compose up jupyter
# Open http://localhost:8888
```

## Adapting for Your Disease

1. **Create a pathway GMT file** with disease-relevant gene sets
2. **Copy an example config** and point to your data
3. **Run the pipeline** — validation gates will tell you if subtypes are meaningful

See the full guide: [Adapting for Your Disease](docs/guides/adapting-for-your-disease.md)

## How It Works

```
VCF Input → Variant Filter → Gene Burden → Pathway Aggregation → [Ancestry Correction] → [Batch Correction] → GMM Clustering → [Sensitivity Analysis] → Validation → Report
```

### 1. Pathway Scoring
Rare damaging variants are aggregated into pathway-level disruption scores:
- Loss-of-function variants weighted higher
- Missense variants weighted by CADD score
- Scores normalized across samples

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
├── src/pathway_subtyping/     # Core Python package
│   ├── pipeline.py            # Main pipeline
│   ├── clustering.py          # Multiple clustering algorithms
│   ├── statistical_rigor.py   # FDR, effect sizes, burden weights
│   ├── simulation.py          # Synthetic data & power analysis
│   ├── validation.py          # Validation gates
│   ├── ancestry.py            # Population stratification correction
│   ├── batch_correction.py    # Batch effect detection & correction
│   ├── sensitivity.py         # Parameter sensitivity analysis
│   └── data_quality.py        # VCF quality checks
├── configs/                   # Example YAML configurations
├── data/
│   ├── pathways/              # Pathway GMT files (6 diseases)
│   └── sample/                # Synthetic test data
├── docs/
│   ├── METHODS.md             # Statistical methods documentation
│   └── guides/                # User guides
├── examples/notebooks/        # Jupyter tutorials
├── tests/                     # Test suite (383 tests)
├── Dockerfile                 # Container support
└── docker-compose.yml         # Easy orchestration
```

## Development

```bash
# Install with dev dependencies (from cloned repo)
pip install -e ".[dev,vcf]"

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

- **[Autism Pathway Framework](https://github.com/topmist-admin/autism-pathway-framework)** — The original autism-focused implementation with SFARI cohort validation

## Contributing

Contributions welcome! Areas where help is needed:
- Additional disease pathway definitions
- Performance optimization for large cohorts
- Documentation and tutorials

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Citation

If you use this framework, please cite:

```
Chauhan R. Pathway Subtyping Framework. GitHub. 2026.
https://github.com/topmist-admin/pathway-subtyping-framework
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
- GitHub: [@topmist-admin](https://github.com/topmist-admin)

---

> **RESEARCH USE ONLY** — This framework is for hypothesis generation. Not for clinical diagnosis or treatment decisions.
