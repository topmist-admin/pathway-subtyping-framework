# Pathway Subtyping Framework

**A Disease-Agnostic Tool for Pathway-Based Molecular Subtype Discovery**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
| *Your disease* | [Adapt it →](docs/guides/adapting-for-your-disease.md) | `your_pathways.gmt` |

## Key Features

| Feature | Description |
|---------|-------------|
| **Pathway Scoring** | Aggregate gene burdens across biological pathways |
| **Subtype Discovery** | GMM clustering with automatic model selection (BIC) |
| **Validation Gates** | Negative controls + bootstrap stability testing |
| **Reproducibility** | Deterministic execution, pinned dependencies, Docker |
| **Config-Driven** | YAML configuration for all parameters |

## Quick Start

### Installation

```bash
git clone https://github.com/topmist-admin/pathway-subtyping-framework
cd pathway-subtyping-framework
pip install -r requirements.txt
pip install -e .
```

### Run Demo

```bash
# Run with example autism data
python -m pathway_subtyping --config configs/example_autism.yaml

# Run with your disease config
python -m pathway_subtyping --config configs/your_disease.yaml
```

### Try in Browser (No Installation)

[Open in Google Colab](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/01_demo_generic.ipynb)

## Adapting for Your Disease

1. **Create a pathway GMT file** with disease-relevant gene sets
2. **Copy an example config** and point to your data
3. **Run the pipeline** — validation gates will tell you if subtypes are meaningful

See the full guide: [Adapting for Your Disease](docs/guides/adapting-for-your-disease.md)

## How It Works

```
VCF Input → Variant Filter → Gene Burden → Pathway Aggregation → GMM Clustering → Validation → Report
```

### 1. Pathway Scoring
Rare damaging variants are aggregated into pathway-level disruption scores:
- Loss-of-function variants weighted higher
- Missense variants weighted by CADD score
- Scores normalized across samples

### 2. Subtype Discovery
Gaussian Mixture Model clustering identifies patient subgroups:
- Automatic cluster selection via BIC
- Configurable cluster range (default: 2-8)

### 3. Validation Gates
Built-in tests prevent overfitting:
- **Label shuffle**: Randomized labels should NOT cluster (ARI < 0.15)
- **Random genes**: Fake pathways should NOT work (ARI < 0.15)
- **Bootstrap**: Clusters should be stable under resampling (ARI > 0.8)

## Data Requirements

| Input | Format | Notes |
|-------|--------|-------|
| Variants | VCF | Annotated with gene symbols, consequences |
| Phenotypes | CSV | Sample IDs + clinical features |
| Pathways | GMT | Gene sets for your disease |

**Your data stays on your infrastructure.** The framework runs locally or in your cloud environment.

## Project Structure

```
pathway-subtyping-framework/
├── src/pathway_subtyping/     # Core Python package
├── configs/                   # Example YAML configurations
│   ├── example_autism.yaml
│   ├── example_schizophrenia.yaml
│   └── example_epilepsy.yaml
├── data/pathways/             # Pathway GMT files
│   ├── autism_pathways.gmt
│   ├── schizophrenia_pathways.gmt
│   └── README.md              # How to create your own
├── docs/guides/               # Documentation
│   ├── adapting-for-your-disease.md
│   ├── pathway-curation-guide.md
│   └── validation-gates.md
├── examples/notebooks/        # Jupyter notebooks
└── tests/                     # Test suite
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
