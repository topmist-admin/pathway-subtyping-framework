# Pathway Subtyping Framework

**A Research Tool for Pathway-Based Molecular Subtyping of Genetically Heterogeneous Diseases**

---

## The Problem

Neuropsychiatric and neurodevelopmental disorders are genetically complex—hundreds of genes, diverse mechanisms, and significant heterogeneity across individuals. Traditional gene-centric analyses often fail to generalize across cohorts.

## Our Approach

The **Pathway Subtyping Framework** shifts the unit of analysis from individual genes to **biological pathways**, enabling:

- Integration of rare variants into pathway-level disruption scores
- Unsupervised clustering to identify biologically coherent subgroups
- Built-in validation gates to prevent overfitting
- Disease-agnostic design adaptable to any condition

## Key Features

| Feature | Description |
|---------|-------------|
| **Pathway Scoring** | Aggregate gene burdens across curated biological pathways |
| **Subtype Discovery** | GMM clustering with automatic cluster selection |
| **Validation Gates** | Negative controls + bootstrap stability testing |
| **Reproducibility** | Deterministic execution with pinned dependencies |
| **Colab-Ready** | Run the demo in your browser—no installation |

## Pre-Built Disease Pathways

- Autism Spectrum Disorder
- Schizophrenia
- Epilepsy
- Intellectual Disability

## What We're Looking For

We're seeking **collaborators with cohort data** (N ≥ 100) to:

1. **Validate** the framework on independent cohorts
2. **Extend** pathway definitions with domain expertise
3. **Co-author** publications on subtype discovery

## Data Requirements

| Input | Format | Notes |
|-------|--------|-------|
| Variants | VCF | Annotated with gene symbols |
| Phenotypes | CSV | Sample IDs + clinical features |
| (Optional) Expression | CSV | For multi-omic integration (v0.3) |

**Your data stays on your infrastructure.** The framework runs locally or in your cloud environment.

## Getting Started

```bash
# 5-minute setup
pip install pathway-subtyping
pathway-subtyping run --config configs/your_disease.yaml
```

Or try the **[Colab notebook](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/01_demo.ipynb)** (no installation required).

## Links

- **GitHub**: https://github.com/topmist-admin/pathway-subtyping-framework
- **PyPI**: https://pypi.org/project/pathway-subtyping/
- **DOI**: https://doi.org/10.5281/zenodo.14765212
- **Documentation**: https://github.com/topmist-admin/pathway-subtyping-framework/tree/main/docs

## Contact

**Rohit Chauhan**
GitHub: [@topmist-admin](https://github.com/topmist-admin)

---

> **RESEARCH USE ONLY** — This framework is for hypothesis generation. Not for clinical diagnosis or treatment decisions.
