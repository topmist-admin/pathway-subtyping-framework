# One-Pager Template

Use this template to create a one-pager for your disease application.

---

# [Disease Name] Pathway Framework

**A Research Tool for Pathway-Based Analysis of Genetic Heterogeneity in [Disease]**

---

## The Problem

[Disease] is genetically complex—[number] of genes, diverse mechanisms, and significant heterogeneity across individuals. Traditional gene-centric analyses often fail to generalize across cohorts.

## Our Approach

The **[Disease] Pathway Framework** shifts the unit of analysis from individual genes to **biological pathways**, enabling:

- Integration of rare variants into pathway-level disruption scores
- Network-based signal refinement using protein-protein interactions
- Unsupervised clustering to identify biologically coherent subgroups
- Built-in validation gates to prevent overfitting

## Key Features

| Feature | Description |
|---------|-------------|
| **Pathway Scoring** | Aggregate gene burdens across [X] biological pathways |
| **Subtype Discovery** | GMM clustering with automatic cluster selection |
| **Validation Gates** | Negative controls + bootstrap stability testing |
| **Reproducibility** | Deterministic execution with pinned dependencies |
| **Colab-Ready** | Run the demo in your browser—no installation |

## What We're Looking For

We're seeking **collaborators with [disease] cohort data** (N ≥ 100) to:

1. **Validate** the framework on independent cohorts
2. **Extend** pathway definitions with domain expertise
3. **Co-author** publications on subtype discovery

## Data Requirements

| Input | Format | Notes |
|-------|--------|-------|
| Variants | VCF | Annotated with gene symbols |
| Phenotypes | CSV | Sample IDs + clinical features |
| (Optional) Expression | CSV | For multi-omic integration |

**Your data stays on your infrastructure.** The framework runs locally or in your cloud environment.

## Getting Started

```bash
# 5-minute setup
git clone https://github.com/[username]/[disease]-pathway-framework
cd [disease]-pathway-framework
pip install -r requirements.txt && pip install -e .
make demo  # Run on synthetic data
```

Or try the **[Colab notebook](https://colab.research.google.com/github/[username]/[disease]-pathway-framework/blob/main/examples/notebooks/01_demo.ipynb)** (no installation required).

## Links

- **GitHub**: https://github.com/[username]/[disease]-pathway-framework
- **DOI**: [Your Zenodo DOI]
- **Preprint**: [Your preprint DOI]
- **Documentation**: https://github.com/[username]/[disease]-pathway-framework/tree/main/docs

## Contact

**[Your Name]**
Email: [your email]
GitHub: [@username](https://github.com/username)

---

> **RESEARCH USE ONLY** — This framework is for hypothesis generation. Not for clinical diagnosis or treatment decisions.

---

# Instructions for Customization

Replace the following placeholders:

| Placeholder | Replace With |
|-------------|--------------|
| `[Disease Name]` | Your disease (e.g., "Schizophrenia", "Epilepsy") |
| `[Disease]` | Same, lowercase where appropriate |
| `[number]` | Approximate number of implicated genes |
| `[X]` | Number of pathways in your GMT file |
| `[username]` | Your GitHub username |
| `[Your Name]` | Your name |
| `[your email]` | Your contact email |
| `[Your Zenodo DOI]` | Your archival DOI (create at zenodo.org) |
| `[Your preprint DOI]` | Preprint DOI if available |

## Tips

1. **Keep it to one page** when printed
2. **Include the Colab link** — lowers barrier to trying
3. **Be specific** about data requirements
4. **State clearly** what you're looking for (validation, collaboration, etc.)
