---
title: 'Pathway Subtyping Framework: Validated Molecular Subtype Discovery from Pathway-Level Genomic Data'
tags:
  - Python
  - bioinformatics
  - genomics
  - molecular subtyping
  - pathway analysis
  - precision medicine
  - clustering
  - autism
  - schizophrenia
  - transcriptomics
authors:
  - name: Rohit Chauhan
    orcid: 0009-0003-9895-4629
    affiliation: 1
affiliations:
  - name: Topmist LLC, Georgia, USA
    index: 1
date: 18 February 2026
bibliography: paper.bib
---

# Summary

The Pathway Subtyping Framework is an open-source Python package for discovering molecular subtypes in genetically heterogeneous diseases. It shifts subtyping from individual genes to biological pathways — curated gene sets representing shared molecular processes — providing biologically motivated dimensionality reduction that improves statistical power and interpretability. The framework accepts gene expression matrices (bulk RNA-seq, microarray) or variant burden scores (VCF), aggregates them into pathway-level scores using ssGSEA, GSVA, or mean-Z methods, clusters patients using Gaussian Mixture Models (GMM), and enforces mandatory validation gates before accepting results. It is disease-agnostic: the same pipeline operates on any condition for which curated pathway definitions exist.

# Statement of Need

Over 100 risk genes contribute to autism spectrum disorder alone, each explaining a small fraction of cases [@Satterstrom2020]. Similar genetic architectures characterize schizophrenia [@Trubetskoy2022], epilepsy, and intellectual disability. This heterogeneity means that patients sharing a single clinical diagnosis may have fundamentally different molecular disruptions — and may respond to different treatments. Gene-level clustering is underpowered for typical cohort sizes, and existing pathway analysis tools focus on differential enrichment between predefined groups rather than unsupervised subtype discovery.

The Pathway Subtyping Framework addresses this gap by combining pathway-level feature engineering with validated unsupervised clustering. Unlike GSEA [@Subramanian2005] or GSVA [@Hanzelmann2013], which compare predefined groups, this framework discovers groups from the data. Unlike NMF-based approaches [@Brunet2004], it leverages curated biological pathway knowledge. And unlike any existing tool, it enforces mandatory validation gates that distinguish real molecular subtypes from statistical artifacts.

# Key Features

**Multi-modal input.** The framework accepts bulk RNA-seq expression matrices, microarray data, and VCF files (rare variant burden scoring). All inputs converge on a shared pathway aggregation engine.

**Three scoring methods.** Pathway-level scores are computed via ssGSEA [@Barbie2009], GSVA [@Hanzelmann2013], or mean-Z (Z-score averaging across member genes).

**Gaussian Mixture Model clustering.** GMM with full covariance, BIC-based model selection over k=2–8, 10 random restarts, and deterministic seeding.

**Five mandatory validation gates:**

1. **Label shuffle** (negative control): Verifies clusters are not random artifacts
2. **Random gene sets** (negative control): Verifies curated pathways drive clustering, not arbitrary gene groupings
3. **Bootstrap stability**: Verifies clusters are robust to resampling
4. **Ancestry independence**: Verifies clusters do not reflect population stratification
5. **Cross-modal concordance**: Verifies subtypes replicate across data modalities

**Correction modules.** ComBat batch correction, ancestry regression (PCA-based), and sensitivity analysis across clustering algorithms, feature subsets, and normalization methods.

**Benchmark suite.** Built-in comparison against PCA+K-means, NMF+K-means, gene-level K-means, and random baseline.

**Disease templates.** Pre-curated GMT pathway files for autism, schizophrenia, epilepsy, intellectual disability, Parkinson's disease, and bipolar disorder.

**Multi-omic integration (v0.3.0).** Single-cell pathway scoring (h5ad support), bulk deconvolution (NNLS), signaling pathway databases (CellPhoneDB, CellChatDB), and multi-omic fusion (concatenate, weighted average, intersection strategies).

# Real-Data Validation

The framework has been validated on three independent postmortem brain transcriptome datasets:

- **GSE28521** [@Voineagu2011]: 32 frontal cortex samples (16 ASD, 16 control). Discovered a GABA-Collapsed subtype (n=9, 100% ASD, Cohen's d=3.21 for GABA signaling).
- **GSE64018** [@Gupta2014]: 24 temporal cortex RNA-seq samples. Independent replication: disease-enriched subtype (83% ASD) with the same top 3 disrupted pathways (GABA, glutamate, cell adhesion).
- **GSE80655** [@Ramaker2017]: 141 schizophrenia + control samples, 3 brain regions. Three subtypes identified (Dopamine-Hyperactive, Neurodevelopment-Activated, Synaptic-Collapsed). First real-data analysis to pass all 3 validation gates (bootstrap ARI=0.923). Cross-disease analysis showed ASD and SCZ pathway sets produce 87% identical subtypes (ARI=0.870).

Across 444 brain samples and 3 independent datasets, the framework demonstrated that molecular subtypes are diagnosis-independent (chi-squared p=0.72) but brain-region-dependent (p=3×10⁻¹⁶), with direct implications for precision drug selection in psychiatric clinical trials.

# Software Architecture

The framework is implemented in Python (≥3.9) with core dependencies limited to NumPy, pandas, scikit-learn, and SciPy. Optional extras include `[vcf]` (pysam for VCF parsing), `[viz]` (Plotly for interactive visualization), and `[sc]` (anndata for single-cell data). It is distributed via PyPI (`pip install pathway-subtyping`), tested with 968 automated tests, and archived on Zenodo (DOI: 10.5281/zenodo.18442426). Thirteen tutorial notebooks (including three real-data GEO validation notebooks) are provided as Google Colab links for zero-install reproducibility.

# Acknowledgements

The GSE28521 dataset was generated by Voineagu et al. (2011), GSE64018 by Gupta et al. (2014), and GSE80655 by Ramaker et al. (2017), all publicly available via NCBI GEO. Pathway definitions are curated from SFARI Gene, KEGG, Reactome, and MSigDB.

# References
