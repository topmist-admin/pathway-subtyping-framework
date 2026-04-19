# Terminology and Glossary

> This document defines key terms used throughout the Pathway Subtyping Framework to ensure consistent usage across all documentation and code.

## Genetic Concepts

### Variant
A genetic difference from the reference genome at a specific position. In this framework, we typically focus on **rare variants** (allele frequency < 1%) that may contribute to disease risk.

### Consequence
The predicted functional effect of a variant on a gene. Common consequences include:
- **Frameshift**: Insertion or deletion that shifts the reading frame
- **Stop-gained (nonsense)**: Creates a premature stop codon
- **Missense**: Changes one amino acid to another
- **Synonymous**: No change to amino acid sequence
- **Splice site**: Affects mRNA splicing

### Loss-of-Function (LoF)
Variants that are predicted to severely disrupt gene function. Includes:
- Frameshift variants
- Stop-gained variants
- Splice donor/acceptor variants
- Start-lost variants

LoF variants are typically given the highest weight in burden calculations.

### Allele Frequency (AF)
The proportion of chromosomes in a population carrying a specific allele. Key thresholds:
- **Ultra-rare**: AF < 0.0001 (< 0.01%)
- **Rare**: AF < 0.01 (< 1%)
- **Common**: AF >= 0.01 (>= 1%)

### gnomAD
The Genome Aggregation Database - a reference database of genetic variation from ~140,000 individuals. Used for:
- Allele frequency filtering
- Gene constraint scores

## Scoring Concepts

### Gene Burden (Burden Score)
A numerical score representing the cumulative impact of variants in a gene for a single individual. Higher burden = more predicted disruption.

**Calculation**: Sum of weighted variant contributions
```
burden(gene, sample) = sum(weight(variant))
```

Where `weight(variant)` considers:
- Consequence type (LoF > missense > other)
- Pathogenicity scores (CADD, REVEL)
- Allele frequency (rarer = higher weight)

**Synonyms**: Gene-level score, disruption score

### CADD Score
Combined Annotation Dependent Depletion - a pathogenicity predictor that integrates multiple annotations. Scores are typically presented as **Phred-scaled**:
- CADD >= 10: Top 10% most deleterious
- CADD >= 20: Top 1% most deleterious
- CADD >= 30: Top 0.1% most deleterious

Used in this framework to weight missense variants.

### Pathway Score (Pathway Disruption Score)
A numerical score representing the aggregate genetic burden across all genes in a biological pathway for a single individual.

**Calculation**: Aggregate of gene burdens with size normalization
```
pathway_score = sum(burden(gene)) / sqrt(pathway_size)
```

**Synonyms**: Pathway disruption, pathway perturbation

### Z-Score
A standardized score indicating how many standard deviations a value is from the population mean.
```
z-score = (value - mean) / std
```

Used to normalize pathway scores across the cohort for comparability.

## Constraint Concepts

### Gene Constraint
A measure of how intolerant a gene is to mutations, inferred from the deficit of observed variants in gnomAD compared to expectation.

### pLI Score
Probability of being Loss-of-function Intolerant. Ranges 0-1.
- pLI >= 0.9: Highly constrained (LoF variants likely pathogenic)
- pLI < 0.1: Tolerant to LoF

### LOEUF Score
Loss-of-function Observed/Expected Upper bound Fraction. Lower = more constrained.
- LOEUF < 0.35: Highly constrained
- LOEUF > 1.0: Tolerant

## Pathway Concepts

### Pathway
A set of genes that function together in a biological process. Sources include:
- **GO (Gene Ontology)**: Hierarchical ontology of biological processes, cellular components, molecular functions
- **Reactome**: Curated pathway database with detailed molecular interactions
- **KEGG**: Pathway maps including metabolism, signaling, disease

### Pathway Size
The number of genes annotated to a pathway. We typically filter:
- Minimum size: 5 genes (too small = unstable)
- Maximum size: 500 genes (too large = uninformative)

### GMT Format
Gene Matrix Transposed - simple format for pathway definitions:
```
PATHWAY_NAME<tab>DESCRIPTION<tab>GENE1<tab>GENE2<tab>...
```

## Clustering Concepts

### Subtype
A genetically coherent subgroup of individuals discovered through unsupervised clustering of pathway profiles. Subtypes represent **hypotheses** about distinct etiological mechanisms.

### GMM (Gaussian Mixture Model)
A probabilistic clustering method that models data as a mixture of Gaussian distributions. Provides soft (probabilistic) cluster assignments.

Key properties:
- Assumes clusters are Gaussian-distributed
- Returns probability of membership in each cluster
- Number of clusters selected by BIC

### Stability
The consistency of clustering results under resampling (bootstrap). Measured as the proportion of sample pairs that are consistently co-clustered. Target: stability >= 0.7.

### BIC (Bayesian Information Criterion)
A model selection criterion balancing fit and complexity. Lower BIC = better model. Used to select optimal number of clusters.

Formula:
```
BIC = -2 * log_likelihood + k * log(n)
```
Where k is the number of parameters and n is the number of samples.

## Validation Concepts

### ARI (Adjusted Rand Index)
A measure of agreement between two clusterings, adjusted for chance. Ranges from -1 to 1:
- ARI = 1: Perfect agreement
- ARI = 0: Agreement by chance
- ARI < 0: Less agreement than expected by chance

Used to compare:
- Discovered clusters vs. ground truth
- Clusters across bootstrap samples

### Negative Control
A test designed to produce a null (negative) result. In validation:
- **Label shuffle**: Shuffling labels should yield ARI ~ 0
- **Random gene sets**: Random pathways should not cluster meaningfully

### Bootstrap
A resampling method where samples are drawn with replacement. Used to:
- Estimate stability of clusters
- Calculate confidence intervals
- Assess robustness of findings

### Validation Gate
A mandatory statistical test that must pass before results are considered trustworthy. The framework implements three gates:
1. Label shuffle (negative control)
2. Random gene sets (negative control)
3. Bootstrap stability (positive control)

## Data Format Terms

### VCF
Variant Call Format - standard file format for genetic variants. Contains:
- Chromosome, position, reference/alternate alleles
- Quality and filter status
- Per-sample genotype information

Required INFO fields for this framework:
- GENE: Gene symbol
- CONSEQUENCE: Variant consequence
- CADD: Pathogenicity score (optional but recommended)

### Phenotype File
CSV file containing sample metadata:
```csv
sample_id,sex,age,cohort
SAMPLE_001,M,25,discovery
SAMPLE_002,F,32,discovery
```

### Configuration File
YAML file specifying pipeline parameters:
```yaml
pipeline:
  seed: 42
clustering:
  n_clusters_range: [2, 8]
```

## Single-Cell Concepts

### scRNA-seq (Single-Cell RNA Sequencing)
High-throughput sequencing of mRNA from individual cells. Produces a cell x gene expression matrix, typically with high sparsity (60-95% zeros) due to dropout.

### AnnData / h5ad
The standard Python data structure and file format for single-cell data, used by Scanpy and many scRNA-seq tools. An h5ad file stores:
- `.X`: Expression matrix (often sparse)
- `.obs`: Cell-level metadata (cell type, sample, etc.)
- `.var`: Gene-level metadata

### Pseudobulk
Aggregating single-cell expression within cell types to produce a cell_types x genes matrix. This matrix is structurally identical to a bulk expression matrix and can be scored using the same methods (ssGSEA, GSVA, mean-Z).

### UMI (Unique Molecular Identifier)
Short barcodes added during library prep to count individual mRNA molecules. UMI counts are the standard measure of expression in droplet-based scRNA-seq (10X Genomics).

### Dropout
The phenomenon where a gene expressed in a cell is not detected due to sampling limitations. Results in excess zeros in scRNA-seq data. The framework handles this via sparse matrix operations.

### Cell Type
A category of cells with shared transcriptomic profiles (e.g., T cells, B cells, neurons). Cell type annotations in `.obs` are used to group cells for pseudobulk aggregation.

## v0.6 Terms (Rigor + Foundation-Model Layers)

### MSV (Molecular State Vector)
The per-sample vector of pathway scores that PSF v0.5 produces. In v0.6 the term is reused for richer variants:
- **Bootstrap MSV** — the distribution of MSVs across bootstrap resamples; used to summarise uncertainty on each pathway (F1).
- **MSV-from-embedding** — the linear head (`MSVFromEmbedding`) that maps a foundation-model's per-cell embedding into a pathway-score vector. Lets F5 perturbation screens predict ΔMSV directly from embedding shifts.

### Conformal (split conformal prediction)
A distribution-free uncertainty quantification method: given any black-box predictor, a held-out calibration set, and a target coverage ``α``, conformal intervals are guaranteed to cover the true value with probability ``α`` — assuming exchangeability between calibration and test data. In PSF F1, pathway-score predictions are wrapped with `ConformalPathwayPredictor`; intervals are finite-sample oracle-tight, not just asymptotic. See [guides/uncertainty.md](guides/uncertainty.md).

### Calibration (probabilistic, not instrument-wise)
Refers to whether predicted probabilities (or prediction intervals) are empirically accurate — e.g., a 90% conformal interval that covers the truth on 90% of held-out cases is "calibrated." `CalibrationReport` is the F1 diagnostic for this. Distinct from the v0.5 QC sense of calibration (which usually means "tuning a threshold").

### Hallmark (as in "hallmark pathways")
The 50 curated gene sets in MSigDB's Hallmark collection (Liberzon et al. 2015) — each represents a well-defined biological process (apoptosis, hypoxia, MYC targets, …). PSF uses them as the default pathway set for validation and real-data acceptance runs because they are well-studied, non-overlapping, and relatively robust. The GMT file ships at `data/pathways/hallmark_200genes.gmt`.

### Cascade (molecular QC sense)
One of 12 manufacturing QC signals (`pathway_subtyping.qc`) that look for *runaway* transcriptional states — where cells have left a "manufacturable" state and drifted toward, e.g., an unrecovered stress response or senescence trajectory. Not to be confused with the generic biological usage ("signaling cascade"). See [api/qc.md](api/qc.md).

### Dosage (QC / variant sense)
Dosage-aware variant scoring: when a variant's functional impact is expected to scale with the number of affected alleles (copy number or zygosity), the QC cascade adjusts its contribution to the pathway score. The AlphaMissense extension (F4) refines gene weights using per-variant pathogenicity probabilities.

### Drift / Stress / Off-target / Tension / Crosstalk
Five of the 12 molecular QC features (`pathway_subtyping.qc`):
- **Drift** — detects when a cell state moves away from its reference distribution across passages.
- **Stress** — detects elevated stress-response pathway activity (glycolysis, UPR, hypoxia, P53).
- **Off-target** — flags unintended gene activation relative to the on-target signature.
- **Tension** — flags pathway pairs that disagree under expected biological constraints.
- **Crosstalk** — flags pathway pairs that co-activate unexpectedly (regulatory coupling that shouldn't be present).
See [api/qc.md](api/qc.md) for the full cascade.

### UCE (Universal Cell Embeddings)
A foundation model from Rosen et al. (*Nature*, 2024) that produces platform-invariant 1280-dim per-cell embeddings. In F2 (`pathway_subtyping.harmonize`) the embedding is the *anchor* that the `CrossPlatformAligner` uses to separate biology from platform drift. Opt-in via `[harmonize]`; a deterministic `FallbackEmbedder` (PCA on pooled input) preserves the API for CI.

### Geneformer
A foundation model from Theodoris et al. (*Nature*, 2023) that tokenises scRNA-seq by median-normalised gene rank and predicts cell-context from masked tokens. In F5 (`pathway_subtyping.perturb`), the `OfficialBackend` wraps the published V2 104M checkpoint for in-silico knockout: removing a gene's token from the sequence simulates loss of function. Apache-2.0 licensed. A `FallbackPerturber` (PCA) preserves the API without the checkpoint. See [guides/perturbation.md](guides/perturbation.md).

### Content-hashed embedding cache
An on-disk cache (`OfficialBackend(cache_dir=...)` in F5; `EmbeddingCache` in F6) keyed by a SHA-256 over `(checkpoint + embedding-config + expression bytes)`. Reruns on the same cohort return in sub-millisecond time instead of a full forward pass through the foundation model. See [guides/embeddings.md](guides/embeddings.md).

### Validation Gate 4 / 5
Extensions of the three v0.5 validation gates. Gate 4 (ancestry independence) asserts subtype assignments are not simply tracking ancestry PCs; Gate 5 (cross-modal consistency) asserts that subtypes agree across independent data modalities (e.g., WES vs RNA-seq) when `per_modality_scores` is supplied. Automatically included by `ValidationGates.run_all()` when their inputs are present.

### ICP (Invariant Causal Prediction)
A causal-inference approach used in F11 (`pathway_subtyping.causal`): test which subset of pathways produces invariant conditional distributions of the outcome across multiple "environments" (cohorts, platforms, tissues). The invariant subset is the causally-grounded feature set; non-invariant subsets are rejected as confounded. PSF's implementation extends standard ICP with a combined mean+variance invariance test.

### Active Sample Selection
F12 (`pathway_subtyping.active`): given a partially-labelled cohort, pick the next few samples to send to the wet lab that maximise expected information gain. Three strategies ship — uncertainty (pick samples with widest conformal intervals), diversity (pick samples farthest from labelled set), and hybrid (convex combination).

## Abbreviations

| Abbreviation | Full Term |
|--------------|-----------|
| AF | Allele Frequency |
| ARI | Adjusted Rand Index |
| BIC | Bayesian Information Criterion |
| CADD | Combined Annotation Dependent Depletion |
| GMM | Gaussian Mixture Model |
| GMT | Gene Matrix Transposed |
| GO | Gene Ontology |
| ICP | Invariant Causal Prediction |
| KEGG | Kyoto Encyclopedia of Genes and Genomes |
| LoF | Loss of Function |
| LOEUF | Loss-of-function Observed/Expected Upper bound Fraction |
| MSV | Molecular State Vector |
| PCA | Principal Component Analysis |
| QC | Quality Control |
| scRNA-seq | Single-Cell RNA Sequencing |
| UCE | Universal Cell Embeddings |
| UMI | Unique Molecular Identifier |
| VCF | Variant Call Format |
| VEP | Variant Effect Predictor |

## See Also

- [Framework Overview](framework_overview.md) - High-level architecture
- [Data Formats](data_formats.md) - Detailed format specifications
- [API Reference](api/index.md) - Module documentation
