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
| KEGG | Kyoto Encyclopedia of Genes and Genomes |
| LoF | Loss of Function |
| LOEUF | Loss-of-function Observed/Expected Upper bound Fraction |
| PCA | Principal Component Analysis |
| QC | Quality Control |
| VCF | Variant Call Format |
| VEP | Variant Effect Predictor |

## See Also

- [Framework Overview](framework_overview.md) - High-level architecture
- [Data Formats](data_formats.md) - Detailed format specifications
- [API Reference](api/index.md) - Module documentation
