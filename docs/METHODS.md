# Statistical Methods Documentation

This document provides detailed statistical methodology for the Pathway Subtyping Framework, suitable for methods sections in peer-reviewed publications.

## Overview

The Pathway Subtyping Framework identifies molecular subtypes in genetically heterogeneous diseases by:
1. Applying variant quality control filters (optional)
2. Computing gene-level burden scores from rare variant data
3. Aggregating burdens into pathway-level scores
4. Correcting for population stratification (optional)
5. Clustering samples based on pathway profiles
6. Validating discoveries through rigorous statistical tests

## Gene Burden Scoring

### Variant Quality Control

Before burden computation, an optional variant QC step removes technical artifacts and common variants that would dilute rare variant signal. Filters are applied sequentially:

1. **QUAL score**: Variants with phred-scaled quality below threshold (default: 30) are removed
2. **Call rate**: Variants with genotype missingness exceeding threshold (default: > 10%) are removed
3. **Hardy-Weinberg equilibrium**: Variants deviating from HWE (chi-squared test, default p < 1e-6) are removed, as extreme deviation often indicates genotyping error
4. **Minor Allele Frequency**: Variants with MAF above threshold (default: 1%) are removed to retain rare variants only

Per-genotype filters (GQ, DP) may also be applied to mask low-confidence individual calls before the variant-level filters above.

Reference: Anderson CA et al. Data quality control in genetic case-control association studies. *Nat Protoc*. 2010;5(9):1564-73.

### Variant Selection

After QC, qualifying variants are selected based on:
- **Minor Allele Frequency (MAF)**: < 1% in gnomAD (configurable threshold)
- **Consequence**: Loss-of-function (LoF) or missense variants
- **Quality**: Standard VCF quality filters (PASS status)

### Burden Weighting

Gene burdens are computed as weighted sums of variant counts:

```
burden_g = Σ_v (w_consequence × w_CADD × allele_count_v)
```

Where:
- `w_consequence`: LoF = 2.0, missense = 1.0 (default)
- `w_CADD`: CADD score normalized by percentile rank
- `allele_count`: Number of alternate alleles (0, 1, or 2)

#### Literature-Based Weight Schemes

The framework supports multiple weighting schemes (see `BurdenWeightScheme`):

| Scheme | LoF Weight | High-Impact Missense | Moderate Missense | VUS |
|--------|------------|---------------------|-------------------|-----|
| DEFAULT | 2.0 | 1.5 | 1.0 | 0.5 |
| GNOMAD_CONSTRAINT | pLI-scaled | CADD-scaled | CADD-scaled | 0.1 |
| ACMG_INSPIRED | 3.0 | 2.0 | 1.5 | 0.5 |
| UNIFORM | 1.0 | 1.0 | 1.0 | 1.0 |

### CADD Score Handling

For variants without CADD annotations, consequence-based defaults are applied:
- LoF variants: CADD = 35 (highly deleterious)
- Missense variants: CADD = 20 (moderate impact)
- Other variants: CADD = 10 (low impact)

## Pathway Score Computation

### Gene-to-Pathway Aggregation

Pathway scores are computed from gene burdens using GMT file definitions:

```
pathway_score_p = normalize(aggregate(burden_g for g in pathway_p))
```

### Aggregation Methods

The framework supports multiple aggregation approaches:

| Method | Formula | Use Case |
|--------|---------|----------|
| MEAN | Σburden / n | Default; simple average |
| MEDIAN | median(burdens) | Robust to outliers |
| SIZE_NORMALIZED | Σburden / √n | Corrects for pathway size bias |
| PCA | PC1(burdens) | Captures dominant signal |

### Normalization

After aggregation, pathway scores are Z-score normalized across samples:

```
z_p,i = (x_p,i - μ_p) / σ_p
```

This ensures equal contribution from each pathway regardless of scale.

## Signaling Pathway Databases

### Cell-Cell Signaling as Pathway Gene Sets

In addition to metabolic/disease pathway databases (Reactome, KEGG), cell-cell signaling interaction databases provide an alternative source of biologically meaningful gene sets. Ligand-receptor interactions are grouped by signaling pathway classification, and the union of all participating gene symbols forms the pathway gene set.

### CellPhoneDB

CellPhoneDB (Efremova et al., 2020) catalogs curated ligand-receptor interactions including multi-subunit complexes. The loading process:

1. **Download** three CSV files: `interaction_input.csv`, `gene_input.csv`, `complex_input.csv`
2. **Gene resolution**: Map UniProt accession IDs to HGNC gene symbols via `gene_input.csv`
3. **Complex resolution**: Decompose multi-subunit complexes (e.g., `integrin_a2b1_complex`) into constituent gene symbols via `complex_input.csv` subunit columns
4. **Pathway grouping**: Group interactions by their `classification` field (e.g., "WNT signaling", "Notch signaling")
5. **Gene set construction**: For each group, collect the union of all ligand and receptor gene symbols

### CellChatDB

CellChatDB (Jin et al., 2021) provides ligand-receptor interactions with explicit `pathway_name` annotations. Since CellChatDB distributes data in R binary format (`.rda`), users export to CSV from R and the framework parses the `pathway_name`, `ligand`, and `receptor` columns directly.

### Pathway Name Prefixing

All signaling pathway gene sets are prefixed with `"Signaling: "` (e.g., `"Signaling: WNT signaling"`) to distinguish them from metabolic or disease pathways when merged into a unified pathway dictionary.

### Merging Multiple Sources

When gene sets from multiple databases share the same pathway name, their gene lists are unioned to maximize coverage:

```
merged_genes(p) = genes_CellPhoneDB(p) ∪ genes_CellChatDB(p)
```

References:
- Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. *Nat Protoc*. 2020;15(4):1484-1506.
- Jin S, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun*. 2021;12(1):1088.

## Single-Cell Pathway Scoring

### Pseudobulk Aggregation

For scRNA-seq data, the framework computes pseudobulk expression profiles by averaging per-cell expression within each cell type:

```
pseudobulk_g,t = mean(expression_g,c)  for all cells c in cell type t
```

Where sparse matrices are handled via row-wise slicing to avoid densification. Cell types with fewer than `min_cells_per_type` (default: 10) cells are excluded.

The resulting pseudobulk matrix (cell_types x genes) is structurally identical to a bulk expression matrix and is scored using the same methods (ssGSEA, GSVA, mean-Z) described in [Pathway Score Computation](#pathway-score-computation).

### Per-Cell Scoring

For cell-level resolution, per-cell mean-Z scoring computes pathway scores for each individual cell:

```
score_p,c = mean(z(expression_g,c))  for g in pathway_p
```

Per-cell scores are then aggregated per cell type:

```
score_p,t = mean(score_p,c)  for all cells c in cell type t
```

Per-cell scoring is performed in chunks (default: 5000 cells) for memory efficiency and supports sparse matrices natively.

### Normalization

Raw UMI counts are log-normalized before scoring:

```
normalized_g,c = log(1 + expression_g,c / total_c × 10000)
```

Where `total_c` is the total UMI count for cell `c`. This follows the standard Seurat/Scanpy normalization procedure.

### Output

All single-cell scoring methods produce a Z-normalized cell_types x pathways matrix, directly compatible with the same clustering, validation, and characterization tools used for VCF and bulk expression data.

Reference: Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*.

## Expression Pathway Scoring

### Overview

For bulk RNA-seq data, pathway scores are computed directly from gene expression without variant calling. Three methods are available, all producing Z-normalized samples x pathways matrices identical in format to VCF-based scores.

### ssGSEA (Recommended)

Single-sample Gene Set Enrichment Analysis (Barbie et al., 2009):

1. For each sample, rank all genes by expression value
2. Walk along the ranked list; for gene at rank *r*:
   - If gene is in the pathway: step up by |r|^α
   - If gene is not in the pathway: step down by 1/(N - n_pathway)
3. Enrichment score = sum of the running statistic

```
ES_p,i = Σ_g∈pathway |rank(g)|^α / Σ_g∈pathway |rank(g)|^α  -  Σ_g∉pathway 1/(N - n_p)
```

Where α = 0.25 (default) controls the weight of the rank. Positive scores indicate pathway enrichment; negative scores indicate depletion.

Reference: Barbie DA et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*.

### GSVA (Simplified)

Gene Set Variation Analysis (Hänzelmann et al., 2013):

1. Compute empirical CDF per gene across all samples
2. Rank CDF values within each sample
3. Compute KS-like statistic: maximum deviation between the distributions of CDF-ranked pathway genes vs non-pathway genes

```
GSVA_p,i = max_j |F_pathway(j) - F_non-pathway(j)|
```

For publication-grade GSVA, precompute scores using the R GSVA package and import with `input_type=LOG2`.

Reference: Hänzelmann S et al. (2013). GSVA: gene set variation analysis for microarray and RNA-Seq data. *BMC Bioinformatics*.

### Mean-Z

Fast Z-score averaging:

1. Z-score normalize each gene across samples: `z_g,i = (x_g,i - μ_g) / σ_g`
2. For each pathway: compute mean of member gene Z-scores: `score_p,i = mean(z_g,i for g in pathway_p)`
3. Z-normalize resulting pathway scores

Best for quick exploration and very large datasets (>10,000 samples).

### Input Preprocessing

Expression matrices are preprocessed based on input type:

| Input Type | Transformation |
|-----------|----------------|
| COUNTS | `log2(x + 1)` |
| TPM/FPKM | `log2(x + 1)` if max > 20, else no-op |
| LOG2 | None (already transformed) |

Genes with all-zero expression, zero variance, or appearing in fewer than 3 samples are filtered.

## Bulk Deconvolution

### Overview

Estimates cell-type proportions from bulk RNA-seq using a single-cell reference profile. Proportions can be combined with pathway scores for cell-type-aware subtype discovery.

### NNLS Deconvolution

Non-negative least squares (NNLS) is used because cell-type proportions are inherently non-negative. For each bulk sample *i*:

```
minimize  ||R^T × p_i - b_i||²
subject to:  p_i ≥ 0
```

Where:
- `R` = reference profile matrix (cell_types × genes)
- `p_i` = proportion vector for sample *i* (length = n_cell_types)
- `b_i` = bulk expression vector for sample *i*

Post-processing: normalize proportions to sum to 1:

```
p_normalized = p / Σ(p)    if Σ(p) > 0
p_normalized = 1/K         otherwise (uniform)
```

This approach is used in CIBERSORT (Newman et al., 2015) and MuSiC (Wang et al., 2019).

### Reference Profile Construction

From single-cell data, a reference profile is built by averaging expression per cell type:

```
reference_g,t = mean(expression_g,c)  for all cells c of type t
```

Cell types with fewer than `min_cells_per_type` (default: 5) cells are excluded.

### Feature Combination

Pathway scores and cell-type proportions are combined for clustering:

```
combined = [(1-w) × Z(pathway_scores) | w × Z(proportions)]
```

Where `w` is the proportion weight (default: 0.5) and `Z()` denotes Z-normalization. Cell-type columns are prefixed with `CELLTYPE:` to distinguish from pathway features.

References:
- Newman AM et al. (2015). Robust enumeration of cell subsets from tissue expression profiles. *Nature Methods*.
- Wang X et al. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. *Nature Communications*.

## Multi-Omic Fusion

### Overview

Fuses pathway scores from multiple data modalities (VCF, expression, single-cell, deconvolution) into a unified feature matrix for integrated subtype discovery.

### Concatenation Strategy (Default)

Column-binds all modality scores with label prefixes:

```
fused = [mod1:pathway1, mod1:pathway2, ..., mod2:pathway1, mod2:pathway2, ...]
```

Handles partial sample overlap via missing data strategies:
- **IMPUTE_ZERO**: Fill missing values with 0
- **IMPUTE_MEAN**: Fill with column mean
- **DROP**: Drop samples missing from any modality

After concatenation, columns are Z-normalized if `renormalize=True`.

### Weighted Average Strategy

For shared pathways only, computes a weighted mean:

```
fused_p = Σ_m (w_m × score_p,m)    where Σ w_m = 1
```

Weights default to uniform (1/n_modalities). Produces a matrix with the same dimensionality as a single modality.

### Intersection Strategy

Restricts to samples and pathways present in all modalities. No imputation needed. Conservative but avoids missing data assumptions.

### Quality Metrics

The fusion result includes a quality report with:
- Sample overlap fraction (shared / total)
- Pathway overlap count
- Warnings for low overlap (<50% samples shared)

---

## Clustering Methodology

### Gaussian Mixture Model (Default)

The primary clustering method uses Gaussian Mixture Models (GMM):

```python
GMM(n_components=k, covariance_type='full', n_init=10, reg_covar=1e-6)
```

Key parameters:
- **n_components**: Number of clusters (determined by model selection)
- **covariance_type**: Full covariance matrices (captures correlations)
- **n_init**: Multiple random initializations for robustness
- **reg_covar**: Regularization for numerical stability

### Alternative Algorithms

For validation, multiple clustering algorithms are supported:

| Algorithm | Characteristics | Best For |
|-----------|-----------------|----------|
| GMM | Soft assignments, BIC selection | Primary analysis |
| K-means | Fast, spherical clusters | Large datasets |
| Hierarchical | Dendogram, no K required | Exploratory |
| Spectral | Nonlinear boundaries | Complex structure |

### Model Selection (Number of Clusters)

The optimal number of clusters is determined by:

**BIC Method** (default):
```
BIC = -2 × log(L) + k × log(n)
```
Select k that minimizes BIC.

**Silhouette Method** (alternative):
```
silhouette = (b - a) / max(a, b)
```
Select k that maximizes mean silhouette score.

## Validation Framework

### Negative Control 1: Label Shuffle Test

Tests whether discovered subtypes are spurious:
1. Shuffle sample labels randomly
2. Recompute pathway scores
3. Recluster and compute ARI vs original labels
4. **Pass criterion**: ARI < 0.15

### Negative Control 2: Random Gene Sets Test

Tests whether pathway structure is meaningful:
1. Replace curated pathways with random gene sets (same sizes)
2. Recompute pathway scores with random genes
3. Recluster and compute ARI vs original labels
4. **Pass criterion**: ARI < 0.15

### Stability Test: Bootstrap Resampling

Tests cluster stability across samples:
1. Generate 100 bootstrap samples (with replacement)
2. Recluster each bootstrap sample
3. Compare to original clustering via ARI
4. **Pass criterion**: Mean ARI ≥ 0.80

### Negative Control 3: Ancestry Independence Test

Tests whether discovered subtypes are confounded by population structure:
1. Compute Kruskal-Wallis H-statistic for each ancestry PC across cluster labels
2. Apply Bonferroni correction for multiple PCs tested
3. **Pass criterion**: No PC significantly associated with clusters after correction

This gate only runs when ancestry principal components are provided. See [Ancestry Correction](#ancestry--population-stratification-correction) below.

### Validation Threshold Calibration

The default validation thresholds (0.15 for null ARI, 0.80 for stability) are fixed values that do not account for sample size or number of clusters. The threshold calibration module replaces these with data-driven values.

#### Motivation

- Under the null hypothesis, ARI distributions narrow as sample size increases (tighter convergence to zero)
- Chance-level ARI increases with more clusters (more possible agreements by chance)
- A single threshold is therefore either too permissive for small samples or too strict for large samples

#### Pre-Computed Lookup Tables

Thresholds are pre-computed via simulation across a grid of configurations:

```
n_samples ∈ {30, 50, 75, 100, 150, 200, 300, 500}
n_clusters ∈ {2, 3, 4, 5, 6, 7, 8}
```

For each (n, k) pair, 500 simulations are run:

**Null ARI table** (95th percentile):
1. Generate random data (n samples, 15 pathways, no structure)
2. Fit GMM with k clusters
3. Compute ARI between random labels and GMM labels
4. Record 95th percentile across 500 simulations

**Stability table** (5th percentile):
1. Generate structured data with k planted subtypes (effect size = 1.5)
2. Bootstrap resample and re-cluster
3. Compute ARI between bootstrap and original labels
4. Record 5th percentile across 500 simulations

#### Interpolation

For (n_samples, n_clusters) values between grid points, bilinear interpolation is used:

```
threshold(n, k) = Σᵢ Σⱼ wᵢⱼ × table(nᵢ, kⱼ)
```

Where weights wᵢⱼ are inverse-distance based on the four nearest grid points.

#### Fallback Simulation

For configurations outside the table range, fresh simulations are run on-the-fly using the same methodology.

### Cross-Validation

K-fold cross-validation for clustering stability:
1. Split data into K folds
2. Train clustering on K-1 folds
3. Predict labels for held-out fold
4. Compare held-out predictions to training labels
5. Report mean and standard deviation of fold ARIs

## Statistical Corrections

### Multiple Testing Correction

All pathway-level comparisons use Benjamini-Hochberg FDR correction:

```
p_adjusted_i = min(1, p_i × m / rank_i)
```

Where m = total tests, rank = ascending p-value rank.

### Permutation-Based P-Values

For effect size significance:
```
p_perm = (# permutations with |effect| ≥ |observed|) / n_permutations
```

Using ≥1000 permutations for robust estimation.

### Confidence Intervals

95% bootstrap confidence intervals for effect sizes:
1. Generate 1000 bootstrap samples
2. Compute effect size for each
3. Report 2.5th and 97.5th percentiles

## Ancestry / Population Stratification Correction

Population stratification—systematic differences in allele frequencies across ancestral groups—can confound pathway-based clustering, producing subtypes that reflect ancestry rather than disease biology. The framework provides optional ancestry correction to mitigate this.

### Ancestry PCA

Principal Component Analysis on the genotype matrix captures population structure:

```
PC_matrix = PCA(standardize(genotype_matrix), n_components=k)
```

Where:
- Genotype matrix is per-variant standardized (zero mean, unit variance)
- Monomorphic variants (zero variance) are excluded
- Default `k = 10` components (adjustable)
- Explained variance per PC is reported for diagnostics

Reference: Price AL et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nature Genetics*.

### Regression-Based Correction (Default)

Ancestry-correlated variance is removed from pathway scores via OLS residualization:

```
score_adjusted_p = score_p - X_ancestry × β_p
```

For each pathway `p`:
1. Fit linear regression: `score_p ~ PC_1 + PC_2 + ... + PC_k`
2. Compute R² (proportion of variance explained by ancestry)
3. Replace pathway score with residuals
4. Re-normalize residuals to zero mean, unit variance

### Confounding Detection

Pathways with high ancestry-explained variance are flagged:

```
confounded = {p : R²_p > threshold}
```

Default threshold: R² > 0.1 (10% of variance explained by ancestry PCs). Flagged pathways remain in the analysis but are reported for transparency.

### Ancestry Independence Testing

Post-clustering, the framework tests whether subtypes are independent of ancestry:

```
H_j = KruskalWallis(PC_j ~ cluster_labels)    for j = 1..k
p_corrected_j = p_j × k                        (Bonferroni)
```

- **Kruskal-Wallis H-test**: Non-parametric test for differences in PC distributions across clusters
- **Bonferroni correction**: Controls family-wise error rate across multiple PCs
- **Pass criterion**: No corrected p-value < significance threshold (default 0.05)

### Stratified Analysis

For datasets with known ancestry groups, per-group clustering provides additional validation:

1. Cluster each ancestry group independently using GMM
2. Compute cross-group centroid concordance via cosine similarity
3. High concordance (> 0.7) indicates subtypes are consistent across populations

### Correction Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| REGRESS_OUT | OLS residualization of ancestry PCs | Default; removes linear ancestry effects |
| COVARIATE_AWARE | Covariate-adjusted residualization | Future: mixed-model extensions |
| STRATIFIED | Per-ancestry-group analysis | When groups are known and large enough |

## Batch Effect Correction

### Batch Effect Detection

Batch effects (e.g., sequencing site, sample processing date) are detected via one-way ANOVA per pathway:

```
F_p = MS_between / MS_within
η²_p = SS_between / SS_total
```

Where η² (eta-squared) measures the proportion of variance explained by batch. Features with FDR-corrected p < 0.05 are flagged as batch-affected.

### ComBat Correction (Default)

The framework implements ComBat empirical Bayes batch correction (Johnson et al., 2007):

1. Estimate batch-specific location (γ) and scale (δ) parameters
2. Shrink estimates toward pooled priors using empirical Bayes
3. Adjust scores: `score_adjusted = (score - γ_batch) / δ_batch × δ_pooled + μ_grand`

ComBat is preferred for multi-site genomic studies as it borrows strength across features.

### Alternative Correction Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| COMBAT | Empirical Bayes location/scale | Default; best for genomic data |
| MEAN_CENTER | Per-batch mean subtraction | Simple; when only location differs |
| STANDARDIZE | Per-batch z-score then re-standardize | When both location and scale differ |

### Post-Correction Validation

After correction, the framework validates:
- **Batch variance reduction**: Mean batch-explained variance should decrease
- **Biological signal preservation**: Variance explained by known biological groups should be maintained
- **Correlation with original**: Corrected scores should remain correlated with originals (signal not destroyed)

## Sensitivity Analysis

### Parameter Robustness Testing

Sensitivity analysis evaluates how stable clustering results are under perturbation of key parameters:

```
stability = mean(ARI(labels_config_i, labels_config_j))  for all pairs i,j
```

### Parameter Axes

| Parameter | Variations | Purpose |
|-----------|------------|---------|
| Clustering algorithm | GMM, K-means, Hierarchical | Algorithm dependence |
| Number of clusters | k = 2..max_k | K sensitivity |
| Feature subset | Leave-one-out per pathway | Feature importance |
| Normalization | Z-score, min-max, robust, rank | Pre-processing sensitivity |

### Robustness Scoring

Overall stability is the mean of per-parameter mean ARIs:

```
overall_stability = mean(mean_ARI_parameter)  for all parameters
```

A result is considered **robust** if overall_stability exceeds a user-defined threshold (default: 0.7).

### Most/Least Sensitive Parameter

The parameter with the lowest mean ARI is reported as most sensitive (results change most when this parameter varies). The parameter with the highest mean ARI is least sensitive.

Reference: Hennig C (2007). Cluster-wise assessment of cluster stability. *Computational Statistics & Data Analysis*.

## Effect Size Metrics

### Cohen's d

For continuous pathway scores between clusters:
```
d = (μ_1 - μ_2) / σ_pooled
```

Interpretation:
- |d| < 0.2: Negligible
- 0.2 ≤ |d| < 0.5: Small
- 0.5 ≤ |d| < 0.8: Medium
- |d| ≥ 0.8: Large

### Adjusted Rand Index (ARI)

For comparing cluster assignments:
```
ARI = (RI - E[RI]) / (max(RI) - E[RI])
```

Range: -1 to 1 (1 = perfect agreement, 0 = random)

### Normalized Mutual Information (NMI)

Alternative clustering comparison metric:
```
NMI = 2 × I(U,V) / (H(U) + H(V))
```

Range: 0 to 1 (1 = perfect agreement)

## Power Analysis

### Type I Error Estimation

Estimate false positive rate under null:
1. Generate random data with no true structure
2. Apply clustering pipeline
3. Measure rate of ARI > threshold

Reported as Type I error rate at threshold = 0.15.

### Power Curves

Power as a function of effect size:
1. Generate data with planted subtypes at various effect sizes
2. Apply clustering pipeline
3. Measure recovery rate (ARI ≥ 0.80)
4. Report power (proportion successful) at each effect size

### Sample Size Recommendations

Based on power analysis:
- **Minimum**: n = 50 per subtype (effect size ≥ 1.5)
- **Recommended**: n = 100 per subtype (effect size ≥ 1.0)
- **Optimal**: n = 200 per subtype (effect size ≥ 0.5)

## Reproducibility

### Random Seed Control

All stochastic operations use controlled random seeds:
```python
np.random.seed(config.seed)
sklearn_model(..., random_state=config.seed)
```

### Metadata Logging

Each analysis logs:
- Software versions (package, Python, sklearn)
- Configuration parameters
- Input file checksums
- Timestamp and execution environment

## Simulation Validation

### Synthetic Data Generation

For validation, synthetic data with known structure is generated:
1. Define n_subtypes with specified proportions
2. Assign subset of pathways to each subtype
3. Elevate pathway scores by effect_size × noise_level
4. Add Gaussian noise
5. Optionally simulate ancestry groups with configurable confounding between ancestry and subtypes

### Ground Truth Recovery

Clustering performance is evaluated against planted labels:
- **ARI**: Primary metric (accounts for label permutation)
- **NMI**: Secondary metric (information-theoretic)
- **Correct K**: Whether optimal K matches true K

## Reporting Standards

### Required Outputs

All analyses report:
1. Number of samples and pathways
2. Optimal K and selection method
3. Per-cluster sample sizes
4. Top contributing pathways per cluster (FDR < 0.05)
5. Validation gate results (pass/fail)
6. Effect sizes with 95% CIs

### Visualization Standards

Recommended plots:
- Heatmap of pathway scores by cluster
- PCA/UMAP projection colored by cluster
- BIC/silhouette curves for K selection
- Bootstrap ARI distribution

## Data Provenance

### External Reference Data

The framework references the following publicly available databases for variant annotation and filtering. **No data from these databases is stored in this repository** — they are referenced at runtime only when users supply their own data.

| Database | Use in Framework | Access | License |
|----------|-----------------|--------|---------|
| gnomAD | MAF filtering threshold (< 1%) | Open (https://gnomad.broadinstitute.org/) | ODC-ODbL |
| CADD | Variant pathogenicity scoring | Open (https://cadd.gs.washington.edu/) | Non-commercial |
| ClinVar | Validation of pathway gene lists | Open (https://www.ncbi.nlm.nih.gov/clinvar/) | Public domain |
| Reactome | Cross-reference for pathway validation | Open (https://reactome.org/) | CC BY 4.0 |

### Data Shipped with This Repository

All data files in this repository are either **computationally generated** (synthetic VCF, phenotypes, test fixtures) or **curated from public literature** (pathway GMT gene lists). No real patient data, proprietary data, or commercially licensed data is included. See [DISCLAIMER.md](../DISCLAIMER.md) for the complete data provenance statement.

### Pathway Gene Lists

The GMT files in `data/pathways/` contain gene symbols curated from:
- **SFARI Gene** (https://gene.sfari.org/) — autism gene scoring, open access
- **KEGG** (https://www.kegg.jp/) — curated biological pathways
- **Reactome** (https://reactome.org/) — peer-reviewed pathway database
- **MSigDB** (https://www.gsea-msigdb.org/) — gene set collections
- **Gene Ontology** (http://geneontology.org/) — standardized gene annotations

Gene symbols (e.g., SHANK3, CHD8, NRXN1) are standard HGNC identifiers used in thousands of published research papers and are not proprietary information.

## References

- Anderson CA et al. (2010). Data quality control in genetic case-control association studies. Nat Protoc.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. J R Stat Soc B.
- Cohen J (1988). Statistical Power Analysis for the Behavioral Sciences.
- Hubert L, Arabie P (1985). Comparing partitions. J Classif.
- McLachlan GJ, Peel D (2000). Finite Mixture Models.
- Patterson N, Price AL, Reich D (2006). Population structure and eigenanalysis. PLoS Genet.
- Price AL, Patterson NJ, Plenge RM, et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. Nat Genet.
- Hennig C (2007). Cluster-wise assessment of cluster stability. Comput Stat Data Anal.
- Johnson WE, Li C, Rabinovic A (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics.
- Leek JT, et al. (2010). Tackling the widespread and critical impact of batch effects in high-throughput data. Nat Rev Genet.
- Strehl A, Ghosh J (2002). Cluster ensembles. J Mach Learn Res.
- Landrum MJ et al. (2020). ClinVar: improvements to accessing data. Nucleic Acids Res.
- Gillespie M et al. (2022). The Reactome Pathway Knowledgebase 2022. Nucleic Acids Res.
- Karczewski KJ et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. Nature.
- Wigginton JE, Cutler DJ, Gravel A (2005). A note on exact tests of Hardy-Weinberg equilibrium. Am J Hum Genet.
- Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. Nat Commun.

---

*Document version: 0.2.4-dev | Last updated: February 2026*
