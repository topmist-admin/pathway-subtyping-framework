# Changelog

All notable changes to the Pathway Subtyping Framework will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

#### Notebooks
- **Notebook 13** (`13_geo_blood_ados.ipynb`): GSE111175 blood transcriptomics with ADOS clinical correlation
  - 141 blood leukocyte samples (28 ASD, 113 control), Illumina BeadChip
  - Pathway subtypes detectable in peripheral tissue
  - Synaptic transmission correlates with ADOS Social Affect (rho=-0.52, FDR p=0.032)
- **Notebook 14** (`14_geo_blood_large_cohort.ipynb`): GSE18123 large blood cohort validation
  - 285 blood samples (72 ASD, 213 control) across two Affymetrix platforms
  - Cross-cohort projection from GSE111175 achieves ARI=0.374 (exceeds 0.3 threshold)
  - Cross-tissue pathway correlation with brain subtypes (rho=0.371)
- **Notebook 15** (`15_geo_scz_replication.ipynb`): GSE53987 SCZ replication on Affymetrix microarray
  - 205 samples across 4 diagnoses (SCZ, BD, MDD, CTL), 3 brain regions
  - Cross-platform projection from GSE80655 RNA-seq: ARI=0.319 (PASS)
  - Cross-disease ARI=0.792 (SCZ/ASD pathway convergence)
  - Multi-diagnosis pooled clustering: k=5, silhouette=0.450
- **Notebook 16** (`16_knowledge_graph_analysis.ipynb`): Cross-disease knowledge graph & drug repurposing
  - 336 genes from 21 unique pathways (15 ASD + 14 SCZ, 8 shared)
  - STRING PPI network: 4,378 edges, 97.9% gene coverage
  - 20 hub genes (11 cross-disease bridges) ranked by betweenness centrality
  - DGIdb drug repurposing: 1,546 unique drugs targeting 44 hub/core genes
  - 6 Louvain communities (all cross-disease), convergence subnetwork (260 genes)
  - Manuscript-ready network figure, drug overlay, pathway crosstalk heatmap

#### Scripts
- **Cytoscape figure generator** (`scripts/generate_cytoscape_figures.py`): Publication-ready network figures via py4cytoscape
  - Hub gene subnetwork (20 hubs + top neighbors, community-colored)
  - Drug-hub target map (bipartite gene-drug network with interaction types)
  - Community overview (336 nodes with community-aware layout)
  - Requires `[graph]` extra + Cytoscape desktop app

#### Dependencies
- New `[graph]` optional extra: `networkx>=3.0`, `py4cytoscape>=1.0.0`
- Updated `requirements.txt` to match `pyproject.toml` (removed stale `pysam` from core)

#### Documentation
- Notebook execution guide with dependency diagram (`docs/notebook-guide.md`)
- Notebook execution registry (`research-results/NOTEBOOK-EXECUTION-REGISTRY.md`)
- Updated JOSS paper with all 6 validation datasets (1,075 samples)

## [0.3.0] - 2026-02-16

### Added

#### Signaling Pathway Databases (#32)
- **Signaling databases module** (`signaling_databases.py`): Load cell-cell signaling databases as pathway gene sets
  - `SignalingDatabase`: Enum for CellPhoneDB and CellChatDB
  - `SignalingInteraction`: Ligand-receptor interaction dataclass with gene resolution
  - `SignalingDatabaseResult`: Result with pathway gene sets, citations (`.to_dict()`, `.format_report()`, `.get_citations()`)
  - `load_cellphonedb()`: Auto-download from GitHub, resolve complexes to genes, group by signaling classification
  - `load_cellchatdb()`: Load user-exported CellChatDB CSV (R binary not Python-readable)
  - `convert_interactions_to_pathways()`: Group interactions into pathway gene sets
  - `merge_signaling_databases()`: Merge gene sets from multiple databases
- Signaling pathway gene sets directly compatible with `score_pathways_from_expression()`, `run_clustering()`, and validation gates

#### Bulk Deconvolution Integration (#31)
- **Deconvolution module** (`deconvolution.py`): Estimate cell-type proportions from bulk RNA-seq
  - `DeconvolutionMethod`: Enum for NNLS deconvolution
  - `DeconvolutionQualityReport`: Reference coverage, proportion validation, gene overlap
  - `DeconvolutionResult`: Cell-type proportion matrix with `.to_dict()`, `.format_report()`, `.get_citations()`
  - `build_reference_profile()`: Aggregate single-cell reference to cell-type mean expression profiles
  - `deconvolve_bulk()`: Main entry point — NNLS deconvolution with quality checks
  - `combine_features()`: Merge pathway scores + cell-type proportions into unified feature matrix
  - `generate_synthetic_bulk()`: Create synthetic bulk from known proportions for testing
- **Multi-omic integration**: Added `ModalityType.DECONVOLUTION` for deconvolution-derived features

#### Cross-Modal Validation Gate (#33)
- **Cross-modal validation module** (`cross_modal_validation.py`): Gate 5 — tests whether subtypes replicate across data modalities
  - `CrossModalPairResult`: Per-pair concordance metrics (ARI, NMI, bidirectional transfer ARI)
  - `CrossModalValidationResult`: Top-level result with gate pass/fail, null calibration, format_report/get_citations
  - `SingleCellCompositionResult`: ANOVA-based test for distinct cell-type compositions across subtypes
  - `cross_modal_concordance()`: Main entry point — independent clustering per modality, pairwise ARI/NMI, transfer validation, permutation null calibration
  - `single_cell_composition_test()`: One-way ANOVA per cell type with Bonferroni correction
  - `generate_synthetic_multimodal_data()`: Planted shared subtype structure for testing and calibration
- **Validation integration**: Gate 5 added to `ValidationGates.run_all()` when `per_modality_scores` is provided (>= 2 modalities)
- **Pipeline integration**: `run_validation_gates()` passes `per_modality_scores` from multi-omic fusion result

#### Multi-Omic Pipeline Integration (#35)
- **Multi-omic fusion module** (`multi_omic.py`): Fuse pathway scores from VCF, bulk RNA-seq, and scRNA-seq modalities
  - `ModalityType`: Enum for VCF, expression, and single-cell modalities
  - `FusionStrategy`: Enum for concatenate, weighted_average, intersection_only fusion strategies
  - `MissingStrategy`: Enum for handling samples missing from some modalities (impute_zero, impute_mean, drop)
  - `ModalityInput`: Dataclass wrapping a single modality's scored output with validation
  - `SampleOverlapStats`: Statistics about sample overlap across modalities
  - `MultiOmicFusionResult`: Result container with fused scores, per-modality scores, overlap stats (`.to_dict()`, `.format_report()`, `.get_citations()`)
  - `MultiOmicQualityReport`: Per-modality quality reports, sample/pathway overlap statistics, warnings
  - `prepare_modality()`: Validate and wrap modality output for fusion
  - `fuse_modalities()`: Main entry point — concatenation (default), weighted average, or intersection-only fusion
  - `compute_sample_overlap()`: Analyze sample overlap across modalities with pairwise statistics
  - `correlation_analysis()`: Cross-modality pathway correlation (Pearson/Spearman)
- **Pipeline integration**: `input_type: "multi_omic"` in YAML config with `multi_omic` section for modality list, fusion strategy, weights
- **Config validation**: New `multi_omic` section validation with modality-type, path, weight, and strategy checks

#### Single-Cell Pathway Scoring (#30)
- **Single-cell module** (`single_cell.py`): Per-cell and pseudobulk pathway scoring from scRNA-seq data
  - `SingleCellScoringMethod`: Enum for per-cell mean_z and pseudobulk (mean_z, ssGSEA, GSVA) methods
  - `SingleCellInputType`: Enum for raw counts, log-normalized, and h5ad input types
  - `SingleCellQualityReport`: QC report with sparsity, cell counts, gene detection, pathway coverage
  - `SingleCellScoringResult`: Result container with cell-type-level and optional per-cell scores (`.to_dict()`, `.format_report()`, `.get_citations()`)
  - `load_single_cell_data()`: Load and validate h5ad or CSV/TSV single-cell data with auto-normalization
  - `score_single_cell_pathways()`: Main scoring entry point; pseudobulk reuses expression.py internals
  - Memory-efficient chunked per-cell scoring for datasets up to 50K cells
  - Sparse matrix support (scipy CSR/CSC) throughout
- **Optional `[sc]` dependency group**: `anndata>=0.9.0`

#### Advanced Visualization (#9)
- **Visualization module** (`visualization.py`): Interactive and publication-quality visualizations
  - `DimReductionMethod`: Enum for PCA, t-SNE, UMAP dimensionality reduction
  - `FigureFormat`: Enum for PNG, SVG, PDF, HTML export formats
  - `ReportConfig`: Configuration dataclass for interactive report generation
  - `VisualizationResult`: Result dataclass with output paths and metadata (with `.to_dict()`)
  - `compute_dim_reduction()`: PCA, t-SNE, or UMAP embedding with metadata (explained variance, KL divergence, etc.)
  - `plot_interactive_scatter()`: Plotly scatter plot of samples in reduced space (PCA/t-SNE/UMAP)
  - `plot_interactive_heatmap()`: Plotly heatmap of mean pathway Z-scores per subtype
  - `plot_cluster_distribution()`: Plotly bar chart of cluster sizes
  - `plot_subtype_trajectories()`: Plotly radar chart of subtype pathway profiles
  - `plot_static_scatter()`: Matplotlib fallback scatter plot (works without Plotly)
  - `export_figure()`: Multi-format export (PNG, SVG, PDF, HTML) for both Plotly and matplotlib figures
  - `create_interactive_report()`: Self-contained HTML report with all charts, summary statistics, and styling
  - `generate_all_figures()`: Convenience function to generate all visualizations at once
- **Optional `[viz]` dependency group**: `plotly>=5.15.0`, `umap-learn>=0.5.0`, `kaleido>=0.2.0`
- **Pipeline integration**: `generate_interactive_report` and `interactive_dim_reduction` config options in `PipelineConfig`
  - YAML config: `output.generate_interactive_report: true` and `output.interactive_dim_reduction: "umap"`
  - Graceful fallback when Plotly not installed (warning + skip)
- 53 new tests covering all visualization functions, dimensionality reduction, export formats, edge cases
- All interactive features degrade gracefully to static matplotlib when `[viz]` extra not installed

#### Pipeline Input Validation
- **Phenotype file validation**: Checks for empty files, missing `sample_id` column, duplicate sample IDs
- **Minimum sample size guard**: Error if samples < 2× max clusters, warning if < 5× max clusters
- **Sample overlap reporting**: Logs overlap between data samples and phenotype samples
- 6 new tests for pipeline input validation

#### CI and Code Quality
- Fixed all black/isort/flake8 lint issues across 16 files
- Added `-m "not network"` to CI test commands to skip flaky network tests
- Removed unused imports in expression.py, threshold_calibration.py, validation_datasets.py
- Closed GitHub issues #7 (cross-cohort), #8 (performance), #26 (variant QC)
- Created `py.typed` PEP 561 marker file

#### CI and Testing
- Added `requires_plotly` and `requires_umap` skip markers in `test_visualization.py` so visualization tests skip gracefully when `[viz]` extras are not installed (fixes CI failures)
- 49 new single-cell tests, 768 total tests passing

### Changed
- Total test count: 971 (up from 716)

## [0.2.3] - 2026-02-14

### Added

#### Cross-Cohort Validation Enhancements (#7)
- **CohortResult.to_dict()**: Serialize cohort results to dictionary
- **CrossCohortResult methods**: Added `.to_dict()`, `.format_report()`, `.get_citations()` for publication-ready output
- **`generate_synthetic_cohort_pair()`**: Convenience function for creating matched synthetic cohort pairs for testing and demos
- **Cross-cohort example script** (`scripts/cross_cohort_example.py`): CLI example with argparse
- **Cross-cohort user guide** (`docs/guides/cross-cohort-validation.md`): Interpretation tables, real-world workflow, common pitfalls
- **Cross-cohort API reference** (`docs/api/cross_cohort.md`): Full API documentation
- 12 new tests for cross-cohort enhancements

#### Performance Optimization (#8)
- **tqdm progress bars**: Added to validation gates (label shuffle, random gene sets, bootstrap stability), simulation analysis (Type I error, power analysis, sample size analysis), expression scoring (ssGSEA, GSVA), and sensitivity analysis (feature LOO)
- **`show_progress` parameter**: Added to `ValidationGates`, `estimate_type_i_error()`, `run_power_analysis()`, `run_sample_size_analysis()`, `score_pathways_from_expression()`, `vary_feature_subset()`, `run_sensitivity_analysis()`
- **Chunked processing**: `PipelineConfig` now accepts `use_chunked_processing` and `chunk_size` options, delegating to existing `compute_gene_burdens_chunked()` for large VCF files
- **Benchmark script** (`scripts/benchmark_performance.py`): Measures wall-clock time and peak memory via `tracemalloc`, reports pass/fail against targets (30 min, 8 GB for 10K samples)
- **Hardware guide** (`docs/guides/performance-and-hardware.md`): Recommended hardware table, memory estimation, chunked processing guide, Colab constraints, performance tips
- 9 new tests for performance parameters

#### Variant Quality Control (#26)
- **Variant QC module** (`variant_qc.py`): Standard genetic variant quality control filters applied before burden computation
  - `VariantQCConfig`: Dataclass with min_qual, min_call_rate, hwe_p_threshold, max_maf, min_gq, min_dp (with `.to_dict()`)
  - `VariantQCResult`: Dataclass with removal counts, per-variant metrics, retention rate (with `.to_dict()`, `.format_report()`, `.get_citations()`)
  - `compute_call_rate()`: Per-variant genotype call rate
  - `compute_maf()`: Minor allele frequency computation for diploid samples
  - `check_hwe()`: Hardy-Weinberg equilibrium chi-squared test per variant
  - `apply_genotype_filters()`: Per-genotype GQ/DP masking
  - `filter_variants()`: Applies QUAL, call rate, HWE, MAF filters in sequence with per-filter removal tracking
- **Pipeline integration**: `variant_qc_enabled` config option runs QC between data loading and burden computation
  - `PipelineConfig`: Added `variant_qc_enabled`, `variant_qc_min_qual`, `variant_qc_min_call_rate`, `variant_qc_hwe_p_threshold`, `variant_qc_max_maf` fields
  - `from_yaml()`: Now parses `variant_qc:` section
  - QC results included in JSON and Markdown reports
- **Config validation** (`config.py`): `_validate_variant_qc_section()` validates all QC parameter ranges
- 40 new tests covering all QC functions, config validation, and package imports

#### Data-Driven Validation Threshold Calibration (#19)
- **Threshold calibration module** (`threshold_calibration.py`): Replaces hard-coded validation thresholds with data-driven values that adjust for sample size and number of clusters
  - `CalibratedThresholds`: Dataclass with null ARI threshold, stability threshold, calibration method, interpolation flag (with `.to_dict()`, `.format_report()`, `.get_citations()`)
  - `CalibrationSimulationResult`: Dataclass for simulation distributions
  - `calibrate_thresholds(n_samples, n_clusters, ...)`: Auto-calibrate thresholds via lookup table, interpolation, or simulation fallback
  - `get_default_thresholds()`: Returns legacy 0.15/0.8 values for backward compatibility
  - `generate_calibration_table()`: Regenerate lookup tables from simulations
  - Pre-computed lookup tables: 56-entry grid (8 sample sizes × 7 cluster counts) with empirically-derived 95th percentile null ARI and 5th percentile stability ARI
  - Bilinear interpolation for intermediate configurations
  - On-the-fly simulation fallback for out-of-range configurations
- **Pipeline integration**: Auto-calibration in `run_validation_gates()` when thresholds are `null` in config
  - `PipelineConfig`: Added `validation_calibrate`, `validation_stability_threshold`, `validation_null_ari_max`, `validation_alpha`, `validation_n_permutations`, `validation_n_bootstrap` fields
  - `from_yaml()`: Now parses `validation:` section (previously ignored)
  - Calibration info included in JSON and Markdown reports
- **Config validation** (`config.py`): `_validate_validation_section()` validates threshold ranges, alpha, iteration counts
- **Threshold calibration tests** (`tests/test_threshold_calibration.py`): 46 tests covering lookup tables, interpolation, simulation, calibration modes, reproducibility
- **Table generation script** (`scripts/generate_calibration_table.py`): CLI script to regenerate lookup tables

### Fixed

#### ClinVar and Reactome Parser Updates
- **ClinVar parser** (`validation_datasets.py`): Handle NCBI's updated `gene_specific_summary.txt` column format
  - New format uses `Alleles_reported_Pathogenic_Likely_pathogenic` (combined column) instead of separate `Number_Pathogenic`/`Number_Likely_Pathogenic` columns
  - Parser auto-detects format and handles both old and new column names
  - Handles `Number_uncertain` (new) vs `Number_Uncertain_Significance` (old) column naming
- **Reactome parser** (`validation_datasets.py`): Handle Reactome's updated GMT layout
  - New format: `Pathway Name\tR-HSA-ID\tGenes` (R-HSA-ID moved to column 1)
  - Old format: `R-HSA-ID\tHomo sapiens: Description\tGenes`
  - Parser now checks species prefix in description field alongside existing name checks

## [0.2.0] - 2026-02-09

### Added

#### Ancestry / Population Stratification Correction
- **Ancestry correction module** (`ancestry.py`): PCA-based population stratification detection and correction
  - `AncestryMethod`: Enum for correction approaches (REGRESS_OUT, COVARIATE_AWARE, STRATIFIED)
  - `compute_ancestry_pcs()`: PCA on genotype matrix with monomorphic variant handling
  - `adjust_pathway_scores()`: Regression residualization to remove ancestry-correlated variance
  - `check_ancestry_independence()`: Kruskal-Wallis test with Bonferroni correction for cluster-ancestry independence
  - `stratified_analysis()`: Per-ancestry-group clustering with cross-group concordance
  - `compute_ancestry_correlation()`: Pearson correlation matrix between pathways and ancestry PCs
  - Dataclasses: `AncestryPCs`, `AncestryAdjustmentResult`, `AncestryStratificationReport` (all with `.to_dict()`, `.format_report()`, `.get_citations()`)
- **Ancestry validation gate** (`validation.py`): 4th validation gate tests cluster-ancestry independence
- **Ancestry simulation** (`simulation.py`): Configurable ancestry confounding for synthetic data
  - `n_ancestry_groups`, `ancestry_effect_size`, `ancestry_confounding` parameters in `SimulationConfig`
  - Simulated ancestry PCs and group labels in `SimulatedData`
- **Pipeline integration** (`pipeline.py`): Optional ancestry correction between pathway scoring and clustering
  - `ancestry_pcs_path`, `ancestry_correction`, `ancestry_n_pcs` in `PipelineConfig`
  - Ancestry section in JSON and Markdown reports
- **Config validation** (`config.py`): Ancestry section validation (method, PCs file, n_pcs)
- **Ancestry test suite** (`tests/test_ancestry.py`): 44 tests covering all functions, edge cases, and end-to-end validation

#### Batch Correction & Sensitivity Analysis
- **Batch correction module** (`batch_correction.py`): ComBat-style batch effect detection and correction
  - `BatchCorrectionMethod`: Enum for correction approaches (COMBAT, MEAN_CENTER, STANDARDIZE)
  - `detect_batch_effects()`: ANOVA-based batch effect detection with eta-squared variance explained
  - `correct_batch_effects()`: ComBat empirical Bayes correction, mean centering, or standardization
  - `validate_batch_correction()`: Post-correction validation of variance reduction and signal preservation
  - Dataclasses: `BatchEffectReport`, `BatchCorrectionResult` (all with `.to_dict()`, `.format_report()`, `.get_citations()`)
- **Sensitivity analysis module** (`sensitivity.py`): Systematic parameter robustness testing
  - `SensitivityParameter`: Enum for parameter axes (CLUSTERING_ALGORITHM, N_CLUSTERS, NORMALIZATION, FEATURE_SUBSET)
  - `vary_clustering_algorithm()`: Compare GMM, K-means, Hierarchical across algorithms
  - `vary_n_clusters()`: Sweep cluster count range with pairwise ARI
  - `vary_feature_subset()`: Leave-one-out pathway sensitivity
  - `vary_normalization()`: Compare z-score, min-max, robust, rank normalization
  - `run_sensitivity_analysis()`: Full sensitivity analysis with robustness scoring
  - Dataclasses: `ParameterVariationResult`, `SensitivityAnalysisResult` (all with `.to_dict()`, `.format_report()`, `.get_citations()`)
- **Batch correction tests** (`tests/test_batch_correction.py`): 34 tests covering detection, correction, validation, edge cases
- **Sensitivity analysis tests** (`tests/test_sensitivity.py`): 27 tests covering all parameter axes, reproducibility, dataclasses

#### Scientific Rigor Modules (Publication Readiness)
- **Statistical rigor module** (`statistical_rigor.py`): Publication-quality statistics
  - `benjamini_hochberg()`: FDR correction for multiple testing
  - `BurdenWeightScheme`: Literature-based variant weighting (DEFAULT, GNOMAD_CONSTRAINT, ACMG_INSPIRED, UNIFORM)
  - `PathwayNormalization`: Multiple aggregation methods (MEAN, MEDIAN, SIZE_NORMALIZED, PCA)
  - `compute_pathway_effect_sizes()`: Cohen's d with bootstrap confidence intervals
  - `compute_pathway_pvalues()`: Permutation-based p-value computation
  - `run_statistical_analysis()`: Comprehensive statistical analysis pipeline
- **Multiple clustering algorithms** (`clustering.py`): Algorithm comparison framework
  - `ClusteringAlgorithm`: GMM, K-means, Hierarchical, Spectral clustering
  - `run_clustering()`: Unified interface for all algorithms
  - `select_n_clusters()`: BIC or silhouette-based model selection
  - `cross_validate_clustering()`: K-fold cross-validation for stability
  - `compare_algorithms()`: Pairwise ARI comparison, consensus labels
- **Simulation framework** (`simulation.py`): Ground truth validation
  - `SimulationConfig`: Configurable synthetic data generation
  - `generate_synthetic_data()`: Planted subtype structure with effect size control
  - `estimate_type_i_error()`: False positive rate estimation under null
  - `run_power_analysis()`: Power curves across effect sizes
  - `run_sample_size_analysis()`: Sample size recommendations for target power
  - `validate_framework()`: Comprehensive framework validation
- **Formal methods documentation** (`docs/METHODS.md`): Statistical methodology for publications

#### Additional Disease Pathways (Week 5)
- **Parkinson's Disease pathways** (`parkinsons_pathways.gmt`): 14 pathways, ~280 genes
  - Alpha-synuclein aggregation, mitochondrial function, autophagy-lysosomal pathway
  - Dopamine metabolism, endolysosomal trafficking, immune/inflammation
  - Sources: Nalls et al. 2019 (Lancet Neurol), Blauwendraat et al. 2020, IPDGC
- **Bipolar Disorder pathways** (`bipolar_pathways.gmt`): 14 pathways, ~290 genes
  - Calcium signaling, circadian rhythm, WNT/GSK3 signaling
  - Glutamate/GABA signaling, HPA stress response, neuroplasticity
  - Sources: Mullins et al. 2021 (Nat Genet), Stahl et al. 2019, BDgene
- **Literature citations** added to autism pathway file header
- Updated pathway documentation with new disease recommendations

#### Real-World Data Support (Week 4)
- **Multi-allelic variant support**: Automatically expands multi-allelic variants (e.g., A→G,T) into separate bi-allelic records
- **Data quality module** (`data_quality.py`): Comprehensive VCF parsing with quality checks
  - `DataQualityReport`: Reports annotation coverage, multi-allelic handling, and data usability
  - `VCFDataQualityError`: User-friendly exceptions with fix suggestions
  - `load_vcf_with_quality_check()`: Robust VCF loading with quality validation
  - `validate_vcf_for_pipeline()`: Pre-flight validation function
- **Graceful handling of missing annotations**: Pipeline continues with warnings instead of failing
- **Enhanced annotation helper** (`scripts/annotate_vcf.py`):
  - Verbose mode with detailed statistics
  - Validation-only mode to check existing VCF
  - Better VEP/ANNOVAR format detection
  - Comprehensive error messages with fix suggestions

#### Performance Module Enhancements (`utils/performance.py`)
- **Gzip support**: `chunked_vcf_reader()` now handles `.vcf.gz` files automatically
- **Multi-allelic expansion**: Chunked reader expands multi-allelic variants with allele-specific genotype counting
- **Consistent genotype parsing**: Uses `parse_genotype()` with `target_allele` parameter for accurate multi-allelic handling
- **CADD defaults in chunked processing**: `compute_gene_burdens_chunked()` now applies consequence-based CADD defaults (35/20/10)
- **Zero-variance pathway filtering**: `parallel_pathway_scores()` removes zero-variance pathways before Z-score normalization with clear error messages

#### Configuration Validation (`config.py`)
- **ConfigValidationError class**: Custom exception with field tracking and actionable fix suggestions
- **Enhanced `validate_config()`**: Now accepts `check_files` parameter to skip file existence checks during testing
- **`validate_gmt_file()` function**: Validates GMT pathway files with detailed error reporting
  - Checks for minimum 3 tab-separated fields per line
  - Validates minimum 2 genes per pathway
  - Reports duplicate pathway names and parsing errors

#### Analytical Reliability Improvements
- **GMM convergence checking**: All GMM fits now verify convergence and log warnings if not converged
- **GMM covariance regularization**: Added `reg_covar=1e-6` to all GMM calls for numerical stability
- **Zero-variance pathway handling**: Automatically detects and removes pathways with zero variance before normalization
- **CADD missing value handling**: Uses consequence-based defaults (35/20/10) instead of silent zeros
- **Consistent genotype parsing**: Unified allele-specific counting between bi-allelic and multi-allelic variants
- **Empty ARI array handling**: Validation gates now handle edge cases where no GMM fits converge

#### Test Suite Expansion
- **CLI test suite** (`tests/test_cli.py`): 20+ tests for command-line interface
  - Version and help display tests
  - Config loading and validation tests
  - Command-line override tests (--output, --seed, --quiet)
  - Error handling and exit code tests
- **Performance module tests** (`tests/test_performance.py`): 25+ tests for performance utilities
  - Chunked VCF reader tests (plain and gzipped)
  - Multi-allelic expansion and genotype parsing tests
  - Gene burden computation tests
  - Parallel pathway scoring tests with zero-variance handling
  - Memory estimation and downsampling tests
  - Progress tracking tests
- **Updated config tests** (`tests/test_config.py`): Tests for new validation functions
- **Updated data quality tests**: Tests for multi-allelic `parse_genotype()` with `target_allele`
- **Statistical rigor tests** (`tests/test_statistical_rigor.py`): 32 tests for FDR, burden weights, effect sizes
- **Clustering tests** (`tests/test_clustering.py`): 26 tests for algorithms, CV, comparison
- **Simulation tests** (`tests/test_simulation.py`): 24 tests for synthetic data, power analysis
- **Ancestry tests** (`tests/test_ancestry.py`): 44 tests for PCA, correction, independence, stratified analysis
- **Batch correction tests** (`tests/test_batch_correction.py`): 34 tests for detection, correction, validation
- **Sensitivity analysis tests** (`tests/test_sensitivity.py`): 27 tests for parameter variation, robustness
- **Total test count**: 347 tests (up from 242)

### Changed
- Pipeline now uses `data_quality` module for VCF loading
- Pipeline reports include data quality section
- Version bumped to 0.2.0-dev
- `parse_genotype()` now takes `target_allele` parameter for consistent multi-allelic handling
- `validate_config()` raises `ConfigValidationError` instead of generic `ValueError`

### Documentation
- Updated troubleshooting guide with comprehensive real-world data section
- Added VCF validation instructions
- Added multi-allelic variant handling explanation
- Added CADD score coverage guidance
- Updated API documentation for config and validation modules

### Other
- PyPI package publishing preparation

## [0.1.0] - 2026-01-29

### Added

#### Core Pipeline
- Complete pathway subtyping pipeline with VCF → clustering → report workflow
- Gene burden computation with LoF/missense weighting and CADD score normalization
- Pathway score aggregation using GMT file definitions
- GMM clustering with automatic cluster selection via BIC
- Cluster labeling based on top contributing pathways

#### Validation Gates
- Negative Control 1: Label shuffle test (ARI should be < 0.15)
- Negative Control 2: Random gene sets test (ARI should be < 0.15)
- Stability Test: Bootstrap resampling (ARI should be >= 0.8)
- Comprehensive validation reporting in JSON and Markdown

#### Disease Support
- Autism pathway definitions (15 pathways, ~200 genes) - Validated
- Schizophrenia pathway template
- Epilepsy pathway template
- Intellectual disability pathway template
- Guide for adapting to new diseases

#### Infrastructure
- `pyproject.toml` for modern Python packaging
- CLI entry points (`psf`, `pathway-subtyping`)
- YAML-based configuration system
- Reproducibility features (seed control, metadata logging)

#### Testing
- 64 unit and integration tests
- Test fixtures for synthetic data
- CI/CD with GitHub Actions
- Multi-OS (Ubuntu, macOS) and multi-Python (3.9-3.12) testing

#### Documentation
- Comprehensive README with quick start guide
- Disease adaptation guide
- Pathway curation guide
- Validation gates explanation
- Getting-started Jupyter notebook tutorial

#### Containerization
- Multi-stage Dockerfile (runtime, dev, jupyter targets)
- docker-compose.yml for easy orchestration
- Pre-configured development environment

#### Sample Data
- Synthetic VCF with 60 samples and 30 variants
- Synthetic phenotypes with 4 planted subtypes
- Test configuration ready to run out of the box

### Technical Details
- Python 3.8+ support
- Dependencies: numpy, pandas, scikit-learn, scipy, pysam, matplotlib, seaborn
- MIT License

## [0.0.1] - 2026-01-29

### Added
- Initial project structure
- Core module scaffolding
- Basic documentation

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| 0.3.0 | 2026-02-16 | Multi-omic integration: single-cell, deconvolution, signaling databases, cross-modal validation, fusion |
| 0.2.3 | 2026-02-14 | Cross-cohort validation, performance optimization, threshold calibration, expression scoring |
| 0.2.0 | 2026-02-09 | Scientific rigor, ancestry/batch correction, benchmarks, sensitivity analysis |
| 0.1.0 | 2026-01-29 | First public release with full pipeline |
| 0.0.1 | 2026-01-29 | Initial project setup |

[Unreleased]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.2.3...v0.3.0
[0.2.3]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.2.0...v0.2.3
[0.2.0]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/topmist-admin/pathway-subtyping-framework/releases/tag/v0.1.0
[0.0.1]: https://github.com/topmist-admin/pathway-subtyping-framework/releases/tag/v0.0.1
