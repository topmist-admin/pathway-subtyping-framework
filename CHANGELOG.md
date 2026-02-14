# Changelog

All notable changes to the Pathway Subtyping Framework will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.3] - 2026-02-14

### Added

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
| 0.2.0 | 2026-02-09 | Scientific rigor, ancestry/batch correction, benchmarks, sensitivity analysis |
| 0.1.0 | 2026-01-29 | First public release with full pipeline |
| 0.0.1 | 2026-01-29 | Initial project setup |

[Unreleased]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/topmist-admin/pathway-subtyping-framework/releases/tag/v0.1.0
[0.0.1]: https://github.com/topmist-admin/pathway-subtyping-framework/releases/tag/v0.0.1
