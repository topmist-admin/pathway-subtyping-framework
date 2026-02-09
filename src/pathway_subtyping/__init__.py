"""
Pathway Subtyping Framework

A disease-agnostic framework for pathway-based molecular subtype discovery
in genetically heterogeneous conditions.
"""

__version__ = "0.2.2"
__author__ = "Rohit Chauhan"

from .ancestry import (
    AncestryAdjustmentResult,
    AncestryMethod,
    AncestryPCs,
    AncestryStratificationReport,
    adjust_pathway_scores,
    check_ancestry_independence,
    compute_ancestry_correlation,
    compute_ancestry_pcs,
    stratified_analysis,
)
from .benchmark import (
    BenchmarkComparisonResult,
    BenchmarkMethod,
    BenchmarkResult,
    BenchmarkSweepResult,
    run_benchmark_comparison,
    run_benchmark_sweep,
    run_single_benchmark,
)
from .batch_correction import (
    BatchCorrectionMethod,
    BatchCorrectionResult,
    BatchEffectReport,
    correct_batch_effects,
    detect_batch_effects,
    validate_batch_correction,
)
from .clustering import (
    AlgorithmComparisonResult,
    ClusteringAlgorithm,
    ClusteringResult,
    CrossValidationResult,
    ModelSelectionResult,
    compare_algorithms,
    cross_validate_clustering,
    run_clustering,
    select_n_clusters,
)
from .cross_cohort import (
    CohortResult,
    CrossCohortResult,
    batch_compare_cohorts,
    compare_cohorts,
    load_cohort_result,
)
from .data_quality import (
    DataQualityReport,
    VCFDataQualityError,
    load_vcf_with_quality_check,
    validate_vcf_for_pipeline,
)
from .pipeline import DemoPipeline, PipelineConfig
from .sensitivity import (
    ParameterVariationResult,
    SensitivityAnalysisResult,
    SensitivityParameter,
    run_sensitivity_analysis,
    vary_clustering_algorithm,
    vary_feature_subset,
    vary_n_clusters,
    vary_normalization,
)
from .simulation import (
    PowerAnalysisResult,
    RecoveryResult,
    SampleSizeAnalysisResult,
    SimulatedData,
    SimulationConfig,
    TypeIErrorResult,
    estimate_type_i_error,
    evaluate_recovery,
    generate_synthetic_data,
    run_power_analysis,
    run_sample_size_analysis,
    validate_framework,
)

# Scientific rigor modules (v0.2)
from .statistical_rigor import (
    BurdenWeights,
    BurdenWeightScheme,
    FDRResult,
    PathwayNormalization,
    StatisticalRigorResult,
    benjamini_hochberg,
    compute_pathway_effect_sizes,
    compute_pathway_pvalues,
    run_statistical_analysis,
)
from .validation import ValidationGates, ValidationGatesResult, ValidationResult

__all__ = [
    # Pipeline
    "DemoPipeline",
    "PipelineConfig",
    # Validation
    "ValidationGates",
    "ValidationGatesResult",
    "ValidationResult",
    # Cross-cohort
    "compare_cohorts",
    "load_cohort_result",
    "batch_compare_cohorts",
    "CohortResult",
    "CrossCohortResult",
    # Data quality
    "DataQualityReport",
    "VCFDataQualityError",
    "load_vcf_with_quality_check",
    "validate_vcf_for_pipeline",
    # Statistical rigor
    "BurdenWeights",
    "BurdenWeightScheme",
    "FDRResult",
    "PathwayNormalization",
    "StatisticalRigorResult",
    "benjamini_hochberg",
    "compute_pathway_pvalues",
    "compute_pathway_effect_sizes",
    "run_statistical_analysis",
    # Clustering
    "ClusteringAlgorithm",
    "ClusteringResult",
    "ModelSelectionResult",
    "CrossValidationResult",
    "AlgorithmComparisonResult",
    "run_clustering",
    "select_n_clusters",
    "cross_validate_clustering",
    "compare_algorithms",
    # Simulation
    "SimulationConfig",
    "SimulatedData",
    "RecoveryResult",
    "TypeIErrorResult",
    "PowerAnalysisResult",
    "SampleSizeAnalysisResult",
    "generate_synthetic_data",
    "evaluate_recovery",
    "estimate_type_i_error",
    "run_power_analysis",
    "run_sample_size_analysis",
    "validate_framework",
    # Benchmark comparison
    "BenchmarkComparisonResult",
    "BenchmarkMethod",
    "BenchmarkResult",
    "BenchmarkSweepResult",
    "run_benchmark_comparison",
    "run_benchmark_sweep",
    "run_single_benchmark",
    # Batch correction
    "BatchCorrectionMethod",
    "BatchCorrectionResult",
    "BatchEffectReport",
    "correct_batch_effects",
    "detect_batch_effects",
    "validate_batch_correction",
    # Sensitivity analysis
    "ParameterVariationResult",
    "SensitivityAnalysisResult",
    "SensitivityParameter",
    "run_sensitivity_analysis",
    "vary_clustering_algorithm",
    "vary_feature_subset",
    "vary_n_clusters",
    "vary_normalization",
    # Ancestry correction
    "AncestryAdjustmentResult",
    "AncestryMethod",
    "AncestryPCs",
    "AncestryStratificationReport",
    "adjust_pathway_scores",
    "check_ancestry_independence",
    "compute_ancestry_correlation",
    "compute_ancestry_pcs",
    "stratified_analysis",
    # Meta
    "__version__",
]
