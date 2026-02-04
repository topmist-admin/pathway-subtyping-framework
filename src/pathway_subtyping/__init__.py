"""
Pathway Subtyping Framework

A disease-agnostic framework for pathway-based molecular subtype discovery
in genetically heterogeneous conditions.
"""

__version__ = "0.2.0-dev"
__author__ = "Rohit Chauhan"

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
from .validation import ValidationGates, ValidationGatesResult, ValidationResult

# Scientific rigor modules (v0.2)
from .statistical_rigor import (
    BurdenWeights,
    BurdenWeightScheme,
    FDRResult,
    PathwayNormalization,
    StatisticalRigorResult,
    benjamini_hochberg,
    compute_pathway_pvalues,
    compute_pathway_effect_sizes,
    run_statistical_analysis,
)
from .clustering import (
    ClusteringAlgorithm,
    ClusteringResult,
    ModelSelectionResult,
    CrossValidationResult,
    AlgorithmComparisonResult,
    run_clustering,
    select_n_clusters,
    cross_validate_clustering,
    compare_algorithms,
)
from .simulation import (
    SimulationConfig,
    SimulatedData,
    RecoveryResult,
    TypeIErrorResult,
    PowerAnalysisResult,
    SampleSizeAnalysisResult,
    generate_synthetic_data,
    evaluate_recovery,
    estimate_type_i_error,
    run_power_analysis,
    run_sample_size_analysis,
    validate_framework,
)

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
    # Meta
    "__version__",
]
