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

__all__ = [
    "DemoPipeline",
    "PipelineConfig",
    "ValidationGates",
    "ValidationGatesResult",
    "ValidationResult",
    "compare_cohorts",
    "load_cohort_result",
    "batch_compare_cohorts",
    "CohortResult",
    "CrossCohortResult",
    # Data quality (v0.2)
    "DataQualityReport",
    "VCFDataQualityError",
    "load_vcf_with_quality_check",
    "validate_vcf_for_pipeline",
    "__version__",
]
