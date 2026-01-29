"""
Pathway Subtyping Framework

A disease-agnostic framework for pathway-based molecular subtype discovery
in genetically heterogeneous conditions.
"""

__version__ = "0.1.0"
__author__ = "Rohit Chauhan"

from .pipeline import DemoPipeline, PipelineConfig
from .validation import ValidationGates, ValidationGatesResult, ValidationResult
from .cross_cohort import (
    compare_cohorts,
    load_cohort_result,
    batch_compare_cohorts,
    CohortResult,
    CrossCohortResult,
)

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
    "__version__",
]
