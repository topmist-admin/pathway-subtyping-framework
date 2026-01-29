"""
Pathway Subtyping Framework

A disease-agnostic framework for pathway-based molecular subtype discovery
in genetically heterogeneous conditions.
"""

__version__ = "0.1.0"
__author__ = "Rohit Chauhan"

from .pipeline import DemoPipeline, PipelineConfig
from .validation import ValidationGates, ValidationGatesResult, ValidationResult

__all__ = [
    "DemoPipeline",
    "PipelineConfig",
    "ValidationGates",
    "ValidationGatesResult",
    "ValidationResult",
    "__version__",
]
