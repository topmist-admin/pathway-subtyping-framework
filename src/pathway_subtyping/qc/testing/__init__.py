"""
Test infrastructure for the Molecular QC Layer.

Provides:
    - ManufacturingSimulator: Synthetic batch generation with injectable defects
    - OrthogonalityMatrix: 12x12 feature confusion matrix runner
    - SeverityTitration: Multi-severity detection threshold calibration
"""

from .simulator import (
    DefectLabel,
    DefectMetadata,
    InjectedBatch,
    ManufacturingSimulator,
    ManufacturingSpec,
)

__all__ = [
    "ManufacturingSimulator",
    "ManufacturingSpec",
    "InjectedBatch",
    "DefectLabel",
    "DefectMetadata",
]
