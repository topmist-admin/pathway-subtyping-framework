"""
Severity Titration Harness for Molecular QC Testing.

Tests each defect at 5 severity levels (0.1, 0.3, 0.5, 0.7, 0.9) to
establish detection thresholds. Validates:
    - Severity >= 0.5: >90% sensitivity (obvious defects)
    - Severity >= 0.3: >70% sensitivity (subtle defects)
    - Severity 0.1: <10% false positive rate (near-normal)

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from .orthogonality import FeatureDetection, FeatureDetector
from .simulator import DefectType, ManufacturingSimulator

logger = logging.getLogger(__name__)

SEVERITY_LEVELS = [0.1, 0.3, 0.5, 0.7, 0.9]


@dataclass
class TitrationCurve:
    """Detection rate at each severity level for one feature x defect pair."""

    feature_id: str
    defect_type: str
    severity_levels: List[float]
    detection_rates: List[float]
    threshold_50: Optional[float] = None  # Severity at 50% detection

    def to_dict(self) -> Dict[str, Any]:
        return {
            "feature_id": self.feature_id,
            "defect_type": self.defect_type,
            "severity_levels": self.severity_levels,
            "detection_rates": self.detection_rates,
            "threshold_50": self.threshold_50,
        }


@dataclass
class TitrationResult:
    """Result of severity titration across all defect types."""

    curves: List[TitrationCurve]
    high_severity_pass: bool = False  # All >=0.5 severity detected at >90%
    medium_severity_pass: bool = False  # All >=0.3 severity detected at >70%
    low_severity_pass: bool = False  # All 0.1 severity FPR <10%
    summary: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "high_severity_pass": self.high_severity_pass,
            "medium_severity_pass": self.medium_severity_pass,
            "low_severity_pass": self.low_severity_pass,
            "summary": self.summary,
            "n_curves": len(self.curves),
        }


class SeverityTitration:
    """Runs severity titration across all defect types for registered detectors.

    Usage::

        titration = SeverityTitration(simulator)
        titration.register_detector("F6_offtarget", my_offtarget_detector)
        result = titration.run(n_trials=10)
        for curve in result.curves:
            print(f"{curve.feature_id} x {curve.defect_type}: {curve.detection_rates}")
    """

    def __init__(
        self,
        simulator: Optional[ManufacturingSimulator] = None,
        seed: Optional[int] = None,
    ):
        self.simulator = simulator or ManufacturingSimulator(seed=seed)
        self._detectors: Dict[str, FeatureDetector] = {}

    def register_detector(self, feature_id: str, detector: FeatureDetector) -> None:
        """Register a detector function for a feature."""
        self._detectors[feature_id] = detector

    def run(
        self,
        n_trials: int = 10,
        n_cells: int = 200,
        severity_levels: Optional[List[float]] = None,
        defect_types: Optional[List[DefectType]] = None,
    ) -> TitrationResult:
        """Run severity titration.

        Args:
            n_trials: Independent trials per severity level.
            n_cells: Cells per synthetic batch.
            severity_levels: Override default levels.
            defect_types: Subset of defects to test. Defaults to all.

        Returns:
            TitrationResult with detection curves.
        """
        levels = severity_levels or SEVERITY_LEVELS
        dt_list = defect_types or list(DefectType)
        curves: List[TitrationCurve] = []

        for fid, detector in self._detectors.items():
            for dt in dt_list:
                rates = []
                for sev in levels:
                    detections = []
                    for _ in range(n_trials):
                        batch = self._generate_at_severity(dt, sev, n_cells)
                        result = detector(batch)
                        detections.append(result.detected)
                    rates.append(np.mean(detections))

                # Estimate threshold at 50% detection
                threshold_50 = None
                for i, rate in enumerate(rates):
                    if rate >= 0.5:
                        threshold_50 = levels[i]
                        break

                curves.append(
                    TitrationCurve(
                        feature_id=fid,
                        defect_type=dt.value,
                        severity_levels=list(levels),
                        detection_rates=rates,
                        threshold_50=threshold_50,
                    )
                )

        # Evaluate pass criteria
        high_results = []
        med_results = []
        low_results = []

        for curve in curves:
            for sev, rate in zip(curve.severity_levels, curve.detection_rates):
                if sev >= 0.5:
                    high_results.append(rate >= 0.9)
                if sev >= 0.3:
                    med_results.append(rate >= 0.7)
                if sev <= 0.1:
                    low_results.append(rate <= 0.1)

        high_pass = all(high_results) if high_results else True
        med_pass = all(med_results) if med_results else True
        low_pass = all(low_results) if low_results else True

        summary = (
            f"Titration: {len(curves)} curves across {len(self._detectors)} features. "
            f"High severity (>=0.5): {'PASS' if high_pass else 'FAIL'}. "
            f"Medium (>=0.3): {'PASS' if med_pass else 'FAIL'}. "
            f"Low (0.1 FPR): {'PASS' if low_pass else 'FAIL'}."
        )

        logger.info("[QC Titration] %s", summary)

        return TitrationResult(
            curves=curves,
            high_severity_pass=high_pass,
            medium_severity_pass=med_pass,
            low_severity_pass=low_pass,
            summary=summary,
        )

    def _generate_at_severity(
        self, defect_type: DefectType, severity: float, n_cells: int
    ) -> "InjectedBatch":
        """Generate a defect batch at a specific severity level."""
        from .orthogonality import OrthogonalityMatrix

        ortho = OrthogonalityMatrix(simulator=self.simulator)
        return ortho._generate_defect_batch(defect_type, severity, n_cells)
