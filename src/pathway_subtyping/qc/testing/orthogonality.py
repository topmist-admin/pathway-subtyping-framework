"""
Feature Orthogonality Matrix for Molecular QC Testing.

Runs the 12x12 confusion matrix: each defect type against all 12 QC features.
Validates that primary detectors fire (>90% sensitivity) and non-detectors
do not (< 10% false positive rate).

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from .simulator import DefectType, InjectedBatch, ManufacturingSimulator

logger = logging.getLogger(__name__)

# Feature IDs for matrix columns
FEATURE_IDS = [
    "F1_cascade",
    "F2_temporal",
    "F3_tension",
    "F4_resolution_gate",
    "F5_drift",
    "F6_offtarget",
    "F7_heterogeneity",
    "F8_dosage",
    "F9_crosstalk",
    "F10_feedback",
    "F11_stress",
    "F12_atlas",
]

# Expected detection matrix from test plan
# YES = primary (must fire), yes = secondary (may fire), - = must NOT fire
EXPECTED_DETECTIONS: Dict[str, Dict[str, str]] = {
    DefectType.STALLED_CASCADE.value: {
        "F1_cascade": "YES", "F3_tension": "yes", "F4_resolution_gate": "yes",
    },
    DefectType.TEMPORAL_DRIFT.value: {
        "F2_temporal": "YES", "F5_drift": "YES", "F3_tension": "yes", "F12_atlas": "yes",
    },
    DefectType.OFFTARGET.value: {
        "F6_offtarget": "YES",
    },
    DefectType.OVERDOSE.value: {
        "F8_dosage": "YES", "F3_tension": "yes",
    },
    DefectType.SUBPOPULATION.value: {
        "F7_heterogeneity": "YES", "F12_atlas": "yes",
    },
    DefectType.CROSSTALK.value: {
        "F9_crosstalk": "YES", "F1_cascade": "yes",
    },
    DefectType.BROKEN_FEEDBACK.value: {
        "F10_feedback": "YES",
    },
    DefectType.STRESS.value: {
        "F11_stress": "YES", "F3_tension": "yes",
    },
    DefectType.ATLAS_DRIFT.value: {
        "F12_atlas": "YES", "F7_heterogeneity": "yes",
    },
}


@dataclass
class FeatureDetection:
    """Result of running one QC feature on one injected batch."""

    feature_id: str
    detected: bool
    score: float = 0.0
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class OrthogonalityResult:
    """Result of the full 12x12 orthogonality test."""

    matrix: pd.DataFrame  # defect_type x feature_id -> detection rate
    primary_sensitivity: Dict[str, float] = field(default_factory=dict)
    false_positive_rate: Dict[str, float] = field(default_factory=dict)
    passed: bool = False
    summary: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "summary": self.summary,
            "primary_sensitivity": self.primary_sensitivity,
            "false_positive_rate": self.false_positive_rate,
        }


# Type alias for feature detector functions
FeatureDetector = Callable[[InjectedBatch], FeatureDetection]


class OrthogonalityMatrix:
    """Runs the 12x12 feature confusion matrix.

    Register QC feature detector functions, then run against all defect types.
    Features that are not yet implemented return stub results.

    Usage::

        ortho = OrthogonalityMatrix(simulator)
        ortho.register_detector("F6_offtarget", my_offtarget_detector)
        result = ortho.run(n_trials=10, severity=0.5)
        print(result.matrix)
    """

    def __init__(
        self,
        simulator: Optional[ManufacturingSimulator] = None,
        seed: Optional[int] = None,
    ):
        self.simulator = simulator or ManufacturingSimulator(seed=seed)
        self._detectors: Dict[str, FeatureDetector] = {}
        self._seed = seed

        # Register stubs for all features
        for fid in FEATURE_IDS:
            self._detectors[fid] = self._make_stub(fid)

    def register_detector(self, feature_id: str, detector: FeatureDetector) -> None:
        """Register a real detector function for a feature."""
        if feature_id not in FEATURE_IDS:
            raise ValueError(f"Unknown feature ID: {feature_id}. Must be one of {FEATURE_IDS}")
        self._detectors[feature_id] = detector
        logger.info("[QC Orthogonality] Registered detector: %s", feature_id)

    def run(
        self,
        n_trials: int = 10,
        severity: float = 0.5,
        n_cells: int = 200,
    ) -> OrthogonalityResult:
        """Run the full orthogonality matrix.

        Args:
            n_trials: Number of independent trials per defect type.
            severity: Default severity for all defect injections.
            n_cells: Cells per synthetic batch.

        Returns:
            OrthogonalityResult with detection rates and pass/fail assessment.
        """
        defect_types = list(DefectType)
        results = {dt.value: {fid: [] for fid in FEATURE_IDS} for dt in defect_types}

        for trial in range(n_trials):
            for dt in defect_types:
                batch = self._generate_defect_batch(dt, severity, n_cells)
                for fid in FEATURE_IDS:
                    detection = self._detectors[fid](batch)
                    results[dt.value][fid].append(detection.detected)

        # Compute detection rates
        matrix_data = {}
        for dt_val, feature_results in results.items():
            matrix_data[dt_val] = {
                fid: np.mean(detections) for fid, detections in feature_results.items()
            }

        matrix = pd.DataFrame(matrix_data).T
        matrix.index.name = "defect_type"
        matrix.columns.name = "feature"

        # Evaluate against expected detections
        primary_sens = {}
        fpr = {}
        for dt_val, expected in EXPECTED_DETECTIONS.items():
            for fid in FEATURE_IDS:
                level = expected.get(fid, "-")
                rate = matrix_data.get(dt_val, {}).get(fid, 0.0)

                if level == "YES":
                    key = f"{dt_val}/{fid}"
                    primary_sens[key] = rate

        # Compute false positive rates (non-detectors)
        for dt_val in [dt.value for dt in defect_types]:
            expected = EXPECTED_DETECTIONS.get(dt_val, {})
            for fid in FEATURE_IDS:
                if fid not in expected:
                    key = f"{dt_val}/{fid}"
                    rate = matrix_data.get(dt_val, {}).get(fid, 0.0)
                    fpr[key] = rate

        # Pass criteria
        all_primary_pass = all(s >= 0.9 for s in primary_sens.values()) if primary_sens else True
        all_fpr_pass = all(r <= 0.1 for r in fpr.values()) if fpr else True
        passed = all_primary_pass and all_fpr_pass

        n_primary = len(primary_sens)
        n_primary_pass = sum(1 for s in primary_sens.values() if s >= 0.9)
        n_fpr = len(fpr)
        n_fpr_pass = sum(1 for r in fpr.values() if r <= 0.1)

        summary = (
            f"Primary detectors: {n_primary_pass}/{n_primary} pass (>=90% sensitivity). "
            f"Non-detectors: {n_fpr_pass}/{n_fpr} pass (<=10% FPR)."
        )

        logger.info("[QC Orthogonality] %s", summary)

        return OrthogonalityResult(
            matrix=matrix,
            primary_sensitivity=primary_sens,
            false_positive_rate=fpr,
            passed=passed,
            summary=summary,
        )

    def _generate_defect_batch(
        self, defect_type: DefectType, severity: float, n_cells: int
    ) -> InjectedBatch:
        """Generate a batch with a specific defect type."""
        batch = self.simulator.generate_healthy_batch(n_cells=n_cells)

        if defect_type == DefectType.STALLED_CASCADE:
            return self.simulator.inject_stalled_cascades(
                batch, pathways=self.simulator.pathways[:2], severity=severity
            )
        elif defect_type == DefectType.TEMPORAL_DRIFT:
            passages = self.simulator.inject_temporal_drift(
                batch, n_passages=5, drift_rate=severity * 0.05
            )
            return passages[-1]
        elif defect_type == DefectType.OFFTARGET:
            return self.simulator.inject_offtarget(
                batch, pathways=self.simulator.pathways[5:7],
                cell_fraction=0.3, activation_level=severity,
            )
        elif defect_type == DefectType.OVERDOSE:
            return self.simulator.inject_overdose(
                batch, pathways=self.simulator.pathways[:2],
                multiplier=1.0 + severity * 4.0,
            )
        elif defect_type == DefectType.SUBPOPULATION:
            return self.simulator.inject_subpopulation(batch, fraction=0.15)
        elif defect_type == DefectType.CROSSTALK:
            return self.simulator.inject_crosstalk(
                batch, self.simulator.pathways[0], self.simulator.pathways[1],
            )
        elif defect_type == DefectType.BROKEN_FEEDBACK:
            return self.simulator.inject_broken_feedback(
                batch, loops=[(self.simulator.pathways[0], self.simulator.pathways[6])],
            )
        elif defect_type == DefectType.STRESS:
            from .simulator import StressorType
            return self.simulator.inject_stress(
                batch, stressor=StressorType.HYPOXIA, severity=severity,
            )
        elif defect_type == DefectType.ATLAS_DRIFT:
            return self.simulator.inject_atlas_drift(
                batch, drift_distance=severity,
            )
        else:
            return batch

    @staticmethod
    def _make_stub(feature_id: str) -> FeatureDetector:
        """Create a stub detector that always returns False (not detected)."""
        def stub(batch: InjectedBatch) -> FeatureDetection:
            return FeatureDetection(feature_id=feature_id, detected=False, score=0.0)
        return stub
