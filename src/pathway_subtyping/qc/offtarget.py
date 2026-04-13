"""
Off-Target Pathway Activation Detection (F6).

Detects unintended pathway activations in manufactured cells by comparing
observed pathway activation against an engineering specification.

Classification:
    - INTENDED: Pathway in spec and activated -> expected
    - TOLERATED: Not in spec but below noise threshold -> acceptable
    - OFF_TARGET: Not in spec and above noise threshold -> investigate
    - EXCLUDED_VIOLATION: Explicitly excluded and activated -> safety flag

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class ActivationClass(Enum):
    """Classification of pathway activation relative to spec."""

    INTENDED = "intended"
    TOLERATED = "tolerated"
    OFF_TARGET = "off_target"
    EXCLUDED_VIOLATION = "excluded_violation"


@dataclass
class PathwayActivationResult:
    """Classification result for a single pathway."""

    pathway: str
    classification: ActivationClass
    mean_activation: float
    threshold: float
    n_cells_active: int = 0
    n_cells_total: int = 0

    @property
    def fraction_active(self) -> float:
        if self.n_cells_total == 0:
            return 0.0
        return self.n_cells_active / self.n_cells_total

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway": self.pathway,
            "classification": self.classification.value,
            "mean_activation": round(self.mean_activation, 4),
            "threshold": round(self.threshold, 4),
            "fraction_active": round(self.fraction_active, 4),
        }


@dataclass
class OffTargetResult:
    """Result of off-target pathway activation analysis."""

    pathway_results: List[PathwayActivationResult]
    n_off_target: int = 0
    n_excluded_violations: int = 0
    safety_flag: bool = False
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_excluded_violations == 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "safety_flag": self.safety_flag,
            "n_off_target": self.n_off_target,
            "n_excluded_violations": self.n_excluded_violations,
            "summary": self.summary,
            "pathways": [r.to_dict() for r in self.pathway_results],
        }

    def get_off_target_pathways(self) -> List[str]:
        """Return names of off-target pathways."""
        return [
            r.pathway
            for r in self.pathway_results
            if r.classification == ActivationClass.OFF_TARGET
        ]

    def get_excluded_violations(self) -> List[str]:
        """Return names of excluded pathways that violated safety spec."""
        return [
            r.pathway
            for r in self.pathway_results
            if r.classification == ActivationClass.EXCLUDED_VIOLATION
        ]


class OffTargetDetector:
    """Detects off-target pathway activations against an engineering spec.

    Usage::

        detector = OffTargetDetector(
            intended_pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            excluded_pathways=["HALLMARK_APOPTOSIS"],
        )
        result = detector.scan(pathway_scores)
    """

    def __init__(
        self,
        intended_pathways: Optional[List[str]] = None,
        excluded_pathways: Optional[List[str]] = None,
        noise_threshold: float = 0.5,
        cell_fraction_threshold: float = 0.1,
        seed: Optional[int] = None,
    ):
        """Initialize the off-target detector.

        Args:
            intended_pathways: Pathways that SHOULD be active.
            excluded_pathways: Pathways that must NOT be active (safety-critical).
            noise_threshold: Z-score threshold above which activation is considered real.
            cell_fraction_threshold: Minimum fraction of cells that must show
                activation for a pathway to be flagged.
            seed: Random seed (unused, kept for API consistency).
        """
        self.intended_pathways = set(intended_pathways or [])
        self.excluded_pathways = set(excluded_pathways or [])
        self.noise_threshold = noise_threshold
        self.cell_fraction_threshold = cell_fraction_threshold

    def calibrate_threshold(
        self,
        control_scores: pd.DataFrame,
        percentile: float = 95.0,
    ) -> float:
        """Calibrate noise threshold from control sample data.

        Args:
            control_scores: Pathway scores from control/reference samples.
            percentile: Percentile of control activation to use as threshold.

        Returns:
            Calibrated noise threshold.
        """
        values = control_scores.values.flatten()
        threshold = float(np.percentile(values, percentile))
        self.noise_threshold = threshold
        logger.info(
            "[QC OffTarget] Calibrated noise threshold: %.4f (%.0fth percentile)",
            threshold,
            percentile,
        )
        return threshold

    def scan(self, pathway_scores: pd.DataFrame) -> OffTargetResult:
        """Scan pathway scores for off-target activations.

        Args:
            pathway_scores: DataFrame (cells x pathways) of pathway activation scores.

        Returns:
            OffTargetResult with per-pathway classification.
        """
        n_cells = len(pathway_scores)
        results: List[PathwayActivationResult] = []
        n_off_target = 0
        n_excluded = 0
        safety_flag = False

        for pathway in pathway_scores.columns:
            scores = pathway_scores[pathway].values
            mean_activation = float(np.mean(scores))
            n_active = int(np.sum(scores > self.noise_threshold))

            if pathway in self.intended_pathways:
                classification = ActivationClass.INTENDED
            elif pathway in self.excluded_pathways:
                fraction_active = n_active / n_cells if n_cells > 0 else 0.0
                if fraction_active > self.cell_fraction_threshold:
                    classification = ActivationClass.EXCLUDED_VIOLATION
                    n_excluded += 1
                    safety_flag = True
                else:
                    classification = ActivationClass.TOLERATED
            else:
                fraction_active = n_active / n_cells if n_cells > 0 else 0.0
                if fraction_active > self.cell_fraction_threshold:
                    classification = ActivationClass.OFF_TARGET
                    n_off_target += 1
                else:
                    classification = ActivationClass.TOLERATED

            results.append(
                PathwayActivationResult(
                    pathway=pathway,
                    classification=classification,
                    mean_activation=mean_activation,
                    threshold=self.noise_threshold,
                    n_cells_active=n_active,
                    n_cells_total=n_cells,
                )
            )

        # Sort by severity: excluded violations first, then off-target
        results.sort(
            key=lambda r: (
                0 if r.classification == ActivationClass.EXCLUDED_VIOLATION else
                1 if r.classification == ActivationClass.OFF_TARGET else 2
            )
        )

        summary = (
            f"Off-target scan: {n_off_target} off-target, "
            f"{n_excluded} excluded violations. "
            f"{'SAFETY FLAG' if safety_flag else 'No safety issues'}."
        )

        logger.info("[QC OffTarget] %s", summary)

        return OffTargetResult(
            pathway_results=results,
            n_off_target=n_off_target,
            n_excluded_violations=n_excluded,
            safety_flag=safety_flag,
            summary=summary,
        )
