"""
Base scenario class for end-to-end QC pipeline tests.

Research use only. Not for clinical decision-making.
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

import pandas as pd

from ..simulator import InjectedBatch, ManufacturingSimulator, ManufacturingSpec

logger = logging.getLogger(__name__)


@dataclass
class TimePointResult:
    """QC results at a single timepoint."""

    timepoint_id: str
    batch: InjectedBatch
    qc_scores: Dict[str, float] = field(default_factory=dict)
    detections: Dict[str, bool] = field(default_factory=dict)
    decision: str = "PENDING"  # PASS / HOLD / REJECT / PENDING

    def to_dict(self) -> Dict[str, Any]:
        return {
            "timepoint": self.timepoint_id,
            "decision": self.decision,
            "n_cells": batch.n_cells if (batch := self.batch) else 0,
            "qc_scores": self.qc_scores,
            "detections": self.detections,
        }


@dataclass
class ScenarioReport:
    """Full scenario QC report."""

    scenario_name: str
    timepoints: List[TimePointResult]
    final_decision: str = "PENDING"
    identified_defects: List[str] = field(default_factory=list)
    remediation_recommendations: List[str] = field(default_factory=list)
    passed_rubric: Dict[str, bool] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "scenario": self.scenario_name,
            "final_decision": self.final_decision,
            "n_timepoints": len(self.timepoints),
            "identified_defects": self.identified_defects,
            "remediation_recommendations": self.remediation_recommendations,
            "rubric": self.passed_rubric,
        }


# Type for QC feature runner functions
FeatureRunner = Callable[[InjectedBatch], Dict[str, Any]]


class BaseScenario(ABC):
    """Abstract base class for end-to-end QC scenarios.

    Subclasses define the timeline, injection points, and expected outcomes.
    The base class provides the execution and report generation framework.
    """

    def __init__(
        self,
        simulator: Optional[ManufacturingSimulator] = None,
        seed: Optional[int] = None,
    ):
        self.simulator = simulator or ManufacturingSimulator(seed=seed)
        self._feature_runners: Dict[str, FeatureRunner] = {}

    def register_feature(self, feature_id: str, runner: FeatureRunner) -> None:
        """Register a QC feature runner."""
        self._feature_runners[feature_id] = runner

    @abstractmethod
    def get_spec(self) -> ManufacturingSpec:
        """Return the manufacturing specification for this scenario."""

    @abstractmethod
    def get_timepoints(self) -> List[str]:
        """Return ordered timepoint identifiers."""

    @abstractmethod
    def generate_timepoint(
        self, timepoint_id: str, previous: Optional[InjectedBatch]
    ) -> InjectedBatch:
        """Generate the batch for a specific timepoint."""

    @abstractmethod
    def get_expected_outcomes(self) -> Dict[str, Dict[str, bool]]:
        """Return expected detection outcomes per timepoint per feature."""

    def run(self) -> ScenarioReport:
        """Execute the full scenario and generate a report."""
        timepoints = self.get_timepoints()
        results: List[TimePointResult] = []
        previous_batch: Optional[InjectedBatch] = None

        for tp_id in timepoints:
            batch = self.generate_timepoint(tp_id, previous_batch)

            # Run all registered features
            qc_scores = {}
            detections = {}
            for fid, runner in self._feature_runners.items():
                try:
                    result = runner(batch)
                    qc_scores[fid] = result.get("score", 0.0)
                    detections[fid] = result.get("detected", False)
                except Exception as e:
                    logger.warning(
                        "[QC Scenario] Feature %s failed at %s: %s", fid, tp_id, e
                    )
                    qc_scores[fid] = 0.0
                    detections[fid] = False

            tp_result = TimePointResult(
                timepoint_id=tp_id,
                batch=batch,
                qc_scores=qc_scores,
                detections=detections,
            )
            results.append(tp_result)
            previous_batch = batch

        return ScenarioReport(
            scenario_name=self.__class__.__name__,
            timepoints=results,
        )
