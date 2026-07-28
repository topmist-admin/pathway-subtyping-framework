"""
Drift Detection via Open-Loop Accumulation (F5).

Detects slow pathway-level drift across passages by measuring cumulative
departure from a known-good baseline.

Formula:
    D(t) = sum_p |CCR_p(t) - CCR_p(t0)| * A_p(t)

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class DriftMeasurement:
    """Drift measurement at a single timepoint."""

    timepoint_id: str
    drift_score: float
    per_pathway_drift: Dict[str, float] = field(default_factory=dict)
    alarm_triggered: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "timepoint": self.timepoint_id,
            "drift_score": round(self.drift_score, 4),
            "alarm_triggered": self.alarm_triggered,
            "top_drivers": dict(
                sorted(self.per_pathway_drift.items(), key=lambda x: -abs(x[1]))[:5]
            ),
        }


@dataclass
class DriftResult:
    """Result of drift detection across timepoints."""

    measurements: List[DriftMeasurement]
    alarm_triggered: bool = False
    alarm_timepoint: Optional[str] = None
    drift_drivers: List[str] = field(default_factory=list)
    summary: str = ""

    @property
    def passed(self) -> bool:
        return not self.alarm_triggered

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "alarm_triggered": self.alarm_triggered,
            "alarm_timepoint": self.alarm_timepoint,
            "drift_drivers": self.drift_drivers,
            "summary": self.summary,
            "measurements": [m.to_dict() for m in self.measurements],
        }


class DriftDetector:
    """Detects pathway-level drift from a baseline.

    Usage::

        detector = DriftDetector(max_drift=0.5)
        detector.establish_baseline(baseline_ccr)
        for passage_ccr, passage_id in passages:
            measurement = detector.measure_drift(passage_ccr, passage_id)
    """

    def __init__(
        self,
        max_drift: float = 0.5,
        per_pathway_max: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None,
    ):
        """Initialize the drift detector.

        Args:
            max_drift: Maximum acceptable cumulative drift score.
            per_pathway_max: Optional per-pathway drift limits.
            seed: Random seed (unused, kept for API consistency).
        """
        self.max_drift = max_drift
        self.per_pathway_max = per_pathway_max or {}
        self._baseline: Optional[Dict[str, float]] = None
        self._measurements: List[DriftMeasurement] = []

    def establish_baseline(self, ccr_dict: Dict[str, float]) -> None:
        """Set the reference CCR profile (known-good state).

        Args:
            ccr_dict: Pathway -> mean cascade completion ratio at baseline.
        """
        self._baseline = dict(ccr_dict)
        self._measurements = []
        logger.info("[QC Drift] Baseline established: %d pathways", len(self._baseline))

    def measure_drift(
        self,
        ccr_dict: Dict[str, float],
        timepoint_id: str,
    ) -> DriftMeasurement:
        """Measure drift from baseline at a timepoint.

        Args:
            ccr_dict: Current pathway -> mean CCR.
            timepoint_id: Identifier for this timepoint.

        Returns:
            DriftMeasurement with per-pathway and aggregate drift.
        """
        if self._baseline is None:
            raise ValueError("Baseline not established. Call establish_baseline() first.")

        per_pathway: Dict[str, float] = {}
        total_drift = 0.0

        for pathway, baseline_ccr in self._baseline.items():
            current_ccr = ccr_dict.get(pathway, baseline_ccr)
            drift = abs(current_ccr - baseline_ccr)
            per_pathway[pathway] = drift
            total_drift += drift

        alarm = total_drift >= self.max_drift
        measurement = DriftMeasurement(
            timepoint_id=timepoint_id,
            drift_score=total_drift,
            per_pathway_drift=per_pathway,
            alarm_triggered=alarm,
        )
        self._measurements.append(measurement)

        if alarm:
            logger.warning(
                "[QC Drift] ALARM at %s: drift=%.3f (threshold=%.3f)",
                timepoint_id,
                total_drift,
                self.max_drift,
            )

        return measurement

    def drift_trajectory(self) -> DriftResult:
        """Get the full drift trajectory across all measured timepoints.

        Returns:
            DriftResult with all measurements and alarm status.
        """
        alarm_triggered = any(m.alarm_triggered for m in self._measurements)
        alarm_tp = None
        for m in self._measurements:
            if m.alarm_triggered:
                alarm_tp = m.timepoint_id
                break

        # Identify top drift drivers from latest measurement
        drivers: List[str] = []
        if self._measurements:
            latest = self._measurements[-1]
            sorted_pw = sorted(latest.per_pathway_drift.items(), key=lambda x: -x[1])
            drivers = [pw for pw, _ in sorted_pw[:3]]

        summary = (
            f"Drift: {len(self._measurements)} timepoints measured. "
            f"{'ALARM at ' + alarm_tp if alarm_triggered else 'No alarm'}. "
            f"{'FAIL' if alarm_triggered else 'PASS'}."
        )

        logger.info("[QC Drift] %s", summary)

        return DriftResult(
            measurements=list(self._measurements),
            alarm_triggered=alarm_triggered,
            alarm_timepoint=alarm_tp,
            drift_drivers=drivers,
            summary=summary,
        )

    def identify_drivers(self) -> List[str]:
        """Identify pathways driving the most drift."""
        if not self._measurements:
            return []
        latest = self._measurements[-1]
        sorted_pw = sorted(latest.per_pathway_drift.items(), key=lambda x: -x[1])
        return [pw for pw, _ in sorted_pw[:5]]
