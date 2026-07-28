"""
Temporal Pathway Tracking (F2).

Tracks pathway activation trajectories across manufacturing timepoints
to detect stalls, reversals, and oscillations.

Trajectory types:
    - RESOLVING: Activation increases monotonically toward target
    - STALLED: Activation plateaus below target for >=2 timepoints
    - REVERSING: Activation peaks then declines
    - OSCILLATING: Activation fluctuates without trend

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class TrajectoryType(Enum):
    """Classification of pathway activation trajectory."""

    RESOLVING = "resolving"
    STALLED = "stalled"
    REVERSING = "reversing"
    OSCILLATING = "oscillating"
    INSUFFICIENT_DATA = "insufficient_data"


@dataclass
class PathwayTrajectory:
    """Trajectory of a single pathway across timepoints."""

    pathway: str
    timepoint_ids: List[str]
    scores: List[float]
    trajectory_type: TrajectoryType = TrajectoryType.INSUFFICIENT_DATA
    expected_target: Optional[float] = None

    @property
    def n_timepoints(self) -> int:
        return len(self.timepoint_ids)

    @property
    def trend(self) -> float:
        """Linear trend slope (positive = increasing)."""
        if len(self.scores) < 2:
            return 0.0
        x = np.arange(len(self.scores), dtype=float)
        coeffs = np.polyfit(x, self.scores, 1)
        return float(coeffs[0])

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway": self.pathway,
            "trajectory_type": self.trajectory_type.value,
            "n_timepoints": self.n_timepoints,
            "scores": [round(s, 4) for s in self.scores],
            "trend": round(self.trend, 4),
            "expected_target": self.expected_target,
        }


@dataclass
class TemporalResult:
    """Result of temporal pathway tracking."""

    trajectories: List[PathwayTrajectory]
    n_stalled: int = 0
    n_reversing: int = 0
    n_oscillating: int = 0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_stalled == 0 and self.n_reversing == 0

    def get_problem_pathways(self) -> List[str]:
        """Return pathways with stalled or reversing trajectories."""
        return [
            t.pathway
            for t in self.trajectories
            if t.trajectory_type in (TrajectoryType.STALLED, TrajectoryType.REVERSING)
        ]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_stalled": self.n_stalled,
            "n_reversing": self.n_reversing,
            "n_oscillating": self.n_oscillating,
            "summary": self.summary,
            "trajectories": [t.to_dict() for t in self.trajectories],
        }


class TemporalTracker:
    """Tracks pathway activation trajectories across timepoints.

    Usage::

        tracker = TemporalTracker()
        tracker.register_timepoint(expr_day0, "day_0")
        tracker.register_timepoint(expr_day3, "day_3")
        tracker.register_timepoint(expr_day7, "day_7")
        result = tracker.compute_all_trajectories()
    """

    def __init__(
        self,
        expected_profile: Optional[Dict[str, float]] = None,
        stall_threshold: float = 0.05,
        min_timepoints_for_stall: int = 2,
        seed: Optional[int] = None,
    ):
        """Initialize the temporal tracker.

        Args:
            expected_profile: Target pathway scores for "done" state.
            stall_threshold: Max change between timepoints to count as plateau.
            min_timepoints_for_stall: Consecutive plateaus before flagging stall.
            seed: Random seed (unused, kept for API consistency).
        """
        self.expected_profile = expected_profile or {}
        self.stall_threshold = stall_threshold
        self.min_timepoints_for_stall = min_timepoints_for_stall
        self._timepoints: List[str] = []
        self._scores: List[pd.DataFrame] = []

    @property
    def n_timepoints(self) -> int:
        return len(self._timepoints)

    def register_timepoint(
        self,
        pathway_scores: pd.DataFrame,
        timepoint_id: str,
    ) -> None:
        """Register a timepoint's pathway scores.

        Args:
            pathway_scores: DataFrame (cells x pathways) at this timepoint.
            timepoint_id: Identifier for this timepoint (e.g., "day_3").
        """
        self._timepoints.append(timepoint_id)
        self._scores.append(pathway_scores)
        logger.info(
            "[QC Temporal] Registered timepoint %s (%d cells, %d pathways)",
            timepoint_id,
            len(pathway_scores),
            pathway_scores.shape[1],
        )

    def compute_trajectory(self, pathway: str) -> PathwayTrajectory:
        """Compute trajectory for a single pathway.

        Args:
            pathway: Pathway column name.

        Returns:
            PathwayTrajectory with scores and classification.
        """
        scores = []
        for df in self._scores:
            if pathway in df.columns:
                scores.append(float(df[pathway].mean()))
            else:
                scores.append(0.0)

        target = self.expected_profile.get(pathway)
        traj_type = self._classify_trajectory(scores, target)

        return PathwayTrajectory(
            pathway=pathway,
            timepoint_ids=list(self._timepoints),
            scores=scores,
            trajectory_type=traj_type,
            expected_target=target,
        )

    def compute_all_trajectories(
        self,
        pathways: Optional[List[str]] = None,
    ) -> TemporalResult:
        """Compute trajectories for all pathways.

        Args:
            pathways: Subset of pathways. If None, uses all pathways
                from the first timepoint.

        Returns:
            TemporalResult with all trajectories classified.
        """
        if not self._scores:
            return TemporalResult(trajectories=[], summary="No timepoints registered.")

        if pathways is None:
            pathways = list(self._scores[0].columns)

        trajectories = [self.compute_trajectory(pw) for pw in pathways]

        n_stalled = sum(1 for t in trajectories if t.trajectory_type == TrajectoryType.STALLED)
        n_reversing = sum(1 for t in trajectories if t.trajectory_type == TrajectoryType.REVERSING)
        n_oscillating = sum(
            1 for t in trajectories if t.trajectory_type == TrajectoryType.OSCILLATING
        )

        summary = (
            f"Temporal: {len(trajectories)} pathways across {self.n_timepoints} timepoints. "
            f"{n_stalled} stalled, {n_reversing} reversing, {n_oscillating} oscillating. "
            f"{'PASS' if n_stalled == 0 and n_reversing == 0 else 'FAIL'}."
        )

        logger.info("[QC Temporal] %s", summary)

        return TemporalResult(
            trajectories=trajectories,
            n_stalled=n_stalled,
            n_reversing=n_reversing,
            n_oscillating=n_oscillating,
            summary=summary,
        )

    def _classify_trajectory(
        self,
        scores: List[float],
        target: Optional[float],
    ) -> TrajectoryType:
        """Classify a trajectory into one of four types.

        Order matters: oscillation is checked before reversal because
        oscillating trajectories also have a peak (and would match reversal).
        """
        if len(scores) < 2:
            return TrajectoryType.INSUFFICIENT_DATA

        diffs = [scores[i + 1] - scores[i] for i in range(len(scores) - 1)]

        # Check for stall: consecutive small changes below target
        consecutive_stalls = 0
        for d in diffs:
            if abs(d) < self.stall_threshold:
                consecutive_stalls += 1
            else:
                consecutive_stalls = 0
        if consecutive_stalls >= self.min_timepoints_for_stall:
            if target is not None and scores[-1] < target:
                return TrajectoryType.STALLED

        # Check for oscillation BEFORE reversal: multiple sign changes in diffs
        sign_changes = sum(1 for i in range(len(diffs) - 1) if diffs[i] * diffs[i + 1] < 0)
        if sign_changes >= 2:
            return TrajectoryType.OSCILLATING

        # Check for reversal: rises then falls (single peak)
        if len(scores) >= 3:
            peak_idx = int(np.argmax(scores))
            if 0 < peak_idx < len(scores) - 1:
                rise = scores[peak_idx] - scores[0]
                fall = scores[peak_idx] - scores[-1]
                if rise > self.stall_threshold and fall > self.stall_threshold:
                    return TrajectoryType.REVERSING

        return TrajectoryType.RESOLVING
