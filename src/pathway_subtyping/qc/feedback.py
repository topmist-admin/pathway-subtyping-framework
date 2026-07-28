"""
Feedback Loop Integrity Monitoring (F10).

Monitors negative feedback loops that prevent runaway pathway activation.
Broken feedback = pathway runs unchecked.

Integrity score:
    FI > 0.3:  Feedback intact (inhibitor responds to activator)
    FI ~ 0:    Feedback decoupled (inhibitor not tracking activator)
    FI < -0.3: Feedback inverted (pathological)

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class FeedbackStatus(Enum):
    """Status of a feedback loop."""

    INTACT = "intact"
    DECOUPLED = "decoupled"
    INVERTED = "inverted"
    UNKNOWN = "unknown"


@dataclass
class FeedbackLoopResult:
    """Result for a single feedback loop."""

    activator: str
    inhibitor: str
    integrity_score: float
    status: FeedbackStatus
    description: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "activator": self.activator,
            "inhibitor": self.inhibitor,
            "integrity_score": round(self.integrity_score, 4),
            "status": self.status.value,
            "description": self.description,
        }


@dataclass
class FeedbackResult:
    """Result of feedback loop monitoring."""

    loops: List[FeedbackLoopResult]
    n_intact: int = 0
    n_decoupled: int = 0
    n_inverted: int = 0
    runaway_risk_score: float = 0.0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_inverted == 0 and self.n_decoupled == 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_intact": self.n_intact,
            "n_decoupled": self.n_decoupled,
            "n_inverted": self.n_inverted,
            "runaway_risk_score": round(self.runaway_risk_score, 4),
            "summary": self.summary,
            "loops": [loop.to_dict() for loop in self.loops],
        }


class FeedbackMonitor:
    """Monitors feedback loop integrity.

    Usage::

        monitor = FeedbackMonitor(
            loops=[("HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_KRAS_SIGNALING_DN")],
        )
        result = monitor.check(pathway_scores)
    """

    def __init__(
        self,
        loops: Optional[List[Tuple[str, str]]] = None,
        intact_threshold: float = 0.3,
        inverted_threshold: float = -0.3,
        pathway_weights: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None,
    ):
        """Initialize the feedback monitor.

        Args:
            loops: List of (activator_pathway, inhibitor_pathway) pairs.
            intact_threshold: Min correlation for "intact" classification.
            inverted_threshold: Correlation below this = "inverted."
            pathway_weights: Weights for runaway risk scoring.
            seed: Random seed (unused, kept for API consistency).
        """
        self.loops = loops or []
        self.intact_threshold = intact_threshold
        self.inverted_threshold = inverted_threshold
        self.pathway_weights = pathway_weights or {}

    def check(self, pathway_scores: pd.DataFrame) -> FeedbackResult:
        """Check all feedback loops.

        Args:
            pathway_scores: DataFrame (cells x pathways).

        Returns:
            FeedbackResult with per-loop integrity scores.
        """
        results: List[FeedbackLoopResult] = []
        n_intact = 0
        n_decoupled = 0
        n_inverted = 0
        risk = 0.0

        for activator, inhibitor in self.loops:
            if activator not in pathway_scores.columns or inhibitor not in pathway_scores.columns:
                results.append(
                    FeedbackLoopResult(
                        activator=activator,
                        inhibitor=inhibitor,
                        integrity_score=0.0,
                        status=FeedbackStatus.UNKNOWN,
                        description="Missing pathway data",
                    )
                )
                continue

            act_vals = pathway_scores[activator].values
            inh_vals = pathway_scores[inhibitor].values

            if len(act_vals) < 3:
                fi = 0.0
            else:
                fi = float(np.corrcoef(act_vals, inh_vals)[0, 1])
                if np.isnan(fi):
                    fi = 0.0

            if fi >= self.intact_threshold:
                status = FeedbackStatus.INTACT
                n_intact += 1
            elif fi <= self.inverted_threshold:
                status = FeedbackStatus.INVERTED
                n_inverted += 1
                w = self.pathway_weights.get(activator, 1.0)
                risk += w * abs(fi)
            else:
                status = FeedbackStatus.DECOUPLED
                n_decoupled += 1
                w = self.pathway_weights.get(activator, 1.0)
                risk += w * 0.5

            results.append(
                FeedbackLoopResult(
                    activator=activator,
                    inhibitor=inhibitor,
                    integrity_score=fi,
                    status=status,
                    description=f"{activator} -> {inhibitor}: {status.value} (FI={fi:.3f})",
                )
            )

        summary = (
            f"Feedback: {len(results)} loops checked. "
            f"{n_intact} intact, {n_decoupled} decoupled, {n_inverted} inverted. "
            f"Runaway risk: {risk:.3f}. "
            f"{'PASS' if n_inverted == 0 and n_decoupled == 0 else 'FAIL'}."
        )

        logger.info("[QC Feedback] %s", summary)

        return FeedbackResult(
            loops=results,
            n_intact=n_intact,
            n_decoupled=n_decoupled,
            n_inverted=n_inverted,
            runaway_risk_score=risk,
            summary=summary,
        )
