"""
Active Tension Scoring (F3).

Quantifies the molecular tension from maintaining multiple incomplete
signaling cascades simultaneously. Each open signaling loop has
bioenergetic cost; cells with high tension are under molecular stress.

Formula:
    T(cell) = sum_p w_p * (1 - CCR_p) * A_p

    where:
        p     = pathway
        w_p   = pathway importance weight (default: 1.0)
        CCR_p = cascade completion ratio (0 = stalled, 1 = complete)
        A_p   = initiator activation level (only count started cascades)

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class TensionResult:
    """Result of tension scoring for a batch."""

    per_cell_tension: np.ndarray
    per_pathway_contribution: Dict[str, float]
    mean_tension: float = 0.0
    median_tension: float = 0.0
    max_tension: float = 0.0
    gate_passed: bool = True
    threshold: float = 0.0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.gate_passed

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.gate_passed,
            "mean_tension": round(self.mean_tension, 4),
            "median_tension": round(self.median_tension, 4),
            "max_tension": round(self.max_tension, 4),
            "threshold": round(self.threshold, 4),
            "per_pathway_contribution": {
                k: round(v, 4) for k, v in self.per_pathway_contribution.items()
            },
            "summary": self.summary,
        }


class TensionScorer:
    """Scores molecular tension from incomplete cascades.

    Usage::

        from pathway_subtyping.qc.cascade import CascadeAnalyzer

        cascade_result = analyzer.analyze(expression, pathways)
        scorer = TensionScorer(pathway_weights={"MAPK": 2.0})
        tension = scorer.score_batch(cascade_result)
    """

    def __init__(
        self,
        pathway_weights: Optional[Dict[str, float]] = None,
        tension_threshold: float = 1.0,
        seed: Optional[int] = None,
    ):
        """Initialize the tension scorer.

        Args:
            pathway_weights: Per-pathway importance weights. Higher weight
                means more tension contribution when that pathway stalls.
            tension_threshold: Maximum acceptable mean tension score.
            seed: Random seed (unused, kept for API consistency).
        """
        self.pathway_weights = pathway_weights or {}
        self.tension_threshold = tension_threshold

    def score_sample(
        self,
        per_cell_ccr: pd.DataFrame,
        upstream_scores: Optional[Dict[str, np.ndarray]] = None,
    ) -> np.ndarray:
        """Compute tension score per cell.

        Args:
            per_cell_ccr: DataFrame (cells x pathways) of cascade completion ratios.
            upstream_scores: Optional dict of pathway -> per-cell upstream activation.
                If None, uses (1 - CCR) directly without activation weighting.

        Returns:
            1D array of per-cell tension scores.
        """
        n_cells = len(per_cell_ccr)
        tension = np.zeros(n_cells)

        for pathway in per_cell_ccr.columns:
            ccr = per_cell_ccr[pathway].values
            w = self.pathway_weights.get(pathway, 1.0)
            incompleteness = 1.0 - ccr

            if upstream_scores is not None and pathway in upstream_scores:
                activation = np.maximum(upstream_scores[pathway], 0.0)
            else:
                activation = np.ones(n_cells)

            tension += w * incompleteness * activation

        return tension

    def score_batch(
        self,
        cascade_result: Any,  # CascadeResult — Any to avoid circular import
    ) -> TensionResult:
        """Compute tension for a full batch from cascade analysis results.

        Args:
            cascade_result: CascadeResult from CascadeAnalyzer.analyze().

        Returns:
            TensionResult with per-cell and aggregate tension metrics.
        """
        if cascade_result.per_cell_ccr is None or cascade_result.per_cell_ccr.empty:
            return TensionResult(
                per_cell_tension=np.array([]),
                per_pathway_contribution={},
                summary="No cascade data available for tension scoring.",
            )

        ccr_df = cascade_result.per_cell_ccr
        n_cells = len(ccr_df)
        tension = np.zeros(n_cells)
        per_pathway: Dict[str, float] = {}

        for pathway in ccr_df.columns:
            ccr = ccr_df[pathway].values
            w = self.pathway_weights.get(pathway, 1.0)
            incompleteness = 1.0 - ccr
            contribution = w * incompleteness
            tension += contribution
            per_pathway[pathway] = float(np.mean(contribution))

        mean_t = float(np.mean(tension))
        median_t = float(np.median(tension))
        max_t = float(np.max(tension)) if len(tension) > 0 else 0.0
        gate_passed = mean_t <= self.tension_threshold

        summary = (
            f"Tension: mean={mean_t:.3f}, median={median_t:.3f}, "
            f"max={max_t:.3f}, threshold={self.tension_threshold:.3f}. "
            f"{'PASS' if gate_passed else 'FAIL'}."
        )

        logger.info("[QC Tension] %s", summary)

        return TensionResult(
            per_cell_tension=tension,
            per_pathway_contribution=per_pathway,
            mean_tension=mean_t,
            median_tension=median_t,
            max_tension=max_t,
            gate_passed=gate_passed,
            threshold=self.tension_threshold,
            summary=summary,
        )

    def tension_gate(
        self,
        tension_scores: np.ndarray,
        threshold: Optional[float] = None,
    ) -> bool:
        """Binary QC gate: PASS if mean tension below threshold.

        Args:
            tension_scores: Per-cell tension scores.
            threshold: Override threshold (default: self.tension_threshold).

        Returns:
            True if batch passes tension gate.
        """
        t = threshold if threshold is not None else self.tension_threshold
        return bool(np.mean(tension_scores) <= t)
