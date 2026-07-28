"""
Pathway Dosage & Stoichiometry Analysis (F8).

Classifies pathway activation into three states relative to therapeutic windows:
    - UNDER: activation < min -> insufficient engineering
    - IN_RANGE: min <= activation <= max -> target achieved
    - OVER: activation > max -> over-activation risk (e.g., cytokine storm)

Also checks stoichiometric ratios between pathway pairs.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class DosageState(Enum):
    """Classification of pathway activation dosage."""

    UNDER = "under"
    IN_RANGE = "in_range"
    OVER = "over"
    NOT_ASSESSED = "not_assessed"


@dataclass
class PathwayDosageResult:
    """Dosage classification for a single pathway."""

    pathway: str
    state: DosageState
    mean_activation: float
    window_min: float
    window_max: float
    n_cells_under: int = 0
    n_cells_in_range: int = 0
    n_cells_over: int = 0
    n_cells_total: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway": self.pathway,
            "state": self.state.value,
            "mean_activation": round(self.mean_activation, 4),
            "window": [round(self.window_min, 4), round(self.window_max, 4)],
            "n_under": self.n_cells_under,
            "n_in_range": self.n_cells_in_range,
            "n_over": self.n_cells_over,
        }


@dataclass
class StoichiometryResult:
    """Result of stoichiometric ratio check for a pathway pair."""

    pathway_a: str
    pathway_b: str
    observed_ratio: float
    expected_ratio: float
    tolerance: float
    within_tolerance: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway_a": self.pathway_a,
            "pathway_b": self.pathway_b,
            "observed_ratio": round(self.observed_ratio, 4),
            "expected_ratio": round(self.expected_ratio, 4),
            "tolerance": round(self.tolerance, 4),
            "within_tolerance": self.within_tolerance,
        }


@dataclass
class DosageAnalysisResult:
    """Result of full dosage and stoichiometry analysis."""

    pathway_dosages: List[PathwayDosageResult]
    stoichiometry_results: List[StoichiometryResult]
    n_under: int = 0
    n_in_range: int = 0
    n_over: int = 0
    n_ratio_violations: int = 0
    toxicity_risk_score: float = 0.0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_over == 0 and self.n_ratio_violations == 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_under": self.n_under,
            "n_in_range": self.n_in_range,
            "n_over": self.n_over,
            "n_ratio_violations": self.n_ratio_violations,
            "toxicity_risk_score": round(self.toxicity_risk_score, 4),
            "summary": self.summary,
            "pathways": [d.to_dict() for d in self.pathway_dosages],
            "stoichiometry": [s.to_dict() for s in self.stoichiometry_results],
        }


class DosageAnalyzer:
    """Analyzes pathway activation dosage against therapeutic windows.

    Usage::

        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8),
                "HALLMARK_APOPTOSIS": (0.0, 0.3),
            },
        )
        result = analyzer.analyze(pathway_scores)
    """

    def __init__(
        self,
        therapeutic_windows: Optional[Dict[str, Tuple[float, float]]] = None,
        pathway_ratios: Optional[List[Dict[str, Any]]] = None,
        toxicity_weights: Optional[Dict[str, float]] = None,
        max_violations: int = 0,
        seed: Optional[int] = None,
    ):
        """Initialize the dosage analyzer.

        Args:
            therapeutic_windows: Per-pathway (min, max) acceptable activation range.
            pathway_ratios: List of dicts with keys: pathway_a, pathway_b,
                expected_ratio, tolerance.
            toxicity_weights: Per-pathway weight for toxicity risk scoring.
                Higher = more dangerous when over-activated.
            max_violations: Max number of pathway violations before gate fails.
            seed: Random seed (unused, kept for API consistency).
        """
        self.therapeutic_windows = therapeutic_windows or {}
        self.pathway_ratios = pathway_ratios or []
        self.toxicity_weights = toxicity_weights or {}
        self.max_violations = max_violations

    def define_therapeutic_window(self, pathway: str, window_min: float, window_max: float) -> None:
        """Define or update a therapeutic window for a pathway."""
        self.therapeutic_windows[pathway] = (window_min, window_max)

    def analyze(self, pathway_scores: pd.DataFrame) -> DosageAnalysisResult:
        """Analyze dosage for all pathways with defined windows.

        Args:
            pathway_scores: DataFrame (cells x pathways) of pathway scores.

        Returns:
            DosageAnalysisResult with per-pathway dosage and stoichiometry.
        """
        n_cells = len(pathway_scores)
        dosages: List[PathwayDosageResult] = []
        n_under_total = 0
        n_in_range_total = 0
        n_over_total = 0
        toxicity_score = 0.0

        for pathway, (w_min, w_max) in self.therapeutic_windows.items():
            if pathway not in pathway_scores.columns:
                continue

            scores = pathway_scores[pathway].values
            mean_act = float(np.mean(scores))

            n_under = int(np.sum(scores < w_min))
            n_in_range = int(np.sum((scores >= w_min) & (scores <= w_max)))
            n_over = int(np.sum(scores > w_max))

            # Classify based on batch mean
            if mean_act < w_min:
                state = DosageState.UNDER
                n_under_total += 1
            elif mean_act > w_max:
                state = DosageState.OVER
                n_over_total += 1
                # Toxicity contribution
                weight = self.toxicity_weights.get(pathway, 1.0)
                overshoot = (mean_act - w_max) / max(w_max, 0.01)
                toxicity_score += weight * overshoot
            else:
                state = DosageState.IN_RANGE
                n_in_range_total += 1

            dosages.append(
                PathwayDosageResult(
                    pathway=pathway,
                    state=state,
                    mean_activation=mean_act,
                    window_min=w_min,
                    window_max=w_max,
                    n_cells_under=n_under,
                    n_cells_in_range=n_in_range,
                    n_cells_over=n_over,
                    n_cells_total=n_cells,
                )
            )

        # Stoichiometry checks
        stoich_results: List[StoichiometryResult] = []
        n_ratio_violations = 0

        for ratio_spec in self.pathway_ratios:
            pw_a = ratio_spec["pathway_a"]
            pw_b = ratio_spec["pathway_b"]
            expected = ratio_spec["expected_ratio"]
            tolerance = ratio_spec["tolerance"]

            if pw_a not in pathway_scores.columns or pw_b not in pathway_scores.columns:
                continue

            mean_a = float(np.mean(pathway_scores[pw_a].values))
            mean_b = float(np.mean(pathway_scores[pw_b].values))
            observed = mean_a / max(abs(mean_b), 1e-6)

            within = abs(observed - expected) <= tolerance

            if not within:
                n_ratio_violations += 1

            stoich_results.append(
                StoichiometryResult(
                    pathway_a=pw_a,
                    pathway_b=pw_b,
                    observed_ratio=observed,
                    expected_ratio=expected,
                    tolerance=tolerance,
                    within_tolerance=within,
                )
            )

        # Sort dosages: OVER first (most dangerous), then UNDER, then IN_RANGE
        dosages.sort(
            key=lambda d: (
                0 if d.state == DosageState.OVER else 1 if d.state == DosageState.UNDER else 2
            )
        )

        summary = (
            f"Dosage: {n_in_range_total} in-range, {n_under_total} under, "
            f"{n_over_total} over. "
            f"Stoichiometry: {n_ratio_violations} violation(s). "
            f"Toxicity risk: {toxicity_score:.3f}. "
            f"{'PASS' if n_over_total == 0 and n_ratio_violations == 0 else 'FAIL'}."
        )

        logger.info("[QC Dosage] %s", summary)

        return DosageAnalysisResult(
            pathway_dosages=dosages,
            stoichiometry_results=stoich_results,
            n_under=n_under_total,
            n_in_range=n_in_range_total,
            n_over=n_over_total,
            n_ratio_violations=n_ratio_violations,
            toxicity_risk_score=toxicity_score,
            summary=summary,
        )
