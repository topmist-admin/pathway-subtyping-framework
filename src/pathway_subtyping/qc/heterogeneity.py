"""
Batch Heterogeneity Profiling (F7).

Detects non-conforming subpopulations within a manufacturing batch by
measuring each cell's distance from a target pathway profile specification.

Metrics:
    - Per-cell conformity score (distance from target profile)
    - Heterogeneity index: 1 - (n_conforming / n_total)
    - Bimodality coefficient for distribution shape
    - Density-based subpopulation detection

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from sklearn.cluster import DBSCAN

logger = logging.getLogger(__name__)


@dataclass
class SubpopulationInfo:
    """Information about a detected non-conforming subpopulation."""

    cluster_id: int
    n_cells: int
    mean_conformity_score: float
    deviating_pathways: List[str] = field(default_factory=list)
    deviation_directions: Dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cluster_id": self.cluster_id,
            "n_cells": self.n_cells,
            "mean_conformity_score": round(self.mean_conformity_score, 4),
            "deviating_pathways": self.deviating_pathways,
            "deviation_directions": self.deviation_directions,
        }


@dataclass
class HeterogeneityResult:
    """Result of batch heterogeneity analysis."""

    heterogeneity_index: float
    bimodality_coefficient: float
    n_conforming: int
    n_nonconforming: int
    n_total: int
    conformity_scores: np.ndarray
    subpopulations: List[SubpopulationInfo] = field(default_factory=list)
    conformity_threshold: float = 0.0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.heterogeneity_index < 0.15  # Default: <15% non-conforming

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "heterogeneity_index": round(self.heterogeneity_index, 4),
            "bimodality_coefficient": round(self.bimodality_coefficient, 4),
            "n_conforming": self.n_conforming,
            "n_nonconforming": self.n_nonconforming,
            "n_total": self.n_total,
            "n_subpopulations": len(self.subpopulations),
            "summary": self.summary,
            "subpopulations": [s.to_dict() for s in self.subpopulations],
        }


class HeterogeneityProfiler:
    """Profiles batch heterogeneity against a target specification.

    Usage::

        profiler = HeterogeneityProfiler(
            target_profile={"HALLMARK_INFLAMMATORY_RESPONSE": 0.6, ...},
            conformity_threshold=2.0,
        )
        result = profiler.profile(pathway_scores)
    """

    def __init__(
        self,
        target_profile: Optional[Dict[str, float]] = None,
        conformity_threshold: float = 2.0,
        heterogeneity_threshold: float = 0.15,
        dbscan_eps: float = 1.5,
        dbscan_min_samples: int = 5,
        seed: Optional[int] = None,
    ):
        """Initialize the heterogeneity profiler.

        Args:
            target_profile: Target pathway activation scores. If None, uses
                batch mean as the specification.
            conformity_threshold: Max distance from target to be "conforming."
            heterogeneity_threshold: Max fraction of non-conforming cells to pass.
            dbscan_eps: DBSCAN epsilon for subpopulation detection.
            dbscan_min_samples: DBSCAN min_samples parameter.
            seed: Random seed (unused, kept for API consistency).
        """
        self.target_profile = target_profile
        self.conformity_threshold = conformity_threshold
        self.heterogeneity_threshold = heterogeneity_threshold
        self.dbscan_eps = dbscan_eps
        self.dbscan_min_samples = dbscan_min_samples

    def set_specification(self, target_profile: Dict[str, float]) -> None:
        """Set the target pathway profile."""
        self.target_profile = target_profile

    def profile(self, pathway_scores: pd.DataFrame) -> HeterogeneityResult:
        """Profile batch heterogeneity.

        Args:
            pathway_scores: DataFrame (cells x pathways) of pathway scores.

        Returns:
            HeterogeneityResult with conformity scores and subpopulation info.
        """
        n_total = len(pathway_scores)

        # Build target vector
        if self.target_profile is not None:
            target = np.array([
                self.target_profile.get(pw, 0.0) for pw in pathway_scores.columns
            ])
        else:
            target = pathway_scores.values.mean(axis=0)

        # Per-cell conformity score (Euclidean distance from target)
        scores_matrix = pathway_scores.values
        distances = np.sqrt(np.sum((scores_matrix - target) ** 2, axis=1))

        # Classify conforming vs non-conforming
        conforming = distances <= self.conformity_threshold
        n_conforming = int(conforming.sum())
        n_nonconforming = n_total - n_conforming

        # Heterogeneity index
        h_index = 1.0 - (n_conforming / n_total) if n_total > 0 else 0.0

        # Bimodality coefficient
        bc = self._bimodality_coefficient(distances)

        # Detect subpopulations among non-conforming cells
        subpopulations: List[SubpopulationInfo] = []
        if n_nonconforming >= self.dbscan_min_samples:
            nonconf_indices = np.where(~conforming)[0]
            nonconf_scores = scores_matrix[nonconf_indices]
            subpopulations = self._detect_subpopulations(
                nonconf_scores, nonconf_indices, pathway_scores.columns.tolist(), target
            )

        summary = (
            f"Heterogeneity: {h_index:.1%} non-conforming "
            f"({n_nonconforming}/{n_total}), "
            f"BC={bc:.3f}, "
            f"{len(subpopulations)} subpopulation(s). "
            f"{'PASS' if h_index < self.heterogeneity_threshold else 'FAIL'}."
        )

        logger.info("[QC Heterogeneity] %s", summary)

        return HeterogeneityResult(
            heterogeneity_index=h_index,
            bimodality_coefficient=bc,
            n_conforming=n_conforming,
            n_nonconforming=n_nonconforming,
            n_total=n_total,
            conformity_scores=distances,
            subpopulations=subpopulations,
            conformity_threshold=self.conformity_threshold,
            summary=summary,
        )

    def _bimodality_coefficient(self, values: np.ndarray) -> float:
        """Compute bimodality coefficient (BC).

        BC > 0.555 suggests bimodal distribution.
        BC = (skewness^2 + 1) / (kurtosis + 3 * (n-1)^2 / ((n-2)(n-3)))
        """
        n = len(values)
        if n < 4:
            return 0.0
        skew = float(scipy_stats.skew(values))
        kurt = float(scipy_stats.kurtosis(values, fisher=True))
        adjustment = 3.0 * ((n - 1) ** 2) / ((n - 2) * (n - 3))
        denom = kurt + adjustment
        if denom == 0:
            return 0.0
        return (skew ** 2 + 1) / denom

    def _detect_subpopulations(
        self,
        nonconf_scores: np.ndarray,
        nonconf_indices: np.ndarray,
        pathway_names: List[str],
        target: np.ndarray,
    ) -> List[SubpopulationInfo]:
        """Detect subpopulations among non-conforming cells using DBSCAN."""
        clustering = DBSCAN(
            eps=self.dbscan_eps, min_samples=self.dbscan_min_samples
        ).fit(nonconf_scores)

        labels = clustering.labels_
        unique_labels = set(labels)
        unique_labels.discard(-1)  # Noise

        subpops: List[SubpopulationInfo] = []
        for label in sorted(unique_labels):
            mask = labels == label
            cluster_scores = nonconf_scores[mask]
            n_cells = int(mask.sum())

            # Mean distance from target for this cluster
            dists = np.sqrt(np.sum((cluster_scores - target) ** 2, axis=1))
            mean_conf = float(np.mean(dists))

            # Find deviating pathways
            cluster_mean = cluster_scores.mean(axis=0)
            deviations = cluster_mean - target
            top_indices = np.argsort(np.abs(deviations))[::-1][:3]
            deviating = [pathway_names[i] for i in top_indices]
            directions = {
                pathway_names[i]: "elevated" if deviations[i] > 0 else "suppressed"
                for i in top_indices
            }

            subpops.append(
                SubpopulationInfo(
                    cluster_id=int(label),
                    n_cells=n_cells,
                    mean_conformity_score=mean_conf,
                    deviating_pathways=deviating,
                    deviation_directions=directions,
                )
            )

        return subpops
