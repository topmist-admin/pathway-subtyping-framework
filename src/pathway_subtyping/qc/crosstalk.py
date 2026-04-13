"""
Pathway Crosstalk & Interference Detection (F9).

Detects resource competition at shared pathway nodes. When multiple pathways
share intermediate transducers, one pathway can starve another.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class InterferenceEdge:
    """Interference relationship between two pathways."""

    pathway_a: str
    pathway_b: str
    competition_score: float
    shared_genes: List[str] = field(default_factory=list)
    dominant: str = ""  # Which pathway dominates shared nodes

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway_a": self.pathway_a,
            "pathway_b": self.pathway_b,
            "competition_score": round(self.competition_score, 4),
            "n_shared_genes": len(self.shared_genes),
            "dominant": self.dominant,
        }


@dataclass
class CrosstalkResult:
    """Result of crosstalk detection."""

    interference_edges: List[InterferenceEdge]
    n_significant: int = 0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_significant == 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_significant": self.n_significant,
            "summary": self.summary,
            "edges": [e.to_dict() for e in self.interference_edges],
        }


class CrosstalkDetector:
    """Detects pathway crosstalk via shared node competition.

    Usage::

        detector = CrosstalkDetector(kg)
        result = detector.detect(expression, pathways=["PW_A", "PW_B"])
    """

    def __init__(
        self,
        kg: Any,
        competition_threshold: float = 0.3,
        seed: Optional[int] = None,
    ):
        self.kg = kg
        self.competition_threshold = competition_threshold

    def detect(
        self,
        expression: pd.DataFrame,
        pathways: List[str],
    ) -> CrosstalkResult:
        """Detect crosstalk between all pathway pairs.

        Args:
            expression: Expression matrix (cells x genes).
            pathways: List of pathway node IDs to check.

        Returns:
            CrosstalkResult with interference edges.
        """
        edges: List[InterferenceEdge] = []
        n_sig = 0

        for i in range(len(pathways)):
            for j in range(i + 1, len(pathways)):
                pw_a, pw_b = pathways[i], pathways[j]
                shared = self.kg.get_shared_genes(pw_a, pw_b)
                if not shared:
                    continue

                score, dominant = self._compute_competition(
                    expression, pw_a, pw_b, shared
                )

                edge = InterferenceEdge(
                    pathway_a=pw_a,
                    pathway_b=pw_b,
                    competition_score=score,
                    shared_genes=shared,
                    dominant=dominant,
                )
                edges.append(edge)

                if abs(score) > self.competition_threshold:
                    n_sig += 1

        summary = (
            f"Crosstalk: {len(edges)} pathway pairs with shared nodes, "
            f"{n_sig} significant. "
            f"{'PASS' if n_sig == 0 else 'FAIL'}."
        )

        logger.info("[QC Crosstalk] %s", summary)

        return CrosstalkResult(
            interference_edges=edges,
            n_significant=n_sig,
            summary=summary,
        )

    def _compute_competition(
        self,
        expression: pd.DataFrame,
        pw_a: str,
        pw_b: str,
        shared: List[str],
    ) -> Tuple[float, str]:
        """Compute competition score between two pathways at shared nodes."""
        genes_a = set(self.kg.get_pathway_genes(pw_a)) - set(shared)
        genes_b = set(self.kg.get_pathway_genes(pw_b)) - set(shared)

        shared_present = [g for g in shared if g in expression.columns]
        a_present = [g for g in genes_a if g in expression.columns]
        b_present = [g for g in genes_b if g in expression.columns]

        if not shared_present or (not a_present and not b_present):
            return 0.0, ""

        shared_mean = expression[shared_present].values.mean()
        a_mean = expression[a_present].values.mean() if a_present else 0.0
        b_mean = expression[b_present].values.mean() if b_present else 0.0

        # Positive = shared nodes favor A; negative = favor B
        score = a_mean - b_mean
        dominant = pw_a if score > 0 else pw_b if score < 0 else ""
        return float(score), dominant
