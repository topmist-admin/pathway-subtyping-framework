"""
KG regression report for v0.6 F3.

Runs a suite of KG-dependent scoring functions against two KnowledgeGraph
instances and flags whichever pathway/feature outputs move by more than a
configurable tolerance. Used to validate that a KG upgrade (v0.5 -> v0.6)
does not silently shift downstream QC feature outputs beyond an acceptable
threshold (roadmap target: < 5%).

The report is data-agnostic: the caller supplies scoring callables and
a benchmark payload (e.g. a TCGA-COAD expression matrix or a pathway-
score DataFrame). The regression harness itself ships no benchmark data;
users with matched cohorts drop in their own.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional

from .builder import KnowledgeGraph

logger = logging.getLogger(__name__)


ScoreFn = Callable[[KnowledgeGraph, Any], Dict[str, float]]


@dataclass
class ScoreDelta:
    """Per-score delta between two KG snapshots.

    Attributes:
        score_name: Arbitrary string identifying the score (e.g. ``F1_cascade``).
        v1_value: Score on the baseline KG.
        v2_value: Score on the new KG.
        absolute_delta: ``v2_value - v1_value``.
        relative_delta: ``(v2 - v1) / max(abs(v1), eps)``.
        flagged: True if ``|relative_delta| > tolerance``.
    """

    score_name: str
    v1_value: float
    v2_value: float
    absolute_delta: float
    relative_delta: float
    flagged: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "score_name": self.score_name,
            "v1_value": self.v1_value,
            "v2_value": self.v2_value,
            "absolute_delta": self.absolute_delta,
            "relative_delta": self.relative_delta,
            "flagged": self.flagged,
        }


@dataclass
class KGRegressionReport:
    """Summary of KG-upgrade impact on one or more scoring functions.

    Attributes:
        tolerance: Relative-delta threshold beyond which a score is flagged.
        deltas: List of ScoreDelta, one per (score_fn output) element.
        passed: True iff nothing is flagged.
    """

    tolerance: float
    deltas: List[ScoreDelta] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not any(d.flagged for d in self.deltas)

    @property
    def flagged_scores(self) -> List[ScoreDelta]:
        return [d for d in self.deltas if d.flagged]

    def summary(self) -> str:
        return (
            f"KGRegressionReport(n={len(self.deltas)}, "
            f"flagged={len(self.flagged_scores)}, "
            f"passed={self.passed}, tolerance={self.tolerance})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "tolerance": self.tolerance,
            "passed": self.passed,
            "n_scores": len(self.deltas),
            "n_flagged": len(self.flagged_scores),
            "deltas": [d.to_dict() for d in self.deltas],
        }


def run_kg_regression(
    v1: KnowledgeGraph,
    v2: KnowledgeGraph,
    score_fns: Iterable[ScoreFn],
    benchmark_inputs: Iterable[Any],
    tolerance: float = 0.05,
    eps: float = 1e-9,
) -> KGRegressionReport:
    """Execute every ``score_fn`` on both KGs for every input; diff outputs.

    Args:
        v1: Baseline knowledge graph.
        v2: New knowledge graph.
        score_fns: Callables with signature ``f(kg, input) -> Dict[str, float]``.
        benchmark_inputs: Iterable of inputs passed to each score_fn.
        tolerance: Relative-delta threshold (roadmap default 0.05 = 5%).
        eps: Floor on |v1_value| when computing relative_delta.

    Returns:
        KGRegressionReport with per-score deltas and pass/fail summary.
    """
    score_fns = list(score_fns)
    benchmark_inputs = list(benchmark_inputs)
    if not score_fns:
        raise ValueError("score_fns must be non-empty")
    if not benchmark_inputs:
        raise ValueError("benchmark_inputs must be non-empty")

    report = KGRegressionReport(tolerance=tolerance)

    for score_idx, score_fn in enumerate(score_fns):
        for input_idx, benchmark_input in enumerate(benchmark_inputs):
            v1_out = score_fn(v1, benchmark_input)
            v2_out = score_fn(v2, benchmark_input)
            if not isinstance(v1_out, dict) or not isinstance(v2_out, dict):
                raise TypeError(
                    f"score_fn index {score_idx} must return Dict[str, float]"
                )
            keys = set(v1_out) | set(v2_out)
            for key in sorted(keys):
                v1_v = float(v1_out.get(key, 0.0))
                v2_v = float(v2_out.get(key, 0.0))
                abs_d = v2_v - v1_v
                rel_d = abs_d / max(abs(v1_v), eps)
                flagged = abs(rel_d) > tolerance
                report.deltas.append(
                    ScoreDelta(
                        score_name=f"fn{score_idx}:input{input_idx}:{key}",
                        v1_value=v1_v,
                        v2_value=v2_v,
                        absolute_delta=abs_d,
                        relative_delta=rel_d,
                        flagged=flagged,
                    )
                )

    logger.info("[KGRegression] %s", report.summary())
    return report
