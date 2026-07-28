"""
Resolution Gates (F4).

Manufacturing release criterion: has this batch completed the pathway
transitions required by the engineering specification?

Integrates signals from both Workstream A (F6/F7/F8) and Workstream B
(F1/F3/F5) into a unified RELEASE / HOLD / REJECT decision.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class ReleaseDecision(Enum):
    """Batch release decision."""

    RELEASE = "release"
    HOLD = "hold"
    REJECT = "reject"


@dataclass
class GateCheck:
    """Result of a single gate criterion check."""

    name: str
    passed: bool
    value: float
    threshold: float
    details: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "passed": self.passed,
            "value": round(self.value, 4),
            "threshold": round(self.threshold, 4),
            "details": self.details,
        }


@dataclass
class ResolutionResult:
    """Result of resolution gate evaluation."""

    decision: ReleaseDecision
    gate_checks: List[GateCheck]
    n_passed: int = 0
    n_failed: int = 0
    batch_pass_rate: float = 0.0
    justification: str = ""

    @property
    def passed(self) -> bool:
        return self.decision == ReleaseDecision.RELEASE

    def to_dict(self) -> Dict[str, Any]:
        return {
            "decision": self.decision.value,
            "passed": self.passed,
            "n_passed": self.n_passed,
            "n_failed": self.n_failed,
            "batch_pass_rate": round(self.batch_pass_rate, 4),
            "justification": self.justification,
            "gates": [g.to_dict() for g in self.gate_checks],
        }


class ResolutionGate:
    """Manufacturing release gate integrating all QC signals.

    Usage::

        gate = ResolutionGate(
            min_cascade_ccr=0.7,
            max_tension=1.0,
            max_heterogeneity=0.15,
            release_threshold=0.9,
        )
        result = gate.evaluate(
            cascade_result=cascade_result,
            tension_result=tension_result,
            heterogeneity_result=het_result,
            dosage_result=dosage_result,
            offtarget_result=offtarget_result,
        )
    """

    def __init__(
        self,
        min_cascade_ccr: float = 0.7,
        max_tension: float = 1.0,
        max_heterogeneity: float = 0.15,
        max_offtarget: int = 0,
        max_dosage_violations: int = 0,
        release_threshold: float = 0.9,
        hold_threshold: float = 0.7,
        seed: Optional[int] = None,
    ):
        self.min_cascade_ccr = min_cascade_ccr
        self.max_tension = max_tension
        self.max_heterogeneity = max_heterogeneity
        self.max_offtarget = max_offtarget
        self.max_dosage_violations = max_dosage_violations
        self.release_threshold = release_threshold
        self.hold_threshold = hold_threshold

    def evaluate(
        self,
        cascade_result: Optional[Any] = None,
        tension_result: Optional[Any] = None,
        heterogeneity_result: Optional[Any] = None,
        dosage_result: Optional[Any] = None,
        offtarget_result: Optional[Any] = None,
    ) -> ResolutionResult:
        """Evaluate all available QC results against gate criteria.

        Pass any subset of results — gates for missing results are skipped.
        """
        checks: List[GateCheck] = []

        # F1: Cascade completion
        if cascade_result is not None:
            ccr = cascade_result.mean_completion_ratio
            checks.append(
                GateCheck(
                    name="cascade_ccr",
                    passed=ccr >= self.min_cascade_ccr,
                    value=ccr,
                    threshold=self.min_cascade_ccr,
                    details=f"{cascade_result.n_incomplete} incomplete pathway(s)",
                )
            )

        # F3: Tension
        if tension_result is not None:
            t = tension_result.mean_tension
            checks.append(
                GateCheck(
                    name="tension",
                    passed=t <= self.max_tension,
                    value=t,
                    threshold=self.max_tension,
                    details=f"Max cell tension: {tension_result.max_tension:.3f}",
                )
            )

        # F6: Off-target
        if offtarget_result is not None:
            n_viol = offtarget_result.n_excluded_violations + offtarget_result.n_off_target
            checks.append(
                GateCheck(
                    name="offtarget",
                    passed=n_viol <= self.max_offtarget,
                    value=float(n_viol),
                    threshold=float(self.max_offtarget),
                    details=f"Safety flag: {offtarget_result.safety_flag}",
                )
            )

        # F7: Heterogeneity
        if heterogeneity_result is not None:
            h = heterogeneity_result.heterogeneity_index
            checks.append(
                GateCheck(
                    name="heterogeneity",
                    passed=h <= self.max_heterogeneity,
                    value=h,
                    threshold=self.max_heterogeneity,
                    details=f"{heterogeneity_result.n_nonconforming} non-conforming cells",
                )
            )

        # F8: Dosage
        if dosage_result is not None:
            n_over = dosage_result.n_over
            checks.append(
                GateCheck(
                    name="dosage",
                    passed=n_over <= self.max_dosage_violations,
                    value=float(n_over),
                    threshold=float(self.max_dosage_violations),
                    details=f"Toxicity risk: {dosage_result.toxicity_risk_score:.3f}",
                )
            )

        n_passed = sum(1 for c in checks if c.passed)
        n_failed = len(checks) - n_passed
        pass_rate = n_passed / len(checks) if checks else 1.0

        if pass_rate >= self.release_threshold:
            decision = ReleaseDecision.RELEASE
        elif pass_rate >= self.hold_threshold:
            decision = ReleaseDecision.HOLD
        else:
            decision = ReleaseDecision.REJECT

        failed_names = [c.name for c in checks if not c.passed]
        justification = (
            f"{decision.value.upper()}: {n_passed}/{len(checks)} gates passed "
            f"({pass_rate:.0%}). "
            + (f"Failed: {', '.join(failed_names)}." if failed_names else "All gates passed.")
        )

        logger.info("[QC Resolution] %s", justification)

        return ResolutionResult(
            decision=decision,
            gate_checks=checks,
            n_passed=n_passed,
            n_failed=n_failed,
            batch_pass_rate=pass_rate,
            justification=justification,
        )
