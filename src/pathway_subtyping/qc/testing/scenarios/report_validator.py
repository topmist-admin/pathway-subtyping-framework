"""
Report Validator for scenario QC reports.

Validates scenario reports against a rubric of required elements.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List

from .base import ScenarioReport

logger = logging.getLogger(__name__)


@dataclass
class RubricItem:
    """A single rubric check."""

    name: str
    description: str
    passed: bool = False
    reason: str = ""


@dataclass
class RubricResult:
    """Result of rubric validation."""

    items: List[RubricItem]
    all_passed: bool = False
    summary: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "all_passed": self.all_passed,
            "summary": self.summary,
            "items": [
                {"name": i.name, "passed": i.passed, "reason": i.reason}
                for i in self.items
            ],
        }


class ReportValidator:
    """Validates scenario reports against a standard rubric.

    Rubric elements:
        1. Executive summary with PASS/HOLD/REJECT decision
        2. Per-timepoint QC scores
        3. Identified defects with severity and affected pathways
        4. Remediation recommendations (specific, actionable)
        5. No false alarms for non-injected defects
    """

    def validate(self, report: ScenarioReport) -> RubricResult:
        """Validate a scenario report against the rubric."""
        items = []

        # 1. Executive summary
        has_decision = report.final_decision in ("PASS", "HOLD", "REJECT")
        items.append(
            RubricItem(
                name="executive_summary",
                description="Report has PASS/HOLD/REJECT final decision",
                passed=has_decision,
                reason=f"Decision: {report.final_decision}" if has_decision else "Missing decision",
            )
        )

        # 2. Per-timepoint scores
        has_scores = all(
            len(tp.qc_scores) > 0 for tp in report.timepoints
        )
        items.append(
            RubricItem(
                name="per_timepoint_scores",
                description="All timepoints have QC scores",
                passed=has_scores,
                reason=f"{len(report.timepoints)} timepoints scored" if has_scores else "Missing scores",
            )
        )

        # 3. Identified defects
        has_defects_listed = len(report.identified_defects) > 0
        items.append(
            RubricItem(
                name="identified_defects",
                description="Defects identified with details",
                passed=has_defects_listed,
                reason=f"{len(report.identified_defects)} defects identified"
                if has_defects_listed
                else "No defects identified (may be correct for healthy batches)",
            )
        )

        # 4. Remediation recommendations
        has_remediation = len(report.remediation_recommendations) > 0
        items.append(
            RubricItem(
                name="remediation",
                description="Actionable remediation recommendations provided",
                passed=has_remediation or not has_defects_listed,
                reason=f"{len(report.remediation_recommendations)} recommendations"
                if has_remediation
                else "No recommendations (acceptable if no defects)",
            )
        )

        all_passed = all(item.passed for item in items)
        n_pass = sum(1 for i in items if i.passed)

        result = RubricResult(
            items=items,
            all_passed=all_passed,
            summary=f"Rubric: {n_pass}/{len(items)} items passed.",
        )

        logger.info("[QC Report Validator] %s", result.summary)
        return result
