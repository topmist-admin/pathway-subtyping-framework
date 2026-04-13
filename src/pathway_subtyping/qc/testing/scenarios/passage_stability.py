"""
Scenario C: Passage Stability Over 20 Passages.

20 timepoints (Passage 1-20), 200 cells per passage.
Gradual drift injected P6-10, feedback breakdown P11-15, compounding failure P16-20.

Research use only. Not for clinical decision-making.
"""

import logging
from typing import Dict, List, Optional

from ..simulator import (
    FeedbackBreakType,
    InjectedBatch,
    ManufacturingSimulator,
    ManufacturingSpec,
)
from .base import BaseScenario

logger = logging.getLogger(__name__)


class PassageStabilityScenario(BaseScenario):
    """20-passage stability test with progressive drift and feedback breakdown."""

    N_CELLS = 200

    def get_spec(self) -> ManufacturingSpec:
        return ManufacturingSpec(
            intended_pathways=[
                "HALLMARK_INFLAMMATORY_RESPONSE",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                "HALLMARK_IL6_JAK_STAT3_SIGNALING",
            ],
            excluded_pathways=[
                "HALLMARK_APOPTOSIS",
            ],
            feedback_loops=[
                ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_APOPTOSIS"),
                ("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
            ],
        )

    def get_timepoints(self) -> List[str]:
        return [f"P{i}" for i in range(1, 21)]

    def generate_timepoint(
        self, timepoint_id: str, previous: Optional[InjectedBatch]
    ) -> InjectedBatch:
        spec = self.get_spec()
        batch = self.simulator.generate_healthy_batch(
            n_cells=self.N_CELLS, spec=spec
        )

        passage = int(timepoint_id[1:])

        if 6 <= passage <= 10:
            # Gradual drift
            drift_severity = (passage - 5) * 0.1
            batch = self.simulator.inject_atlas_drift(
                batch, drift_distance=drift_severity
            )

        elif 11 <= passage <= 15:
            # Drift + feedback breakdown
            drift_severity = min(1.0, (passage - 5) * 0.1)
            batch = self.simulator.inject_atlas_drift(
                batch, drift_distance=drift_severity
            )
            batch = self.simulator.inject_broken_feedback(
                batch,
                loops=[("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_APOPTOSIS")],
                break_type=FeedbackBreakType.DECOUPLED,
            )

        elif passage > 15:
            # Compounding failure
            batch = self.simulator.inject_atlas_drift(
                batch, drift_distance=0.8
            )
            batch = self.simulator.inject_broken_feedback(
                batch,
                loops=[("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_APOPTOSIS")],
                break_type=FeedbackBreakType.INVERTED,
            )
            batch = self.simulator.inject_subpopulation(
                batch, fraction=0.1 + (passage - 15) * 0.05
            )

        return batch

    def get_expected_outcomes(self) -> Dict[str, Dict[str, bool]]:
        outcomes = {}
        for p in range(1, 21):
            key = f"P{p}"
            if p <= 5:
                outcomes[key] = {}
            elif p <= 10:
                outcomes[key] = {"F5_drift": True}
            elif p <= 15:
                outcomes[key] = {"F5_drift": True, "F10_feedback": True}
            else:
                outcomes[key] = {
                    "F5_drift": True,
                    "F10_feedback": True,
                    "F7_heterogeneity": True,
                    "F12_atlas": True,
                }
        return outcomes
