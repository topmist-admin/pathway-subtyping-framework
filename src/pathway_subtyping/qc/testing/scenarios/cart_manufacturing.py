"""
Scenario A: CAR T Manufacturing Run.

5 timepoints (Day 0, 3, 7, 10, 14), 500 cells per timepoint.
Nutrient depletion injected at Day 7, recovery by Day 14.

Research use only. Not for clinical decision-making.
"""

import logging
from typing import Dict, List, Optional

from ..simulator import (
    InjectedBatch,
    ManufacturingSimulator,
    ManufacturingSpec,
    StressorType,
)
from .base import BaseScenario

logger = logging.getLogger(__name__)


class CARTManufacturingScenario(BaseScenario):
    """CAR T cell manufacturing run with nutrient depletion at Day 7."""

    N_CELLS = 500

    def get_spec(self) -> ManufacturingSpec:
        return ManufacturingSpec(
            intended_pathways=[
                "HALLMARK_INFLAMMATORY_RESPONSE",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_MTORC1_SIGNALING",
            ],
            excluded_pathways=[
                "HALLMARK_APOPTOSIS",
                "HALLMARK_P53_PATHWAY",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            ],
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8),
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB": (0.3, 0.7),
                "HALLMARK_IL6_JAK_STAT3_SIGNALING": (0.2, 0.6),
                "HALLMARK_INTERFERON_GAMMA_RESPONSE": (0.3, 0.8),
                "HALLMARK_MTORC1_SIGNALING": (0.3, 0.7),
            },
        )

    def get_timepoints(self) -> List[str]:
        return ["day_0", "day_3", "day_7", "day_10", "day_14"]

    def generate_timepoint(
        self, timepoint_id: str, previous: Optional[InjectedBatch]
    ) -> InjectedBatch:
        spec = self.get_spec()
        batch = self.simulator.generate_healthy_batch(n_cells=self.N_CELLS, spec=spec)

        if timepoint_id == "day_7":
            # Inject nutrient depletion stress
            batch = self.simulator.inject_stress(
                batch, stressor=StressorType.NUTRIENT_DEPLETION, severity=0.5
            )
        elif timepoint_id == "day_10":
            # Still recovering from Day 7 stress
            batch = self.simulator.inject_stress(
                batch, stressor=StressorType.NUTRIENT_DEPLETION, severity=0.3
            )

        return batch

    def get_expected_outcomes(self) -> Dict[str, Dict[str, bool]]:
        return {
            "day_0": {},
            "day_3": {},
            "day_7": {"F11_stress": True},
            "day_10": {"F11_stress": True},
            "day_14": {},
        }
