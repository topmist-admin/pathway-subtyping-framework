"""
Curated Biological Rules for Autism Variant Interpretation.

Seven rules (R1-R7) encoding domain knowledge about autism-specific
gene-pathway-variant relationships. Ported from autism-pathway-framework
Module 09.

Evidence sources:
    - Samocha et al. 2014 (constraint scores)
    - Willsey et al. 2013 (prenatal cortical expression)
    - Sanders et al. 2015 (autism de novo variants)
    - Voineagu et al. 2011 (transcriptomic convergence)
    - Cotney et al. 2015 (CHD8 ChIP-seq targets)
    - Sugathan et al. 2014 (CHD8 knockdown transcriptomics)
    - Koopmans et al. 2019 (SynGO database)
    - Satterstrom et al. 2020 (cell-type expression)
    - Rubenstein & Merzenich 2003 (E/I imbalance)

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

from .conditions import Condition

logger = logging.getLogger(__name__)


class ConclusionType(Enum):
    """Types of conclusions a rule can produce."""

    PATHWAY_DISRUPTION = "pathway_disruption"
    PATHWAY_CONVERGENCE = "pathway_convergence"
    SUBTYPE_INDICATOR = "subtype_indicator"
    EFFECT_MODIFIER = "effect_modifier"
    THERAPEUTIC_HYPOTHESIS = "therapeutic_hypothesis"


@dataclass
class Conclusion:
    """Conclusion produced when a rule fires."""

    type: str
    attributes: Dict[str, Any] = field(default_factory=dict)
    confidence_modifier: float = 1.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.type,
            "attributes": self.attributes,
            "confidence_modifier": self.confidence_modifier,
        }


@dataclass
class Rule:
    """A biological interpretation rule."""

    id: str
    name: str
    description: str
    conditions: List[Condition]
    conclusion: Conclusion
    base_confidence: float = 0.8
    evidence_sources: List[str] = field(default_factory=list)
    logic: str = "AND"  # "AND" or "OR"
    priority: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "conditions": [c.to_dict() for c in self.conditions],
            "conclusion": self.conclusion.to_dict(),
            "base_confidence": self.base_confidence,
            "evidence_sources": self.evidence_sources,
            "logic": self.logic,
            "priority": self.priority,
        }


class BiologicalRules:
    """Factory class providing curated autism biological rules."""

    @staticmethod
    def R1_constrained_lof_developing_cortex() -> Rule:
        """pLI>=0.9 + prenatal cortex expression + LoF -> haploinsufficiency."""
        return Rule(
            id="R1",
            name="Constrained LoF in Developing Cortex",
            description=(
                "Loss-of-function variant in a constrained gene (pLI>=0.9) "
                "that is expressed in prenatal cortex indicates high-confidence "
                "haploinsufficiency-driven pathway disruption."
            ),
            conditions=[
                Condition(predicate="has_lof_variant"),
                Condition(predicate="is_constrained", arguments={"pli_threshold": 0.9}),
                Condition(predicate="prenatally_expressed"),
            ],
            conclusion=Conclusion(
                type=ConclusionType.PATHWAY_DISRUPTION.value,
                attributes={"mechanism": "haploinsufficiency", "tissue": "prenatal_cortex"},
            ),
            base_confidence=0.90,
            evidence_sources=[
                "Samocha et al. 2014 (pLI constraint)",
                "Willsey et al. 2013 (prenatal cortical expression)",
            ],
            priority=10,
        )

    @staticmethod
    def R2_pathway_convergence() -> Rule:
        """>=2 independent hits in same pathway -> convergence."""
        return Rule(
            id="R2",
            name="Pathway Convergence",
            description=(
                "Two or more independent damaging variants in the same "
                "pathway indicate pathway-level convergence signal."
            ),
            conditions=[
                Condition(predicate="has_multiple_hits", arguments={"min_hits": 2}),
                Condition(predicate="hits_are_independent"),
            ],
            conclusion=Conclusion(
                type=ConclusionType.PATHWAY_CONVERGENCE.value,
                attributes={"mechanism": "convergent_disruption"},
            ),
            base_confidence=0.85,
            evidence_sources=[
                "Voineagu et al. 2011 (transcriptomic convergence)",
                "Sanders et al. 2015 (de novo convergence)",
            ],
            priority=8,
        )

    @staticmethod
    def R3_chd8_cascade() -> Rule:
        """CHD8 disruption -> chromatin remodeling cascade."""
        return Rule(
            id="R3",
            name="CHD8 Chromatin Cascade",
            description=(
                "Disruption of CHD8, a master chromatin regulator, triggers "
                "a cascade affecting hundreds of downstream targets."
            ),
            conditions=[
                Condition(predicate="has_damaging_variant"),
                Condition(predicate="is_chd8_target", arguments={"gene_is_chd8": True}),
            ],
            conclusion=Conclusion(
                type=ConclusionType.SUBTYPE_INDICATOR.value,
                attributes={"subtype": "chromatin_remodeling", "master_regulator": "CHD8"},
            ),
            base_confidence=0.92,
            evidence_sources=[
                "Cotney et al. 2015 (CHD8 ChIP-seq)",
                "Sugathan et al. 2014 (CHD8 knockdown)",
            ],
            priority=9,
        )

    @staticmethod
    def R3b_chd8_target_cascade() -> Rule:
        """CHD8 regulatory target hit -> indirect cascade."""
        return Rule(
            id="R3b",
            name="CHD8 Target Disruption",
            description=(
                "Damaging variant in a CHD8 regulatory target gene indicates "
                "indirect disruption of the chromatin remodeling cascade."
            ),
            conditions=[
                Condition(predicate="has_damaging_variant"),
                Condition(predicate="is_chd8_target"),
            ],
            conclusion=Conclusion(
                type=ConclusionType.PATHWAY_DISRUPTION.value,
                attributes={"mechanism": "chd8_target_disruption"},
            ),
            base_confidence=0.75,
            evidence_sources=["Cotney et al. 2015 (CHD8 ChIP-seq targets)"],
            priority=7,
        )

    @staticmethod
    def R4_synaptic_excitatory() -> Rule:
        """Synaptic gene + LoF + excitatory neuron expression -> E/I imbalance."""
        return Rule(
            id="R4",
            name="Synaptic Excitatory Disruption",
            description=(
                "LoF in a synaptic gene expressed in excitatory neurons "
                "indicates excitatory/inhibitory imbalance subtype."
            ),
            conditions=[
                Condition(predicate="has_lof_variant"),
                Condition(predicate="is_synaptic_gene"),
                Condition(
                    predicate="cell_type_enriched", arguments={"cell_type": "excitatory_neuron"}
                ),
            ],
            conclusion=Conclusion(
                type=ConclusionType.SUBTYPE_INDICATOR.value,
                attributes={
                    "subtype": "excitatory_inhibitory_imbalance",
                    "direction": "excitatory",
                },
            ),
            base_confidence=0.82,
            evidence_sources=[
                "Koopmans et al. 2019 (SynGO)",
                "Rubenstein & Merzenich 2003 (E/I imbalance)",
                "Satterstrom et al. 2020 (cell-type expression)",
            ],
            priority=7,
        )

    @staticmethod
    def R5_compensatory_paralog() -> Rule:
        """LoF + intact expressed paralog -> effect modifier."""
        return Rule(
            id="R5",
            name="Paralog Compensation",
            description=(
                "Damaging variant in a gene with an intact, expressed paralog "
                "may have reduced penetrance due to compensation."
            ),
            conditions=[
                Condition(predicate="has_damaging_variant"),
                Condition(predicate="has_paralog"),
            ],
            conclusion=Conclusion(
                type=ConclusionType.EFFECT_MODIFIER.value,
                attributes={"modifier": "paralog_compensation", "direction": "protective"},
                confidence_modifier=0.7,
            ),
            base_confidence=0.70,
            evidence_sources=["Diss & Bhatt 2020 (paralog compensation)"],
            priority=5,
        )

    @staticmethod
    def R6_drug_pathway_target() -> Rule:
        """Disrupted pathway + druggable gene -> therapeutic hypothesis."""
        return Rule(
            id="R6",
            name="Therapeutic Pathway Targeting",
            description=(
                "Druggable gene in a disrupted pathway with mechanism alignment "
                "generates a therapeutic repurposing hypothesis."
            ),
            conditions=[
                Condition(predicate="pathway_disrupted"),
                Condition(predicate="drug_targets"),
                Condition(predicate="mechanism_alignment"),
            ],
            conclusion=Conclusion(
                type=ConclusionType.THERAPEUTIC_HYPOTHESIS.value,
                attributes={"action": "drug_repurposing"},
            ),
            base_confidence=0.60,
            evidence_sources=[],
            priority=4,
        )

    @staticmethod
    def R7_sfari_high_confidence() -> Rule:
        """High-confidence SFARI gene (score=1) + damaging variant."""
        return Rule(
            id="R7",
            name="SFARI High-Confidence Gene",
            description=(
                "Damaging variant in a SFARI confidence level 1 gene "
                "indicates high-confidence ASD gene involvement."
            ),
            conditions=[
                Condition(predicate="has_damaging_variant"),
                Condition(predicate="is_high_confidence_sfari"),
            ],
            conclusion=Conclusion(
                type=ConclusionType.PATHWAY_DISRUPTION.value,
                attributes={"mechanism": "sfari_high_confidence"},
            ),
            base_confidence=0.88,
            evidence_sources=["SFARI Gene Database (https://gene.sfari.org)"],
            priority=9,
        )

    @classmethod
    def get_all_rules(cls) -> List[Rule]:
        """Return all 8 rules (R1-R7 + R3b)."""
        return [
            cls.R1_constrained_lof_developing_cortex(),
            cls.R2_pathway_convergence(),
            cls.R3_chd8_cascade(),
            cls.R3b_chd8_target_cascade(),
            cls.R4_synaptic_excitatory(),
            cls.R5_compensatory_paralog(),
            cls.R6_drug_pathway_target(),
            cls.R7_sfari_high_confidence(),
        ]

    @classmethod
    def get_core_rules(cls) -> List[Rule]:
        """Return core rules R1-R6 (without R7)."""
        return [
            cls.R1_constrained_lof_developing_cortex(),
            cls.R2_pathway_convergence(),
            cls.R3_chd8_cascade(),
            cls.R3b_chd8_target_cascade(),
            cls.R4_synaptic_excitatory(),
            cls.R5_compensatory_paralog(),
            cls.R6_drug_pathway_target(),
        ]

    @classmethod
    def get_rule_by_id(cls, rule_id: str) -> Optional[Rule]:
        """Get a specific rule by ID."""
        for rule in cls.get_all_rules():
            if rule.id == rule_id:
                return rule
        return None

    @classmethod
    def get_rules_by_conclusion_type(cls, conclusion_type: ConclusionType) -> List[Rule]:
        """Get rules producing a specific conclusion type."""
        return [r for r in cls.get_all_rules() if r.conclusion.type == conclusion_type.value]
