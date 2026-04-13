"""
Explanation generator for autism rule firings.

Produces human-readable reasoning chains from fired rules,
with clinical summaries and cross-individual comparison.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger(__name__)


@dataclass
class ReasoningStep:
    """A single step in a reasoning chain."""

    step_number: int
    description: str
    evidence: Dict[str, Any] = field(default_factory=dict)
    rule_id: Optional[str] = None
    confidence: float = 1.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "step": self.step_number,
            "description": self.description,
            "rule_id": self.rule_id,
            "confidence": round(self.confidence, 4),
        }


@dataclass
class ReasoningChain:
    """Complete reasoning chain for one individual."""

    individual_id: str
    fired_rules: List[Any] = field(default_factory=list)
    pathway_conclusions: Dict[str, float] = field(default_factory=dict)
    subtype_indicators: List[str] = field(default_factory=list)
    therapeutic_hypotheses: List[Dict[str, Any]] = field(default_factory=list)
    explanation_text: str = ""
    steps: List[ReasoningStep] = field(default_factory=list)
    timestamp: datetime = field(default_factory=datetime.now)

    @property
    def n_rules_fired(self) -> int:
        return len(self.fired_rules)

    @property
    def average_confidence(self) -> float:
        if not self.fired_rules:
            return 0.0
        return sum(fr.confidence for fr in self.fired_rules) / len(self.fired_rules)

    @property
    def genes_affected(self) -> Set[str]:
        return {fr.gene for fr in self.fired_rules if fr.gene} - {None}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "individual_id": self.individual_id,
            "n_rules_fired": self.n_rules_fired,
            "average_confidence": round(self.average_confidence, 4),
            "genes_affected": list(self.genes_affected),
            "pathway_conclusions": self.pathway_conclusions,
            "subtype_indicators": self.subtype_indicators,
            "n_therapeutic_hypotheses": len(self.therapeutic_hypotheses),
            "explanation": self.explanation_text,
            "steps": [s.to_dict() for s in self.steps],
        }


class ExplanationGenerator:
    """Generates human-readable explanations from fired rules.

    Usage::

        gen = ExplanationGenerator()
        chain = gen.generate_reasoning_chain("sample_001", fired_rules)
        print(chain.explanation_text)
    """

    DISCLAIMER = (
        "DISCLAIMER: These interpretations are computational hypotheses "
        "for research purposes only. They are NOT clinical diagnoses or "
        "treatment recommendations. All findings require independent "
        "validation by qualified researchers and clinicians."
    )

    def __init__(
        self,
        include_evidence: bool = True,
        include_disclaimers: bool = True,
        verbose: bool = False,
    ):
        self.include_evidence = include_evidence
        self.include_disclaimers = include_disclaimers
        self.verbose = verbose

    def generate_reasoning_chain(
        self,
        individual_id: str,
        fired_rules: List[Any],
    ) -> ReasoningChain:
        """Generate a full reasoning chain from fired rules.

        Args:
            individual_id: Sample identifier.
            fired_rules: List of FiredRule instances.

        Returns:
            ReasoningChain with structured explanation.
        """
        pathway_conclusions = self._extract_pathway_conclusions(fired_rules)
        subtype_indicators = self._extract_subtype_indicators(fired_rules)
        therapeutic = self._extract_therapeutic_hypotheses(fired_rules)
        steps = self._generate_steps(fired_rules)

        chain = ReasoningChain(
            individual_id=individual_id,
            fired_rules=fired_rules,
            pathway_conclusions=pathway_conclusions,
            subtype_indicators=subtype_indicators,
            therapeutic_hypotheses=therapeutic,
            steps=steps,
        )

        chain.explanation_text = self._generate_explanation_text(chain)
        return chain

    def _extract_pathway_conclusions(self, fired_rules: List[Any]) -> Dict[str, float]:
        conclusions: Dict[str, float] = {}
        for fr in fired_rules:
            ct = fr.rule.conclusion.type
            if ct in ("pathway_disruption", "pathway_convergence"):
                gene = fr.gene or "unknown"
                key = f"{ct}:{gene}"
                conclusions[key] = max(conclusions.get(key, 0.0), fr.confidence)
        return conclusions

    def _extract_subtype_indicators(self, fired_rules: List[Any]) -> List[str]:
        indicators = []
        for fr in fired_rules:
            if fr.rule.conclusion.type == "subtype_indicator":
                subtype = fr.rule.conclusion.attributes.get("subtype", "unknown")
                if subtype not in indicators:
                    indicators.append(subtype)
        return indicators

    def _extract_therapeutic_hypotheses(self, fired_rules: List[Any]) -> List[Dict[str, Any]]:
        hypotheses = []
        for fr in fired_rules:
            if fr.rule.conclusion.type == "therapeutic_hypothesis":
                hypotheses.append({
                    "rule_id": fr.rule.id,
                    "gene": fr.gene,
                    "confidence": fr.confidence,
                    "bindings": fr.bindings,
                })
        return hypotheses

    def _generate_steps(self, fired_rules: List[Any]) -> List[ReasoningStep]:
        steps = []
        for i, fr in enumerate(fired_rules, 1):
            steps.append(ReasoningStep(
                step_number=i,
                description=fr.explanation,
                evidence=fr.evidence,
                rule_id=fr.rule.id,
                confidence=fr.confidence,
            ))
        return steps

    def _generate_explanation_text(self, chain: ReasoningChain) -> str:
        lines = [f"Reasoning chain for individual: {chain.individual_id}"]
        lines.append(f"Rules fired: {chain.n_rules_fired}")
        lines.append(f"Average confidence: {chain.average_confidence:.2f}")
        lines.append(f"Genes affected: {', '.join(chain.genes_affected) or 'none'}")
        lines.append("")

        for step in chain.steps:
            lines.append(f"  Step {step.step_number} [{step.rule_id}]: {step.description}")

        if chain.subtype_indicators:
            lines.append(f"\nSubtype indicators: {', '.join(chain.subtype_indicators)}")

        if chain.therapeutic_hypotheses:
            lines.append(f"\nTherapeutic hypotheses: {len(chain.therapeutic_hypotheses)}")

        if self.include_disclaimers:
            lines.append(f"\n{self.DISCLAIMER}")

        return "\n".join(lines)

    def generate_clinical_summary(self, chain: ReasoningChain) -> str:
        """Generate a brief clinical-style summary."""
        genes = ", ".join(sorted(chain.genes_affected)[:5]) or "no specific genes"
        n_rules = chain.n_rules_fired
        subtypes = ", ".join(chain.subtype_indicators) or "no subtype indicators"

        return (
            f"Individual {chain.individual_id}: {n_rules} biological rules fired "
            f"(avg confidence {chain.average_confidence:.2f}). "
            f"Affected genes: {genes}. Subtype indicators: {subtypes}. "
            f"[Research hypothesis only — not a clinical diagnosis.]"
        )
