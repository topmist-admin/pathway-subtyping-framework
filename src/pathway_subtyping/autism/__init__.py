"""
Autism-Specific Interpretation Layer.

Optional module providing autism-only biological interpretation on top of
the disease-agnostic PSF pipeline. Implements Mohit's standing rule:
tool is disease-agnostic, biological interpretation is autism-only.

Submodules:
    - rules: Curated biological rules (R1-R7) for autism variant interpretation
    - neurosymbolic: Combines GNN predictions with symbolic rule conclusions
    - therapeutic: Drug repurposing hypothesis generation with safety scoring

Install: pip install pathway-subtyping[autism]

Research use only. Not for clinical decision-making.
"""

from .rules.biological_rules import BiologicalRules, Conclusion, ConclusionType, Rule
from .rules.conditions import Condition, ConditionResult
from .rules.engine import (
    AutismRuleEngine,
    BiologicalContext,
    FiredRule,
    IndividualData,
)
from .rules.explanation import ExplanationGenerator, ReasoningChain, ReasoningStep

__all__ = [
    # Rules
    "BiologicalRules",
    "Rule",
    "Conclusion",
    "ConclusionType",
    "Condition",
    "ConditionResult",
    # Engine
    "AutismRuleEngine",
    "BiologicalContext",
    "IndividualData",
    "FiredRule",
    # Explanation
    "ExplanationGenerator",
    "ReasoningChain",
    "ReasoningStep",
]
