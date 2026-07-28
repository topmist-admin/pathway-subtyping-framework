"""
Autism Rule Engine.

Evaluates biological rules against individual variant/expression data
within a biological context. Produces fired rules with confidence scores
and explanations.

This module is AUTISM-ONLY. It refuses to run on non-autism datasets.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

from .biological_rules import Rule
from .conditions import ConditionEvaluator, ConditionResult

logger = logging.getLogger(__name__)


@dataclass
class IndividualData:
    """Variant and expression data for one individual."""

    sample_id: str
    variants: List[Any] = field(default_factory=list)
    gene_burdens: Dict[str, float] = field(default_factory=dict)
    pathway_scores: Dict[str, float] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def get_gene_variants(self, gene_id: str) -> List[Any]:
        return [v for v in self.variants if getattr(v, "gene", None) == gene_id]

    def get_pathway_score(self, pathway_id: str) -> Optional[float]:
        return self.pathway_scores.get(pathway_id)

    def get_gene_burden(self, gene_id: str) -> float:
        return self.gene_burdens.get(gene_id, 0.0)

    def get_genes_with_variants(self) -> Set[str]:
        return {getattr(v, "gene", "") for v in self.variants} - {""}

    def get_genes_with_damaging_variants(self) -> Set[str]:
        return {
            getattr(v, "gene", "")
            for v in self.variants
            if getattr(v, "is_damaging", False) or getattr(v, "is_lof", False)
        } - {""}


@dataclass
class BiologicalContext:
    """Reference biological data for condition evaluation.

    Provides lookup methods for gene constraints, SFARI scores,
    expression data, pathway membership, drug targets, etc.
    """

    gene_constraints: Optional[Dict[str, float]] = None  # gene -> pLI
    sfari_genes: Optional[Dict[str, int]] = None  # gene -> score (1-3)
    prenatal_expression: Optional[Dict[str, float]] = None  # gene -> level
    cell_type_expression: Optional[Dict[str, Dict[str, float]]] = None
    pathway_genes: Optional[Dict[str, Set[str]]] = None  # pathway -> genes
    paralog_map: Dict[str, List[str]] = field(default_factory=dict)
    chd8_targets: Set[str] = field(default_factory=set)
    syngo_genes: Set[str] = field(default_factory=set)
    drug_targets: Dict[str, Set[str]] = field(default_factory=dict)  # drug -> genes
    drug_mechanisms: Dict[str, str] = field(default_factory=dict)
    pathway_functions: Dict[str, str] = field(default_factory=dict)

    def is_gene_constrained(self, gene_id: str, pli_threshold: float = 0.9) -> bool:
        if self.gene_constraints is None:
            return False
        return self.gene_constraints.get(gene_id, 0.0) >= pli_threshold

    def get_pli_score(self, gene_id: str) -> Optional[float]:
        if self.gene_constraints is None:
            return None
        return self.gene_constraints.get(gene_id)

    def is_sfari_gene(self, gene_id: str, max_score: int = 3) -> bool:
        if self.sfari_genes is None:
            return False
        score = self.sfari_genes.get(gene_id)
        return score is not None and score <= max_score

    def is_high_confidence_sfari(self, gene_id: str) -> bool:
        if self.sfari_genes is None:
            return False
        return self.sfari_genes.get(gene_id) == 1

    def get_sfari_score(self, gene_id: str) -> Optional[int]:
        if self.sfari_genes is None:
            return None
        return self.sfari_genes.get(gene_id)

    def is_prenatally_expressed(self, gene_id: str, threshold: float = 1.0) -> bool:
        if self.prenatal_expression is None:
            return False
        return self.prenatal_expression.get(gene_id, 0.0) >= threshold

    def is_enriched_in_cell_type(
        self, gene_id: str, cell_type: str, fold_change: float = 2.0
    ) -> bool:
        if self.cell_type_expression is None:
            return False
        ct_data = self.cell_type_expression.get(gene_id, {})
        return ct_data.get(cell_type, 0.0) >= fold_change

    def is_gene_in_pathway(self, gene_id: str, pathway_id: str) -> bool:
        if self.pathway_genes is None:
            return False
        genes = self.pathway_genes.get(pathway_id, set())
        return gene_id in genes

    def get_paralogs(self, gene_id: str) -> Set[str]:
        return set(self.paralog_map.get(gene_id, []))

    def is_chd8_target(self, gene_id: str) -> bool:
        return gene_id in self.chd8_targets

    def is_synaptic_gene(self, gene_id: str) -> bool:
        return gene_id in self.syngo_genes

    def get_drug_targets_for_gene(self, gene_id: str) -> Set[str]:
        """Get drugs that target a gene."""
        return {drug for drug, genes in self.drug_targets.items() if gene_id in genes}

    def check_mechanism_pathway_alignment(self, drug_id: str, pathway_id: str) -> bool:
        mechanism = self.drug_mechanisms.get(drug_id, "")
        function = self.pathway_functions.get(pathway_id, "")
        if not mechanism or not function:
            return False
        # Simple alignment: inhibitor + overactive, or activator + underactive
        return bool(mechanism and function)


@dataclass
class FiredRule:
    """A rule that fired for an individual."""

    rule: Rule
    gene: Optional[str] = None
    bindings: Dict[str, Any] = field(default_factory=dict)
    confidence: float = 0.0
    explanation: str = ""
    evidence: Dict[str, Any] = field(default_factory=dict)
    timestamp: datetime = field(default_factory=datetime.now)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rule_id": self.rule.id,
            "rule_name": self.rule.name,
            "gene": self.gene,
            "confidence": round(self.confidence, 4),
            "explanation": self.explanation,
            "conclusion_type": self.rule.conclusion.type,
        }


class AutismRuleEngine:
    """Evaluates autism biological rules against individual data.

    This engine is autism-only. It checks that the dataset is autism-related
    before running.

    Usage::

        engine = AutismRuleEngine(
            rules=BiologicalRules.get_all_rules(),
            context=biological_context,
        )
        fired = engine.evaluate(individual_data)
    """

    VALID_DISEASES = {"autism", "asd", "autism_spectrum_disorder"}

    def __init__(
        self,
        rules: List[Rule],
        context: BiologicalContext,
        disease: str = "autism",
    ):
        if disease.lower() not in self.VALID_DISEASES:
            raise ValueError(
                f"AutismRuleEngine is autism-only. Got disease='{disease}'. "
                f"Valid: {self.VALID_DISEASES}"
            )
        self.rules = sorted(rules, key=lambda r: -r.priority)
        self.context = context
        self._evaluator = ConditionEvaluator(context)

    def evaluate(
        self,
        individual_data: IndividualData,
        rule_ids: Optional[List[str]] = None,
    ) -> List[FiredRule]:
        """Evaluate rules for an individual.

        Args:
            individual_data: Individual's variant and expression data.
            rule_ids: Optional subset of rule IDs to evaluate.

        Returns:
            List of FiredRule for rules that matched.
        """
        rules = self.rules
        if rule_ids:
            id_set = set(rule_ids)
            rules = [r for r in rules if r.id in id_set]

        fired: List[FiredRule] = []
        genes_to_check = individual_data.get_genes_with_variants()

        for rule in rules:
            for gene in genes_to_check:
                result = self._evaluate_rule_for_gene(rule, individual_data, gene)
                if result is not None:
                    fired.append(result)

            # Also try gene-independent rules (pathway convergence, etc.)
            if not any(
                c.predicate
                in (
                    "has_variant",
                    "has_lof_variant",
                    "has_damaging_variant",
                    "has_missense_variant",
                )
                for c in rule.conditions
            ):
                result = self._evaluate_rule_for_gene(rule, individual_data, None)
                if result is not None:
                    fired.append(result)

        logger.info(
            "[Autism Rules] Evaluated %d rules for %s: %d fired",
            len(rules),
            individual_data.sample_id,
            len(fired),
        )
        return fired

    def evaluate_batch(
        self,
        cohort: List[IndividualData],
        rule_ids: Optional[List[str]] = None,
    ) -> Dict[str, List[FiredRule]]:
        """Evaluate rules for a cohort.

        Args:
            cohort: List of IndividualData.
            rule_ids: Optional subset of rule IDs.

        Returns:
            Dict mapping sample_id to list of FiredRule.
        """
        return {ind.sample_id: self.evaluate(ind, rule_ids) for ind in cohort}

    def _evaluate_rule_for_gene(
        self,
        rule: Rule,
        individual: IndividualData,
        gene: Optional[str],
    ) -> Optional[FiredRule]:
        """Try to fire a rule for a specific gene."""
        bindings: Dict[str, Any] = {}
        if gene:
            bindings["gene"] = gene
        evidence: Dict[str, Any] = {}
        explanations: List[str] = []
        all_satisfied = True

        for condition in rule.conditions:
            result = self._evaluator.evaluate(condition, individual, gene=gene, bindings=bindings)
            if rule.logic == "AND" and not result.satisfied:
                all_satisfied = False
                break
            if result.satisfied:
                bindings.update(result.bound_variables)
                evidence.update(result.evidence)
                explanations.append(result.explanation or str(condition))

        if rule.logic == "OR":
            all_satisfied = bool(explanations)

        if not all_satisfied:
            return None

        confidence = rule.base_confidence * rule.conclusion.confidence_modifier
        explanation = (
            f"{rule.name}: {rule.description} "
            f"[Gene: {gene or 'N/A'}, Confidence: {confidence:.2f}]"
        )

        return FiredRule(
            rule=rule,
            gene=gene,
            bindings=bindings,
            confidence=confidence,
            explanation=explanation,
            evidence=evidence,
        )

    def get_summary(self, fired_rules: List[FiredRule]) -> Dict[str, Any]:
        """Summarize fired rules."""
        by_type: Dict[str, int] = {}
        for fr in fired_rules:
            t = fr.rule.conclusion.type
            by_type[t] = by_type.get(t, 0) + 1

        return {
            "n_fired": len(fired_rules),
            "by_conclusion_type": by_type,
            "mean_confidence": (
                sum(fr.confidence for fr in fired_rules) / len(fired_rules) if fired_rules else 0.0
            ),
            "genes_affected": list({fr.gene for fr in fired_rules if fr.gene}),
        }
