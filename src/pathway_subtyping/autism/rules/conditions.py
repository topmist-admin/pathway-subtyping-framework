"""
Condition predicates for autism biological rules.

Each condition tests a biological property of a gene/variant/sample.
The ConditionEvaluator checks conditions against a BiologicalContext.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


class ConditionType(Enum):
    """Types of condition predicates."""

    # Variant conditions
    HAS_VARIANT = "has_variant"
    HAS_LOF_VARIANT = "has_lof_variant"
    HAS_MISSENSE_VARIANT = "has_missense_variant"
    HAS_DAMAGING_VARIANT = "has_damaging_variant"
    HAS_MULTIPLE_HITS = "has_multiple_hits"
    # Gene conditions
    IS_CONSTRAINED = "is_constrained"
    IS_SFARI_GENE = "is_sfari_gene"
    IS_HIGH_CONFIDENCE_SFARI = "is_high_confidence_sfari"
    HAS_PARALOG = "has_paralog"
    IS_CHD8_TARGET = "is_chd8_target"
    IS_SYNAPTIC_GENE = "is_synaptic_gene"
    # Expression conditions
    PRENATALLY_EXPRESSED = "prenatally_expressed"
    CELL_TYPE_ENRICHED = "cell_type_enriched"
    # Pathway conditions
    GENE_IN_PATHWAY = "gene_in_pathway"
    PATHWAY_DISRUPTED = "pathway_disrupted"
    HITS_ARE_INDEPENDENT = "hits_are_independent"
    # Drug conditions
    DRUG_TARGETS = "drug_targets"
    MECHANISM_ALIGNMENT = "mechanism_alignment"


@dataclass
class Condition:
    """A condition predicate to test."""

    predicate: str
    arguments: Dict[str, Any] = field(default_factory=dict)
    negated: bool = False

    def __str__(self) -> str:
        neg = "NOT " if self.negated else ""
        args = f"({self.arguments})" if self.arguments else ""
        return f"{neg}{self.predicate}{args}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "predicate": self.predicate,
            "arguments": self.arguments,
            "negated": self.negated,
        }


@dataclass
class ConditionResult:
    """Result of evaluating a condition."""

    satisfied: bool
    evidence: Dict[str, Any] = field(default_factory=dict)
    bound_variables: Dict[str, Any] = field(default_factory=dict)
    explanation: str = ""


class ConditionEvaluator:
    """Evaluates conditions against a biological context and individual data.

    Dispatches each predicate to its evaluator method.
    """

    def __init__(self, biological_context: Any):
        """Initialize with a BiologicalContext providing reference data."""
        self.ctx = biological_context
        self._evaluators = {
            "has_variant": self._eval_has_variant,
            "has_lof_variant": self._eval_has_lof_variant,
            "has_missense_variant": self._eval_has_missense_variant,
            "has_damaging_variant": self._eval_has_damaging_variant,
            "has_multiple_hits": self._eval_has_multiple_hits,
            "is_constrained": self._eval_is_constrained,
            "is_sfari_gene": self._eval_is_sfari_gene,
            "is_high_confidence_sfari": self._eval_is_high_confidence_sfari,
            "has_paralog": self._eval_has_paralog,
            "is_chd8_target": self._eval_is_chd8_target,
            "is_synaptic_gene": self._eval_is_synaptic_gene,
            "prenatally_expressed": self._eval_prenatally_expressed,
            "cell_type_enriched": self._eval_cell_type_enriched,
            "gene_in_pathway": self._eval_gene_in_pathway,
            "pathway_disrupted": self._eval_pathway_disrupted,
            "hits_are_independent": self._eval_hits_are_independent,
            "drug_targets": self._eval_drug_targets,
            "mechanism_alignment": self._eval_mechanism_alignment,
        }

    def evaluate(
        self,
        condition: Condition,
        individual_data: Any,
        gene: Optional[str] = None,
        bindings: Optional[Dict[str, Any]] = None,
    ) -> ConditionResult:
        """Evaluate a condition for an individual.

        Args:
            condition: Condition to evaluate.
            individual_data: IndividualData instance.
            gene: Optional gene context for gene-level conditions.
            bindings: Variable bindings from prior conditions.

        Returns:
            ConditionResult with satisfaction status and evidence.
        """
        evaluator = self._evaluators.get(condition.predicate)
        if evaluator is None:
            logger.warning(
                "[Autism Rules] Unknown predicate: %s", condition.predicate
            )
            return ConditionResult(satisfied=False, explanation=f"Unknown: {condition.predicate}")

        result = evaluator(condition.arguments, individual_data, gene, bindings or {})

        if condition.negated:
            result.satisfied = not result.satisfied

        return result

    # ── Variant evaluators ────────────────────────────────────────────

    def _eval_has_variant(self, args, individual, gene, bindings):
        has = gene in individual.get_genes_with_variants() if gene else bool(individual.variants)
        return ConditionResult(satisfied=has, explanation="Has variant" if has else "No variant")

    def _eval_has_lof_variant(self, args, individual, gene, bindings):
        if gene:
            variants = individual.get_gene_variants(gene)
            has = any(
                getattr(v, "consequence", "") in ("frameshift", "stop_gained", "splice_donor", "splice_acceptor")
                or getattr(v, "is_lof", False)
                for v in variants
            )
        else:
            has = bool(individual.get_genes_with_damaging_variants())
        return ConditionResult(satisfied=has, explanation="Has LoF variant" if has else "No LoF")

    def _eval_has_missense_variant(self, args, individual, gene, bindings):
        if gene:
            variants = individual.get_gene_variants(gene)
            has = any(getattr(v, "consequence", "") == "missense" for v in variants)
        else:
            has = False
        return ConditionResult(satisfied=has)

    def _eval_has_damaging_variant(self, args, individual, gene, bindings):
        if gene:
            has = gene in individual.get_genes_with_damaging_variants()
        else:
            has = bool(individual.get_genes_with_damaging_variants())
        return ConditionResult(satisfied=has)

    def _eval_has_multiple_hits(self, args, individual, gene, bindings):
        min_hits = args.get("min_hits", 2)
        n = len(individual.get_genes_with_damaging_variants())
        return ConditionResult(
            satisfied=n >= min_hits,
            evidence={"n_hits": n, "min_required": min_hits},
        )

    # ── Gene evaluators ───────────────────────────────────────────────

    def _eval_is_constrained(self, args, individual, gene, bindings):
        threshold = args.get("pli_threshold", 0.9)
        if gene:
            is_c = self.ctx.is_gene_constrained(gene, threshold)
            pli = self.ctx.get_pli_score(gene)
            return ConditionResult(
                satisfied=is_c,
                evidence={"pli": pli, "threshold": threshold},
            )
        return ConditionResult(satisfied=False)

    def _eval_is_sfari_gene(self, args, individual, gene, bindings):
        if gene:
            return ConditionResult(satisfied=self.ctx.is_sfari_gene(gene))
        return ConditionResult(satisfied=False)

    def _eval_is_high_confidence_sfari(self, args, individual, gene, bindings):
        if gene:
            return ConditionResult(satisfied=self.ctx.is_high_confidence_sfari(gene))
        return ConditionResult(satisfied=False)

    def _eval_has_paralog(self, args, individual, gene, bindings):
        if gene:
            paralogs = self.ctx.get_paralogs(gene)
            return ConditionResult(
                satisfied=bool(paralogs),
                evidence={"paralogs": list(paralogs)[:3]},
            )
        return ConditionResult(satisfied=False)

    def _eval_is_chd8_target(self, args, individual, gene, bindings):
        gene_is_chd8 = args.get("gene_is_chd8", False)
        if gene_is_chd8:
            return ConditionResult(satisfied=(gene == "CHD8"))
        if gene:
            return ConditionResult(satisfied=self.ctx.is_chd8_target(gene))
        return ConditionResult(satisfied=False)

    def _eval_is_synaptic_gene(self, args, individual, gene, bindings):
        if gene:
            return ConditionResult(satisfied=self.ctx.is_synaptic_gene(gene))
        return ConditionResult(satisfied=False)

    # ── Expression evaluators ─────────────────────────────────────────

    def _eval_prenatally_expressed(self, args, individual, gene, bindings):
        if gene:
            return ConditionResult(satisfied=self.ctx.is_prenatally_expressed(gene))
        return ConditionResult(satisfied=False)

    def _eval_cell_type_enriched(self, args, individual, gene, bindings):
        cell_type = args.get("cell_type", "")
        if gene and cell_type:
            return ConditionResult(
                satisfied=self.ctx.is_enriched_in_cell_type(gene, cell_type)
            )
        return ConditionResult(satisfied=False)

    # ── Pathway evaluators ────────────────────────────────────────────

    def _eval_gene_in_pathway(self, args, individual, gene, bindings):
        pathway = args.get("pathway_id", bindings.get("pathway_id", ""))
        if gene and pathway:
            return ConditionResult(satisfied=self.ctx.is_gene_in_pathway(gene, pathway))
        return ConditionResult(satisfied=False)

    def _eval_pathway_disrupted(self, args, individual, gene, bindings):
        threshold = args.get("score_threshold", 2.0)
        pathway = bindings.get("pathway_id", "")
        if pathway:
            score = individual.get_pathway_score(pathway)
            if score is not None:
                return ConditionResult(
                    satisfied=abs(score) >= threshold,
                    evidence={"pathway_score": score, "threshold": threshold},
                )
        # Check any pathway
        for pw, score in individual.pathway_scores.items():
            if abs(score) >= threshold:
                return ConditionResult(
                    satisfied=True,
                    evidence={"pathway_id": pw, "score": score},
                    bound_variables={"pathway_id": pw},
                )
        return ConditionResult(satisfied=False)

    def _eval_hits_are_independent(self, args, individual, gene, bindings):
        genes = individual.get_genes_with_damaging_variants()
        return ConditionResult(satisfied=len(genes) >= 2)

    # ── Drug evaluators ───────────────────────────────────────────────

    def _eval_drug_targets(self, args, individual, gene, bindings):
        if gene:
            drugs = self.ctx.get_drug_targets_for_gene(gene)
            return ConditionResult(
                satisfied=bool(drugs),
                evidence={"drugs": list(drugs)[:3]},
            )
        return ConditionResult(satisfied=False)

    def _eval_mechanism_alignment(self, args, individual, gene, bindings):
        pathway = bindings.get("pathway_id", "")
        drug = bindings.get("drug_id", "")
        if pathway and drug:
            return ConditionResult(
                satisfied=self.ctx.check_mechanism_pathway_alignment(drug, pathway)
            )
        return ConditionResult(satisfied=False)
