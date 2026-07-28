"""
Therapeutic Hypothesis Ranking.

Ranks drug repurposing hypotheses by multi-criteria evidence scoring
with diversity constraints to ensure coverage across pathways.

All hypotheses carry requires_validation=True — these are computational
hypotheses, NOT clinical recommendations.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional

from .drug_mapping import DrugCandidate, PathwayDrugMapper
from .evidence import EvidenceScore, EvidenceScorer

logger = logging.getLogger(__name__)


@dataclass
class TherapeuticHypothesis:
    """A ranked therapeutic repurposing hypothesis.

    requires_validation is ALWAYS True — enforced in __post_init__.
    """

    drug_id: str
    drug_name: str
    target_pathway: str
    target_genes: List[str]
    mechanism: str
    score: float
    evidence: EvidenceScore
    explanation: str
    confidence: float
    rank: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)
    requires_validation: bool = field(default=True, init=False)

    def __post_init__(self) -> None:
        self.requires_validation = True  # Cannot be overridden

    @property
    def is_high_evidence(self) -> bool:
        return self.evidence.overall >= 0.7

    @property
    def has_safety_concerns(self) -> bool:
        return self.evidence.has_critical_safety_flags

    def summary(self) -> str:
        return (
            f"#{self.rank} {self.drug_name} -> {self.target_pathway} "
            f"(score={self.score:.2f}, {self.evidence.level.value}). "
            f"{'⚠ SAFETY' if self.has_safety_concerns else ''}"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rank": self.rank,
            "drug_id": self.drug_id,
            "drug_name": self.drug_name,
            "target_pathway": self.target_pathway,
            "target_genes": self.target_genes,
            "mechanism": self.mechanism,
            "score": round(self.score, 4),
            "confidence": round(self.confidence, 4),
            "evidence": self.evidence.to_dict(),
            "explanation": self.explanation,
            "requires_validation": self.requires_validation,
        }


@dataclass
class RankingConfig:
    """Configuration for hypothesis ranking."""

    weight_evidence: float = 0.4
    weight_pathway_score: float = 0.3
    weight_drug_relevance: float = 0.2
    weight_mechanism_match: float = 0.1
    min_evidence_score: float = 0.2
    max_hypotheses: int = 50
    exclude_critical_safety: bool = False
    max_per_pathway: int = 5
    max_per_mechanism: int = 3


@dataclass
class RankingResult:
    """Result of therapeutic hypothesis ranking."""

    hypotheses: List[TherapeuticHypothesis]
    pathways_analyzed: List[str] = field(default_factory=list)
    drugs_considered: int = 0
    config: RankingConfig = field(default_factory=RankingConfig)
    timestamp: datetime = field(default_factory=datetime.now)

    @property
    def top_hypotheses(self) -> List[TherapeuticHypothesis]:
        return self.hypotheses[:10]

    @property
    def high_evidence_count(self) -> int:
        return sum(1 for h in self.hypotheses if h.is_high_evidence)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_hypotheses": len(self.hypotheses),
            "n_high_evidence": self.high_evidence_count,
            "pathways_analyzed": self.pathways_analyzed,
            "drugs_considered": self.drugs_considered,
            "hypotheses": [h.to_dict() for h in self.hypotheses],
        }

    def summary(self) -> str:
        return (
            f"Ranking: {len(self.hypotheses)} hypotheses across "
            f"{len(self.pathways_analyzed)} pathways. "
            f"{self.high_evidence_count} high-evidence. "
            f"[Computational hypotheses only — not clinical recommendations.]"
        )


class HypothesisRanker:
    """Ranks therapeutic hypotheses with diversity constraints.

    Usage::

        ranker = HypothesisRanker(drug_mapper, evidence_scorer)
        result = ranker.rank(
            pathway_scores={"MAPK": 2.5, "PI3K": -1.8},
            pathway_genes={"MAPK": ["RAS", "RAF"], "PI3K": ["PIK3CA"]},
        )
    """

    def __init__(
        self,
        drug_mapper: Optional[PathwayDrugMapper] = None,
        evidence_scorer: Optional[EvidenceScorer] = None,
        config: Optional[RankingConfig] = None,
    ):
        self.drug_mapper = drug_mapper
        self.evidence_scorer = evidence_scorer or EvidenceScorer()
        self.config = config or RankingConfig()

    def rank(
        self,
        pathway_scores: Dict[str, float],
        pathway_genes: Optional[Dict[str, List[str]]] = None,
        disrupted_genes: Optional[List[str]] = None,
        fired_rules: Optional[List[Any]] = None,
    ) -> RankingResult:
        """Rank therapeutic hypotheses.

        Args:
            pathway_scores: Pathway -> z-score.
            pathway_genes: Pathway -> gene list.
            disrupted_genes: Genes with damaging variants.
            fired_rules: Optional fired rules for confidence boosting.

        Returns:
            RankingResult with ranked hypotheses.
        """
        pw_genes = pathway_genes or {}
        all_hypotheses: List[TherapeuticHypothesis] = []
        drugs_seen: set = set()

        # Get drug candidates for each disrupted pathway
        pathways_analyzed = []
        for pw, score in pathway_scores.items():
            if abs(score) < 1.5:
                continue
            pathways_analyzed.append(pw)

            candidates: List[DrugCandidate] = []
            if self.drug_mapper:
                candidates = self.drug_mapper.map_pathway(pw, pw_genes.get(pw, []), disrupted_genes)

            for drug in candidates:
                drugs_seen.add(drug.drug_id)

                evidence = self.evidence_scorer.score(
                    drug_id=drug.drug_id,
                    target_pathway=pw,
                    target_genes=drug.target_genes,
                    mechanism=drug.mechanism,
                    pathway_genes=pw_genes.get(pw, []),
                    disrupted_genes=disrupted_genes,
                )

                if evidence.overall < self.config.min_evidence_score:
                    continue

                if self.config.exclude_critical_safety and evidence.has_critical_safety_flags:
                    continue

                combined = (
                    self.config.weight_evidence * evidence.overall
                    + self.config.weight_pathway_score * min(abs(score) / 5.0, 1.0)
                    + self.config.weight_drug_relevance * drug.asd_relevance_score
                    + self.config.weight_mechanism_match * (0.7 if drug.mechanism else 0.3)
                )

                explanation = (
                    f"{drug.drug_name} ({drug.mechanism or 'unknown'}) targets "
                    f"{', '.join(drug.target_genes[:3])} in {pw} "
                    f"(pathway z={score:.2f}). "
                    f"Evidence: {evidence.level.value} ({evidence.overall:.2f}). "
                    f"[Computational hypothesis — requires validation.]"
                )

                all_hypotheses.append(
                    TherapeuticHypothesis(
                        drug_id=drug.drug_id,
                        drug_name=drug.drug_name,
                        target_pathway=pw,
                        target_genes=drug.target_genes,
                        mechanism=drug.mechanism,
                        score=combined,
                        evidence=evidence,
                        explanation=explanation,
                        confidence=evidence.confidence,
                    )
                )

        # Sort by score, apply diversity, assign ranks
        all_hypotheses.sort(key=lambda h: -h.score)
        all_hypotheses = self._apply_diversity(all_hypotheses)
        for i, h in enumerate(all_hypotheses):
            h.rank = i + 1

        # Trim to max
        all_hypotheses = all_hypotheses[: self.config.max_hypotheses]

        logger.info(
            "[Autism Therapeutic] Ranked %d hypotheses across %d pathways",
            len(all_hypotheses),
            len(pathways_analyzed),
        )

        return RankingResult(
            hypotheses=all_hypotheses,
            pathways_analyzed=pathways_analyzed,
            drugs_considered=len(drugs_seen),
        )

    def _apply_diversity(
        self, hypotheses: List[TherapeuticHypothesis]
    ) -> List[TherapeuticHypothesis]:
        """Enforce per-pathway and per-mechanism caps."""
        pw_counts: Dict[str, int] = {}
        mech_counts: Dict[str, int] = {}
        filtered: List[TherapeuticHypothesis] = []

        for h in hypotheses:
            pw_c = pw_counts.get(h.target_pathway, 0)
            mech_c = mech_counts.get(h.mechanism, 0)

            if pw_c >= self.config.max_per_pathway:
                continue
            if mech_c >= self.config.max_per_mechanism:
                continue

            filtered.append(h)
            pw_counts[h.target_pathway] = pw_c + 1
            mech_counts[h.mechanism] = mech_c + 1

        return filtered
