"""
Evidence Scoring for Therapeutic Hypotheses.

Multi-criteria scoring: biological plausibility, mechanistic alignment,
literature support, clinical evidence, safety flags.

Critical safety features for pediatric autism population.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class EvidenceLevel(Enum):
    """Evidence strength classification."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INSUFFICIENT = "insufficient"


class SafetyFlag(Enum):
    """Safety concern categories — critical for pediatric ASD population."""

    CNS_EFFECTS = "cns_effects"
    DEVELOPMENTAL_CONCERNS = "developmental_concerns"
    DRUG_INTERACTIONS = "drug_interactions"
    CONTRAINDICATED_PEDIATRIC = "contraindicated_pediatric"
    BLACK_BOX_WARNING = "black_box_warning"
    OFF_LABEL_USE = "off_label_use"
    IMMUNOSUPPRESSION = "immunosuppression"
    HEPATOTOXICITY = "hepatotoxicity"
    CARDIOTOXICITY = "cardiotoxicity"
    TERATOGENIC = "teratogenic"
    WITHDRAWAL_RISK = "withdrawal_risk"


@dataclass
class EvidenceScore:
    """Multi-criteria evidence score for a therapeutic hypothesis."""

    biological_plausibility: float = 0.0
    mechanistic_alignment: float = 0.0
    literature_support: float = 0.0
    clinical_evidence: float = 0.0
    safety_flags: List[str] = field(default_factory=list)
    overall: float = 0.0
    confidence: float = 0.0
    level: EvidenceLevel = EvidenceLevel.INSUFFICIENT
    explanation: str = ""

    @property
    def has_critical_safety_flags(self) -> bool:
        critical = {
            SafetyFlag.CONTRAINDICATED_PEDIATRIC.value,
            SafetyFlag.TERATOGENIC.value,
            SafetyFlag.BLACK_BOX_WARNING.value,
        }
        return bool(set(self.safety_flags) & critical)

    @property
    def safety_summary(self) -> str:
        if not self.safety_flags:
            return "No safety concerns identified"
        return f"{len(self.safety_flags)} flag(s): {', '.join(self.safety_flags)}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "biological_plausibility": round(self.biological_plausibility, 4),
            "mechanistic_alignment": round(self.mechanistic_alignment, 4),
            "literature_support": round(self.literature_support, 4),
            "clinical_evidence": round(self.clinical_evidence, 4),
            "safety_flags": self.safety_flags,
            "has_critical_safety": self.has_critical_safety_flags,
            "overall": round(self.overall, 4),
            "confidence": round(self.confidence, 4),
            "level": self.level.value,
            "explanation": self.explanation,
        }


@dataclass
class EvidenceScorerConfig:
    """Weights and thresholds for evidence scoring."""

    weight_biological: float = 0.30
    weight_mechanistic: float = 0.25
    weight_literature: float = 0.25
    weight_clinical: float = 0.20
    threshold_high: float = 0.70
    threshold_moderate: float = 0.50
    threshold_low: float = 0.30
    safety_penalty_per_flag: float = 0.05
    max_safety_penalty: float = 0.30


class EvidenceScorer:
    """Scores therapeutic hypotheses on multiple evidence criteria.

    Usage::

        scorer = EvidenceScorer()
        score = scorer.score(
            drug_id="SORAFENIB",
            target_pathway="MAPK",
            target_genes=["RAF", "MEK"],
            mechanism="inhibitor",
        )
    """

    def __init__(
        self,
        config: Optional[EvidenceScorerConfig] = None,
        literature_db: Optional[Dict[str, float]] = None,
        clinical_trials_db: Optional[Dict[str, float]] = None,
        safety_db: Optional[Dict[str, List[str]]] = None,
    ):
        self.config = config or EvidenceScorerConfig()
        self.literature_db = literature_db or {}
        self.clinical_trials_db = clinical_trials_db or {}
        self.safety_db = safety_db or {}

    def score(
        self,
        drug_id: str,
        target_pathway: str,
        target_genes: List[str],
        mechanism: str = "",
        pathway_genes: Optional[List[str]] = None,
        disrupted_genes: Optional[List[str]] = None,
        rule_confidence: Optional[float] = None,
    ) -> EvidenceScore:
        """Score a therapeutic hypothesis.

        Args:
            drug_id: Drug identifier.
            target_pathway: Target pathway.
            target_genes: Genes targeted by the drug.
            mechanism: Drug mechanism (inhibitor, agonist, etc.).
            pathway_genes: All genes in the target pathway.
            disrupted_genes: Genes with damaging variants.
            rule_confidence: Confidence from symbolic rule that generated this hypothesis.

        Returns:
            EvidenceScore with multi-criteria assessment.
        """
        bio = self._score_biological_plausibility(
            target_genes, pathway_genes or [], disrupted_genes or []
        )
        mech = self._score_mechanistic_alignment(mechanism, target_pathway)
        lit = self._score_literature(drug_id)
        clin = self._score_clinical(drug_id)
        flags = self._get_safety_flags(drug_id, mechanism)

        # Weighted overall
        cfg = self.config
        overall = (
            cfg.weight_biological * bio
            + cfg.weight_mechanistic * mech
            + cfg.weight_literature * lit
            + cfg.weight_clinical * clin
        )

        # Safety penalty
        penalty = min(
            len(flags) * cfg.safety_penalty_per_flag,
            cfg.max_safety_penalty,
        )
        overall = max(0.0, overall - penalty)

        # Confidence
        confidence = overall
        if rule_confidence is not None:
            confidence = (confidence + rule_confidence) / 2.0

        # Level
        if overall >= cfg.threshold_high:
            level = EvidenceLevel.HIGH
        elif overall >= cfg.threshold_moderate:
            level = EvidenceLevel.MODERATE
        elif overall >= cfg.threshold_low:
            level = EvidenceLevel.LOW
        else:
            level = EvidenceLevel.INSUFFICIENT

        explanation = (
            f"Drug {drug_id} targeting {target_pathway} via {mechanism or 'unknown mechanism'}. "
            f"Bio={bio:.2f}, Mech={mech:.2f}, Lit={lit:.2f}, Clin={clin:.2f}. "
            f"Overall={overall:.2f} ({level.value}). "
            f"{len(flags)} safety flag(s)."
        )

        return EvidenceScore(
            biological_plausibility=bio,
            mechanistic_alignment=mech,
            literature_support=lit,
            clinical_evidence=clin,
            safety_flags=flags,
            overall=overall,
            confidence=confidence,
            level=level,
            explanation=explanation,
        )

    def _score_biological_plausibility(
        self, target_genes: List[str], pathway_genes: List[str], disrupted_genes: List[str]
    ) -> float:
        if not target_genes:
            return 0.0
        overlap = set(target_genes) & set(disrupted_genes)
        if overlap:
            return min(1.0, 0.5 + 0.1 * len(overlap))
        pw_overlap = set(target_genes) & set(pathway_genes)
        if pw_overlap:
            return min(1.0, 0.3 + 0.05 * len(pw_overlap))
        return 0.1

    def _score_mechanistic_alignment(self, mechanism: str, pathway: str) -> float:
        """Heuristic prior on drug mechanism. ⚠️ Does NOT use ``pathway``.

        Despite "alignment", this scores the mechanism STRING alone against a
        fixed lookup (inhibitor 0.7 / agonist 0.6 / modulator 0.5 / other 0.3 /
        missing 0.2); no property of ``pathway`` enters. The constants are an
        unsourced prior, not a fitted or cited quantity. Treat the output as a
        mechanism-type prior and not as evidence that a drug mechanism matches a
        particular pathway -- ranking built on it inherits that limit.
        """
        if not mechanism:
            return 0.2
        mechanism_lower = mechanism.lower()
        if mechanism_lower in ("inhibitor", "antagonist", "blocker"):
            return 0.7
        if mechanism_lower in ("agonist", "activator", "enhancer"):
            return 0.6
        if mechanism_lower in ("modulator", "stabilizer"):
            return 0.5
        return 0.3

    def _score_literature(self, drug_id: str) -> float:
        return self.literature_db.get(drug_id, 0.0)

    def _score_clinical(self, drug_id: str) -> float:
        return self.clinical_trials_db.get(drug_id, 0.0)

    def _get_safety_flags(self, drug_id: str, mechanism: str) -> List[str]:
        flags = list(self.safety_db.get(drug_id, []))
        # Auto-flag CNS effects for CNS-active drugs
        if mechanism.lower() in ("agonist", "antagonist", "blocker") and not flags:
            flags.append(SafetyFlag.CNS_EFFECTS.value)
        return flags

    def add_literature_evidence(self, drug_id: str, score: float) -> None:
        self.literature_db[drug_id] = score

    def add_clinical_evidence(self, drug_id: str, score: float) -> None:
        self.clinical_trials_db[drug_id] = score

    def add_safety_flags(self, drug_id: str, flags: List[str]) -> None:
        self.safety_db[drug_id] = flags
