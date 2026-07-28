"""
Environmental Stress Fingerprinting (F11).

Maps QC failures to fixable manufacturing conditions by matching observed
pathway signatures to known stressor profiles.

Stressors: hypoxia, nutrient depletion, oxidative, mechanical, thermal, pH.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class StressSignature:
    """Reference signature for a known stressor."""

    name: str
    activated_pathways: List[str] = field(default_factory=list)
    suppressed_pathways: List[str] = field(default_factory=list)
    remediation: str = ""


@dataclass
class StressMatch:
    """Match result for a single stressor."""

    stressor: str
    similarity: float
    confidence: float
    remediation: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stressor": self.stressor,
            "similarity": round(self.similarity, 4),
            "confidence": round(self.confidence, 4),
            "remediation": self.remediation,
        }


@dataclass
class StressResult:
    """Result of stress fingerprinting."""

    matches: List[StressMatch]
    primary_stressor: Optional[str] = None
    confidence: float = 0.0
    remediation: str = ""
    summary: str = ""

    @property
    def stress_detected(self) -> bool:
        return self.primary_stressor is not None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stress_detected": self.stress_detected,
            "primary_stressor": self.primary_stressor,
            "confidence": round(self.confidence, 4),
            "remediation": self.remediation,
            "summary": self.summary,
            "matches": [m.to_dict() for m in self.matches],
        }


class StressFingerprinter:
    """Fingerprints environmental stress from pathway activation patterns.

    Usage::

        fp = StressFingerprinter()
        fp.load_default_signatures()
        result = fp.fingerprint(pathway_scores)
    """

    def __init__(
        self,
        signatures: Optional[List[StressSignature]] = None,
        min_confidence: float = 0.3,
        seed: Optional[int] = None,
    ):
        self.signatures = signatures or []
        self.min_confidence = min_confidence

    def load_default_signatures(self) -> None:
        """Load the 6 default stress signatures."""
        self.signatures = [
            StressSignature(
                name="hypoxia",
                activated_pathways=[
                    "HALLMARK_HYPOXIA",
                    "HALLMARK_GLYCOLYSIS",
                    "HALLMARK_ANGIOGENESIS",
                ],
                suppressed_pathways=[
                    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                    "HALLMARK_FATTY_ACID_METABOLISM",
                ],
                remediation="Check dissolved oxygen in bioreactor; verify gas exchange rate",
            ),
            StressSignature(
                name="nutrient_depletion",
                activated_pathways=["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_P53_PATHWAY"],
                suppressed_pathways=["HALLMARK_MTORC1_SIGNALING", "HALLMARK_ADIPOGENESIS"],
                remediation="Verify media replenishment schedule; check glucose/glutamine levels",
            ),
            StressSignature(
                name="oxidative",
                activated_pathways=["HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS"],
                suppressed_pathways=["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
                remediation="Check antioxidant supplementation; verify ROS scavenging capacity",
            ),
            StressSignature(
                name="mechanical",
                activated_pathways=[
                    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "HALLMARK_MYOGENESIS",
                ],
                suppressed_pathways=["HALLMARK_ADIPOGENESIS"],
                remediation="Reduce bioreactor agitation speed; check shear stress parameters",
            ),
            StressSignature(
                name="thermal",
                activated_pathways=["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_P53_PATHWAY"],
                suppressed_pathways=["HALLMARK_MTORC1_SIGNALING"],
                remediation="Verify incubator temperature stability; check for thermal cycling",
            ),
            StressSignature(
                name="ph",
                activated_pathways=["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
                suppressed_pathways=[
                    "HALLMARK_MTORC1_SIGNALING",
                    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                ],
                remediation="Check media pH buffering capacity; verify CO2 levels",
            ),
        ]

    def fingerprint(
        self,
        pathway_scores: pd.DataFrame,
        reference_scores: Optional[pd.DataFrame] = None,
    ) -> StressResult:
        """Match observed pathway pattern to known stress signatures.

        Args:
            pathway_scores: DataFrame (cells x pathways) from the batch.
            reference_scores: Optional control/reference for comparison.

        Returns:
            StressResult with ranked stressor matches.
        """
        if not self.signatures:
            return StressResult(matches=[], summary="No stress signatures loaded.")

        batch_mean = pathway_scores.mean()
        if reference_scores is not None:
            ref_mean = reference_scores.mean()
            diff = batch_mean - ref_mean
        else:
            diff = batch_mean

        matches: List[StressMatch] = []
        for sig in self.signatures:
            sim = self._compute_similarity(diff, sig)
            matches.append(
                StressMatch(
                    stressor=sig.name,
                    similarity=sim,
                    confidence=max(0.0, sim),
                    remediation=sig.remediation,
                )
            )

        matches.sort(key=lambda m: -m.similarity)

        primary = None
        confidence = 0.0
        remediation = ""
        if matches and matches[0].similarity >= self.min_confidence:
            primary = matches[0].stressor
            confidence = matches[0].confidence
            remediation = matches[0].remediation

        summary = (
            f"Stress fingerprint: {primary or 'none detected'}"
            + (f" (confidence={confidence:.2f})" if primary else "")
            + f". {len(matches)} signatures tested."
        )

        logger.info("[QC Stress] %s", summary)

        return StressResult(
            matches=matches,
            primary_stressor=primary,
            confidence=confidence,
            remediation=remediation,
            summary=summary,
        )

    def _compute_similarity(
        self,
        diff: pd.Series,
        sig: StressSignature,
    ) -> float:
        """Cosine-like similarity between observed changes and signature."""
        score = 0.0
        n = 0

        for pw in sig.activated_pathways:
            if pw in diff.index:
                score += diff[pw]
                n += 1

        for pw in sig.suppressed_pathways:
            if pw in diff.index:
                score -= diff[pw]
                n += 1

        return score / max(n, 1)
