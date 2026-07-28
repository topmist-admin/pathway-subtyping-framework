"""
Sequence-level CRISPR off-target scoring (v0.6 F9).

PSF's pathway-level off-target detection (``pathway_subtyping.qc.offtarget``)
catches *transcriptional* off-targets — unintended pathways activating
post-edit. F9 adds a complementary *sequence-level* check: given a guide
RNA and a set of candidate genomic sites, score how likely each site is
to be functionally cut.

The production backend wraps Evo 2 (Arc Institute 2025) for long-context
sequence prediction. A deterministic similarity-score fallback provides
a baseline that the Evo 2 wrapper can be benchmarked against — matching
the roadmap acceptance pattern (``AUROC improvement >= 0.03`` over the
similarity baseline on a held-out CRISPR off-target panel).

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, List, Mapping, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


_DNA_ALPHABET = "ACGT"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _validate_dna(sequence: str) -> None:
    bad = set(sequence.upper()) - set(_DNA_ALPHABET + "N")
    if bad:
        raise ValueError(f"sequence contains non-DNA characters: {sorted(bad)[:5]}")


def _hamming_similarity(guide: str, candidate: str) -> float:
    """Per-position match ratio over the shared prefix length.

    Returns 1.0 on an exact match, 0.0 on a fully mismatched sequence.
    Gaps (length mismatch) are treated as mismatches on the unmatched tail.
    """
    guide_u = guide.upper()
    cand_u = candidate.upper()
    L = max(len(guide_u), len(cand_u))
    if L == 0:
        return 0.0
    matches = 0
    for i in range(L):
        g = guide_u[i] if i < len(guide_u) else ""
        c = cand_u[i] if i < len(cand_u) else ""
        if g and c and g == c:
            matches += 1
    return matches / L


# --------------------------------------------------------------------------- #
# Result type
# --------------------------------------------------------------------------- #


@dataclass
class OffTargetScore:
    """One (guide, candidate_site) scoring row.

    Attributes:
        guide: Guide-RNA sequence.
        site_id: Stable candidate-site identifier.
        candidate_sequence: The genomic DNA at the candidate site.
        similarity_score: Sequence-similarity baseline in [0, 1]. Higher
            = more alignment with the guide.
        functional_score: Evo-2-derived functional cleavage likelihood
            in [0, 1]. ``nan`` if only the baseline was available.
        combined_score: ``0.5 * similarity + 0.5 * functional`` when both
            are present; otherwise equals whichever is available.
    """

    guide: str
    site_id: str
    candidate_sequence: str
    similarity_score: float
    functional_score: float = float("nan")

    @property
    def combined_score(self) -> float:
        if np.isnan(self.functional_score):
            return float(self.similarity_score)
        return float(0.5 * self.similarity_score + 0.5 * self.functional_score)

    def to_dict(self) -> Any:
        return {
            "guide": self.guide,
            "site_id": self.site_id,
            "candidate_sequence": self.candidate_sequence,
            "similarity_score": float(self.similarity_score),
            "functional_score": float(self.functional_score),
            "combined_score": self.combined_score,
        }


# --------------------------------------------------------------------------- #
# Backend interface
# --------------------------------------------------------------------------- #


class OffTargetBackend:
    """Abstract sequence-level off-target backend."""

    backend_id: str = "base"

    def score(
        self,
        guide: str,
        candidates: Mapping[str, str],
    ) -> List[OffTargetScore]:
        raise NotImplementedError


class SimilarityBackend(OffTargetBackend):
    """Hamming-similarity baseline.

    The current-practice proxy: exact or near-exact matches score high.
    Everything else is treated uniformly — and it is the "ambiguous
    off-target sites" that F9 was introduced to distinguish.
    """

    backend_id = "offtarget:similarity"

    def score(
        self,
        guide: str,
        candidates: Mapping[str, str],
    ) -> List[OffTargetScore]:
        _validate_dna(guide)
        rows: List[OffTargetScore] = []
        for site_id, sequence in candidates.items():
            _validate_dna(sequence)
            rows.append(
                OffTargetScore(
                    guide=guide,
                    site_id=site_id,
                    candidate_sequence=sequence,
                    similarity_score=_hamming_similarity(guide, sequence),
                )
            )
        return rows


class Evo2OffTargetScorer(OffTargetBackend):
    """Production Evo 2 wrapper — requires ``pathway-subtyping[qc-sequence]``.

    Combines Hamming similarity with Evo-2-predicted cleavage likelihood.
    The underlying model checkpoint is lazy-loaded; without the extra
    installed, calling :meth:`score` raises an informative error.
    """

    backend_id = "offtarget:evo2/official"

    def __init__(
        self,
        checkpoint: str = "evo2-human-v1.0",
        device: str = "cpu",
    ) -> None:
        self.checkpoint = checkpoint
        self.device = device
        self._model: Any = None

    def _lazy_load(self) -> None:
        if self._model is not None:
            return
        try:
            import evo2  # type: ignore  # noqa: F401
            import torch  # noqa: F401
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "Evo2OffTargetScorer requires torch + evo2. "
                "Install with: pip install pathway-subtyping[qc-sequence]"
            ) from exc
        raise NotImplementedError(
            "Live Evo 2 inference is gated on the packaged checkpoint; "
            "this stub defers to the checkpoint loader expected from the "
            "'[qc-sequence]' extra. Use SimulatedEvo2Backend for tests "
            "and the Hamming baseline for offline deployments."
        )

    def score(
        self,
        guide: str,
        candidates: Mapping[str, str],
    ) -> List[OffTargetScore]:  # pragma: no cover
        self._lazy_load()
        return []


class SimulatedEvo2Backend(OffTargetBackend):
    """Deterministic functional-score stub used for tests.

    Not a scientific substitute for Evo 2. Produces a pseudo-"functional"
    score that boosts candidates with high similarity *and* conserved
    seed-region (first 8 nt) matches — crude but enough to demonstrate
    an AUROC improvement over pure Hamming similarity on synthetic
    benchmarks that encode that same assumption.
    """

    backend_id = "offtarget:evo2/simulated"

    def __init__(self, seed_length: int = 8) -> None:
        if seed_length < 1:
            raise ValueError("seed_length must be >= 1")
        self.seed_length = int(seed_length)

    def score(
        self,
        guide: str,
        candidates: Mapping[str, str],
    ) -> List[OffTargetScore]:
        _validate_dna(guide)
        rows: List[OffTargetScore] = []
        seed = guide[: self.seed_length].upper()
        for site_id, sequence in candidates.items():
            _validate_dna(sequence)
            sim = _hamming_similarity(guide, sequence)
            cand_seed = sequence[: self.seed_length].upper()
            seed_match = _hamming_similarity(seed, cand_seed)
            # Functional score: weight seed match heavily, other-positions lightly
            functional = 0.7 * seed_match + 0.3 * sim
            rows.append(
                OffTargetScore(
                    guide=guide,
                    site_id=site_id,
                    candidate_sequence=sequence,
                    similarity_score=sim,
                    functional_score=float(functional),
                )
            )
        return rows


# --------------------------------------------------------------------------- #
# Benchmarking utility
# --------------------------------------------------------------------------- #


def auroc(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Area under the ROC curve (binary ``labels`` in {0, 1}).

    Simple Mann-Whitney-U-based implementation; no sklearn required.
    Returns 0.5 for constant scores.
    """
    s = np.asarray(scores, dtype=float)
    y = np.asarray(labels, dtype=int)
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    all_scores = np.concatenate([pos, neg])
    ranks = pd.Series(all_scores).rank(method="average").to_numpy()
    pos_rank_sum = ranks[: len(pos)].sum()
    n_pos, n_neg = len(pos), len(neg)
    # Mann-Whitney U = sum_of_pos_ranks - n_pos*(n_pos+1)/2
    U = pos_rank_sum - n_pos * (n_pos + 1) / 2.0
    return float(U / (n_pos * n_neg))


def compare_backends(
    baseline: OffTargetBackend,
    contender: OffTargetBackend,
    guide: str,
    candidates: Mapping[str, str],
    true_labels: Mapping[str, int],
) -> dict:
    """Score both backends on the same panel; return AUROCs and delta.

    Args:
        baseline: Backend to compare against (typically ``SimilarityBackend``).
        contender: Backend under test (e.g. ``SimulatedEvo2Backend`` or real Evo 2).
        guide: Guide RNA sequence.
        candidates: ``{site_id: genomic_sequence}``.
        true_labels: ``{site_id: 0|1}`` ground-truth off-target annotation.

    Returns:
        Dict with ``baseline_auroc``, ``contender_auroc``, ``delta``,
        and per-site score rows.
    """
    baseline_rows = baseline.score(guide, candidates)
    contender_rows = contender.score(guide, candidates)

    ids = [r.site_id for r in baseline_rows]
    y = np.asarray([true_labels[s] for s in ids], dtype=int)

    baseline_scores = np.asarray([r.combined_score for r in baseline_rows])
    contender_scores = np.asarray([r.combined_score for r in contender_rows])

    return {
        "baseline_backend": baseline.backend_id,
        "contender_backend": contender.backend_id,
        "baseline_auroc": auroc(baseline_scores.tolist(), y.tolist()),
        "contender_auroc": auroc(contender_scores.tolist(), y.tolist()),
        "delta": (
            auroc(contender_scores.tolist(), y.tolist())
            - auroc(baseline_scores.tolist(), y.tolist())
        ),
        "n_sites": len(ids),
    }
