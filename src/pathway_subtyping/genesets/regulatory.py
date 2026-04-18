"""
Regulatory-evidence gene-set expansion (v0.6 F7).

Given a user-supplied seed gene set, ``RegulatoryGeneSetExpander``
suggests candidate additions that share regulatory signal with the
seed. Two backends are supported:

    - ``BorzoiBackend`` — production path wrapping the Borzoi
      sequence-to-expression model (Linder et al. 2024). Requires the
      optional ``[genesets]`` extra (torch + borzoi + human genome
      sequence). Stub raises an informative ImportError until the
      checkpoint is packaged.
    - ``CoexpressionBackend`` — deterministic fallback using a supplied
      expression matrix to compute correlation-based regulatory
      similarity. Suitable for CI and for users without Borzoi weights.

Both backends return ranked :class:`ExpansionCandidate` objects carrying
a co-regulation score plus a flag for whether the candidate is already
a known member of any tracked curated database. This is a
*human-in-the-loop* tool — callers review and decide, never auto-merge.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Result types
# --------------------------------------------------------------------------- #

@dataclass
class ExpansionCandidate:
    """A single suggested gene with its supporting score.

    Attributes:
        gene: Candidate gene symbol.
        score: Co-regulation score in [0, 1]; higher = stronger signal
            shared with the seed set.
        in_curated_db: True iff the candidate is already listed in any
            curated database supplied via ``known_members`` — the review
            UI can deprioritise these because they add no novelty.
        source_genes: Up to ``source_k`` seed genes that contributed most
            to this candidate's score. Useful for explaining *why* a
            gene was suggested.
    """

    gene: str
    score: float
    in_curated_db: bool = False
    source_genes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gene": self.gene,
            "score": float(self.score),
            "in_curated_db": bool(self.in_curated_db),
            "source_genes": list(self.source_genes),
        }


@dataclass
class ExpansionResult:
    """Full expansion report.

    Attributes:
        seed_genes: Genes the caller supplied as the seed set.
        candidates: Ranked list (descending score) of up to ``top_n``
            suggestions, excluding genes already in the seed.
        backend: Backend identifier string.
    """

    seed_genes: List[str]
    candidates: List[ExpansionCandidate]
    backend: str

    def as_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([c.to_dict() for c in self.candidates])

    def recall_against(self, ground_truth: Iterable[str]) -> float:
        """Fraction of ``ground_truth`` members present in the suggestion list."""
        gt = set(ground_truth) - set(self.seed_genes)
        if not gt:
            return float("nan")
        suggested = {c.gene for c in self.candidates}
        return float(len(suggested & gt)) / float(len(gt))

    def top(self, n: int) -> List[ExpansionCandidate]:
        return self.candidates[:n]


# --------------------------------------------------------------------------- #
# Backend interface
# --------------------------------------------------------------------------- #

class RegulatoryBackend:
    """Abstract interface for regulatory-similarity backends.

    Subclasses expose a ``similarity_matrix(genes) -> (n_genes, n_genes)``
    method. Entries should be in [0, 1] with 1.0 on the diagonal. The
    high-level :class:`RegulatoryGeneSetExpander` consumes the matrix
    and derives candidates.
    """

    backend_id: str = "base"

    def similarity_matrix(self, genes: Sequence[str]) -> pd.DataFrame:
        raise NotImplementedError


class BorzoiBackend(RegulatoryBackend):
    """Production Borzoi wrapper — requires ``pathway-subtyping[genesets]``."""

    backend_id = "borzoi:official"

    def __init__(self, checkpoint: str = "borzoi-human-v1.0") -> None:
        self.checkpoint = checkpoint

    def _lazy_load(self) -> None:
        try:
            import torch  # noqa: F401
            import borzoi  # type: ignore  # noqa: F401
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "BorzoiBackend requires torch + borzoi + the human genome "
                "sequence. Install with: pip install pathway-subtyping[genesets]"
            ) from exc
        raise NotImplementedError(
            "Live Borzoi inference is gated on the packaged checkpoint; this "
            "stub defers to the checkpoint loader expected from the "
            "'[genesets]' extra. Use CoexpressionBackend for local testing."
        )

    def similarity_matrix(self, genes: Sequence[str]) -> pd.DataFrame:  # pragma: no cover
        self._lazy_load()
        return pd.DataFrame(np.eye(len(genes)), index=list(genes), columns=list(genes))


class CoexpressionBackend(RegulatoryBackend):
    """Deterministic fallback using an expression matrix.

    Computes a Pearson-correlation-based regulatory-similarity matrix
    scaled to [0, 1]. Requires a DataFrame of shape ``(n_samples,
    n_genes)`` at construction.
    """

    backend_id = "regulatory:coexpression"

    def __init__(self, expression: pd.DataFrame) -> None:
        if not isinstance(expression, pd.DataFrame):
            raise TypeError("CoexpressionBackend requires a pandas DataFrame")
        if len(expression) < 3:
            raise ValueError("need >= 3 samples for a meaningful correlation")
        self.expression = expression

    def similarity_matrix(self, genes: Sequence[str]) -> pd.DataFrame:
        genes = list(genes)
        present = [g for g in genes if g in self.expression.columns]
        missing = [g for g in genes if g not in self.expression.columns]
        if missing:
            logger.warning(
                "[CoexpressionBackend] %d gene(s) missing from expression "
                "matrix; they will receive zero similarity: %s",
                len(missing), missing[:5],
            )
        sub = self.expression[present]
        corr = sub.corr(method="pearson")
        # Scale Pearson r in [-1, 1] to similarity in [0, 1]
        sim = (corr + 1.0) / 2.0
        sim = sim.fillna(0.0)
        # Add missing genes as zero-similarity rows/columns
        full = pd.DataFrame(0.0, index=genes, columns=genes)
        full.loc[present, present] = sim.values
        # Diagonal = 1.0 (self-similarity)
        np.fill_diagonal(full.values, 1.0)
        return full


# --------------------------------------------------------------------------- #
# High-level expander
# --------------------------------------------------------------------------- #

class RegulatoryGeneSetExpander:
    """Suggest candidate additions to a seed gene set.

    Usage::

        backend = CoexpressionBackend(expression_df)
        expander = RegulatoryGeneSetExpander(backend=backend)
        result = expander.expand(
            seed_genes=["MYC", "MYCN"],
            candidate_genes=expression_df.columns.tolist(),
            top_n=20,
        )
        print(result.as_dataframe())
    """

    def __init__(
        self,
        backend: RegulatoryBackend,
        source_k: int = 3,
    ) -> None:
        self.backend = backend
        self.source_k = int(source_k)

    def expand(
        self,
        seed_genes: Sequence[str],
        candidate_genes: Sequence[str],
        top_n: int = 20,
        known_members: Optional[Mapping[str, Iterable[str]]] = None,
        min_score: float = 0.0,
    ) -> ExpansionResult:
        """Return the top-N candidate additions to ``seed_genes``.

        Args:
            seed_genes: The user's seed set.
            candidate_genes: Pool of genes from which to draw suggestions.
                Must be a superset of ``seed_genes`` (seed genes are
                filtered out of the final candidate list automatically).
            top_n: Maximum number of candidates to return.
            known_members: Optional mapping ``{database_name: iterable_of_genes}``.
                Candidates appearing in any listed database are flagged.
            min_score: Exclude candidates whose score is below this
                threshold.

        Returns:
            ExpansionResult with ranked candidates (descending score).
        """
        seed_genes = list(seed_genes)
        candidate_genes = list(candidate_genes)
        if not seed_genes:
            raise ValueError("seed_genes must be non-empty")
        if not candidate_genes:
            raise ValueError("candidate_genes must be non-empty")

        universe = list(dict.fromkeys(seed_genes + candidate_genes))
        sim = self.backend.similarity_matrix(universe)
        # Rows: candidates; Cols: seed similarity
        seed_cols = [g for g in seed_genes if g in sim.columns]
        if not seed_cols:
            raise ValueError(
                "none of the seed genes are present in the similarity matrix"
            )

        cand_rows = [g for g in candidate_genes if g in sim.index and g not in set(seed_genes)]
        if not cand_rows:
            return ExpansionResult(
                seed_genes=seed_genes, candidates=[], backend=self.backend.backend_id
            )

        # Score = mean similarity to all seed genes
        score = sim.loc[cand_rows, seed_cols].mean(axis=1).sort_values(ascending=False)

        # Flag membership in any known DB
        flagged: Set[str] = set()
        if known_members:
            for _db, members in known_members.items():
                flagged.update(members)

        candidates: List[ExpansionCandidate] = []
        for gene, s in score.items():
            if s < min_score:
                continue
            # Explain score: top-k source genes
            contributions = sim.loc[gene, seed_cols].sort_values(ascending=False)
            sources = contributions.head(self.source_k).index.tolist()
            candidates.append(
                ExpansionCandidate(
                    gene=str(gene),
                    score=float(s),
                    in_curated_db=gene in flagged,
                    source_genes=sources,
                )
            )
            if len(candidates) >= top_n:
                break

        logger.info(
            "[RegulatoryGeneSetExpander] backend=%s seed=%d candidates=%d",
            self.backend.backend_id, len(seed_genes), len(candidates),
        )
        return ExpansionResult(
            seed_genes=seed_genes,
            candidates=candidates,
            backend=self.backend.backend_id,
        )
