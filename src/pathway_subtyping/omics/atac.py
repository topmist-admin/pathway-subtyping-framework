"""
ATAC-seq pathway scoring.

Given a peak-by-sample accessibility matrix and a peak→gene mapping,
compute pathway-level activity as the sum / mean of accessibility at
peaks whose linked genes are members of each pathway. The result is
a sample-by-pathway score matrix that uses the same API as RNA-side
scoring so fusion logic can treat them interchangeably.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


AggregationMethod = str  # "mean", "sum", "max"
_VALID_AGG = {"mean", "sum", "max"}


@dataclass
class ATACScorer:
    """Score pathways from an ATAC peak-accessibility matrix.

    Attributes:
        peak_to_gene: Mapping from peak ID → linked gene symbol. Peaks
            may link to multiple genes if the upstream pipeline produced
            a multi-gene peak; in that case pass a mapping to a list
            (see :meth:`score`).
        aggregation: How to combine per-peak accessibility into a per-gene
            score before averaging across genes in a pathway. ``'mean'``
            (default) is robust; ``'sum'`` emphasises accessibility mass;
            ``'max'`` caps outliers.
    """

    peak_to_gene: Mapping[str, str]
    aggregation: AggregationMethod = "mean"

    def __post_init__(self) -> None:
        if self.aggregation not in _VALID_AGG:
            raise ValueError(
                f"aggregation must be one of {sorted(_VALID_AGG)}"
            )

    def _gene_scores(self, accessibility: pd.DataFrame) -> pd.DataFrame:
        """Collapse peak-level accessibility to gene-level scores."""
        if not isinstance(accessibility, pd.DataFrame):
            raise TypeError("accessibility must be a pandas DataFrame")
        # Build peak → gene index; skip peaks absent from matrix
        present = [p for p in self.peak_to_gene if p in accessibility.columns]
        if not present:
            return pd.DataFrame(index=accessibility.index)
        grouped: Dict[str, List[str]] = {}
        for p in present:
            g = self.peak_to_gene[p]
            grouped.setdefault(g, []).append(p)
        out = {}
        for gene, peaks in grouped.items():
            block = accessibility[peaks]
            if self.aggregation == "mean":
                out[gene] = block.mean(axis=1)
            elif self.aggregation == "sum":
                out[gene] = block.sum(axis=1)
            else:
                out[gene] = block.max(axis=1)
        return pd.DataFrame(out, index=accessibility.index)

    # ------------------------------------------------------------ score ---
    def score(
        self,
        accessibility: pd.DataFrame,
        pathways: Mapping[str, Iterable[str]],
        min_genes_per_pathway: int = 2,
    ) -> pd.DataFrame:
        """Return a sample-by-pathway accessibility score matrix.

        Args:
            accessibility: Samples × peaks matrix (log-normalized /
                TF-IDF-normalised — caller's responsibility).
            pathways: Mapping pathway → list of gene symbols.
            min_genes_per_pathway: Skip pathways with fewer genes after
                peak→gene mapping.

        Returns:
            DataFrame (samples × pathways). Z-normalised per pathway so
            the output is comparable to RNA-side scores.
        """
        gene_scores = self._gene_scores(accessibility)
        if gene_scores.empty:
            logger.warning(
                "[ATACScorer] no peaks in accessibility matched peak_to_gene"
            )
            return pd.DataFrame(index=accessibility.index)

        available = set(gene_scores.columns)
        records: Dict[str, pd.Series] = {}
        for pw, genes in pathways.items():
            members = [g for g in genes if g in available]
            if len(members) < min_genes_per_pathway:
                continue
            block = gene_scores[members]
            records[pw] = block.mean(axis=1)
        if not records:
            return pd.DataFrame(index=accessibility.index)

        scored = pd.DataFrame(records, index=accessibility.index)
        # Z-normalise per pathway for comparability with RNA scoring.
        mean = scored.mean(axis=0)
        std = scored.std(axis=0).replace(0.0, 1.0)
        return (scored - mean) / std


# --------------------------------------------------------------------------- #
# Convenience entry point
# --------------------------------------------------------------------------- #

def score_atac_pathways(
    accessibility: pd.DataFrame,
    peak_to_gene: Mapping[str, str],
    pathways: Mapping[str, Iterable[str]],
    aggregation: AggregationMethod = "mean",
    min_genes_per_pathway: int = 2,
) -> pd.DataFrame:
    """One-shot helper mirroring ``score_pathways_from_expression``."""
    return ATACScorer(
        peak_to_gene=peak_to_gene, aggregation=aggregation
    ).score(
        accessibility=accessibility,
        pathways=pathways,
        min_genes_per_pathway=min_genes_per_pathway,
    )
