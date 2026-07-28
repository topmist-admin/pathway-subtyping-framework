"""
Proteomics pathway scoring.

Proteomics matrices typically have fewer features than scRNA-seq (3k–5k
proteins vs 20k+ transcripts) and different dynamic range. We score
pathways as the mean-Z of protein-level abundance for pathway members,
matching the API shape of the RNA and ATAC scorers so fusion logic can
treat all three modalities uniformly.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, Mapping, Optional

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ProteomicsScorer:
    """Score pathways from a protein-abundance matrix.

    Attributes:
        gene_to_protein: Optional mapping from gene symbol → protein
            accession. When provided, pathway membership (expressed in
            gene symbols) is translated into protein accessions before
            scoring. When ``None``, pathway members are assumed to be
            protein-matrix column names directly.
    """

    gene_to_protein: Optional[Mapping[str, str]] = None

    # ---------------------------------------------------------- score ---
    def score(
        self,
        abundance: pd.DataFrame,
        pathways: Mapping[str, Iterable[str]],
        min_proteins_per_pathway: int = 2,
    ) -> pd.DataFrame:
        """Return a sample-by-pathway protein-score matrix.

        Args:
            abundance: Samples × proteins (log / VST-normalised).
            pathways: Mapping pathway → list of gene symbols (or direct
                protein IDs if ``gene_to_protein`` is None).
            min_proteins_per_pathway: Skip pathways with fewer members
                after the gene → protein translation.

        Returns:
            DataFrame (samples × pathways), Z-normalised per pathway.
        """
        if not isinstance(abundance, pd.DataFrame):
            raise TypeError("abundance must be a pandas DataFrame")

        available = set(abundance.columns)

        def _translate(members: Iterable[str]) -> list:
            if self.gene_to_protein is None:
                return [m for m in members if m in available]
            return [
                self.gene_to_protein[m]
                for m in members
                if m in self.gene_to_protein and self.gene_to_protein[m] in available
            ]

        records = {}
        for pw, members in pathways.items():
            translated = _translate(members)
            if len(translated) < min_proteins_per_pathway:
                continue
            block = abundance[translated]
            records[pw] = block.mean(axis=1)
        if not records:
            return pd.DataFrame(index=abundance.index)

        scored = pd.DataFrame(records, index=abundance.index)
        mean = scored.mean(axis=0)
        std = scored.std(axis=0).replace(0.0, 1.0)
        return (scored - mean) / std


# --------------------------------------------------------------------------- #
# Convenience entry point
# --------------------------------------------------------------------------- #


def score_proteomics_pathways(
    abundance: pd.DataFrame,
    pathways: Mapping[str, Iterable[str]],
    gene_to_protein: Optional[Mapping[str, str]] = None,
    min_proteins_per_pathway: int = 2,
) -> pd.DataFrame:
    """One-shot helper mirroring ``score_pathways_from_expression``."""
    return ProteomicsScorer(gene_to_protein=gene_to_protein).score(
        abundance=abundance,
        pathways=pathways,
        min_proteins_per_pathway=min_proteins_per_pathway,
    )
