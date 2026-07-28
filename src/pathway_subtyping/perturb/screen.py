"""
Batch perturbation screen.

Takes an expression matrix, a panel of genes, and a fitted
``GeneformerPerturber`` + ``MSVFromEmbedding`` head; runs every
perturbation, and returns a ranked impact table.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

from .geneformer_wrapper import (
    GeneformerPerturber,
    PerturbationMode,
    PerturbationResult,
)
from .msv_from_embedding import MSVFromEmbedding

logger = logging.getLogger(__name__)


@dataclass
class ScreenResult:
    """Outcome of a full perturbation screen.

    Attributes:
        gene_panel: The genes that were actually perturbed (after filtering).
        mode: Knockout or overexpress (uniform across the panel).
        delta_msv_by_gene: DataFrame (n_genes, n_pathways) — mean delta-MSV
            for each gene, averaged across cells.
        l2_by_gene: Series indexed by gene — L2 norm of per-gene delta-MSV row.
        per_gene_results: Raw ``PerturbationResult`` per gene (not included
            in ``to_dict`` because the embeddings are large).
    """

    gene_panel: List[str]
    mode: PerturbationMode
    delta_msv_by_gene: pd.DataFrame
    l2_by_gene: pd.Series
    per_gene_results: Dict[str, PerturbationResult] = field(default_factory=dict)

    def rank(self, top_n: Optional[int] = None) -> pd.DataFrame:
        """Return a gene-by-impact ranking DataFrame (descending L2)."""
        ranked = self.l2_by_gene.sort_values(ascending=False)
        if top_n is not None:
            ranked = ranked.head(top_n)
        return ranked.to_frame(name="l2_delta_msv")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_genes": len(self.gene_panel),
            "mode": self.mode.value,
            "top5_impact": self.rank(5).to_dict()["l2_delta_msv"],
            "mean_impact": float(self.l2_by_gene.mean()),
        }


class PerturbationScreen:
    """Run a panel of single-gene perturbations and rank by MSV impact.

    Usage::

        head = MSVFromEmbedding().fit(baseline_emb, reference_pathway_scores)
        perturber = GeneformerPerturber()  # or OfficialBackend()
        screen = PerturbationScreen(perturber, head)
        result = screen.run(expression_df, gene_panel=["MYC", "MECP2", "TP53"])
        print(result.rank(10))
    """

    def __init__(
        self,
        perturber: GeneformerPerturber,
        head: MSVFromEmbedding,
    ) -> None:
        self.perturber = perturber
        self.head = head

    # ----------------------------------------------------------------- run ---
    def run(
        self,
        expression: pd.DataFrame,
        gene_panel: Iterable[str],
        mode: PerturbationMode = PerturbationMode.KNOCKOUT,
        skip_missing: bool = True,
    ) -> ScreenResult:
        """Perturb each gene in ``gene_panel`` and compute per-gene MSV delta.

        Args:
            expression: Cell x gene expression DataFrame.
            gene_panel: Genes to perturb one at a time.
            mode: Perturbation mode for all genes.
            skip_missing: If True, genes absent from ``expression.columns``
                are silently dropped. If False, raises KeyError on the
                first missing gene.

        Returns:
            ScreenResult with per-gene delta-MSV and L2 impact.
        """
        mode = PerturbationMode(mode) if isinstance(mode, str) else mode
        gene_panel = list(gene_panel)
        available_genes = set(expression.columns)

        valid_genes: List[str] = []
        for g in gene_panel:
            if g in available_genes:
                valid_genes.append(g)
            else:
                if skip_missing:
                    logger.warning(
                        "[PerturbationScreen] gene %r not in expression columns; skipping",
                        g,
                    )
                else:
                    raise KeyError(f"gene {g!r} not in expression columns")

        if not valid_genes:
            raise ValueError("no valid genes in panel after filtering")

        # Baseline embedding is computed once and reused across perturbations.
        baseline_embedding = self.perturber.embed(expression)
        baseline_msv = self.head.transform(baseline_embedding)
        pathway_names = list(baseline_msv.columns)

        rows: List[np.ndarray] = []
        per_gene: Dict[str, PerturbationResult] = {}
        for gene in valid_genes:
            result = self.perturber.perturb(expression, gene, mode)
            per_gene[gene] = result
            perturbed_msv = self.head.transform(result.perturbed_embedding)
            gene_delta = (perturbed_msv - baseline_msv).mean(axis=0).to_numpy()
            rows.append(gene_delta)

        delta_msv_by_gene = pd.DataFrame(
            np.vstack(rows),
            index=valid_genes,
            columns=pathway_names,
        )
        l2_by_gene = pd.Series(
            np.linalg.norm(delta_msv_by_gene.to_numpy(), axis=1),
            index=valid_genes,
            name="l2_delta_msv",
        )

        logger.info(
            "[PerturbationScreen] screened %d genes (%s); top = %s",
            len(valid_genes),
            mode.value,
            l2_by_gene.idxmax() if len(l2_by_gene) else None,
        )
        return ScreenResult(
            gene_panel=valid_genes,
            mode=mode,
            delta_msv_by_gene=delta_msv_by_gene,
            l2_by_gene=l2_by_gene,
            per_gene_results=per_gene,
        )
