"""
Incomplete Cascade Detection (F1).

Detects pathways where upstream initiators are active but downstream
effectors fail to respond — a partial activation indicating engineering
failure or molecular tension.

Uses KG topology to partition pathway genes into ordered layers
(upstream/intermediate/downstream), scores each layer independently,
then flags pathways where upstream > threshold but downstream < threshold.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class LayerScores:
    """Activation scores for each layer of a pathway cascade."""

    pathway: str
    upstream_score: float = 0.0
    intermediate_score: float = 0.0
    downstream_score: float = 0.0
    upstream_genes: List[str] = field(default_factory=list)
    intermediate_genes: List[str] = field(default_factory=list)
    downstream_genes: List[str] = field(default_factory=list)

    @property
    def completion_ratio(self) -> float:
        """Cascade completion ratio: 0.0 = fully stalled, 1.0 = fully resolved.

        Ratio of downstream activation to upstream activation.
        If upstream is inactive (score <= 0), cascade was never initiated — return 1.0.
        """
        if self.upstream_score <= 0:
            return 1.0
        ratio = self.downstream_score / max(self.upstream_score, 1e-10)
        return float(np.clip(ratio, 0.0, 1.0))

    @property
    def is_incomplete(self) -> bool:
        """True if upstream is active but downstream is not responding."""
        return self.upstream_score > 0.0 and self.completion_ratio < 0.5

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway": self.pathway,
            "upstream_score": round(self.upstream_score, 4),
            "intermediate_score": round(self.intermediate_score, 4),
            "downstream_score": round(self.downstream_score, 4),
            "completion_ratio": round(self.completion_ratio, 4),
            "is_incomplete": self.is_incomplete,
            "n_upstream": len(self.upstream_genes),
            "n_intermediate": len(self.intermediate_genes),
            "n_downstream": len(self.downstream_genes),
        }


@dataclass
class CascadeResult:
    """Result of cascade analysis across all pathways for a batch."""

    per_pathway: List[LayerScores]
    per_cell_ccr: Optional[pd.DataFrame] = None  # cells x pathways CCR matrix
    n_incomplete: int = 0
    mean_completion_ratio: float = 0.0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_incomplete == 0

    def get_incomplete_pathways(self) -> List[str]:
        """Return names of pathways with incomplete cascades."""
        return [ls.pathway for ls in self.per_pathway if ls.is_incomplete]

    def get_completion_ratios(self) -> Dict[str, float]:
        """Return pathway -> CCR mapping."""
        return {ls.pathway: ls.completion_ratio for ls in self.per_pathway}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_incomplete": self.n_incomplete,
            "mean_completion_ratio": round(self.mean_completion_ratio, 4),
            "summary": self.summary,
            "pathways": [ls.to_dict() for ls in self.per_pathway],
        }


class CascadeAnalyzer:
    """Detects incomplete signaling cascades using KG topology.

    Partitions pathway genes into upstream (initiators), intermediate
    (transducers), and downstream (effectors) layers using KG directed
    edges, then scores each layer to detect partial activations.

    Usage::

        from pathway_subtyping.knowledge_graph import KnowledgeGraph

        analyzer = CascadeAnalyzer(kg)
        result = analyzer.analyze(pathway_scores, pathways=["MAPK_PATHWAY"])
    """

    def __init__(
        self,
        kg: Any,  # KnowledgeGraph — Any to avoid hard import at module level
        activation_threshold: float = 0.3,
        completion_threshold: float = 0.5,
        seed: Optional[int] = None,
    ):
        """Initialize the cascade analyzer.

        Args:
            kg: KnowledgeGraph instance with pathway topology.
            activation_threshold: Minimum upstream score to consider a
                cascade as "initiated."
            completion_threshold: Minimum CCR for a cascade to be
                considered "complete."
            seed: Random seed (unused, kept for API consistency).
        """
        self.kg = kg
        self.activation_threshold = activation_threshold
        self.completion_threshold = completion_threshold
        self._layer_cache: Dict[str, Dict[str, List[str]]] = {}

    def define_cascade_layers(self, pathway: str) -> Dict[str, List[str]]:
        """Partition pathway genes into ordered layers using KG topology.

        Args:
            pathway: Pathway node ID in the KG.

        Returns:
            Dict with keys 'upstream', 'intermediate', 'downstream'.
        """
        if pathway in self._layer_cache:
            return self._layer_cache[pathway]

        layers = self.kg.partition_pathway_genes(pathway)
        self._layer_cache[pathway] = layers

        logger.info(
            "[QC Cascade] Layers for %s: %d upstream, %d intermediate, %d downstream",
            pathway,
            len(layers["upstream"]),
            len(layers["intermediate"]),
            len(layers["downstream"]),
        )
        return layers

    def score_layer_activation(
        self,
        expression: pd.DataFrame,
        layers: Dict[str, List[str]],
        gene_weights: Optional[pd.DataFrame] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Score each layer's activation across all cells.

        Computes mean z-scored expression per layer per cell, optionally
        down-weighted by a per-cell per-gene weight matrix (e.g., from
        the AlphaMissense F4 extension).

        Args:
            expression: Expression matrix (cells x genes).
            layers: Dict with 'upstream', 'intermediate', 'downstream' gene lists.
            gene_weights: Optional (cells x genes) DataFrame of per-cell,
                per-gene multiplicative weights applied to z-scored
                values before averaging. Missing cells or genes default
                to 1.0. When None, behaviour is identical to v0.5.

        Returns:
            Tuple of (upstream_scores, intermediate_scores, downstream_scores),
            each a 1D array of length n_cells.
        """
        n_cells = len(expression)

        def _layer_mean(gene_list: List[str]) -> np.ndarray:
            present = [g for g in gene_list if g in expression.columns]
            if not present:
                return np.zeros(n_cells)
            vals = expression[present].values
            # Z-score per gene, then mean across genes
            means = vals.mean(axis=0, keepdims=True)
            stds = vals.std(axis=0, keepdims=True)
            stds[stds == 0] = 1.0
            z = (vals - means) / stds
            if gene_weights is not None:
                w = gene_weights.reindex(index=expression.index, columns=present).fillna(1.0).values
                z = z * w
            return z.mean(axis=1)

        up = _layer_mean(layers["upstream"])
        mid = _layer_mean(layers["intermediate"])
        down = _layer_mean(layers["downstream"])
        return up, mid, down

    def analyze_pathway(
        self,
        expression: pd.DataFrame,
        pathway: str,
        gene_weights: Optional[pd.DataFrame] = None,
    ) -> LayerScores:
        """Analyze a single pathway for incomplete cascades.

        Args:
            expression: Expression matrix (cells x genes).
            pathway: Pathway node ID.
            gene_weights: Optional per-cell per-gene down-weighting matrix
                (see :meth:`score_layer_activation`). F4 integration hook
                for AlphaMissense-modulated cascade scoring.

        Returns:
            LayerScores with per-layer activation and completion ratio.
        """
        layers = self.define_cascade_layers(pathway)
        up, mid, down = self.score_layer_activation(expression, layers, gene_weights=gene_weights)

        return LayerScores(
            pathway=pathway,
            upstream_score=float(np.mean(up)),
            intermediate_score=float(np.mean(mid)),
            downstream_score=float(np.mean(down)),
            upstream_genes=layers["upstream"],
            intermediate_genes=layers["intermediate"],
            downstream_genes=layers["downstream"],
        )

    def analyze(
        self,
        expression: pd.DataFrame,
        pathways: Optional[List[str]] = None,
        gene_weights: Optional[pd.DataFrame] = None,
    ) -> CascadeResult:
        """Analyze all pathways for incomplete cascades.

        Args:
            expression: Expression matrix (cells x genes).
            pathways: List of pathway node IDs. If None, analyzes all
                pathways in the KG.
            gene_weights: Optional per-cell per-gene multiplicative
                weights applied during layer-score computation. Typically
                produced by
                :meth:`pathway_subtyping.qc.alphamissense.AlphaMissenseScorer.weights_from_carriers`.
                Passing None reproduces the variant-naive baseline.

        Returns:
            CascadeResult with per-pathway and per-cell cascade metrics.
        """
        if pathways is None:
            from pathway_subtyping.knowledge_graph.schema import NodeType

            pathways = self.kg.get_nodes_by_type(NodeType.PATHWAY)

        per_pathway: List[LayerScores] = []
        ccr_data: Dict[str, np.ndarray] = {}

        for pw in pathways:
            layers = self.define_cascade_layers(pw)
            if not layers["upstream"] and not layers["downstream"]:
                continue

            up, mid, down = self.score_layer_activation(
                expression, layers, gene_weights=gene_weights
            )

            ls = LayerScores(
                pathway=pw,
                upstream_score=float(np.mean(up)),
                intermediate_score=float(np.mean(mid)),
                downstream_score=float(np.mean(down)),
                upstream_genes=layers["upstream"],
                intermediate_genes=layers["intermediate"],
                downstream_genes=layers["downstream"],
            )
            per_pathway.append(ls)

            # Per-cell CCR
            with np.errstate(divide="ignore", invalid="ignore"):
                cell_ccr = np.where(
                    up > 0,
                    np.clip(down / np.maximum(up, 1e-10), 0.0, 1.0),
                    1.0,
                )
            ccr_data[pw] = cell_ccr

        n_incomplete = sum(1 for ls in per_pathway if ls.is_incomplete)
        ccrs = [ls.completion_ratio for ls in per_pathway]
        mean_ccr = float(np.mean(ccrs)) if ccrs else 1.0

        per_cell_ccr = None
        if ccr_data:
            per_cell_ccr = pd.DataFrame(ccr_data, index=expression.index)

        summary = (
            f"Cascade analysis: {len(per_pathway)} pathways, "
            f"{n_incomplete} incomplete, "
            f"mean CCR={mean_ccr:.3f}. "
            f"{'PASS' if n_incomplete == 0 else 'FAIL'}."
        )

        logger.info("[QC Cascade] %s", summary)

        return CascadeResult(
            per_pathway=per_pathway,
            per_cell_ccr=per_cell_ccr,
            n_incomplete=n_incomplete,
            mean_completion_ratio=mean_ccr,
            summary=summary,
        )
