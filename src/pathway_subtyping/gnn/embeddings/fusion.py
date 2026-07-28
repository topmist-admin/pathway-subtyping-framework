"""
Multi-source Embedding Fusion.

Combines embeddings from multiple sources (KG embeddings, expression
features, foundation models) into unified node feature vectors.

Pure numpy — no PyTorch required for fusion operations.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np
from sklearn.decomposition import PCA

from ..kg_embeddings import NodeEmbeddings

logger = logging.getLogger(__name__)


class FusionMethod(Enum):
    """Methods for fusing multiple embedding sources."""

    CONCAT = "concat"
    WEIGHTED_SUM = "weighted_sum"
    AVERAGE = "average"
    PCA = "pca"


@dataclass
class FusionConfig:
    """Configuration for embedding fusion."""

    method: FusionMethod = FusionMethod.CONCAT
    weights: Optional[Dict[str, float]] = None
    output_dim: Optional[int] = None  # For PCA
    normalize: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "method": self.method.value,
            "output_dim": self.output_dim,
            "normalize": self.normalize,
        }


class EmbeddingFusion:
    """Fuses embeddings from multiple sources into unified vectors.

    Usage::

        fusion = EmbeddingFusion(FusionConfig(method=FusionMethod.CONCAT))
        result = fusion.fuse({
            "kg": kg_embeddings,
            "expression": expr_embeddings,
        })
    """

    def __init__(self, config: Optional[FusionConfig] = None):
        self.config = config or FusionConfig()

    def fuse(
        self,
        embeddings: Dict[str, NodeEmbeddings],
        node_ids: Optional[List[str]] = None,
    ) -> NodeEmbeddings:
        """Fuse multiple embedding sources.

        Args:
            embeddings: Dict mapping source_name -> NodeEmbeddings.
            node_ids: Optional ordered list of node IDs. If None, uses
                the intersection of all sources.

        Returns:
            Fused NodeEmbeddings.
        """
        if not embeddings:
            raise ValueError("No embeddings to fuse.")

        if len(embeddings) == 1:
            return list(embeddings.values())[0]

        # Determine common node IDs
        if node_ids is None:
            id_sets = [set(emb.node_ids) for emb in embeddings.values()]
            common = id_sets[0]
            for s in id_sets[1:]:
                common = common & s
            node_ids = sorted(common)

        if not node_ids:
            raise ValueError("No common node IDs across embedding sources.")

        # Collect aligned matrices
        matrices: Dict[str, np.ndarray] = {}
        for source_name, emb in embeddings.items():
            id_to_idx = {nid: i for i, nid in enumerate(emb.node_ids)}
            indices = [id_to_idx[nid] for nid in node_ids if nid in id_to_idx]
            if len(indices) == len(node_ids):
                matrices[source_name] = emb.embeddings[indices]

        if not matrices:
            raise ValueError("Could not align embeddings across sources.")

        # Fuse
        if self.config.method == FusionMethod.CONCAT:
            fused = self._fuse_concat(matrices)
        elif self.config.method == FusionMethod.WEIGHTED_SUM:
            fused = self._fuse_weighted_sum(matrices)
        elif self.config.method == FusionMethod.AVERAGE:
            fused = self._fuse_average(matrices)
        elif self.config.method == FusionMethod.PCA:
            fused = self._fuse_pca(matrices)
        else:
            fused = self._fuse_concat(matrices)

        if self.config.normalize:
            norms = np.linalg.norm(fused, axis=1, keepdims=True)
            fused = fused / np.maximum(norms, 1e-10)

        return NodeEmbeddings(
            node_ids=list(node_ids),
            embeddings=fused,
            source=f"fused({','.join(matrices.keys())})",
        )

    def _fuse_concat(self, matrices: Dict[str, np.ndarray]) -> np.ndarray:
        return np.concatenate(list(matrices.values()), axis=1)

    def _fuse_weighted_sum(self, matrices: Dict[str, np.ndarray]) -> np.ndarray:
        weights = self.config.weights or {}
        # Pad to same dimension
        dims = [m.shape[1] for m in matrices.values()]
        max_dim = max(dims)
        padded = []
        for name, m in matrices.items():
            w = weights.get(name, 1.0)
            if m.shape[1] < max_dim:
                m = np.pad(m, ((0, 0), (0, max_dim - m.shape[1])))
            padded.append(m * w)
        return sum(padded) / sum(weights.get(n, 1.0) for n in matrices)

    def _fuse_average(self, matrices: Dict[str, np.ndarray]) -> np.ndarray:
        dims = [m.shape[1] for m in matrices.values()]
        max_dim = max(dims)
        padded = []
        for m in matrices.values():
            if m.shape[1] < max_dim:
                m = np.pad(m, ((0, 0), (0, max_dim - m.shape[1])))
            padded.append(m)
        return np.mean(padded, axis=0)

    def _fuse_pca(self, matrices: Dict[str, np.ndarray]) -> np.ndarray:
        concatenated = np.concatenate(list(matrices.values()), axis=1)
        n_components = self.config.output_dim or min(64, concatenated.shape[1])
        n_components = min(n_components, concatenated.shape[0], concatenated.shape[1])
        pca = PCA(n_components=n_components)
        return pca.fit_transform(concatenated)
