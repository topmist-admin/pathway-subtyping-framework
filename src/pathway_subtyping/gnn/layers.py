"""
GNN Layers: EdgeTypeTransform, MessagePassingLayer, HierarchicalAggregator,
BioPriorWeighting.

Requires PyTorch. Ported from autism-pathway-framework Module 06.

Research use only. Not for clinical decision-making.
"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

if HAS_TORCH:

    class EdgeTypeTransform(nn.Module):
        """Per-edge-type linear projections."""

        def __init__(self, in_dim: int, out_dim: int, edge_types: List[str]):
            super().__init__()
            self.transforms = nn.ModuleDict({et: nn.Linear(in_dim, out_dim) for et in edge_types})
            self.edge_types = edge_types

        def forward(self, x: torch.Tensor, edge_type: str) -> torch.Tensor:
            return self.transforms[edge_type](x)

        def forward_all(
            self,
            x: torch.Tensor,
            edge_type_indices: torch.Tensor,
            edge_type_names: List[str],
        ) -> torch.Tensor:
            out = torch.zeros_like(x)
            for i, et_name in enumerate(edge_type_names):
                if et_name in self.transforms:
                    mask = edge_type_indices == i
                    if mask.any():
                        out[mask] = self.transforms[et_name](x[mask])
            return out

    class MessagePassingLayer(nn.Module):
        """Edge-type-aware message passing with residual connections."""

        def __init__(
            self,
            in_dim: int,
            out_dim: int,
            edge_types: List[str],
            aggregation: str = "mean",
            dropout: float = 0.1,
            residual: bool = True,
        ):
            super().__init__()
            self.edge_transform = EdgeTypeTransform(in_dim, out_dim, edge_types)
            self.update_mlp = nn.Sequential(
                nn.Linear(out_dim, out_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
            )
            self.layer_norm = nn.LayerNorm(out_dim)
            self.aggregation = aggregation
            self.residual = residual and (in_dim == out_dim)

        def forward(
            self,
            x: torch.Tensor,
            edge_index: torch.Tensor,
            edge_type: torch.Tensor,
            edge_type_names: List[str],
            edge_weight: Optional[torch.Tensor] = None,
        ) -> torch.Tensor:
            n_nodes = x.size(0)
            src, dst = edge_index[0], edge_index[1]

            # Transform source node features by edge type
            src_features = x[src]
            transformed = self.edge_transform.forward_all(src_features, edge_type, edge_type_names)

            if edge_weight is not None:
                transformed = transformed * edge_weight.unsqueeze(-1)

            # Aggregate messages
            out = torch.zeros(n_nodes, transformed.size(-1), device=x.device)
            if self.aggregation == "sum":
                out.scatter_add_(0, dst.unsqueeze(-1).expand_as(transformed), transformed)
            else:  # mean
                counts = torch.zeros(n_nodes, 1, device=x.device)
                out.scatter_add_(0, dst.unsqueeze(-1).expand_as(transformed), transformed)
                counts.scatter_add_(
                    0, dst.unsqueeze(-1), torch.ones_like(dst.unsqueeze(-1).float())
                )
                out = out / counts.clamp(min=1)

            out = self.update_mlp(out)
            out = self.layer_norm(out)

            if self.residual:
                out = out + x

            return out

    class HierarchicalAggregator(nn.Module):
        """Propagates hierarchical information (GO/pathway hierarchy)."""

        def __init__(self, hidden_dim: int, num_heads: int = 4):
            super().__init__()
            self.attention = nn.MultiheadAttention(hidden_dim, num_heads, batch_first=True)
            self.norm = nn.LayerNorm(hidden_dim)

        def forward(
            self,
            x: torch.Tensor,
            hierarchy_edges: Optional[torch.Tensor] = None,
        ) -> torch.Tensor:
            if hierarchy_edges is None or hierarchy_edges.numel() == 0:
                return x
            # Self-attention over node features
            x_unsqueezed = x.unsqueeze(0)
            attn_out, _ = self.attention(x_unsqueezed, x_unsqueezed, x_unsqueezed)
            return self.norm(x + attn_out.squeeze(0))

    class BioPriorWeighting(nn.Module):
        """Learned weighting of biological prior information."""

        def __init__(
            self,
            hidden_dim: int,
            prior_types: List[str],
            combination: str = "multiplicative",
        ):
            super().__init__()
            self.prior_transforms = nn.ModuleDict(
                {
                    pt: nn.Sequential(
                        nn.Linear(1, hidden_dim),
                        nn.Sigmoid(),
                    )
                    for pt in prior_types
                }
            )
            self.combination = combination
            if combination == "learned":
                self.combiner = nn.Linear(len(prior_types) * hidden_dim, hidden_dim)

        def forward(
            self,
            x: torch.Tensor,
            bio_priors: Optional[Dict[str, torch.Tensor]] = None,
        ) -> torch.Tensor:
            if bio_priors is None:
                return x

            weights = []
            for pt, transform in self.prior_transforms.items():
                if pt in bio_priors:
                    prior = bio_priors[pt]
                    if prior.dim() == 1:
                        prior = prior.unsqueeze(-1)
                    w = transform(prior)
                    weights.append(w)

            if not weights:
                return x

            if self.combination == "multiplicative":
                combined = weights[0]
                for w in weights[1:]:
                    combined = combined * w
                return x * combined
            elif self.combination == "learned":
                cat = torch.cat(weights, dim=-1)
                combined = torch.sigmoid(self.combiner(cat))
                return x * combined
            else:  # additive
                combined = sum(weights) / len(weights)
                return x + combined

else:
    # Stubs when torch is not available
    class EdgeTypeTransform:
        pass

    class MessagePassingLayer:
        pass

    class HierarchicalAggregator:
        pass

    class BioPriorWeighting:
        pass
