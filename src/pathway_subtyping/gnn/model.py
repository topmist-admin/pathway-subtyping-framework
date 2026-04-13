"""
OntologyAwareGNN: Heterogeneous GNN with biological attention (#21).

Multi-task GNN supporting gene risk classification and link prediction,
with per-edge-type message passing and biological prior weighting.

Requires PyTorch. Ported from autism-pathway-framework Module 06.

EXPERIMENTAL/BETA — all GNN features are experimental.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


@dataclass
class GNNOutput:
    """Output from OntologyAwareGNN forward pass."""

    node_embeddings: Dict[str, Any] = field(default_factory=dict)
    gene_logits: Optional[Any] = None
    link_scores: Optional[Any] = None
    attention_weights: Optional[Dict[str, Any]] = None
    loss: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "has_gene_logits": self.gene_logits is not None,
            "has_link_scores": self.link_scores is not None,
            "loss": self.loss,
        }


if HAS_TORCH:
    from .attention import BiologicalAttention
    from .config import GNNConfig
    from .layers import (
        BioPriorWeighting,
        HierarchicalAggregator,
        MessagePassingLayer,
    )

    class OntologyAwareGNN(nn.Module):
        """Heterogeneous GNN with biological attention and multi-task heads.

        Architecture::

            Input features (per node type)
            → InputProjection (heterogeneous → common dim)
            → N × MessagePassingLayer (edge-type aware)
            → HierarchicalAggregator (GO/pathway hierarchy)
            → BiologicalAttention (prior-informed)
            → BioPriorWeighting
            → OutputProjection
            → TaskHeads (classification / link prediction)

        Usage::

            config = GNNConfig(input_dim=256, hidden_dim=256, output_dim=128)
            model = OntologyAwareGNN(config)
            output = model(node_features, edge_index, edge_type)
        """

        def __init__(
            self,
            config: GNNConfig,
            input_dims: Optional[Dict[str, int]] = None,
        ):
            super().__init__()
            self.config = config

            # Input projections: one per node type
            in_dims = input_dims or {nt: config.input_dim for nt in config.node_types}
            self.input_projections = nn.ModuleDict({
                nt: nn.Linear(dim, config.hidden_dim)
                for nt, dim in in_dims.items()
            })

            # Message passing layers
            self.mp_layers = nn.ModuleList([
                MessagePassingLayer(
                    config.hidden_dim,
                    config.hidden_dim,
                    config.edge_types,
                    aggregation=config.aggregation,
                    dropout=config.dropout,
                    residual=config.use_residual,
                )
                for _ in range(config.num_layers)
            ])

            # Hierarchical aggregator
            self.hierarchy_agg = HierarchicalAggregator(
                config.hidden_dim, num_heads=config.num_heads // 2 or 1
            )

            # Biological attention
            self.bio_attention = BiologicalAttention(
                config.hidden_dim,
                num_heads=config.num_heads,
                dropout=config.dropout,
                prior_types=config.prior_types,
            )

            # Bio prior weighting
            self.bio_weighting = BioPriorWeighting(
                config.hidden_dim,
                config.prior_types,
            )

            # Output projection
            self.output_proj = nn.Sequential(
                nn.Linear(config.hidden_dim, config.output_dim),
                nn.ReLU(),
                nn.Dropout(config.dropout),
            )

            # Task heads
            self.task_heads = nn.ModuleDict()
            if "gene_classification" in config.task_heads:
                self.task_heads["gene_classification"] = nn.Sequential(
                    nn.Linear(config.output_dim, config.output_dim // 2),
                    nn.ReLU(),
                    nn.Dropout(config.dropout),
                    nn.Linear(config.output_dim // 2, 4),  # 4 SFARI risk classes
                )
            if "link_prediction" in config.task_heads:
                self.task_heads["link_prediction"] = nn.Bilinear(
                    config.output_dim, config.output_dim, 1
                )

        def forward(
            self,
            node_features: torch.Tensor,
            edge_index: torch.Tensor,
            edge_type: torch.Tensor,
            edge_type_names: Optional[List[str]] = None,
            node_type_ids: Optional[Dict[str, torch.Tensor]] = None,
            bio_priors: Optional[Dict[str, torch.Tensor]] = None,
            hierarchy_edges: Optional[torch.Tensor] = None,
            labels: Optional[torch.Tensor] = None,
        ) -> GNNOutput:
            """Forward pass.

            Args:
                node_features: (n_nodes, input_dim) or pre-projected features.
                edge_index: (2, n_edges) source/target indices.
                edge_type: (n_edges,) integer edge type indices.
                edge_type_names: Ordered edge type names matching indices.
                node_type_ids: Dict mapping node_type -> indices tensor.
                bio_priors: Dict of prior_name -> (n_nodes,) values.
                hierarchy_edges: (2, n_hier_edges) for hierarchical aggregation.
                labels: Optional labels for loss computation.

            Returns:
                GNNOutput with embeddings, logits, and optional loss.
            """
            et_names = edge_type_names or self.config.edge_types
            h = node_features

            # Project if needed (assumes pre-projected if no node_type_ids)
            if node_type_ids and h.dim() == 2:
                projected = torch.zeros(
                    h.size(0), self.config.hidden_dim, device=h.device
                )
                for nt, indices in node_type_ids.items():
                    if nt in self.input_projections:
                        projected[indices] = self.input_projections[nt](h[indices])
                h = projected
            elif h.size(-1) != self.config.hidden_dim:
                # Simple projection for uniform input
                default_proj = list(self.input_projections.values())[0]
                h = default_proj(h)

            # Message passing
            for mp_layer in self.mp_layers:
                h = mp_layer(h, edge_index, edge_type, et_names)

            # Hierarchical aggregation
            h = self.hierarchy_agg(h, hierarchy_edges)

            # Biological attention (self-attention)
            h_attn, attn_weights = self.bio_attention(h, h, h, bio_priors)
            h = h + h_attn  # Residual

            # Bio prior weighting
            h = self.bio_weighting(h, bio_priors)

            # Output projection
            embeddings = self.output_proj(h)

            # Task heads
            gene_logits = None
            loss = None

            if "gene_classification" in self.task_heads:
                gene_logits = self.task_heads["gene_classification"](embeddings)
                if labels is not None:
                    loss = F.cross_entropy(gene_logits, labels)

            return GNNOutput(
                node_embeddings={"all": embeddings},
                gene_logits=gene_logits,
                attention_weights={"biological": attn_weights},
                loss=float(loss.item()) if loss is not None else None,
            )

    class GNNTrainer:
        """Training loop for OntologyAwareGNN."""

        def __init__(
            self,
            model: OntologyAwareGNN,
            learning_rate: float = 1e-3,
            weight_decay: float = 1e-5,
            device: str = "cpu",
        ):
            self.model = model.to(device)
            self.device = device
            self.optimizer = torch.optim.AdamW(
                model.parameters(), lr=learning_rate, weight_decay=weight_decay
            )

        def train_step(
            self,
            node_features: torch.Tensor,
            edge_index: torch.Tensor,
            edge_type: torch.Tensor,
            labels: torch.Tensor,
            edge_type_names: Optional[List[str]] = None,
            bio_priors: Optional[Dict[str, torch.Tensor]] = None,
        ) -> Dict[str, float]:
            """Single training step.

            Returns:
                Dict with 'loss' and optional metrics.
            """
            self.model.train()
            self.optimizer.zero_grad()

            output = self.model(
                node_features.to(self.device),
                edge_index.to(self.device),
                edge_type.to(self.device),
                edge_type_names=edge_type_names,
                bio_priors=bio_priors,
                labels=labels.to(self.device),
            )

            if output.loss is not None:
                loss_tensor = F.cross_entropy(
                    output.gene_logits, labels.to(self.device)
                )
                loss_tensor.backward()
                self.optimizer.step()
                return {"loss": float(loss_tensor.item())}

            return {"loss": 0.0}

        @torch.no_grad()
        def evaluate(
            self,
            node_features: torch.Tensor,
            edge_index: torch.Tensor,
            edge_type: torch.Tensor,
            labels: torch.Tensor,
            edge_type_names: Optional[List[str]] = None,
        ) -> Dict[str, float]:
            """Evaluate model on validation data."""
            self.model.eval()
            output = self.model(
                node_features.to(self.device),
                edge_index.to(self.device),
                edge_type.to(self.device),
                edge_type_names=edge_type_names,
                labels=labels.to(self.device),
            )

            metrics = {"loss": output.loss or 0.0}
            if output.gene_logits is not None:
                preds = output.gene_logits.argmax(dim=-1)
                acc = float((preds == labels.to(self.device)).float().mean())
                metrics["accuracy"] = acc

            return metrics

else:

    class OntologyAwareGNN:
        """Stub when PyTorch is not available."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise ImportError(
                "OntologyAwareGNN requires PyTorch. "
                "Install with: pip install pathway-subtyping[gnn]"
            )

    class GNNTrainer:
        """Stub when PyTorch is not available."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise ImportError(
                "GNNTrainer requires PyTorch. "
                "Install with: pip install pathway-subtyping[gnn]"
            )
