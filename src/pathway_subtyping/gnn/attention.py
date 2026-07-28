"""
Biological Attention Mechanisms for GNN Interpretability (#23).

Multi-head attention with biological prior injection, edge-type-specific
attention, and pathway co-attention for gene-pathway relationships.

Requires PyTorch. Ported from autism-pathway-framework Module 06.

Research use only. Not for clinical decision-making.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

if HAS_TORCH:

    class BiologicalAttention(nn.Module):
        """Multi-head attention with biological prior bias injection.

        Incorporates pLI constraint scores, expression levels, and
        SFARI scores as attention biases to produce biologically-informed
        attention weights.
        """

        def __init__(
            self,
            hidden_dim: int,
            num_heads: int = 8,
            dropout: float = 0.1,
            prior_types: Optional[List[str]] = None,
        ):
            super().__init__()
            self.hidden_dim = hidden_dim
            self.num_heads = num_heads
            self.head_dim = hidden_dim // num_heads

            self.q_proj = nn.Linear(hidden_dim, hidden_dim)
            self.k_proj = nn.Linear(hidden_dim, hidden_dim)
            self.v_proj = nn.Linear(hidden_dim, hidden_dim)
            self.out_proj = nn.Linear(hidden_dim, hidden_dim)
            self.dropout = nn.Dropout(dropout)

            # Biological prior bias projections
            self.prior_types = prior_types or []
            if self.prior_types:
                self.prior_bias = nn.ModuleDict(
                    {pt: nn.Linear(1, num_heads) for pt in self.prior_types}
                )

        def forward(
            self,
            query: torch.Tensor,
            key: torch.Tensor,
            value: torch.Tensor,
            bio_priors: Optional[Dict[str, torch.Tensor]] = None,
            attention_mask: Optional[torch.Tensor] = None,
        ) -> Tuple[torch.Tensor, torch.Tensor]:
            """Forward pass with optional biological prior injection.

            Args:
                query: (n, hidden_dim)
                key: (m, hidden_dim)
                value: (m, hidden_dim)
                bio_priors: Dict of prior_name -> (m, 1) prior values
                attention_mask: Optional mask

            Returns:
                (output, attention_weights)
            """
            n = query.size(0)
            m = key.size(0)

            q = self.q_proj(query).view(n, self.num_heads, self.head_dim)
            k = self.k_proj(key).view(m, self.num_heads, self.head_dim)
            v = self.v_proj(value).view(m, self.num_heads, self.head_dim)

            # Scaled dot-product attention
            scale = self.head_dim**0.5
            attn = torch.einsum("nhd,mhd->hnm", q, k) / scale

            # Inject biological priors as attention bias
            if bio_priors and self.prior_types:
                for pt in self.prior_types:
                    if pt in bio_priors and pt in self.prior_bias:
                        prior = bio_priors[pt]
                        if prior.dim() == 1:
                            prior = prior.unsqueeze(-1)
                        bias = self.prior_bias[pt](prior)  # (m, num_heads)
                        attn = attn + bias.permute(1, 0).unsqueeze(1)

            if attention_mask is not None:
                attn = attn.masked_fill(~attention_mask, float("-inf"))

            weights = F.softmax(attn, dim=-1)
            weights = self.dropout(weights)

            out = torch.einsum("hnm,mhd->nhd", weights, v)
            out = out.reshape(n, self.hidden_dim)
            out = self.out_proj(out)

            # Average attention weights across heads for interpretability
            avg_weights = weights.mean(dim=0)

            return out, avg_weights

    class PathwayCoAttention(nn.Module):
        """Bi-directional attention between gene and pathway features.

        Reveals which gene-pathway relationships are most important.
        """

        def __init__(
            self,
            gene_dim: int,
            pathway_dim: int,
            num_heads: int = 4,
            dropout: float = 0.1,
        ):
            super().__init__()
            self.gene_to_pathway = nn.MultiheadAttention(
                gene_dim,
                num_heads,
                kdim=pathway_dim,
                vdim=pathway_dim,
                batch_first=True,
                dropout=dropout,
            )
            self.pathway_to_gene = nn.MultiheadAttention(
                pathway_dim,
                num_heads,
                kdim=gene_dim,
                vdim=gene_dim,
                batch_first=True,
                dropout=dropout,
            )

        def forward(
            self,
            gene_features: torch.Tensor,
            pathway_features: torch.Tensor,
        ) -> Tuple[torch.Tensor, torch.Tensor]:
            """Bi-directional cross-attention.

            Args:
                gene_features: (n_genes, gene_dim)
                pathway_features: (n_pathways, pathway_dim)

            Returns:
                (gene_attended, pathway_attended)
            """
            g = gene_features.unsqueeze(0)
            p = pathway_features.unsqueeze(0)

            gene_attended, _ = self.gene_to_pathway(g, p, p)
            pathway_attended, _ = self.pathway_to_gene(p, g, g)

            return gene_attended.squeeze(0), pathway_attended.squeeze(0)

else:

    class BiologicalAttention:
        pass

    class PathwayCoAttention:
        pass
