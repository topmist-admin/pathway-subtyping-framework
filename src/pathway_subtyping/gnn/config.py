"""
GNN Configuration.

Research use only. Not for clinical decision-making.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class GNNConfig:
    """Configuration for the OntologyAwareGNN model."""

    input_dim: int = 256
    hidden_dim: int = 256
    output_dim: int = 128
    num_layers: int = 3
    num_heads: int = 8
    dropout: float = 0.1
    use_residual: bool = True
    use_layer_norm: bool = True
    aggregation: str = "mean"
    activation: str = "relu"
    edge_types: List[str] = field(default_factory=lambda: [
        "gene_interacts_gene",
        "gene_regulates_gene",
        "gene_in_pathway",
        "gene_has_go",
    ])
    node_types: List[str] = field(default_factory=lambda: [
        "gene", "pathway", "go_term",
    ])
    prior_types: List[str] = field(default_factory=lambda: [
        "pli", "expression", "sfari_score",
    ])
    task_heads: List[str] = field(default_factory=lambda: [
        "gene_classification",
    ])

    def to_dict(self) -> dict:
        return {
            "input_dim": self.input_dim,
            "hidden_dim": self.hidden_dim,
            "output_dim": self.output_dim,
            "num_layers": self.num_layers,
            "num_heads": self.num_heads,
            "dropout": self.dropout,
            "edge_types": self.edge_types,
            "node_types": self.node_types,
        }
