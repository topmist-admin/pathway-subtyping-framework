"""
GNN and Graph Embedding Layer for the Pathway Subtyping Framework.

Provides graph neural network models and knowledge graph embeddings
for gene risk classification, patient-level subtyping, and
interpretable subnetwork discovery.

Submodules:
    - kg_embeddings: TransE/RotatE knowledge graph embeddings (numpy, no torch)
    - config: GNN model configuration
    - layers: Message passing and aggregation layers (requires torch)
    - attention: Biological attention with prior injection (requires torch)
    - model: OntologyAwareGNN and GNNEmbedder (requires torch)
    - embeddings: Pretrained embedding extractors and fusion

Install: pip install pathway-subtyping[gnn]

All GNN features are EXPERIMENTAL/BETA.

Research use only. Not for clinical decision-making.
"""

# KG embeddings are pure numpy — always available
from .kg_embeddings import (
    BaseEmbeddingModel,
    EmbeddingTrainer,
    NodeEmbeddings,
    RotatEModel,
    TrainingConfig,
    TrainingHistory,
    TransEModel,
)

__all__ = [
    # KG Embeddings (numpy, always available)
    "BaseEmbeddingModel",
    "TransEModel",
    "RotatEModel",
    "NodeEmbeddings",
    "TrainingHistory",
    "TrainingConfig",
    "EmbeddingTrainer",
]

# Torch-dependent modules — optional import
try:
    from .config import GNNConfig
    from .model import GNNOutput, OntologyAwareGNN

    __all__.extend([
        "GNNConfig",
        "OntologyAwareGNN",
        "GNNOutput",
    ])
except ImportError:
    pass
