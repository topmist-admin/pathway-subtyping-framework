"""
Pretrained embedding extractors and fusion for GNN node features.

Supports multiple embedding sources:
    - KG embeddings (TransE/RotatE from kg_embeddings.py)
    - Expression-based features (from existing expression.py)
    - Fusion: concatenation, weighted sum, PCA

Foundation model extractors (Geneformer, ESM-2, PubMedBERT) require
the transformers library and are optional.
"""

from .fusion import EmbeddingFusion, FusionConfig, FusionMethod

__all__ = [
    "EmbeddingFusion",
    "FusionConfig",
    "FusionMethod",
]
