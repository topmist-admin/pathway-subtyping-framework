# GNN & Graph Embeddings API

> **Module**: `pathway_subtyping.gnn`
> **Install**: `pip install pathway-subtyping[gnn]`
> **Status**: EXPERIMENTAL/BETA

Graph neural network models and knowledge graph embeddings for gene risk classification, patient-level subtyping, and interpretable subnetwork discovery.

**Key reference:** Pelletier et al. (Briefings in Bioinformatics, 2024) — GNNs rarely improve over simpler methods on transcriptomic data when biological networks are imperfect. All GNN features benchmark against non-GNN baselines.

---

## Quick Example

```python
# KG embeddings (pure numpy — no torch required)
from pathway_subtyping.gnn import TransEModel, EmbeddingTrainer, TrainingConfig

triples = [("RAS", "activates", "RAF"), ("RAF", "activates", "MEK"), ...]
trainer = EmbeddingTrainer(TrainingConfig(model_type="transe", embedding_dim=64, epochs=100))
model, history = trainer.train(triples)
embeddings = model.get_node_embeddings()
print(embeddings.most_similar("RAS", k=5))
```

---

## KG Embeddings (numpy — always available)

### `TransEModel`

Translational embedding: `h + r ≈ t`. Bordes et al. (NeurIPS 2013).

```python
model = TransEModel(
    embedding_dim: int = 128,
    margin: float = 1.0,
    norm: int = 1,           # 1 = L1, 2 = L2
    learning_rate: float = 0.01,
    seed: Optional[int] = None,
)

history = model.train(triples, epochs=100, batch_size=256, n_negative=5)
embeddings = model.get_node_embeddings()  -> NodeEmbeddings
score = model.predict_link("RAS", "activates", "RAF")  -> float
```

### `RotatEModel`

Relational rotation in complex space. Sun et al. (ICLR 2019).

```python
model = RotatEModel(
    embedding_dim: int = 128,
    margin: float = 6.0,
    learning_rate: float = 0.001,
    seed: Optional[int] = None,
)
```

Same API as TransE: `train()`, `get_node_embeddings()`, `predict_link()`.

### `EmbeddingTrainer`

Convenience wrapper for model creation and training.

```python
config = TrainingConfig(
    model_type="transe",  # or "rotate"
    embedding_dim=128, epochs=100, batch_size=256, n_negative=5, seed=42,
)
trainer = EmbeddingTrainer(config)
model, history = trainer.train(triples)
```

### `NodeEmbeddings`

Container for learned node embedding vectors.

```python
embeddings.n_nodes    # Number of embedded nodes
embeddings.dim        # Embedding dimensionality
embeddings.get("RAS") -> np.ndarray          # Single node vector
embeddings.get_batch(["RAS", "RAF"]) -> np.ndarray  # Multiple nodes
embeddings.most_similar("RAS", k=10) -> List[Tuple[str, float]]
```

---

## Embedding Fusion (numpy — always available)

### `EmbeddingFusion`

Combines embeddings from multiple sources (KG embeddings, expression features, foundation models) into unified node feature vectors.

```python
from pathway_subtyping.gnn.embeddings import EmbeddingFusion, FusionConfig, FusionMethod

fusion = EmbeddingFusion(FusionConfig(
    method=FusionMethod.CONCAT,  # or WEIGHTED_SUM, AVERAGE, PCA
    output_dim=64,               # For PCA only
    normalize=True,
))

fused = fusion.fuse({
    "kg": kg_embeddings,         # NodeEmbeddings from TransE/RotatE
    "expression": expr_embeddings,  # NodeEmbeddings from expression data
})
```

| Method | Description |
|--------|-------------|
| `CONCAT` | Column-bind all sources (preserves all info) |
| `WEIGHTED_SUM` | Weighted combination with per-source weights |
| `AVERAGE` | Simple mean across sources |
| `PCA` | Concatenate then PCA to `output_dim` |

---

## GNN Model (requires PyTorch)

### `GNNConfig`

```python
from pathway_subtyping.gnn import GNNConfig

config = GNNConfig(
    input_dim=256,       # Node feature dimension
    hidden_dim=256,      # Hidden layer dimension
    output_dim=128,      # Output embedding dimension
    num_layers=3,        # Message passing layers
    num_heads=8,         # Attention heads
    dropout=0.1,
    edge_types=["gene_interacts_gene", "gene_regulates_gene"],
    prior_types=["pli", "expression", "sfari_score"],
    task_heads=["gene_classification"],
)
```

### `OntologyAwareGNN`

Heterogeneous GNN with edge-type-aware message passing, hierarchical aggregation, and biological attention.

```python
from pathway_subtyping.gnn import OntologyAwareGNN, GNNConfig
import torch

config = GNNConfig(input_dim=64, hidden_dim=64, output_dim=32, num_layers=2)
model = OntologyAwareGNN(config)

output = model(
    node_features: torch.Tensor,            # (n_nodes, input_dim)
    edge_index: torch.Tensor,               # (2, n_edges)
    edge_type: torch.Tensor,                # (n_edges,) integer indices
    edge_type_names: List[str],             # Ordered type names
    bio_priors: Dict[str, torch.Tensor],    # Optional: pLI, expression, etc.
    labels: torch.Tensor,                   # Optional: for loss computation
) -> GNNOutput

output.node_embeddings["all"]    # (n_nodes, output_dim) tensor
output.gene_logits               # (n_nodes, n_classes) if gene_classification head
output.attention_weights         # Interpretable attention maps
output.loss                      # Cross-entropy loss if labels provided
```

**Architecture:**
```
Input → InputProjection (per node type) → N × MessagePassingLayer
→ HierarchicalAggregator → BiologicalAttention → BioPriorWeighting
→ OutputProjection → TaskHeads
```

### `GNNTrainer`

```python
from pathway_subtyping.gnn.model import GNNTrainer

trainer = GNNTrainer(model, learning_rate=1e-3, device="cpu")
metrics = trainer.train_step(node_features, edge_index, edge_type, labels)
eval_metrics = trainer.evaluate(node_features, edge_index, edge_type, labels)
print(f"Loss: {eval_metrics['loss']:.4f}, Accuracy: {eval_metrics['accuracy']:.2%}")
```

---

## Attention Modules (requires PyTorch)

### `BiologicalAttention`

Multi-head attention with biological prior bias injection. Priors (pLI, expression, SFARI scores) are projected into attention bias terms before softmax.

```python
from pathway_subtyping.gnn.attention import BiologicalAttention

attn = BiologicalAttention(
    hidden_dim=128, num_heads=8,
    prior_types=["pli", "expression"],
)
output, weights = attn(query, key, value, bio_priors={"pli": pli_scores})
# weights: (n_query, n_key) — interpretable attention map
```

### `PathwayCoAttention`

Bi-directional cross-attention between gene and pathway features.

```python
from pathway_subtyping.gnn.attention import PathwayCoAttention

coattn = PathwayCoAttention(gene_dim=128, pathway_dim=128, num_heads=4)
gene_attended, pathway_attended = coattn(gene_features, pathway_features)
```

---

## Layer Components (requires PyTorch)

| Class | Description |
|-------|-------------|
| `EdgeTypeTransform` | Per-edge-type linear projections |
| `MessagePassingLayer` | Edge-type-aware aggregation with residuals + LayerNorm |
| `HierarchicalAggregator` | Multi-head self-attention for GO/pathway hierarchy |
| `BioPriorWeighting` | Learned combination of biological prior signals |

---

## Graceful Degradation

All torch-dependent classes check for PyTorch at import time:

```python
# Without torch installed:
from pathway_subtyping.gnn import TransEModel      # Works (numpy)
from pathway_subtyping.gnn import OntologyAwareGNN  # ImportError with helpful message
```

---

EXPERIMENTAL/BETA — benchmark against non-GNN baselines before drawing conclusions. Research use only.
