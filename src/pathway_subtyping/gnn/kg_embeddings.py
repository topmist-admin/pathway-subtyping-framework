"""
Knowledge Graph Embeddings: TransE and RotatE.

Pure numpy implementation — no PyTorch required. Learns low-dimensional
vector representations of KG entities and relations for link prediction
and gene similarity computation.

Ported from autism-pathway-framework Module 04.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class NodeEmbeddings:
    """Embedding vectors for a set of nodes."""

    node_ids: List[str]
    embeddings: np.ndarray  # (n_nodes, embedding_dim)
    source: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def dim(self) -> int:
        return self.embeddings.shape[1] if len(self.embeddings.shape) > 1 else 0

    @property
    def n_nodes(self) -> int:
        return len(self.node_ids)

    def get(self, node_id: str) -> Optional[np.ndarray]:
        try:
            idx = self.node_ids.index(node_id)
            return self.embeddings[idx]
        except ValueError:
            return None

    def get_batch(self, node_ids: List[str]) -> np.ndarray:
        indices = [self.node_ids.index(nid) for nid in node_ids if nid in self.node_ids]
        return self.embeddings[indices]

    def most_similar(self, node_id: str, k: int = 10) -> List[Tuple[str, float]]:
        emb = self.get(node_id)
        if emb is None:
            return []
        sims = self.embeddings @ emb / (
            np.linalg.norm(self.embeddings, axis=1) * np.linalg.norm(emb) + 1e-10
        )
        top_idx = np.argsort(-sims)[:k + 1]
        return [
            (self.node_ids[i], float(sims[i]))
            for i in top_idx
            if self.node_ids[i] != node_id
        ][:k]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_nodes": self.n_nodes,
            "dim": self.dim,
            "source": self.source,
        }


@dataclass
class TrainingHistory:
    """Training history for embedding models."""

    epochs: List[int] = field(default_factory=list)
    losses: List[float] = field(default_factory=list)

    def add(self, epoch: int, loss: float) -> None:
        self.epochs.append(epoch)
        self.losses.append(loss)

    @property
    def best_loss(self) -> float:
        return min(self.losses) if self.losses else float("inf")


@dataclass
class TrainingConfig:
    """Configuration for KG embedding training."""

    model_type: str = "transe"  # "transe" or "rotate"
    embedding_dim: int = 128
    learning_rate: float = 0.01
    margin: float = 1.0
    norm: int = 1
    epochs: int = 100
    batch_size: int = 256
    n_negative: int = 5
    seed: Optional[int] = None


class BaseEmbeddingModel:
    """Base class for KG embedding models.

    Subclasses implement _score_triple, _compute_loss, _update_embeddings.
    """

    def __init__(
        self,
        embedding_dim: int = 128,
        learning_rate: float = 0.01,
        seed: Optional[int] = None,
    ):
        self.embedding_dim = embedding_dim
        self.learning_rate = learning_rate
        self._rng = np.random.default_rng(seed)
        self._entity_embeddings: Optional[np.ndarray] = None
        self._relation_embeddings: Optional[np.ndarray] = None
        self._entity_to_idx: Dict[str, int] = {}
        self._relation_to_idx: Dict[str, int] = {}
        self._idx_to_entity: Dict[int, str] = {}

    def train(
        self,
        triples: List[Tuple[str, str, str]],
        epochs: int = 100,
        batch_size: int = 256,
        n_negative: int = 5,
    ) -> TrainingHistory:
        """Train embeddings on (head, relation, tail) triples.

        Args:
            triples: List of (head_id, relation_type, tail_id).
            epochs: Training epochs.
            batch_size: Batch size.
            n_negative: Negative samples per positive.

        Returns:
            TrainingHistory with loss per epoch.
        """
        # Build indices
        entities = set()
        relations = set()
        for h, r, t in triples:
            entities.update([h, t])
            relations.add(r)

        self._entity_to_idx = {e: i for i, e in enumerate(sorted(entities))}
        self._idx_to_entity = {i: e for e, i in self._entity_to_idx.items()}
        self._relation_to_idx = {r: i for i, r in enumerate(sorted(relations))}

        n_entities = len(self._entity_to_idx)
        n_relations = len(self._relation_to_idx)

        # Initialize embeddings
        self._entity_embeddings = self._rng.normal(
            0, 0.1, (n_entities, self.embedding_dim)
        )
        self._relation_embeddings = self._rng.normal(
            0, 0.1, (n_relations, self.embedding_dim)
        )
        # Normalize entities
        norms = np.linalg.norm(self._entity_embeddings, axis=1, keepdims=True)
        self._entity_embeddings /= np.maximum(norms, 1e-10)

        # Convert triples to indices
        idx_triples = np.array([
            [self._entity_to_idx[h], self._relation_to_idx[r], self._entity_to_idx[t]]
            for h, r, t in triples
        ])

        history = TrainingHistory()

        for epoch in range(epochs):
            self._rng.shuffle(idx_triples)
            epoch_loss = 0.0
            n_batches = 0

            for start in range(0, len(idx_triples), batch_size):
                batch = idx_triples[start:start + batch_size]
                # Generate negative samples
                neg_batch = batch.copy()
                for i in range(len(neg_batch)):
                    for _ in range(n_negative):
                        if self._rng.random() < 0.5:
                            neg_batch[i, 0] = self._rng.integers(0, n_entities)
                        else:
                            neg_batch[i, 2] = self._rng.integers(0, n_entities)

                loss = self._compute_loss(batch, neg_batch)
                self._update_embeddings(batch, neg_batch)
                epoch_loss += loss
                n_batches += 1

            avg_loss = epoch_loss / max(n_batches, 1)
            history.add(epoch, avg_loss)

            if epoch % 20 == 0:
                logger.info(
                    "[GNN KG Embeddings] Epoch %d/%d, loss=%.4f", epoch, epochs, avg_loss
                )

        return history

    def get_node_embeddings(self) -> NodeEmbeddings:
        """Return trained entity embeddings."""
        if self._entity_embeddings is None:
            raise ValueError("Model not trained. Call train() first.")
        return NodeEmbeddings(
            node_ids=list(self._entity_to_idx.keys()),
            embeddings=self._entity_embeddings.copy(),
            source=self.__class__.__name__,
        )

    def predict_link(self, head: str, relation: str, tail: str) -> Optional[float]:
        """Score a triple."""
        if head not in self._entity_to_idx or tail not in self._entity_to_idx:
            return None
        if relation not in self._relation_to_idx:
            return None
        h_idx = self._entity_to_idx[head]
        r_idx = self._relation_to_idx[relation]
        t_idx = self._entity_to_idx[tail]
        return float(self._score_triple(h_idx, r_idx, t_idx))

    def _score_triple(self, h: int, r: int, t: int) -> float:
        raise NotImplementedError

    def _compute_loss(self, pos: np.ndarray, neg: np.ndarray) -> float:
        raise NotImplementedError

    def _update_embeddings(self, pos: np.ndarray, neg: np.ndarray) -> None:
        raise NotImplementedError


class TransEModel(BaseEmbeddingModel):
    """TransE: h + r ≈ t in embedding space.

    Bordes et al. "Translating Embeddings for Modeling Multi-relational Data"
    (NeurIPS 2013).
    """

    def __init__(
        self,
        embedding_dim: int = 128,
        margin: float = 1.0,
        norm: int = 1,
        learning_rate: float = 0.01,
        seed: Optional[int] = None,
    ):
        super().__init__(embedding_dim, learning_rate, seed)
        self.margin = margin
        self.norm = norm

    def _distance(self, h: np.ndarray, r: np.ndarray, t: np.ndarray) -> np.ndarray:
        if self.norm == 1:
            return np.sum(np.abs(h + r - t), axis=-1)
        return np.sqrt(np.sum((h + r - t) ** 2, axis=-1))

    def _score_triple(self, h: int, r: int, t: int) -> float:
        d = self._distance(
            self._entity_embeddings[h],
            self._relation_embeddings[r],
            self._entity_embeddings[t],
        )
        return -float(d)

    def _compute_loss(self, pos: np.ndarray, neg: np.ndarray) -> float:
        pos_d = self._distance(
            self._entity_embeddings[pos[:, 0]],
            self._relation_embeddings[pos[:, 1]],
            self._entity_embeddings[pos[:, 2]],
        )
        neg_d = self._distance(
            self._entity_embeddings[neg[:, 0]],
            self._relation_embeddings[neg[:, 1]],
            self._entity_embeddings[neg[:, 2]],
        )
        loss = np.maximum(0, self.margin + pos_d - neg_d)
        return float(np.mean(loss))

    def _update_embeddings(self, pos: np.ndarray, neg: np.ndarray) -> None:
        lr = self.learning_rate
        for i in range(len(pos)):
            h_p, r_p, t_p = pos[i]
            h_n, r_n, t_n = neg[i]

            pos_d = self._distance(
                self._entity_embeddings[h_p],
                self._relation_embeddings[r_p],
                self._entity_embeddings[t_p],
            )
            neg_d = self._distance(
                self._entity_embeddings[h_n],
                self._relation_embeddings[r_n],
                self._entity_embeddings[t_n],
            )

            if self.margin + pos_d - neg_d > 0:
                grad = self._entity_embeddings[h_p] + self._relation_embeddings[r_p] - self._entity_embeddings[t_p]
                sign = np.sign(grad) if self.norm == 1 else grad

                self._entity_embeddings[h_p] -= lr * sign
                self._relation_embeddings[r_p] -= lr * sign
                self._entity_embeddings[t_p] += lr * sign


class RotatEModel(BaseEmbeddingModel):
    """RotatE: Relational rotation in complex space.

    Sun et al. "RotatE: Knowledge Graph Embedding by Relational Rotation"
    (ICLR 2019).
    """

    def __init__(
        self,
        embedding_dim: int = 128,
        margin: float = 6.0,
        learning_rate: float = 0.001,
        seed: Optional[int] = None,
    ):
        super().__init__(embedding_dim, learning_rate, seed)
        self.margin = margin

    def _to_complex(self, emb: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        half = emb.shape[-1] // 2
        return emb[..., :half], emb[..., half:]

    def _score_triple(self, h: int, r: int, t: int) -> float:
        h_re, h_im = self._to_complex(self._entity_embeddings[h])
        r_re, r_im = self._to_complex(self._relation_embeddings[r])
        t_re, t_im = self._to_complex(self._entity_embeddings[t])

        # Rotation: h * r
        hr_re = h_re * r_re - h_im * r_im
        hr_im = h_re * r_im + h_im * r_re

        dist = np.sqrt(np.sum((hr_re - t_re) ** 2 + (hr_im - t_im) ** 2))
        return -float(dist)

    def _compute_loss(self, pos: np.ndarray, neg: np.ndarray) -> float:
        pos_scores = np.array([
            self._score_triple(h, r, t) for h, r, t in pos
        ])
        neg_scores = np.array([
            self._score_triple(h, r, t) for h, r, t in neg
        ])
        loss = np.maximum(0, self.margin - pos_scores + neg_scores)
        return float(np.mean(loss))

    def _update_embeddings(self, pos: np.ndarray, neg: np.ndarray) -> None:
        lr = self.learning_rate
        for i in range(min(len(pos), 50)):  # Cap gradient updates per batch
            h_p, r_p, t_p = pos[i]
            h_re, h_im = self._to_complex(self._entity_embeddings[h_p])
            r_re, r_im = self._to_complex(self._relation_embeddings[r_p])
            t_re, t_im = self._to_complex(self._entity_embeddings[t_p])

            hr_re = h_re * r_re - h_im * r_im
            hr_im = h_re * r_im + h_im * r_re

            diff_re = hr_re - t_re
            diff_im = hr_im - t_im

            half = self.embedding_dim // 2
            self._entity_embeddings[h_p, :half] -= lr * diff_re * r_re
            self._entity_embeddings[h_p, half:] -= lr * diff_im * r_im
            self._entity_embeddings[t_p, :half] += lr * diff_re
            self._entity_embeddings[t_p, half:] += lr * diff_im


class EmbeddingTrainer:
    """Convenience trainer for KG embedding models.

    Usage::

        trainer = EmbeddingTrainer(TrainingConfig(model_type="transe", epochs=50))
        model, history = trainer.train(triples)
    """

    def __init__(self, config: Optional[TrainingConfig] = None):
        self.config = config or TrainingConfig()

    def create_model(self) -> BaseEmbeddingModel:
        if self.config.model_type == "transe":
            return TransEModel(
                embedding_dim=self.config.embedding_dim,
                margin=self.config.margin,
                norm=self.config.norm,
                learning_rate=self.config.learning_rate,
                seed=self.config.seed,
            )
        elif self.config.model_type == "rotate":
            return RotatEModel(
                embedding_dim=self.config.embedding_dim,
                margin=self.config.margin,
                learning_rate=self.config.learning_rate,
                seed=self.config.seed,
            )
        else:
            raise ValueError(f"Unknown model type: {self.config.model_type}")

    def train(
        self,
        triples: List[Tuple[str, str, str]],
    ) -> Tuple[BaseEmbeddingModel, TrainingHistory]:
        """Train a KG embedding model.

        Args:
            triples: List of (head, relation, tail) string triples.

        Returns:
            Tuple of (trained model, training history).
        """
        model = self.create_model()
        history = model.train(
            triples,
            epochs=self.config.epochs,
            batch_size=self.config.batch_size,
            n_negative=self.config.n_negative,
        )
        return model, history
