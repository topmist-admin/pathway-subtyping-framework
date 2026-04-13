"""
Tests for GNN Layer (#21, #22, #23).

KG embedding tests (TransE, RotatE) are pure numpy — always run.
GNN model/attention tests require torch — skipped if unavailable.

EXPERIMENTAL/BETA features.
"""

import numpy as np
import pytest

from pathway_subtyping.gnn.kg_embeddings import (
    BaseEmbeddingModel,
    EmbeddingTrainer,
    NodeEmbeddings,
    RotatEModel,
    TrainingConfig,
    TransEModel,
)

# Check torch availability
try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


# ══════════════════════════════════════════════════════════════════════
# KG EMBEDDINGS (pure numpy — always run)
# ══════════════════════════════════════════════════════════════════════

class TestNodeEmbeddings:

    def test_basic_operations(self):
        emb = NodeEmbeddings(
            node_ids=["A", "B", "C"],
            embeddings=np.array([[1, 0], [0, 1], [1, 1]], dtype=float),
            source="test",
        )
        assert emb.n_nodes == 3
        assert emb.dim == 2
        assert emb.get("A") is not None
        assert emb.get("NONEXISTENT") is None

    def test_most_similar(self):
        emb = NodeEmbeddings(
            node_ids=["A", "B", "C"],
            embeddings=np.array([[1, 0], [0.9, 0.1], [0, 1]], dtype=float),
        )
        similar = emb.most_similar("A", k=2)
        assert len(similar) == 2
        assert similar[0][0] == "B"  # B is most similar to A

    def test_get_batch(self):
        emb = NodeEmbeddings(
            node_ids=["A", "B", "C"],
            embeddings=np.array([[1, 0], [0, 1], [1, 1]], dtype=float),
        )
        batch = emb.get_batch(["A", "C"])
        assert batch.shape == (2, 2)

    def test_to_dict(self):
        emb = NodeEmbeddings(
            node_ids=["A"], embeddings=np.array([[1.0, 2.0]]), source="test"
        )
        d = emb.to_dict()
        assert d["n_nodes"] == 1
        assert d["dim"] == 2


class TestTransE:

    @pytest.fixture
    def triples(self):
        return [
            ("RAS", "activates", "RAF"),
            ("RAF", "activates", "MEK"),
            ("MEK", "activates", "ERK"),
            ("ERK", "activates", "TF1"),
            ("RAS", "interacts", "SCAFFOLD"),
            ("PIK3CA", "activates", "AKT"),
        ]

    def test_train(self, triples):
        model = TransEModel(embedding_dim=32, seed=42)
        history = model.train(triples, epochs=20, batch_size=4)
        assert len(history.losses) == 20
        assert history.best_loss < float("inf")

    def test_get_embeddings(self, triples):
        model = TransEModel(embedding_dim=16, seed=42)
        model.train(triples, epochs=10)
        emb = model.get_node_embeddings()
        assert emb.n_nodes == 8  # 8 unique entities
        assert emb.dim == 16

    def test_predict_link(self, triples):
        model = TransEModel(embedding_dim=16, seed=42)
        model.train(triples, epochs=20)
        score = model.predict_link("RAS", "activates", "RAF")
        assert score is not None
        assert isinstance(score, float)

    def test_predict_missing_entity(self, triples):
        model = TransEModel(embedding_dim=16, seed=42)
        model.train(triples, epochs=5)
        assert model.predict_link("NONEXISTENT", "activates", "RAF") is None

    def test_l2_norm(self, triples):
        model = TransEModel(embedding_dim=16, norm=2, seed=42)
        history = model.train(triples, epochs=10)
        assert len(history.losses) == 10


class TestRotatE:

    @pytest.fixture
    def triples(self):
        return [
            ("A", "r1", "B"),
            ("B", "r1", "C"),
            ("C", "r2", "D"),
            ("A", "r2", "D"),
        ]

    def test_train(self, triples):
        model = RotatEModel(embedding_dim=32, seed=42)
        history = model.train(triples, epochs=15, batch_size=4)
        assert len(history.losses) == 15

    def test_get_embeddings(self, triples):
        model = RotatEModel(embedding_dim=16, seed=42)
        model.train(triples, epochs=10)
        emb = model.get_node_embeddings()
        assert emb.n_nodes == 4
        assert emb.dim == 16

    def test_predict_link(self, triples):
        model = RotatEModel(embedding_dim=16, seed=42)
        model.train(triples, epochs=15)
        score = model.predict_link("A", "r1", "B")
        assert score is not None


class TestEmbeddingTrainer:

    def test_transe_trainer(self):
        config = TrainingConfig(model_type="transe", embedding_dim=16, epochs=10, seed=42)
        trainer = EmbeddingTrainer(config)
        triples = [("A", "r", "B"), ("B", "r", "C"), ("A", "r", "C")]
        model, history = trainer.train(triples)
        assert isinstance(model, TransEModel)
        assert len(history.losses) == 10

    def test_rotate_trainer(self):
        config = TrainingConfig(model_type="rotate", embedding_dim=16, epochs=10, seed=42)
        trainer = EmbeddingTrainer(config)
        triples = [("A", "r", "B"), ("B", "r", "C")]
        model, history = trainer.train(triples)
        assert isinstance(model, RotatEModel)

    def test_invalid_model_type(self):
        config = TrainingConfig(model_type="invalid")
        trainer = EmbeddingTrainer(config)
        with pytest.raises(ValueError, match="Unknown model type"):
            trainer.create_model()


# ══════════════════════════════════════════════════════════════════════
# EMBEDDING FUSION (numpy — always run)
# ══════════════════════════════════════════════════════════════════════

class TestEmbeddingFusion:

    def test_concat_fusion(self):
        from pathway_subtyping.gnn.embeddings.fusion import (
            EmbeddingFusion,
            FusionConfig,
            FusionMethod,
        )
        emb_a = NodeEmbeddings(
            node_ids=["G1", "G2"], embeddings=np.array([[1, 0], [0, 1]], dtype=float)
        )
        emb_b = NodeEmbeddings(
            node_ids=["G1", "G2"], embeddings=np.array([[0.5], [0.8]], dtype=float)
        )
        fusion = EmbeddingFusion(FusionConfig(method=FusionMethod.CONCAT, normalize=False))
        result = fusion.fuse({"a": emb_a, "b": emb_b})
        assert result.dim == 3  # 2 + 1

    def test_pca_fusion(self):
        from pathway_subtyping.gnn.embeddings.fusion import (
            EmbeddingFusion,
            FusionConfig,
            FusionMethod,
        )
        rng = np.random.default_rng(42)
        emb_a = NodeEmbeddings(
            node_ids=[f"G{i}" for i in range(20)],
            embeddings=rng.normal(0, 1, (20, 10)),
        )
        emb_b = NodeEmbeddings(
            node_ids=[f"G{i}" for i in range(20)],
            embeddings=rng.normal(0, 1, (20, 8)),
        )
        fusion = EmbeddingFusion(
            FusionConfig(method=FusionMethod.PCA, output_dim=5, normalize=False)
        )
        result = fusion.fuse({"a": emb_a, "b": emb_b})
        assert result.dim == 5

    def test_single_source(self):
        from pathway_subtyping.gnn.embeddings.fusion import EmbeddingFusion
        emb = NodeEmbeddings(
            node_ids=["G1"], embeddings=np.array([[1.0, 2.0]])
        )
        fusion = EmbeddingFusion()
        result = fusion.fuse({"only": emb})
        assert result is emb  # Returns same object

    def test_no_common_nodes(self):
        from pathway_subtyping.gnn.embeddings.fusion import EmbeddingFusion
        emb_a = NodeEmbeddings(node_ids=["A"], embeddings=np.array([[1.0]]))
        emb_b = NodeEmbeddings(node_ids=["B"], embeddings=np.array([[2.0]]))
        fusion = EmbeddingFusion()
        with pytest.raises(ValueError, match="No common"):
            fusion.fuse({"a": emb_a, "b": emb_b})

    def test_weighted_sum(self):
        from pathway_subtyping.gnn.embeddings.fusion import (
            EmbeddingFusion,
            FusionConfig,
            FusionMethod,
        )
        emb_a = NodeEmbeddings(
            node_ids=["G1"], embeddings=np.array([[2.0, 0.0]])
        )
        emb_b = NodeEmbeddings(
            node_ids=["G1"], embeddings=np.array([[0.0, 2.0]])
        )
        fusion = EmbeddingFusion(FusionConfig(
            method=FusionMethod.WEIGHTED_SUM,
            weights={"a": 0.7, "b": 0.3},
            normalize=False,
        ))
        result = fusion.fuse({"a": emb_a, "b": emb_b})
        assert result.dim == 2


# ══════════════════════════════════════════════════════════════════════
# GNN CONFIG (always available)
# ══════════════════════════════════════════════════════════════════════

class TestGNNConfig:

    def test_default_config(self):
        from pathway_subtyping.gnn.config import GNNConfig
        config = GNNConfig()
        assert config.input_dim == 256
        assert config.num_layers == 3
        assert len(config.edge_types) > 0

    def test_custom_config(self):
        from pathway_subtyping.gnn.config import GNNConfig
        config = GNNConfig(input_dim=64, hidden_dim=64, num_layers=2)
        assert config.input_dim == 64
        assert config.num_layers == 2

    def test_to_dict(self):
        from pathway_subtyping.gnn.config import GNNConfig
        d = GNNConfig().to_dict()
        assert "input_dim" in d
        assert "edge_types" in d


# ══════════════════════════════════════════════════════════════════════
# GNN MODEL + LAYERS + ATTENTION (requires torch)
# ══════════════════════════════════════════════════════════════════════

@pytest.mark.skipif(not HAS_TORCH, reason="PyTorch not installed")
class TestGNNModel:

    def test_ontology_gnn_forward(self):
        from pathway_subtyping.gnn.config import GNNConfig
        from pathway_subtyping.gnn.model import OntologyAwareGNN

        config = GNNConfig(
            input_dim=16, hidden_dim=16, output_dim=8,
            num_layers=2, num_heads=2,
            edge_types=["gene_interacts_gene", "gene_regulates_gene"],
            prior_types=[],
        )
        model = OntologyAwareGNN(config)

        n_nodes = 10
        n_edges = 15
        x = torch.randn(n_nodes, 16)
        edge_index = torch.randint(0, n_nodes, (2, n_edges))
        edge_type = torch.randint(0, 2, (n_edges,))

        output = model(x, edge_index, edge_type)
        assert output.node_embeddings is not None
        assert "all" in output.node_embeddings

    def test_with_labels(self):
        from pathway_subtyping.gnn.config import GNNConfig
        from pathway_subtyping.gnn.model import OntologyAwareGNN

        config = GNNConfig(
            input_dim=16, hidden_dim=16, output_dim=8,
            num_layers=1, num_heads=2,
            edge_types=["gene_interacts_gene"],
            prior_types=[],
            task_heads=["gene_classification"],
        )
        model = OntologyAwareGNN(config)

        n = 10
        x = torch.randn(n, 16)
        ei = torch.randint(0, n, (2, 20))
        et = torch.zeros(20, dtype=torch.long)
        labels = torch.randint(0, 4, (n,))

        output = model(x, ei, et, labels=labels)
        assert output.loss is not None
        assert output.gene_logits is not None
        assert output.gene_logits.shape == (n, 4)

    def test_gnn_output_to_dict(self):
        from pathway_subtyping.gnn.model import GNNOutput
        out = GNNOutput(loss=0.5)
        d = out.to_dict()
        assert d["loss"] == 0.5


@pytest.mark.skipif(not HAS_TORCH, reason="PyTorch not installed")
class TestGNNLayers:

    def test_message_passing(self):
        from pathway_subtyping.gnn.layers import MessagePassingLayer
        layer = MessagePassingLayer(16, 16, ["r1", "r2"])
        x = torch.randn(10, 16)
        ei = torch.randint(0, 10, (2, 20))
        et = torch.randint(0, 2, (20,))
        out = layer(x, ei, et, ["r1", "r2"])
        assert out.shape == (10, 16)

    def test_hierarchical_aggregator(self):
        from pathway_subtyping.gnn.layers import HierarchicalAggregator
        agg = HierarchicalAggregator(16, num_heads=2)
        x = torch.randn(10, 16)
        out = agg(x)
        assert out.shape == (10, 16)

    def test_bio_prior_weighting(self):
        from pathway_subtyping.gnn.layers import BioPriorWeighting
        bpw = BioPriorWeighting(16, ["pli", "expression"])
        x = torch.randn(10, 16)
        priors = {"pli": torch.rand(10), "expression": torch.rand(10)}
        out = bpw(x, priors)
        assert out.shape == (10, 16)


@pytest.mark.skipif(not HAS_TORCH, reason="PyTorch not installed")
class TestAttention:

    def test_biological_attention(self):
        from pathway_subtyping.gnn.attention import BiologicalAttention
        attn = BiologicalAttention(16, num_heads=2, prior_types=["pli"])
        q = torch.randn(10, 16)
        k = torch.randn(10, 16)
        v = torch.randn(10, 16)
        priors = {"pli": torch.rand(10)}
        out, weights = attn(q, k, v, bio_priors=priors)
        assert out.shape == (10, 16)
        assert weights.shape == (10, 10)

    def test_attention_without_priors(self):
        from pathway_subtyping.gnn.attention import BiologicalAttention
        attn = BiologicalAttention(16, num_heads=2)
        q = torch.randn(5, 16)
        out, weights = attn(q, q, q)
        assert out.shape == (5, 16)

    def test_pathway_co_attention(self):
        from pathway_subtyping.gnn.attention import PathwayCoAttention
        coattn = PathwayCoAttention(16, 16, num_heads=2)
        genes = torch.randn(20, 16)
        pathways = torch.randn(5, 16)
        g_out, p_out = coattn(genes, pathways)
        assert g_out.shape == (20, 16)
        assert p_out.shape == (5, 16)


@pytest.mark.skipif(not HAS_TORCH, reason="PyTorch not installed")
class TestGNNTrainer:

    def test_train_step(self):
        from pathway_subtyping.gnn.config import GNNConfig
        from pathway_subtyping.gnn.model import GNNTrainer, OntologyAwareGNN

        config = GNNConfig(
            input_dim=8, hidden_dim=8, output_dim=4,
            num_layers=1, num_heads=2,
            edge_types=["r1"],
            prior_types=[],
            task_heads=["gene_classification"],
        )
        model = OntologyAwareGNN(config)
        trainer = GNNTrainer(model, learning_rate=1e-3)

        n = 10
        x = torch.randn(n, 8)
        ei = torch.randint(0, n, (2, 15))
        et = torch.zeros(15, dtype=torch.long)
        labels = torch.randint(0, 4, (n,))

        metrics = trainer.train_step(x, ei, et, labels, edge_type_names=["r1"])
        assert "loss" in metrics
        assert metrics["loss"] > 0

    def test_evaluate(self):
        from pathway_subtyping.gnn.config import GNNConfig
        from pathway_subtyping.gnn.model import GNNTrainer, OntologyAwareGNN

        config = GNNConfig(
            input_dim=8, hidden_dim=8, output_dim=4,
            num_layers=1, num_heads=2,
            edge_types=["r1"],
            prior_types=[],
            task_heads=["gene_classification"],
        )
        model = OntologyAwareGNN(config)
        trainer = GNNTrainer(model)

        n = 10
        x = torch.randn(n, 8)
        ei = torch.randint(0, n, (2, 15))
        et = torch.zeros(15, dtype=torch.long)
        labels = torch.randint(0, 4, (n,))

        metrics = trainer.evaluate(x, ei, et, labels, edge_type_names=["r1"])
        assert "loss" in metrics
        assert "accuracy" in metrics
