"""Tests for the network_propagation module."""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.network_propagation import (
    NetworkPropagator,
    PropagationConfig,
    PropagationMethod,
    PropagationResult,
    propagate_network_scores,
)

# --- Enum and config tests ---


class TestPropagationMethod:
    def test_enum_values(self):
        assert PropagationMethod.RANDOM_WALK.value == "random_walk"
        assert PropagationMethod.HEAT_DIFFUSION.value == "heat_diffusion"
        assert PropagationMethod.INSULATED.value == "insulated"
        assert PropagationMethod.PAGERANK.value == "pagerank"


class TestPropagationConfig:
    def test_default_config(self):
        config = PropagationConfig()
        assert config.method == PropagationMethod.RANDOM_WALK
        assert config.restart_prob == 0.5
        assert config.n_iterations == 100
        assert config.seed is None

    def test_custom_config(self):
        config = PropagationConfig(
            method=PropagationMethod.HEAT_DIFFUSION,
            diffusion_time=0.5,
            seed=42,
        )
        assert config.method == PropagationMethod.HEAT_DIFFUSION
        assert config.diffusion_time == 0.5
        assert config.seed == 42


# --- PropagationResult tests ---


class TestPropagationResult:
    def test_to_dict(self):
        result = PropagationResult(
            gene_scores={"A": 1.0, "B": 0.5, "C": 0.1},
            n_iterations=10,
            converged=True,
            method="random_walk",
            metadata={"alpha": 0.5},
        )
        d = result.to_dict()
        assert d["method"] == "random_walk"
        assert d["n_iterations"] == 10
        assert d["converged"] is True
        assert d["n_genes_scored"] == 3
        assert "A" in d["top_genes"]

    def test_format_report(self):
        result = PropagationResult(
            gene_scores={"A": 1.0, "B": 0.5},
            n_iterations=5,
            converged=True,
            method="random_walk",
        )
        report = result.format_report()
        assert "Network Propagation Result" in report
        assert "random_walk" in report
        assert "A:" in report

    def test_get_citations_rwr(self):
        result = PropagationResult(
            gene_scores={}, n_iterations=1, converged=True, method="random_walk"
        )
        citations = result.get_citations()
        assert len(citations) >= 1
        assert any("Kohler" in c for c in citations)

    def test_get_citations_heat(self):
        result = PropagationResult(
            gene_scores={}, n_iterations=1, converged=True, method="heat_diffusion"
        )
        citations = result.get_citations()
        assert len(citations) >= 1
        assert any("Vandin" in c for c in citations)


# --- NetworkPropagator tests ---


class TestNetworkPropagator:
    """Test propagation on synthetic networks."""

    @pytest.fixture
    def linear_chain(self):
        """Build a linear chain: A-B-C-D-E."""
        edges = [
            ("A", "B", 1.0),
            ("B", "C", 1.0),
            ("C", "D", 1.0),
            ("D", "E", 1.0),
        ]
        propagator = NetworkPropagator()
        propagator.build_network_from_edges(edges)
        return propagator

    @pytest.fixture
    def clique(self):
        """Build a 4-node clique."""
        edges = [
            ("A", "B", 1.0),
            ("A", "C", 1.0),
            ("A", "D", 1.0),
            ("B", "C", 1.0),
            ("B", "D", 1.0),
            ("C", "D", 1.0),
        ]
        propagator = NetworkPropagator()
        propagator.build_network_from_edges(edges)
        return propagator

    def test_build_network_from_edges(self):
        edges = [("A", "B", 1.0), ("B", "C", 1.0)]
        propagator = NetworkPropagator()
        propagator.build_network_from_edges(edges)
        stats = propagator.get_network_stats()
        assert stats["n_nodes"] == 3
        assert stats["n_edges"] == 2

    def test_propagate_without_build_raises(self):
        propagator = NetworkPropagator()
        with pytest.raises(RuntimeError, match="Network not built"):
            propagator.propagate({"A": 1.0})

    def test_random_walk_propagation(self, linear_chain):
        """Signal should propagate from A outward."""
        result = linear_chain.propagate({"A": 1.0})
        assert result.method == "random_walk"
        assert result.converged
        assert "A" in result.gene_scores
        assert "B" in result.gene_scores

    def test_signal_decay_with_distance(self, linear_chain):
        """Scores should decrease with distance from seed."""
        result = linear_chain.propagate({"A": 1.0})
        scores = result.gene_scores
        # A should have highest score, then B, C, D, E
        assert scores.get("A", 0) >= scores.get("B", 0)
        assert scores.get("B", 0) >= scores.get("C", 0)
        assert scores.get("C", 0) >= scores.get("D", 0)

    def test_heat_diffusion_propagation(self, linear_chain):
        config = PropagationConfig(method=PropagationMethod.HEAT_DIFFUSION)
        propagator = NetworkPropagator(config)
        propagator.build_network_from_edges(
            [("A", "B", 1.0), ("B", "C", 1.0), ("C", "D", 1.0), ("D", "E", 1.0)]
        )
        result = propagator.propagate({"A": 1.0})
        assert result.method == "heat_diffusion"
        assert result.converged
        assert len(result.gene_scores) > 0

    def test_insulated_diffusion_propagation(self, linear_chain):
        config = PropagationConfig(method=PropagationMethod.INSULATED)
        propagator = NetworkPropagator(config)
        propagator.build_network_from_edges([("A", "B", 1.0), ("B", "C", 1.0)])
        result = propagator.propagate({"A": 1.0})
        assert result.method == "insulated"
        assert result.converged

    def test_pagerank_propagation(self, linear_chain):
        config = PropagationConfig(method=PropagationMethod.PAGERANK)
        propagator = NetworkPropagator(config)
        propagator.build_network_from_edges([("A", "B", 1.0), ("B", "C", 1.0)])
        result = propagator.propagate({"A": 1.0})
        assert result.method == "pagerank"

    def test_convergence(self, clique):
        """Small clique should converge quickly."""
        result = clique.propagate({"A": 1.0})
        assert result.converged
        assert result.n_iterations < 100

    def test_get_network_stats(self, linear_chain):
        stats = linear_chain.get_network_stats()
        assert stats["n_nodes"] == 5
        assert stats["n_edges"] == 4
        assert stats["avg_degree"] > 0
        assert stats["max_degree"] >= 2
        assert stats["density"] > 0

    def test_propagate_score_matrix(self, linear_chain):
        """Test batch propagation via DataFrame."""
        df = pd.DataFrame(
            {"A": [1.0, 0.0], "B": [0.0, 1.0], "C": [0.0, 0.0]},
            index=["sample1", "sample2"],
        )
        result = linear_chain.propagate_score_matrix(df)
        assert isinstance(result, pd.DataFrame)
        assert list(result.index) == ["sample1", "sample2"]
        # Sample 1 seeded A, so A should have signal
        assert result.loc["sample1", "A"] > 0
        # Sample 2 seeded B, so B should have signal
        assert result.loc["sample2", "B"] > 0

    def test_multiple_seeds(self, linear_chain):
        """Multiple seeds should both propagate."""
        result = linear_chain.propagate({"A": 1.0, "E": 1.0})
        scores = result.gene_scores
        # Middle node C should receive signal from both ends
        assert "C" in scores
        assert scores["C"] > 0

    def test_weighted_edges(self):
        """Higher weight edges should propagate more signal."""
        edges_strong = [("A", "B", 10.0), ("B", "C", 1.0)]
        edges_weak = [("A", "B", 1.0), ("B", "C", 10.0)]

        prop_strong = NetworkPropagator(PropagationConfig(normalize_edges=False))
        prop_strong.build_network_from_edges(edges_strong)
        result_strong = prop_strong.propagate({"A": 1.0})

        prop_weak = NetworkPropagator(PropagationConfig(normalize_edges=False))
        prop_weak.build_network_from_edges(edges_weak)
        result_weak = prop_weak.propagate({"A": 1.0})

        # With strong A-B edge, B should get at least as much signal. The two
        # steady-state B scores are near-identical here (both ~0.102), so compare
        # with a floating-point tolerance: without it the last-digit difference
        # (~3e-17) flips the >= across numpy/BLAS builds (passes on macOS, fails on
        # Linux). The tolerance still catches a real regression where the strong
        # edge gives B meaningfully less signal.
        assert result_strong.gene_scores.get("B", 0) >= result_weak.gene_scores.get("B", 0) - 1e-9


# --- Convenience function tests ---


class TestConvenienceFunction:
    def test_propagate_network_scores(self):
        edges = [("A", "B", 1.0), ("B", "C", 1.0)]
        result = propagate_network_scores(
            gene_scores={"A": 1.0},
            edges=edges,
            method=PropagationMethod.RANDOM_WALK,
        )
        assert isinstance(result, PropagationResult)
        assert result.converged
        assert "A" in result.gene_scores

    def test_convenience_heat_diffusion(self):
        edges = [("X", "Y", 1.0), ("Y", "Z", 1.0)]
        result = propagate_network_scores(
            gene_scores={"X": 1.0},
            edges=edges,
            method=PropagationMethod.HEAT_DIFFUSION,
            diffusion_time=0.5,
        )
        assert result.method == "heat_diffusion"
        assert len(result.gene_scores) > 0
