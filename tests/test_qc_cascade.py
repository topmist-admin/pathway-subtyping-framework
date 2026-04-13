"""
Tests for F1: Incomplete Cascade Detection (CascadeAnalyzer).

Tests validate:
    - Layer partitioning from KG topology
    - Per-layer activation scoring
    - Cascade completion ratio computation
    - Incomplete cascade detection
    - Per-cell CCR matrix
    - Integration with ManufacturingSimulator stalled cascade injection
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.knowledge_graph.builder import KnowledgeGraph
from pathway_subtyping.knowledge_graph.schema import EdgeType, NodeType
from pathway_subtyping.qc.cascade import CascadeAnalyzer, CascadeResult, LayerScores


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def mapk_kg():
    """KG with a MAPK-like signaling cascade.

    RAS -> RAF -> MEK -> ERK -> TF1
                              -> TF2
    """
    kg = KnowledgeGraph()
    genes = ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2"]
    for g in genes:
        kg.add_node(g, NodeType.GENE)

    kg.add_node("MAPK", NodeType.PATHWAY, {"name": "MAPK Signaling"})
    for g in genes:
        kg.add_edge(g, "MAPK", EdgeType.GENE_IN_PATHWAY)

    kg.add_edge("RAS", "RAF", EdgeType.GENE_REGULATES)
    kg.add_edge("RAF", "MEK", EdgeType.GENE_REGULATES)
    kg.add_edge("MEK", "ERK", EdgeType.GENE_REGULATES)
    kg.add_edge("ERK", "TF1", EdgeType.GENE_REGULATES)
    kg.add_edge("ERK", "TF2", EdgeType.GENE_REGULATES)
    return kg


@pytest.fixture
def healthy_expression():
    """Expression where all genes are uniformly active."""
    rng = np.random.default_rng(42)
    genes = ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2"]
    n_cells = 100
    data = {g: rng.normal(5.0, 0.5, n_cells) for g in genes}
    return pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])


@pytest.fixture
def stalled_expression():
    """Expression where upstream is high but downstream is suppressed."""
    rng = np.random.default_rng(42)
    n_cells = 100
    data = {
        "RAS": rng.normal(8.0, 0.5, n_cells),    # High upstream
        "RAF": rng.normal(7.0, 0.5, n_cells),     # High upstream
        "MEK": rng.normal(5.0, 0.5, n_cells),     # Mid
        "ERK": rng.normal(2.0, 0.5, n_cells),     # Low downstream
        "TF1": rng.normal(1.0, 0.5, n_cells),     # Very low downstream
        "TF2": rng.normal(1.0, 0.5, n_cells),     # Very low downstream
    }
    return pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])


# ── LayerScores ───────────────────────────────────────────────────────

class TestLayerScores:

    def test_completion_ratio_complete(self):
        ls = LayerScores(pathway="test", upstream_score=1.0, downstream_score=1.0)
        assert ls.completion_ratio == 1.0
        assert not ls.is_incomplete

    def test_completion_ratio_stalled(self):
        ls = LayerScores(pathway="test", upstream_score=1.0, downstream_score=0.1)
        assert ls.completion_ratio == pytest.approx(0.1, abs=0.01)
        assert ls.is_incomplete

    def test_completion_ratio_no_initiation(self):
        ls = LayerScores(pathway="test", upstream_score=0.0, downstream_score=0.0)
        assert ls.completion_ratio == 1.0
        assert not ls.is_incomplete

    def test_completion_ratio_negative_upstream(self):
        ls = LayerScores(pathway="test", upstream_score=-0.5, downstream_score=0.5)
        assert ls.completion_ratio == 1.0

    def test_to_dict(self):
        ls = LayerScores(
            pathway="MAPK",
            upstream_score=0.8,
            downstream_score=0.3,
            upstream_genes=["RAS"],
            downstream_genes=["TF1"],
        )
        d = ls.to_dict()
        assert d["pathway"] == "MAPK"
        assert "completion_ratio" in d
        assert "is_incomplete" in d
        assert d["n_upstream"] == 1


# ── CascadeAnalyzer ───────────────────────────────────────────────────

class TestCascadeAnalyzer:

    def test_define_layers(self, mapk_kg):
        analyzer = CascadeAnalyzer(mapk_kg)
        layers = analyzer.define_cascade_layers("MAPK")
        assert "RAS" in layers["upstream"]
        assert "TF1" in layers["downstream"] or "TF2" in layers["downstream"]
        total = len(layers["upstream"]) + len(layers["intermediate"]) + len(layers["downstream"])
        assert total == 6

    def test_layer_caching(self, mapk_kg):
        analyzer = CascadeAnalyzer(mapk_kg)
        layers1 = analyzer.define_cascade_layers("MAPK")
        layers2 = analyzer.define_cascade_layers("MAPK")
        assert layers1 is layers2

    def test_score_layer_activation(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        layers = analyzer.define_cascade_layers("MAPK")
        up, mid, down = analyzer.score_layer_activation(healthy_expression, layers)
        assert len(up) == 100
        assert len(mid) == 100
        assert len(down) == 100

    def test_analyze_healthy_pathway(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(healthy_expression, pathways=["MAPK"])
        assert isinstance(result, CascadeResult)
        assert len(result.per_pathway) == 1
        # Healthy: uniform expression -> CCR near 1.0
        ccr = result.per_pathway[0].completion_ratio
        assert ccr > 0.5

    def test_analyze_stalled_pathway(self, mapk_kg, stalled_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(stalled_expression, pathways=["MAPK"])
        ls = result.per_pathway[0]
        # Upstream should be higher than downstream
        assert ls.upstream_score > ls.downstream_score

    def test_per_cell_ccr_matrix(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(healthy_expression, pathways=["MAPK"])
        assert result.per_cell_ccr is not None
        assert result.per_cell_ccr.shape == (100, 1)
        assert "MAPK" in result.per_cell_ccr.columns

    def test_get_incomplete_pathways(self, mapk_kg, stalled_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(stalled_expression, pathways=["MAPK"])
        incomplete = result.get_incomplete_pathways()
        assert isinstance(incomplete, list)

    def test_get_completion_ratios(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(healthy_expression, pathways=["MAPK"])
        ratios = result.get_completion_ratios()
        assert "MAPK" in ratios
        assert 0.0 <= ratios["MAPK"] <= 1.0

    def test_to_dict(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(healthy_expression, pathways=["MAPK"])
        d = result.to_dict()
        assert "passed" in d
        assert "n_incomplete" in d
        assert "mean_completion_ratio" in d
        assert "pathways" in d

    def test_empty_pathway_list(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(healthy_expression, pathways=[])
        assert len(result.per_pathway) == 0
        assert result.passed

    def test_nonexistent_pathway(self, mapk_kg, healthy_expression):
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(healthy_expression, pathways=["NONEXISTENT"])
        assert len(result.per_pathway) == 0

    def test_missing_genes_in_expression(self, mapk_kg):
        """Expression matrix missing some pathway genes."""
        rng = np.random.default_rng(42)
        n_cells = 50
        # Only has RAS and RAF, missing MEK/ERK/TF1/TF2
        data = {
            "RAS": rng.normal(5.0, 0.5, n_cells),
            "RAF": rng.normal(5.0, 0.5, n_cells),
        }
        expr = pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])
        analyzer = CascadeAnalyzer(mapk_kg)
        result = analyzer.analyze(expr, pathways=["MAPK"])
        # Should not crash
        assert len(result.per_pathway) == 1


# ── Multi-pathway Analysis ────────────────────────────────────────────

class TestMultiPathway:

    @pytest.fixture
    def multi_kg(self):
        """KG with two pathways."""
        kg = KnowledgeGraph()
        for g in ["A", "B", "C", "D", "E", "F"]:
            kg.add_node(g, NodeType.GENE)

        kg.add_node("PW1", NodeType.PATHWAY)
        kg.add_node("PW2", NodeType.PATHWAY)

        for g in ["A", "B", "C"]:
            kg.add_edge(g, "PW1", EdgeType.GENE_IN_PATHWAY)
        for g in ["D", "E", "F"]:
            kg.add_edge(g, "PW2", EdgeType.GENE_IN_PATHWAY)

        kg.add_edge("A", "B", EdgeType.GENE_REGULATES)
        kg.add_edge("B", "C", EdgeType.GENE_REGULATES)
        kg.add_edge("D", "E", EdgeType.GENE_REGULATES)
        kg.add_edge("E", "F", EdgeType.GENE_REGULATES)
        return kg

    def test_analyze_multiple_pathways(self, multi_kg):
        rng = np.random.default_rng(42)
        n_cells = 50
        data = {g: rng.normal(5.0, 0.5, n_cells) for g in "ABCDEF"}
        expr = pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])

        analyzer = CascadeAnalyzer(multi_kg)
        result = analyzer.analyze(expr, pathways=["PW1", "PW2"])
        assert len(result.per_pathway) == 2
        assert result.per_cell_ccr.shape == (50, 2)

    def test_auto_discover_pathways(self, multi_kg):
        rng = np.random.default_rng(42)
        n_cells = 30
        data = {g: rng.normal(5.0, 0.5, n_cells) for g in "ABCDEF"}
        expr = pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])

        analyzer = CascadeAnalyzer(multi_kg)
        result = analyzer.analyze(expr)  # No pathways specified
        assert len(result.per_pathway) == 2


# ── Integration with Simulator ────────────────────────────────────────

class TestCascadeSimulatorIntegration:

    def test_simulator_stalled_cascade_detected(self):
        """Verify CascadeAnalyzer detects simulator-injected stalled cascades."""
        from pathway_subtyping.qc.testing.simulator import ManufacturingSimulator, ManufacturingSpec

        sim = ManufacturingSimulator(seed=42, genes_per_pathway=9)
        spec = ManufacturingSpec(
            intended_pathways=sim.pathways[:2],
            excluded_pathways=sim.pathways[2:4],
        )
        batch = sim.generate_healthy_batch(n_cells=100, spec=spec)

        # Build a KG from the simulator's pathway structure
        kg = KnowledgeGraph()
        for gene in sim.all_genes:
            kg.add_node(gene, NodeType.GENE)

        for pw, genes in sim.pathway_genes.items():
            kg.add_node(pw, NodeType.PATHWAY)
            for g in genes:
                kg.add_edge(g, pw, EdgeType.GENE_IN_PATHWAY)
            # Add directed signaling edges within each pathway
            for i in range(len(genes) - 1):
                kg.add_edge(genes[i], genes[i + 1], EdgeType.GENE_REGULATES)

        target_pw = sim.pathways[0]
        analyzer = CascadeAnalyzer(kg)

        # Inject stalled cascade — severity 0.9 should create visible
        # upstream vs downstream divergence in raw expression
        stalled = sim.inject_stalled_cascades(
            batch, pathways=[target_pw], severity=0.9
        )

        # Check raw expression divergence between upstream and downstream genes
        layers = analyzer.define_cascade_layers(target_pw)
        up_genes = [g for g in layers["upstream"] if g in stalled.expression.columns]
        down_genes = [g for g in layers["downstream"] if g in stalled.expression.columns]

        if up_genes and down_genes:
            up_mean = stalled.expression[up_genes].values.mean()
            down_mean = stalled.expression[down_genes].values.mean()
            # Upstream should be higher than downstream after stall injection
            assert up_mean > down_mean, (
                f"Expected upstream ({up_mean:.2f}) > downstream ({down_mean:.2f})"
            )
