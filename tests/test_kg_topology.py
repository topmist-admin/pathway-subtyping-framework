"""
Tests for KG topology-aware methods (v0.5 QC layer support).

Tests cover:
    - Pathway gene extraction
    - Pathway subgraph extraction
    - Directed edge filtering
    - In/out degree computation
    - Pathway gene partitioning (upstream/intermediate/downstream)
    - Topological sorting
    - Cascade path finding
    - GENE_REGULATES edge type
"""

import pytest

from pathway_subtyping.knowledge_graph.builder import KnowledgeGraph
from pathway_subtyping.knowledge_graph.schema import EdgeType, NodeType


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def signaling_kg():
    """Build a small KG with a directed signaling cascade.

    Pathway structure (MAPK-like):
        RAS -> RAF -> MEK -> ERK -> TF1
                              |-> TF2
    Plus PPI edges (symmetric):
        RAS -- SCAFFOLD
        ERK -- PHOSPHATASE
    """
    kg = KnowledgeGraph()

    # Genes
    genes = ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2", "SCAFFOLD", "PHOSPHATASE"]
    for g in genes:
        kg.add_node(g, NodeType.GENE)

    # Pathway
    kg.add_node("MAPK_PATHWAY", NodeType.PATHWAY, {"name": "MAPK Signaling"})
    for g in ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2"]:
        kg.add_edge(g, "MAPK_PATHWAY", EdgeType.GENE_IN_PATHWAY)

    # Directed signaling edges (GENE_REGULATES)
    kg.add_edge("RAS", "RAF", EdgeType.GENE_REGULATES, weight=0.9)
    kg.add_edge("RAF", "MEK", EdgeType.GENE_REGULATES, weight=0.9)
    kg.add_edge("MEK", "ERK", EdgeType.GENE_REGULATES, weight=0.9)
    kg.add_edge("ERK", "TF1", EdgeType.GENE_REGULATES, weight=0.8)
    kg.add_edge("ERK", "TF2", EdgeType.GENE_REGULATES, weight=0.7)

    # Symmetric PPI edges
    kg.add_edge("RAS", "SCAFFOLD", EdgeType.GENE_INTERACTS, weight=0.5)
    kg.add_edge("ERK", "PHOSPHATASE", EdgeType.GENE_INTERACTS, weight=0.6)

    return kg


@pytest.fixture
def ppi_only_kg():
    """KG with only PPI edges (no directed signaling) for fallback testing."""
    kg = KnowledgeGraph()
    genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E", "GENE_F"]
    for g in genes:
        kg.add_node(g, NodeType.GENE)

    kg.add_node("PPI_PATHWAY", NodeType.PATHWAY)
    for g in genes:
        kg.add_edge(g, "PPI_PATHWAY", EdgeType.GENE_IN_PATHWAY)

    # Hub-spoke PPI topology: GENE_C is the hub
    kg.add_edge("GENE_A", "GENE_C", EdgeType.GENE_INTERACTS)
    kg.add_edge("GENE_B", "GENE_C", EdgeType.GENE_INTERACTS)
    kg.add_edge("GENE_C", "GENE_D", EdgeType.GENE_INTERACTS)
    kg.add_edge("GENE_C", "GENE_E", EdgeType.GENE_INTERACTS)
    kg.add_edge("GENE_D", "GENE_F", EdgeType.GENE_INTERACTS)

    return kg


# ── Pathway Gene Extraction ──────────────────────────────────────────

class TestGetPathwayGenes:

    def test_returns_pathway_genes(self, signaling_kg):
        genes = signaling_kg.get_pathway_genes("MAPK_PATHWAY")
        assert set(genes) == {"RAS", "RAF", "MEK", "ERK", "TF1", "TF2"}

    def test_excludes_non_pathway_genes(self, signaling_kg):
        genes = signaling_kg.get_pathway_genes("MAPK_PATHWAY")
        assert "SCAFFOLD" not in genes
        assert "PHOSPHATASE" not in genes

    def test_empty_for_nonexistent_pathway(self, signaling_kg):
        genes = signaling_kg.get_pathway_genes("NONEXISTENT")
        assert genes == []


# ── Pathway Subgraph ─────────────────────────────────────────────────

class TestGetPathwaySubgraph:

    def test_subgraph_contains_only_pathway_genes(self, signaling_kg):
        sub = signaling_kg.get_pathway_subgraph("MAPK_PATHWAY")
        assert sub.n_nodes == 6  # RAS, RAF, MEK, ERK, TF1, TF2
        assert sub.has_node("RAS")
        assert not sub.has_node("SCAFFOLD")

    def test_subgraph_preserves_directed_edges(self, signaling_kg):
        sub = signaling_kg.get_pathway_subgraph("MAPK_PATHWAY")
        assert sub.has_edge("RAS", "RAF", EdgeType.GENE_REGULATES)
        assert sub.has_edge("ERK", "TF1", EdgeType.GENE_REGULATES)

    def test_subgraph_excludes_external_edges(self, signaling_kg):
        sub = signaling_kg.get_pathway_subgraph("MAPK_PATHWAY")
        # RAS-SCAFFOLD PPI edge should not be in subgraph
        assert not sub.has_node("SCAFFOLD")


# ── Directed Edge Filtering ──────────────────────────────────────────

class TestDirectedEdges:

    def test_returns_only_directed_edges(self, signaling_kg):
        directed = signaling_kg.get_directed_edges_in_pathway(
            "MAPK_PATHWAY", exclude_symmetric=True
        )
        sources_targets = [(s, t) for s, t, _ in directed]
        # Should include GENE_REGULATES but not GENE_INTERACTS
        assert ("RAS", "RAF") in sources_targets
        assert ("MEK", "ERK") in sources_targets
        # No PPI edges (even though RAS-SCAFFOLD exists, SCAFFOLD is not in pathway)
        assert len(directed) == 5  # 5 GENE_REGULATES edges

    def test_includes_symmetric_when_requested(self, signaling_kg):
        all_edges = signaling_kg.get_directed_edges_in_pathway(
            "MAPK_PATHWAY", exclude_symmetric=False
        )
        # Should include GENE_REGULATES edges (within pathway genes only)
        assert len(all_edges) >= 5


# ── In/Out Degree ─────────────────────────────────────────────────────

class TestDegree:

    def test_in_degree_upstream(self, signaling_kg):
        gene_set = set(signaling_kg.get_pathway_genes("MAPK_PATHWAY"))
        # RAS has no incoming edges from pathway genes
        assert signaling_kg.get_in_degree("RAS", within_nodes=gene_set) == 0

    def test_in_degree_downstream(self, signaling_kg):
        gene_set = set(signaling_kg.get_pathway_genes("MAPK_PATHWAY"))
        # TF1 has 1 incoming edge (from ERK)
        assert signaling_kg.get_in_degree("TF1", within_nodes=gene_set) >= 1

    def test_out_degree_upstream(self, signaling_kg):
        gene_set = set(signaling_kg.get_pathway_genes("MAPK_PATHWAY"))
        # RAS has 1 outgoing GENE_REGULATES edge to RAF
        assert signaling_kg.get_out_degree("RAS", within_nodes=gene_set) >= 1

    def test_out_degree_downstream(self, signaling_kg):
        gene_set = set(signaling_kg.get_pathway_genes("MAPK_PATHWAY"))
        # TF1 has no outgoing edges to pathway genes
        assert signaling_kg.get_out_degree("TF1", EdgeType.GENE_REGULATES, within_nodes=gene_set) == 0

    def test_degree_nonexistent_node(self, signaling_kg):
        assert signaling_kg.get_in_degree("NONEXISTENT") == 0
        assert signaling_kg.get_out_degree("NONEXISTENT") == 0


# ── Pathway Gene Partitioning ─────────────────────────────────────────

class TestPartitionPathwayGenes:

    def test_directed_partitioning(self, signaling_kg):
        layers = signaling_kg.partition_pathway_genes("MAPK_PATHWAY")
        assert "upstream" in layers
        assert "intermediate" in layers
        assert "downstream" in layers

        # RAS is upstream (no incoming directed edges within pathway)
        assert "RAS" in layers["upstream"]

        # TF1, TF2 are downstream (no outgoing directed edges)
        downstream_set = set(layers["downstream"])
        assert "TF1" in downstream_set or "TF2" in downstream_set

        # RAF, MEK, ERK are intermediate
        intermediate_set = set(layers["intermediate"])
        assert "RAF" in intermediate_set or "MEK" in intermediate_set

    def test_all_genes_accounted_for(self, signaling_kg):
        layers = signaling_kg.partition_pathway_genes("MAPK_PATHWAY")
        all_genes = set(layers["upstream"] + layers["intermediate"] + layers["downstream"])
        expected = set(signaling_kg.get_pathway_genes("MAPK_PATHWAY"))
        assert all_genes == expected

    def test_ppi_fallback_partitioning(self, ppi_only_kg):
        """When no directed edges exist, falls back to degree-based thirds."""
        layers = ppi_only_kg.partition_pathway_genes("PPI_PATHWAY")
        all_genes = layers["upstream"] + layers["intermediate"] + layers["downstream"]
        assert len(all_genes) == 6

    def test_empty_pathway(self, signaling_kg):
        layers = signaling_kg.partition_pathway_genes("NONEXISTENT")
        assert layers == {"upstream": [], "intermediate": [], "downstream": []}


# ── Topological Sort ──────────────────────────────────────────────────

class TestTopologicalSort:

    def test_topological_order(self, signaling_kg):
        order = signaling_kg.topological_sort_pathway("MAPK_PATHWAY")
        assert len(order) == 6

        # RAS must come before RAF, RAF before MEK, etc.
        idx = {g: i for i, g in enumerate(order)}
        assert idx["RAS"] < idx["RAF"]
        assert idx["RAF"] < idx["MEK"]
        assert idx["MEK"] < idx["ERK"]

    def test_empty_pathway(self, signaling_kg):
        order = signaling_kg.topological_sort_pathway("NONEXISTENT")
        assert order == []


# ── Cascade Path Finding ─────────────────────────────────────────────

class TestCascadePaths:

    def test_finds_direct_path(self, signaling_kg):
        paths = signaling_kg.find_cascade_paths("RAS", "TF1", "MAPK_PATHWAY")
        assert len(paths) >= 1
        # Path should be RAS -> RAF -> MEK -> ERK -> TF1
        assert any(len(p) == 5 for p in paths)

    def test_no_path_for_disconnected(self, signaling_kg):
        # TF1 -> RAS has no path (wrong direction)
        paths = signaling_kg.find_cascade_paths("TF1", "RAS", "MAPK_PATHWAY")
        assert len(paths) == 0

    def test_nonexistent_nodes(self, signaling_kg):
        paths = signaling_kg.find_cascade_paths("NONEXISTENT", "RAS")
        assert paths == []


# ── GENE_REGULATES Edge Type ──────────────────────────────────────────

class TestGeneRegulatesEdge:

    def test_gene_regulates_is_directed(self):
        from pathway_subtyping.knowledge_graph.schema import EDGE_TYPE_METADATA
        meta = EDGE_TYPE_METADATA.get(EdgeType.GENE_REGULATES, {})
        assert meta.get("symmetric") is False

    def test_gene_regulates_in_schema(self):
        from pathway_subtyping.knowledge_graph.schema import GraphSchema
        schema = GraphSchema.default_schema()
        assert EdgeType.GENE_REGULATES in schema.edge_types

    def test_can_add_gene_regulates_edge(self):
        kg = KnowledgeGraph()
        kg.add_node("A", NodeType.GENE)
        kg.add_node("B", NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_REGULATES, weight=0.9)
        assert kg.has_edge("A", "B", EdgeType.GENE_REGULATES)


# ── Shared Genes & Crosstalk ─────────────────────────────────────────

class TestCrosstalk:

    @pytest.fixture
    def crosstalk_kg(self):
        """KG with two pathways sharing some genes."""
        kg = KnowledgeGraph()
        for g in ["A", "B", "C", "D", "E", "F"]:
            kg.add_node(g, NodeType.GENE)

        kg.add_node("PW1", NodeType.PATHWAY)
        kg.add_node("PW2", NodeType.PATHWAY)

        # PW1 has A, B, C, D
        for g in ["A", "B", "C", "D"]:
            kg.add_edge(g, "PW1", EdgeType.GENE_IN_PATHWAY)

        # PW2 has C, D, E, F (shares C, D)
        for g in ["C", "D", "E", "F"]:
            kg.add_edge(g, "PW2", EdgeType.GENE_IN_PATHWAY)

        # Inter-pathway PPI edge
        kg.add_edge("B", "E", EdgeType.GENE_INTERACTS, weight=0.7)

        return kg

    def test_shared_genes(self, crosstalk_kg):
        shared = crosstalk_kg.get_shared_genes("PW1", "PW2")
        assert set(shared) == {"C", "D"}

    def test_no_shared_genes(self, crosstalk_kg):
        shared = crosstalk_kg.get_shared_genes("PW1", "NONEXISTENT")
        assert shared == []

    def test_crosstalk_metrics(self, crosstalk_kg):
        ct = crosstalk_kg.get_pathway_crosstalk("PW1", "PW2")
        assert ct["n_shared"] == 2
        assert ct["pathway_a_size"] == 4
        assert ct["pathway_b_size"] == 4
        assert ct["jaccard_index"] == 2 / 6  # 2 shared / 6 union
        assert ct["n_inter_edges"] >= 1  # B-E edge


# ── KnowledgeGraphBuilder Signaling Edges ─────────────────────────────

class TestBuilderSignalingEdges:

    def test_add_signaling_edges(self):
        from pathway_subtyping.knowledge_graph.builder import KnowledgeGraphBuilder
        builder = KnowledgeGraphBuilder()
        builder.add_genes(["RAS", "RAF", "MEK"])
        builder.add_signaling_edges([("RAS", "RAF"), ("RAF", "MEK")])
        kg = builder.build()
        assert kg.has_edge("RAS", "RAF", EdgeType.GENE_REGULATES)
        assert kg.has_edge("RAF", "MEK", EdgeType.GENE_REGULATES)
        # Directed — reverse should not exist
        assert not kg.has_edge("RAF", "RAS", EdgeType.GENE_REGULATES)

    def test_add_signaling_edges_from_dict(self):
        from pathway_subtyping.knowledge_graph.builder import KnowledgeGraphBuilder
        builder = KnowledgeGraphBuilder()
        builder.add_genes(["A", "B", "C", "D"])
        builder.add_signaling_edges_from_dict({
            "MAPK": [("A", "B"), ("B", "C")],
            "PI3K": [("C", "D")],
        })
        kg = builder.build()
        assert kg.has_edge("A", "B", EdgeType.GENE_REGULATES)
        assert kg.has_edge("C", "D", EdgeType.GENE_REGULATES)

    def test_signaling_edges_with_pathways(self):
        from pathway_subtyping.knowledge_graph.builder import KnowledgeGraphBuilder
        builder = KnowledgeGraphBuilder()
        builder.add_pathways_from_dict({"MAPK": ["RAS", "RAF", "MEK", "ERK"]})
        builder.add_signaling_edges([
            ("RAS", "RAF"), ("RAF", "MEK"), ("MEK", "ERK"),
        ])
        kg = builder.build()

        # Should partition correctly
        layers = kg.partition_pathway_genes("MAPK")
        assert "RAS" in layers["upstream"]
        assert "ERK" in layers["downstream"]

    def test_clear_resets_signaling(self):
        from pathway_subtyping.knowledge_graph.builder import KnowledgeGraphBuilder
        builder = KnowledgeGraphBuilder()
        builder.add_signaling_edges([("A", "B")])
        builder.clear()
        assert len(builder._signaling_edges) == 0
