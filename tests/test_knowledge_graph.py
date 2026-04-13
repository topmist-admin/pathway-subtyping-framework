"""Tests for the knowledge_graph module."""

import json
import os
import tempfile

import numpy as np
import pytest

nx = pytest.importorskip("networkx")

from pathway_subtyping.knowledge_graph import (  # noqa: E402
    EDGE_TYPE_METADATA,
    Edge,
    EdgeType,
    GraphSchema,
    KnowledgeGraph,
    KnowledgeGraphBuilder,
    KnowledgeGraphStats,
    Node,
    NodeFeatures,
    NodeMapping,
    NodeType,
    PPINetwork,
    create_node_mapping,
    load_ppi_from_file,
    to_adjacency_matrix,
    to_csv,
    to_edge_list,
    to_neo4j_cypher,
)

# --- Schema tests ---


class TestNodeType:
    def test_enum_values(self):
        assert NodeType.GENE.value == "gene"
        assert NodeType.PATHWAY.value == "pathway"
        assert NodeType.GO_TERM.value == "go_term"
        assert NodeType.DRUG.value == "drug"

    def test_all_types_present(self):
        expected = {
            "gene",
            "pathway",
            "go_term",
            "cell_type",
            "drug",
            "protein",
            "variant",
            "phenotype",
        }
        assert {nt.value for nt in NodeType} == expected


class TestEdgeType:
    def test_enum_values(self):
        assert EdgeType.GENE_INTERACTS.value == "gene_interacts_gene"
        assert EdgeType.GENE_IN_PATHWAY.value == "gene_in_pathway"
        assert EdgeType.DRUG_TARGETS.value == "drug_targets_gene"

    def test_count(self):
        assert len(EdgeType) == 17


class TestGraphSchema:
    def test_default_schema(self):
        schema = GraphSchema.default_schema()
        assert len(schema.node_types) == 8
        assert len(schema.edge_types) == 17
        assert len(schema.valid_edges) == 17

    def test_is_valid_edge(self):
        schema = GraphSchema.default_schema()
        assert schema.is_valid_edge(NodeType.GENE, NodeType.GENE, EdgeType.GENE_INTERACTS)
        assert schema.is_valid_edge(NodeType.GENE, NodeType.PATHWAY, EdgeType.GENE_IN_PATHWAY)
        # Invalid: gene->gene with drug_targets
        assert not schema.is_valid_edge(NodeType.GENE, NodeType.GENE, EdgeType.DRUG_TARGETS)


class TestNode:
    def test_to_dict(self):
        node = Node(id="BRCA1", node_type=NodeType.GENE, attributes={"symbol": "BRCA1"})
        d = node.to_dict()
        assert d["id"] == "BRCA1"
        assert d["node_type"] == "gene"
        assert d["attributes"]["symbol"] == "BRCA1"

    def test_equality(self):
        n1 = Node(id="BRCA1", node_type=NodeType.GENE)
        n2 = Node(id="BRCA1", node_type=NodeType.GENE)
        n3 = Node(id="BRCA1", node_type=NodeType.PROTEIN)
        assert n1 == n2
        assert n1 != n3

    def test_hash(self):
        n1 = Node(id="BRCA1", node_type=NodeType.GENE)
        n2 = Node(id="BRCA1", node_type=NodeType.GENE)
        assert hash(n1) == hash(n2)
        assert len({n1, n2}) == 1


class TestEdge:
    def test_to_dict(self):
        edge = Edge(source="A", target="B", edge_type=EdgeType.GENE_INTERACTS, weight=0.9)
        d = edge.to_dict()
        assert d["source"] == "A"
        assert d["target"] == "B"
        assert d["edge_type"] == "gene_interacts_gene"
        assert d["weight"] == 0.9

    def test_equality(self):
        e1 = Edge(source="A", target="B", edge_type=EdgeType.GENE_INTERACTS)
        e2 = Edge(source="A", target="B", edge_type=EdgeType.GENE_INTERACTS)
        e3 = Edge(source="A", target="B", edge_type=EdgeType.GENE_COEXPRESSED)
        assert e1 == e2
        assert e1 != e3


class TestNodeFeatures:
    def test_to_dict(self):
        nf = NodeFeatures(
            node_id="BRCA1",
            node_type=NodeType.GENE,
            features=[1.0, 2.0, 3.0],
            feature_names=["f1", "f2", "f3"],
        )
        d = nf.to_dict()
        assert d["node_id"] == "BRCA1"
        assert d["n_features"] == 3


class TestEdgeTypeMetadata:
    def test_metadata_has_entries(self):
        assert EdgeType.GENE_INTERACTS in EDGE_TYPE_METADATA
        meta = EDGE_TYPE_METADATA[EdgeType.GENE_INTERACTS]
        assert meta["symmetric"] is True
        assert meta["source"] == "STRING"

    def test_no_autism_specific_references(self):
        for edge_type, meta in EDGE_TYPE_METADATA.items():
            assert "autism" not in meta.get("source", "").lower()
            assert "autism" not in meta.get("description", "").lower()


# --- KnowledgeGraph tests ---


class TestKnowledgeGraph:
    def test_create_empty(self):
        kg = KnowledgeGraph()
        assert kg.n_nodes == 0
        assert kg.n_edges == 0

    def test_add_node(self):
        kg = KnowledgeGraph()
        kg.add_node("BRCA1", NodeType.GENE)
        assert kg.n_nodes == 1
        assert kg.has_node("BRCA1")
        assert kg.get_node_type("BRCA1") == NodeType.GENE

    def test_add_node_with_attributes(self):
        kg = KnowledgeGraph()
        kg.add_node("BRCA1", NodeType.GENE, {"symbol": "BRCA1", "chr": "17"})
        node = kg.get_node("BRCA1")
        assert node is not None
        assert node.attributes["symbol"] == "BRCA1"
        assert node.attributes["chr"] == "17"

    def test_add_multiple_nodes(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C"], NodeType.GENE)
        assert kg.n_nodes == 3

    def test_add_edge(self):
        kg = KnowledgeGraph()
        kg.add_node("A", NodeType.GENE)
        kg.add_node("B", NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS, weight=0.8)
        assert kg.n_edges == 1
        assert kg.has_edge("A", "B")
        assert kg.has_edge("A", "B", EdgeType.GENE_INTERACTS)

    def test_add_edge_with_weight(self):
        kg = KnowledgeGraph()
        kg.add_node("A", NodeType.GENE)
        kg.add_node("B", NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS, weight=0.95)
        edges = kg.get_edges_by_type(EdgeType.GENE_INTERACTS)
        assert len(edges) == 1
        assert edges[0][2] == 0.95

    def test_add_edge_missing_node_raises(self):
        kg = KnowledgeGraph()
        kg.add_node("A", NodeType.GENE)
        with pytest.raises(ValueError, match="Target node not found"):
            kg.add_edge("A", "B", EdgeType.GENE_INTERACTS)

    def test_add_invalid_edge_type_raises(self):
        kg = KnowledgeGraph()
        kg.add_node("A", NodeType.GENE)
        kg.add_node("B", NodeType.GENE)
        with pytest.raises(ValueError, match="Invalid edge type"):
            kg.add_edge("A", "B", EdgeType.DRUG_TARGETS)

    def test_add_edges_batch(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C", "D"], NodeType.GENE)
        added = kg.add_edges(
            [("A", "B"), ("B", "C"), ("C", "D")],
            EdgeType.GENE_INTERACTS,
        )
        assert added == 3
        assert kg.n_edges == 3

    def test_add_edges_skip_missing(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B"], NodeType.GENE)
        added = kg.add_edges(
            [("A", "B"), ("A", "Z")],
            EdgeType.GENE_INTERACTS,
            skip_missing_nodes=True,
        )
        assert added == 1

    def test_get_neighbors(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C"], NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS)
        kg.add_edge("A", "C", EdgeType.GENE_INTERACTS)
        neighbors = kg.get_neighbors("A", direction="out")
        assert set(neighbors) == {"B", "C"}

    def test_get_neighbors_by_edge_type(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C"], NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS)
        kg.add_edge("A", "C", EdgeType.GENE_COEXPRESSED)
        ppi = kg.get_neighbors("A", edge_type=EdgeType.GENE_INTERACTS)
        assert ppi == ["B"]

    def test_get_nodes_by_type(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B"], NodeType.GENE)
        kg.add_node("P1", NodeType.PATHWAY)
        genes = kg.get_nodes_by_type(NodeType.GENE)
        pathways = kg.get_nodes_by_type(NodeType.PATHWAY)
        assert set(genes) == {"A", "B"}
        assert pathways == ["P1"]

    def test_get_stats(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C"], NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS)
        kg.add_edge("B", "C", EdgeType.GENE_INTERACTS)

        stats = kg.get_stats()
        assert stats.n_nodes == 3
        assert stats.n_edges == 2
        assert "gene" in stats.node_type_counts
        assert stats.node_type_counts["gene"] == 3
        assert stats.avg_degree > 0

        # Test to_dict
        d = stats.to_dict()
        assert d["n_nodes"] == 3
        assert d["n_edges"] == 2

        # Test format_report
        report = stats.format_report()
        assert "Knowledge Graph Statistics" in report
        assert "gene: 3" in report

    def test_subgraph(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C", "D"], NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS)
        kg.add_edge("B", "C", EdgeType.GENE_INTERACTS)
        kg.add_edge("C", "D", EdgeType.GENE_INTERACTS)

        sub = kg.subgraph(["A", "B", "C"])
        assert sub.n_nodes == 3
        assert sub.n_edges == 2
        assert sub.has_node("A")
        assert not sub.has_node("D")

    def test_save_and_load_gpickle(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B"], NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS, weight=0.8)

        with tempfile.NamedTemporaryFile(suffix=".gpickle", delete=False) as f:
            path = f.name

        try:
            kg.save(path)
            loaded = KnowledgeGraph.load(path)
            assert loaded.n_nodes == 2
            assert loaded.n_edges == 1
            assert loaded.has_edge("A", "B")
        finally:
            os.unlink(path)

    def test_save_and_load_json(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B"], NodeType.GENE)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS, weight=0.5)

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name

        try:
            kg.save(path)
            loaded = KnowledgeGraph.load(path)
            assert loaded.n_nodes == 2
            assert loaded.n_edges == 1
        finally:
            os.unlink(path)

    def test_get_citations(self):
        kg = KnowledgeGraph()
        citations = kg.get_citations()
        assert len(citations) >= 2
        assert any("STRING" in c for c in citations)
        assert any("Gene Ontology" in c for c in citations)


# --- KnowledgeGraphBuilder tests ---


class TestKnowledgeGraphBuilder:
    def test_add_genes(self):
        builder = KnowledgeGraphBuilder()
        builder.add_genes(["BRCA1", "TP53", "EGFR"])
        kg = builder.build()
        assert kg.n_nodes == 3
        genes = kg.get_nodes_by_type(NodeType.GENE)
        assert set(genes) == {"BRCA1", "TP53", "EGFR"}

    def test_add_pathways_from_dict(self):
        builder = KnowledgeGraphBuilder()
        pathways = {
            "P1": ["A", "B", "C"],
            "P2": ["B", "C", "D"],
        }
        builder.add_pathways_from_dict(
            pathways, pathway_names={"P1": "Pathway One", "P2": "Pathway Two"}
        )
        kg = builder.build()

        # 4 genes + 2 pathways
        assert kg.n_nodes == 6
        pathway_nodes = kg.get_nodes_by_type(NodeType.PATHWAY)
        assert set(pathway_nodes) == {"P1", "P2"}

        # Gene-pathway edges
        p1_node = kg.get_node("P1")
        assert p1_node is not None
        assert p1_node.attributes["name"] == "Pathway One"

    def test_add_pathways_duck_typed(self):
        class MockDB:
            pathways = {"P1": {"X", "Y"}}
            pathway_names = {"P1": "Test Pathway"}

        builder = KnowledgeGraphBuilder()
        builder.add_pathways(MockDB())
        kg = builder.build()
        assert kg.has_node("P1")
        assert kg.has_node("X")
        assert kg.has_node("Y")

    def test_add_go_terms(self):
        builder = KnowledgeGraphBuilder()
        builder.add_genes(["A", "B"])
        builder.add_go_terms(
            go_terms={
                "GO:0001": {"name": "cell cycle", "namespace": "biological_process"},
            },
            annotations={"A": {"GO:0001"}},
        )
        kg = builder.build()

        go_nodes = kg.get_nodes_by_type(NodeType.GO_TERM)
        assert "GO:0001" in go_nodes
        assert kg.has_edge("A", "GO:0001", EdgeType.GENE_HAS_GO)

    def test_add_go_hierarchy(self):
        builder = KnowledgeGraphBuilder()
        builder.add_go_terms(
            go_terms={
                "GO:0001": {"name": "parent", "namespace": "bp"},
                "GO:0002": {"name": "child", "namespace": "bp"},
            },
            annotations={},
        )
        builder.add_go_hierarchy([("GO:0002", "GO:0001", "is_a")])
        kg = builder.build()
        assert kg.has_edge("GO:0002", "GO:0001", EdgeType.GO_IS_A)

    def test_add_ppi(self):
        ppi = PPINetwork(
            interactions=[("A", "B", 900), ("B", "C", 800)],
            source="STRING",
            score_threshold=700,
        )
        builder = KnowledgeGraphBuilder()
        builder.add_ppi(ppi)
        kg = builder.build()

        assert kg.has_node("A")
        assert kg.has_node("B")
        assert kg.has_node("C")
        ppi_edges = kg.get_edges_by_type(EdgeType.GENE_INTERACTS)
        assert len(ppi_edges) == 2

    def test_builder_chaining(self):
        result = (
            KnowledgeGraphBuilder()
            .add_genes(["A", "B"])
            .add_pathways_from_dict({"P1": ["A", "B"]})
            .build()
        )
        assert result.n_nodes == 3  # 2 genes + 1 pathway

    def test_builder_clear(self):
        builder = KnowledgeGraphBuilder()
        builder.add_genes(["A", "B", "C"])
        builder.clear()
        kg = builder.build()
        assert kg.n_nodes == 0


# --- PPINetwork tests ---


class TestPPINetwork:
    def test_create(self):
        ppi = PPINetwork(
            interactions=[("A", "B", 900), ("C", "D", 500)],
            source="STRING",
        )
        assert len(ppi) == 2

    def test_filter_by_score(self):
        ppi = PPINetwork(
            interactions=[("A", "B", 900), ("C", "D", 500), ("E", "F", 200)],
        )
        filtered = ppi.filter_by_score(600)
        assert len(filtered) == 1
        assert filtered.score_threshold == 600

    def test_to_dict(self):
        ppi = PPINetwork(
            interactions=[("A", "B", 900)],
            source="STRING",
            score_threshold=400,
        )
        d = ppi.to_dict()
        assert d["n_interactions"] == 1
        assert d["source"] == "STRING"
        assert d["score_threshold"] == 400


# --- load_ppi_from_file test ---


class TestLoadPPIFromFile:
    def test_load_string_format(self):
        content = (
            "protein1 protein2 combined_score\n"
            "9606.ENSP001 9606.ENSP002 900\n"
            "9606.ENSP002 9606.ENSP003 500\n"
            "9606.ENSP003 9606.ENSP004 300\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            path = f.name

        try:
            ppi = load_ppi_from_file(path, min_score=400)
            assert len(ppi) == 2  # 300 filtered out
            assert ppi.source == "STRING"
            assert ppi.score_threshold == 400
        finally:
            os.unlink(path)


# --- Exporter tests ---


class TestExporters:
    @pytest.fixture
    def sample_kg(self):
        kg = KnowledgeGraph()
        kg.add_nodes(["A", "B", "C"], NodeType.GENE)
        kg.add_node("P1", NodeType.PATHWAY)
        kg.add_edge("A", "B", EdgeType.GENE_INTERACTS, weight=0.9)
        kg.add_edge("B", "C", EdgeType.GENE_INTERACTS, weight=0.7)
        kg.add_edge("A", "P1", EdgeType.GENE_IN_PATHWAY)
        return kg

    def test_create_node_mapping(self, sample_kg):
        mapping = create_node_mapping(sample_kg)
        assert len(mapping) == 4
        for node_id in ["A", "B", "C", "P1"]:
            idx = mapping.get_idx(node_id)
            assert mapping.get_node(idx) == node_id

    def test_node_mapping_to_dict(self, sample_kg):
        mapping = create_node_mapping(sample_kg)
        d = mapping.to_dict()
        assert d["n_nodes"] == 4

    def test_node_mapping_save_load(self, sample_kg):
        mapping = create_node_mapping(sample_kg)
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name

        try:
            mapping.save(path)
            loaded = NodeMapping.load(path)
            assert len(loaded) == len(mapping)
            for node_id in ["A", "B", "C", "P1"]:
                assert loaded.get_idx(node_id) == mapping.get_idx(node_id)
        finally:
            os.unlink(path)

    def test_to_csv(self, sample_kg):
        with tempfile.TemporaryDirectory() as tmpdir:
            files = to_csv(sample_kg, tmpdir)
            assert "nodes" in files
            assert "edges" in files

            with open(files["nodes"]) as f:
                lines = f.readlines()
            assert len(lines) == 5  # header + 4 nodes

            with open(files["edges"]) as f:
                lines = f.readlines()
            assert len(lines) == 4  # header + 3 edges

    def test_sparse_adjacency_matrix(self, sample_kg):
        adj, mapping = to_adjacency_matrix(sample_kg, use_sparse=True)
        assert adj.shape == (4, 4)
        # A->B edge should be present
        a_idx = mapping.get_idx("A")
        b_idx = mapping.get_idx("B")
        assert adj[a_idx, b_idx] == 0.9

    def test_dense_adjacency_matrix(self, sample_kg):
        adj, mapping = to_adjacency_matrix(sample_kg, use_sparse=False)
        assert isinstance(adj, np.ndarray)
        assert adj.shape == (4, 4)

    def test_adjacency_filter_by_edge_type(self, sample_kg):
        adj, mapping = to_adjacency_matrix(sample_kg, edge_types=[EdgeType.GENE_INTERACTS])
        # Only gene-gene PPI edges should be included
        a_idx = mapping.get_idx("A")
        p_idx = mapping.get_idx("P1")
        assert adj[a_idx, p_idx] == 0.0

    def test_edge_list(self, sample_kg):
        edges = to_edge_list(sample_kg)
        assert len(edges) == 3

    def test_edge_list_to_file(self, sample_kg):
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            path = f.name

        try:
            edges = to_edge_list(sample_kg, output_path=path)
            assert len(edges) == 3
            with open(path) as f:
                lines = f.readlines()
            assert len(lines) == 3
        finally:
            os.unlink(path)

    def test_neo4j_cypher(self, sample_kg):
        with tempfile.NamedTemporaryFile(suffix=".cypher", delete=False) as f:
            path = f.name

        try:
            to_neo4j_cypher(sample_kg, path)
            with open(path) as f:
                content = f.read()
            assert "CREATE" in content
            assert "MATCH" in content
            assert "INDEX" in content
        finally:
            os.unlink(path)


# --- Integration test ---


class TestIntegration:
    def test_full_build_pipeline(self):
        """Test the complete build → stats → export pipeline."""
        builder = KnowledgeGraphBuilder()

        builder.add_genes(["BRCA1", "TP53", "EGFR", "MYC", "PTEN"])
        builder.add_pathways_from_dict(
            {
                "DNA_repair": ["BRCA1", "TP53"],
                "cell_cycle": ["TP53", "MYC", "PTEN"],
            },
            pathway_names={
                "DNA_repair": "DNA Repair Pathway",
                "cell_cycle": "Cell Cycle Regulation",
            },
        )

        ppi = PPINetwork(
            interactions=[
                ("BRCA1", "TP53", 950),
                ("TP53", "MYC", 800),
                ("EGFR", "MYC", 700),
            ]
        )
        builder.add_ppi(ppi)

        kg = builder.build()

        # Verify structure
        assert kg.n_nodes == 7  # 5 genes + 2 pathways
        genes = kg.get_nodes_by_type(NodeType.GENE)
        assert len(genes) == 5
        pathways = kg.get_nodes_by_type(NodeType.PATHWAY)
        assert len(pathways) == 2

        # Verify stats
        stats = kg.get_stats()
        assert stats.n_nodes == 7
        assert stats.n_edges > 0
        assert stats.n_connected_components >= 1

        # Verify subgraph
        sub = kg.subgraph(["BRCA1", "TP53"])
        assert sub.n_nodes == 2

        # Verify export
        adj, mapping = to_adjacency_matrix(kg, edge_types=[EdgeType.GENE_INTERACTS])
        assert adj.shape[0] == 7
        brca1_idx = mapping.get_idx("BRCA1")
        tp53_idx = mapping.get_idx("TP53")
        assert adj[brca1_idx, tp53_idx] > 0
