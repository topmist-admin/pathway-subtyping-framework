"""
Tests for KG advanced methods: topology-weighted scoring,
hierarchical pathway queries, cross-omics entity resolution.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.knowledge_graph.builder import KnowledgeGraph
from pathway_subtyping.knowledge_graph.schema import EdgeType, NodeType


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def rich_kg():
    """KG with pathways, hierarchy, PPI, drugs, and proteins."""
    kg = KnowledgeGraph()

    # Genes
    for g in ["RAS", "RAF", "MEK", "ERK", "TF1", "AKT", "PIK3CA"]:
        kg.add_node(g, NodeType.GENE)

    # Pathways with hierarchy
    kg.add_node("SIGNALING", NodeType.PATHWAY, {"name": "Signaling"})
    kg.add_node("MAPK", NodeType.PATHWAY, {"name": "MAPK Cascade"})
    kg.add_node("PI3K", NodeType.PATHWAY, {"name": "PI3K/AKT"})

    # SIGNALING contains MAPK and PI3K
    kg.add_edge("SIGNALING", "MAPK", EdgeType.PATHWAY_CONTAINS)
    kg.add_edge("SIGNALING", "PI3K", EdgeType.PATHWAY_CONTAINS)

    # MAPK genes
    for g in ["RAS", "RAF", "MEK", "ERK", "TF1"]:
        kg.add_edge(g, "MAPK", EdgeType.GENE_IN_PATHWAY)
    # PI3K genes (RAS is shared)
    for g in ["RAS", "PIK3CA", "AKT"]:
        kg.add_edge(g, "PI3K", EdgeType.GENE_IN_PATHWAY)
    # SIGNALING direct genes (top-level)
    kg.add_edge("RAS", "SIGNALING", EdgeType.GENE_IN_PATHWAY)

    # Directed signaling
    kg.add_edge("RAS", "RAF", EdgeType.GENE_REGULATES)
    kg.add_edge("RAF", "MEK", EdgeType.GENE_REGULATES)
    kg.add_edge("MEK", "ERK", EdgeType.GENE_REGULATES)
    kg.add_edge("ERK", "TF1", EdgeType.GENE_REGULATES)
    kg.add_edge("RAS", "PIK3CA", EdgeType.GENE_REGULATES)
    kg.add_edge("PIK3CA", "AKT", EdgeType.GENE_REGULATES)

    # PPI
    kg.add_edge("RAF", "MEK", EdgeType.GENE_INTERACTS, weight=0.8)
    kg.add_edge("ERK", "AKT", EdgeType.GENE_INTERACTS, weight=0.5)

    # Proteins
    kg.add_node("KRAS_PROTEIN", NodeType.PROTEIN)
    kg.add_node("BRAF_PROTEIN", NodeType.PROTEIN)
    kg.add_edge("RAS", "KRAS_PROTEIN", EdgeType.GENE_ENCODES)
    kg.add_edge("RAF", "BRAF_PROTEIN", EdgeType.GENE_ENCODES)

    # Drugs
    kg.add_node("SORAFENIB", NodeType.DRUG)
    kg.add_node("TRAMETINIB", NodeType.DRUG)
    kg.add_edge("SORAFENIB", "RAF", EdgeType.DRUG_TARGETS)
    kg.add_edge("TRAMETINIB", "MEK", EdgeType.DRUG_TARGETS)

    return kg


@pytest.fixture
def expression(rich_kg):
    rng = np.random.default_rng(42)
    genes = ["RAS", "RAF", "MEK", "ERK", "TF1", "AKT", "PIK3CA"]
    n = 50
    data = {g: rng.normal(5.0, 1.0, n) for g in genes}
    # Make ERK a hub — higher expression
    data["ERK"] = rng.normal(8.0, 0.5, n)
    return pd.DataFrame(data, index=[f"cell_{i}" for i in range(n)])


# ══════════════════════════════════════════════════════════════════════
# TOPOLOGY-WEIGHTED SCORING
# ══════════════════════════════════════════════════════════════════════

class TestCentrality:

    def test_degree_centrality(self, rich_kg):
        c = rich_kg.compute_centrality("MAPK", method="degree")
        assert len(c) == 5
        assert all(0.0 <= v <= 1.0 for v in c.values())

    def test_betweenness_centrality(self, rich_kg):
        c = rich_kg.compute_centrality("MAPK", method="betweenness")
        assert "MEK" in c  # MEK is on the main path — should have centrality

    def test_pagerank_centrality(self, rich_kg):
        c = rich_kg.compute_centrality("MAPK", method="pagerank")
        assert len(c) == 5
        assert all(v > 0 for v in c.values())

    def test_global_centrality(self, rich_kg):
        c = rich_kg.compute_centrality(method="degree")
        # Should include all gene nodes
        assert "RAS" in c
        assert "AKT" in c

    def test_empty_pathway(self, rich_kg):
        c = rich_kg.compute_centrality("NONEXISTENT", method="degree")
        assert c == {}

    def test_invalid_method_raises(self, rich_kg):
        with pytest.raises(ValueError, match="Unknown centrality"):
            rich_kg.compute_centrality("MAPK", method="invalid")


class TestTopologyWeightedScoring:

    def test_weighted_score_shape(self, rich_kg, expression):
        scores = rich_kg.topology_weighted_pathway_score(
            expression, "MAPK", centrality_method="degree"
        )
        assert len(scores) == 50

    def test_weighted_differs_from_uniform(self, rich_kg, expression):
        weighted = rich_kg.topology_weighted_pathway_score(
            expression, "MAPK", centrality_method="degree"
        )
        # Uniform mean for comparison
        genes = [g for g in rich_kg.get_pathway_genes("MAPK") if g in expression.columns]
        uniform = expression[genes].values.mean(axis=1)
        # They should differ (hub genes weighted differently)
        assert not np.allclose(weighted, uniform, atol=0.01)

    def test_missing_pathway(self, rich_kg, expression):
        scores = rich_kg.topology_weighted_pathway_score(
            expression, "NONEXISTENT"
        )
        assert np.allclose(scores, 0.0)


# ══════════════════════════════════════════════════════════════════════
# HIERARCHICAL PATHWAY QUERIES
# ══════════════════════════════════════════════════════════════════════

class TestHierarchicalQueries:

    def test_get_child_pathways(self, rich_kg):
        children = rich_kg.get_child_pathways("SIGNALING")
        assert set(children) == {"MAPK", "PI3K"}

    def test_get_parent_pathways(self, rich_kg):
        parents = rich_kg.get_parent_pathways("MAPK")
        assert "SIGNALING" in parents

    def test_leaf_has_no_children(self, rich_kg):
        children = rich_kg.get_child_pathways("MAPK")
        assert children == []

    def test_root_has_no_parents(self, rich_kg):
        parents = rich_kg.get_parent_pathways("SIGNALING")
        assert parents == []

    def test_get_pathway_hierarchy(self, rich_kg):
        tree = rich_kg.get_pathway_hierarchy("SIGNALING")
        assert tree["id"] == "SIGNALING"
        assert tree["name"] == "Signaling"
        assert len(tree["children"]) == 2
        child_ids = {c["id"] for c in tree["children"]}
        assert child_ids == {"MAPK", "PI3K"}

    def test_hierarchy_includes_genes(self, rich_kg):
        tree = rich_kg.get_pathway_hierarchy("SIGNALING")
        mapk_child = next(c for c in tree["children"] if c["id"] == "MAPK")
        assert mapk_child["n_genes"] == 5

    def test_get_all_descendant_genes(self, rich_kg):
        all_genes = rich_kg.get_all_descendant_genes("SIGNALING")
        # SIGNALING's direct genes + MAPK genes + PI3K genes (deduplicated)
        assert "RAS" in all_genes
        assert "TF1" in all_genes
        assert "AKT" in all_genes
        # RAS is in both MAPK and PI3K — should appear once
        assert len(all_genes) == len(set(all_genes))

    def test_leaf_descendant_genes(self, rich_kg):
        genes = rich_kg.get_all_descendant_genes("MAPK")
        assert set(genes) == {"RAS", "RAF", "MEK", "ERK", "TF1"}


# ══════════════════════════════════════════════════════════════════════
# CROSS-OMICS ENTITY RESOLUTION
# ══════════════════════════════════════════════════════════════════════

class TestEntityResolution:

    def test_gene_to_protein(self, rich_kg):
        proteins = rich_kg.resolve_gene_to_protein("RAS")
        assert "KRAS_PROTEIN" in proteins

    def test_protein_to_gene(self, rich_kg):
        genes = rich_kg.resolve_protein_to_gene("KRAS_PROTEIN")
        assert "RAS" in genes

    def test_no_protein(self, rich_kg):
        proteins = rich_kg.resolve_gene_to_protein("TF1")
        assert proteins == []

    def test_entity_chain_gene_to_drug(self, rich_kg):
        # RAS -> (GENE_REGULATES) -> RAF <- (DRUG_TARGETS) <- SORAFENIB
        # But DRUG_TARGETS is drug -> gene, so from gene perspective we need in-edges
        # resolve_entity_chain follows out-edges, so gene -> drug won't find it
        # Let's test from drug to gene instead
        paths = rich_kg.resolve_entity_chain("SORAFENIB", NodeType.GENE, max_hops=2)
        # SORAFENIB -> RAF (via DRUG_TARGETS)
        assert len(paths) >= 1
        assert any("RAF" in p for p in paths)

    def test_entity_chain_nonexistent(self, rich_kg):
        paths = rich_kg.resolve_entity_chain("NONEXISTENT", NodeType.DRUG)
        assert paths == []

    def test_drug_targets_in_pathway(self, rich_kg):
        gene_drugs = rich_kg.get_drug_targets_in_pathway("MAPK")
        assert "RAF" in gene_drugs
        assert "SORAFENIB" in gene_drugs["RAF"]
        assert "MEK" in gene_drugs
        assert "TRAMETINIB" in gene_drugs["MEK"]

    def test_no_drugs_in_pathway(self, rich_kg):
        gene_drugs = rich_kg.get_drug_targets_in_pathway("PI3K")
        # PI3K has RAS, PIK3CA, AKT — none have drug targets in our fixture
        # (SORAFENIB targets RAF which is in MAPK, not PI3K)
        assert "PIK3CA" not in gene_drugs or len(gene_drugs.get("PIK3CA", [])) == 0
