"""
Tests for the v0.6 F3 KG refresh infrastructure.

Covers:
    - pathway_subtyping.knowledge_graph.sources — versioned manifest
    - pathway_subtyping.knowledge_graph.diff — structural KG diff
    - pathway_subtyping.knowledge_graph.regression — threshold-flagging
      regression harness
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from pathway_subtyping.knowledge_graph import (
    KG_SOURCES,
    KG_SOURCES_V05,
    KG_SOURCES_V06,
    EdgeType,
    KGDiff,
    KGRegressionReport,
    KGSource,
    KnowledgeGraph,
    KnowledgeGraphBuilder,
    NodeType,
    ScoreDelta,
    diff_kgs,
    get_source,
    list_sources,
    manifest_digest,
    run_kg_regression,
)

# --------------------------------------------------------------------------- #
# Source manifest
# --------------------------------------------------------------------------- #


class TestKGSources:

    def test_v06_registry_has_three_sources(self):
        assert set(KG_SOURCES_V06.keys()) == {"omnipath", "signor", "reactome"}

    def test_v05_and_v06_differ(self):
        assert KG_SOURCES_V05 != KG_SOURCES_V06
        for src_id in KG_SOURCES_V05:
            v5 = KG_SOURCES_V05[src_id]
            v6 = KG_SOURCES_V06[src_id]
            assert v5.release != v6.release

    def test_default_alias_is_v06(self):
        assert KG_SOURCES is KG_SOURCES_V06

    def test_sha256_validates_on_construction(self):
        with pytest.raises(ValueError, match="64"):
            KGSource(
                source_id="test",
                release="1.0",
                url="https://example.com",
                sha256="shorthash",
                release_date="2026-01-01",
                citation="N/A",
            )

    def test_list_sources_returns_all(self):
        srcs = list_sources()
        assert len(srcs) == 3
        assert all(isinstance(s, KGSource) for s in srcs)

    def test_get_source_hit(self):
        src = get_source("omnipath")
        assert src.source_id == "omnipath"
        assert src.release == "2025"

    def test_get_source_miss(self):
        with pytest.raises(KeyError, match="no pinned source"):
            get_source("nonexistent")

    def test_manifest_digest_stable(self):
        d1 = manifest_digest()
        d2 = manifest_digest()
        assert d1 == d2
        assert len(d1) == 64

    def test_manifest_digest_distinguishes_v05_v06(self):
        assert manifest_digest(KG_SOURCES_V05) != manifest_digest(KG_SOURCES_V06)

    def test_verify_archive_missing_file(self, tmp_path):
        src = get_source("omnipath")
        assert src.verify_archive(tmp_path / "not_there.tar.gz") is False

    def test_verify_archive_hash_mismatch(self, tmp_path):
        src = get_source("omnipath")
        p = tmp_path / "bundle.tar.gz"
        p.write_bytes(b"not the pinned release")
        assert src.verify_archive(p) is False

    def test_verify_archive_match_on_crafted_file(self, tmp_path):
        # Build a KGSource whose sha256 matches the content we write.
        content = b"hello world"
        sha = hashlib.sha256(content).hexdigest()
        src = KGSource(
            source_id="custom",
            release="1",
            url="https://example.com",
            sha256=sha,
            release_date="2026-01-01",
            citation="N/A",
        )
        p = tmp_path / "bundle"
        p.write_bytes(content)
        assert src.verify_archive(p) is True


# --------------------------------------------------------------------------- #
# KG diff
# --------------------------------------------------------------------------- #


def _gene_kg(edges, extra_genes=()):
    kg = KnowledgeGraph()
    genes = {e[0] for e in edges} | {e[1] for e in edges} | set(extra_genes)
    for g in sorted(genes):
        kg.add_node(g, NodeType.GENE)
    for src, tgt, etype in edges:
        kg.add_edge(src, tgt, etype, weight=1.0)
    return kg


class TestKGDiff:

    def test_identical_kgs_produce_empty_diff(self):
        edges = [("A", "B", EdgeType.GENE_INTERACTS)]
        kg1 = _gene_kg(edges)
        kg2 = _gene_kg(edges)
        diff = diff_kgs(kg1, kg2)
        assert diff.is_identical()
        assert diff.summary["nodes_added_total"] == 0
        assert diff.summary["edges_added_total"] == 0

    def test_nodes_added(self):
        kg1 = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)])
        kg2 = _gene_kg(
            [("A", "B", EdgeType.GENE_INTERACTS)],
            extra_genes=["C", "D"],
        )
        diff = diff_kgs(kg1, kg2)
        assert diff.nodes_added["gene"] == ["C", "D"]
        assert "gene" not in diff.nodes_removed

    def test_nodes_removed(self):
        kg1 = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)], extra_genes=["X"])
        kg2 = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)])
        diff = diff_kgs(kg1, kg2)
        assert diff.nodes_removed["gene"] == ["X"]

    def test_edges_added_and_removed(self):
        kg1 = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)])
        kg2 = _gene_kg(
            [
                ("A", "B", EdgeType.GENE_INTERACTS),
                ("A", "C", EdgeType.GENE_INTERACTS),
            ]
        )
        diff = diff_kgs(kg1, kg2)
        assert ("A", "C") in diff.edges_added["gene_interacts_gene"]
        assert "gene_interacts_gene" not in diff.edges_removed

    def test_direction_flip_detected(self):
        kg1 = _gene_kg([("A", "B", EdgeType.GENE_REGULATES)])
        kg2 = _gene_kg([("B", "A", EdgeType.GENE_REGULATES)])
        diff = diff_kgs(kg1, kg2)
        assert ("A", "B", "gene_regulates_gene") in diff.edges_flipped
        # Non-flipped directed edges should still be in added/removed
        assert ("A", "B") in diff.edges_removed["gene_regulates_gene"]
        assert ("B", "A") in diff.edges_added["gene_regulates_gene"]

    def test_diff_to_dict(self):
        kg1 = _gene_kg([("A", "B", EdgeType.GENE_REGULATES)])
        kg2 = _gene_kg([("A", "C", EdgeType.GENE_REGULATES)])
        diff = diff_kgs(kg1, kg2)
        d = diff.to_dict()
        assert "nodes_added" in d
        assert "summary" in d
        assert d["summary"]["edges_added_total"] == 1
        assert d["summary"]["edges_removed_total"] == 1


# --------------------------------------------------------------------------- #
# KG regression
# --------------------------------------------------------------------------- #


class TestKGRegression:

    def _edge_count_score(self):
        def score(kg: KnowledgeGraph, _input) -> dict:
            return {
                "n_edges_gene_interacts": len(kg.get_edges_by_type(EdgeType.GENE_INTERACTS)),
                "n_edges_gene_regulates": len(kg.get_edges_by_type(EdgeType.GENE_REGULATES)),
            }

        return score

    def test_identical_kgs_no_deltas_flagged(self):
        edges = [("A", "B", EdgeType.GENE_INTERACTS)]
        kg1 = _gene_kg(edges)
        kg2 = _gene_kg(edges)
        report = run_kg_regression(
            kg1,
            kg2,
            score_fns=[self._edge_count_score()],
            benchmark_inputs=["placeholder"],
            tolerance=0.05,
        )
        assert report.passed
        assert all(d.absolute_delta == 0 for d in report.deltas)

    def test_small_change_under_tolerance(self):
        # KG2 grows edge count by 4 / 100 = 4% — under 5% tolerance
        base_edges = [(f"G{i}", f"G{i+1}", EdgeType.GENE_INTERACTS) for i in range(100)]
        extra_edges = [(f"G{i}", f"G{i+2}", EdgeType.GENE_INTERACTS) for i in range(4)]
        kg1 = _gene_kg(base_edges)
        kg2 = _gene_kg(base_edges + extra_edges)
        report = run_kg_regression(
            kg1,
            kg2,
            score_fns=[self._edge_count_score()],
            benchmark_inputs=["placeholder"],
            tolerance=0.05,
        )
        assert (
            report.passed
        ), f"expected passage; flagged={[d.score_name for d in report.flagged_scores]}"

    def test_large_change_flagged(self):
        base_edges = [(f"G{i}", f"G{i+1}", EdgeType.GENE_INTERACTS) for i in range(100)]
        extra_edges = [(f"G{i}", f"G{i+2}", EdgeType.GENE_INTERACTS) for i in range(20)]
        kg1 = _gene_kg(base_edges)
        kg2 = _gene_kg(base_edges + extra_edges)
        report = run_kg_regression(
            kg1,
            kg2,
            score_fns=[self._edge_count_score()],
            benchmark_inputs=["placeholder"],
            tolerance=0.05,
        )
        assert not report.passed
        flagged = [d.score_name for d in report.flagged_scores]
        assert any("n_edges_gene_interacts" in name for name in flagged)

    def test_rejects_empty_inputs(self):
        kg = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)])
        with pytest.raises(ValueError, match="score_fns"):
            run_kg_regression(kg, kg, score_fns=[], benchmark_inputs=["x"])
        with pytest.raises(ValueError, match="benchmark_inputs"):
            run_kg_regression(kg, kg, score_fns=[self._edge_count_score()], benchmark_inputs=[])

    def test_report_to_dict_structure(self):
        kg = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)])
        report = run_kg_regression(
            kg,
            kg,
            score_fns=[self._edge_count_score()],
            benchmark_inputs=["placeholder"],
        )
        d = report.to_dict()
        assert d["passed"] is True
        assert d["n_scores"] == 2  # n_edges_gene_interacts + n_edges_gene_regulates
        assert all("score_name" in x for x in d["deltas"])

    def test_score_fn_non_dict_rejected(self):
        def bad_score_fn(_kg, _input):
            return "not a dict"

        kg = _gene_kg([("A", "B", EdgeType.GENE_INTERACTS)])
        with pytest.raises(TypeError):
            run_kg_regression(
                kg,
                kg,
                score_fns=[bad_score_fn],  # type: ignore[list-item]
                benchmark_inputs=["x"],
            )
