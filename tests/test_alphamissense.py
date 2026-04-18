"""
Tests for v0.6 F4 AlphaMissense-modulated cascade scoring.

Covers:
    - AlphaMissenseScorer table validation and lookup semantics
    - weights_from_carriers: default 1.0, correct down-weighting,
      damage_floor, most-damaging retention across duplicate variants
    - CascadeAnalyzer integration: gene_weights=None is bit-identical
      to the variant-naive baseline; passing AM-derived weights for a
      carrier cohort produces a meaningful score shift vs non-carriers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.knowledge_graph.builder import KnowledgeGraph
from pathway_subtyping.knowledge_graph.schema import EdgeType, NodeType
from pathway_subtyping.qc.alphamissense import AlphaMissenseScorer
from pathway_subtyping.qc.cascade import CascadeAnalyzer


# --------------------------------------------------------------------------- #
# Shared fixtures
# --------------------------------------------------------------------------- #

@pytest.fixture
def mapk_kg():
    """Minimal RAS->RAF->MEK->ERK->TF cascade KG."""
    kg = KnowledgeGraph()
    for g in ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2"]:
        kg.add_node(g, NodeType.GENE)
    kg.add_node("MAPK", NodeType.PATHWAY, {"name": "MAPK"})
    for g in ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2"]:
        kg.add_edge(g, "MAPK", EdgeType.GENE_IN_PATHWAY)
    for a, b in [("RAS", "RAF"), ("RAF", "MEK"), ("MEK", "ERK"),
                 ("ERK", "TF1"), ("ERK", "TF2")]:
        kg.add_edge(a, b, EdgeType.GENE_REGULATES)
    return kg


@pytest.fixture
def healthy_expression():
    rng = np.random.default_rng(42)
    genes = ["RAS", "RAF", "MEK", "ERK", "TF1", "TF2"]
    n_cells = 30
    return pd.DataFrame(
        {g: rng.normal(5.0, 0.5, n_cells) for g in genes},
        index=[f"cell_{i}" for i in range(n_cells)],
    )


@pytest.fixture
def score_table():
    """Small AlphaMissense score table: one benign, two pathogenic."""
    return pd.DataFrame({
        "variant_id": ["RAF:p.V600E", "MEK:p.P124S", "TF1:p.S50A"],
        "am_score":   [0.95, 0.80, 0.10],
    })


# --------------------------------------------------------------------------- #
# Scorer basics
# --------------------------------------------------------------------------- #

class TestAlphaMissenseScorer:

    def test_from_table_loads(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        assert scorer.lookup("RAF:p.V600E") == pytest.approx(0.95)
        assert scorer.lookup("MEK:p.P124S") == pytest.approx(0.80)
        assert np.isnan(scorer.lookup("UNKNOWN"))

    def test_empty_scorer(self):
        scorer = AlphaMissenseScorer.empty()
        assert np.isnan(scorer.lookup("anything"))
        assert scorer.summary()["n_variants"] == 0

    def test_missing_column_raises(self):
        with pytest.raises(ValueError, match="variant_id"):
            AlphaMissenseScorer.from_table(
                pd.DataFrame({"foo": [], "am_score": []})
            )

    def test_score_clipped_to_0_1(self):
        scorer = AlphaMissenseScorer.from_table(pd.DataFrame({
            "variant_id": ["A", "B", "C"],
            "am_score":   [-0.2, 1.5, 0.5],
        }))
        assert scorer.lookup("A") == 0.0
        assert scorer.lookup("B") == 1.0
        assert scorer.lookup("C") == pytest.approx(0.5)

    def test_lookup_many(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        got = scorer.lookup_many(["RAF:p.V600E", "UNKNOWN"])
        assert got[0] == pytest.approx(0.95)
        assert np.isnan(got[1])

    def test_summary_buckets(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        s = scorer.summary()
        assert s["n_variants"] == 3
        # Default AM thresholds: likely_pathogenic >= 0.564
        assert s["n_likely_pathogenic"] == 2
        assert s["n_likely_benign"] == 1


# --------------------------------------------------------------------------- #
# Weight construction
# --------------------------------------------------------------------------- #

class TestWeightsFromCarriers:

    def test_default_weight_is_one(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        weights = scorer.weights_from_carriers(
            carriers=pd.DataFrame(
                columns=["cell_id", "gene", "variant_id"]
            ),
            cells=["c1", "c2"],
            genes=["RAF", "MEK"],
        )
        assert (weights == 1.0).all().all()

    def test_carrier_downweighted(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        carriers = pd.DataFrame([
            {"cell_id": "c1", "gene": "RAF", "variant_id": "RAF:p.V600E"},
        ])
        weights = scorer.weights_from_carriers(
            carriers=carriers,
            cells=["c1", "c2"],
            genes=["RAF", "MEK"],
        )
        assert weights.at["c1", "RAF"] == pytest.approx(1.0 - 0.95)
        assert weights.at["c2", "RAF"] == 1.0
        assert weights.at["c1", "MEK"] == 1.0

    def test_damage_floor(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        carriers = pd.DataFrame([
            {"cell_id": "c1", "gene": "RAF", "variant_id": "RAF:p.V600E"},
        ])
        weights = scorer.weights_from_carriers(
            carriers=carriers,
            cells=["c1"],
            genes=["RAF"],
            damage_floor=0.25,
        )
        assert weights.at["c1", "RAF"] == pytest.approx(0.25)

    def test_multiple_variants_retain_most_damaging(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        carriers = pd.DataFrame([
            {"cell_id": "c1", "gene": "RAF", "variant_id": "TF1:p.S50A"},       # benign, 0.10 -> weight 0.90
            {"cell_id": "c1", "gene": "RAF", "variant_id": "RAF:p.V600E"},      # pathogenic, 0.95 -> weight 0.05
        ])
        weights = scorer.weights_from_carriers(
            carriers=carriers, cells=["c1"], genes=["RAF"],
        )
        assert weights.at["c1", "RAF"] == pytest.approx(0.05)

    def test_unknown_variant_skipped(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        carriers = pd.DataFrame([
            {"cell_id": "c1", "gene": "RAF", "variant_id": "UNKNOWN"},
        ])
        weights = scorer.weights_from_carriers(
            carriers=carriers, cells=["c1"], genes=["RAF"],
        )
        assert weights.at["c1", "RAF"] == 1.0

    def test_missing_cell_or_gene_skipped(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        carriers = pd.DataFrame([
            {"cell_id": "c_missing", "gene": "RAF", "variant_id": "RAF:p.V600E"},
            {"cell_id": "c1", "gene": "NOT_A_GENE", "variant_id": "RAF:p.V600E"},
        ])
        weights = scorer.weights_from_carriers(
            carriers=carriers, cells=["c1"], genes=["RAF"],
        )
        assert (weights == 1.0).all().all()

    def test_carriers_missing_column_raises(self):
        scorer = AlphaMissenseScorer.empty()
        with pytest.raises(ValueError, match="variant_id"):
            scorer.weights_from_carriers(
                carriers=pd.DataFrame({"cell_id": [], "gene": []}),
                cells=["c1"], genes=["RAF"],
            )

    def test_damage_floor_bounds(self, score_table):
        scorer = AlphaMissenseScorer.from_table(score_table)
        with pytest.raises(ValueError, match="damage_floor"):
            scorer.weights_from_carriers(
                carriers=pd.DataFrame(columns=["cell_id", "gene", "variant_id"]),
                cells=["c1"], genes=["g"], damage_floor=1.5,
            )


# --------------------------------------------------------------------------- #
# CascadeAnalyzer integration
# --------------------------------------------------------------------------- #

class TestCascadeIntegration:

    def test_none_weights_bit_identical_to_baseline(
        self, mapk_kg, healthy_expression
    ):
        analyzer = CascadeAnalyzer(mapk_kg)
        r_baseline = analyzer.analyze(healthy_expression, pathways=["MAPK"])
        r_none = analyzer.analyze(
            healthy_expression, pathways=["MAPK"], gene_weights=None
        )
        # Bit-identical score per pathway
        assert len(r_baseline.per_pathway) == len(r_none.per_pathway)
        for a, b in zip(r_baseline.per_pathway, r_none.per_pathway):
            assert a.upstream_score == b.upstream_score
            assert a.intermediate_score == b.intermediate_score
            assert a.downstream_score == b.downstream_score

    def test_carrier_downweighting_shifts_scores(
        self, mapk_kg, healthy_expression, score_table
    ):
        """Knocking down ERK via carriers should depress downstream score."""
        # Put a strong signal into ERK so zeroing its weight has effect
        expr = healthy_expression.copy()
        expr["ERK"] = 10.0 + 0.01 * np.arange(len(expr))  # strongly high

        scorer = AlphaMissenseScorer.from_table(pd.DataFrame({
            "variant_id": ["ERK:p.K71R"],
            "am_score":   [0.99],
        }))
        carriers = pd.DataFrame([
            {"cell_id": cid, "gene": "ERK", "variant_id": "ERK:p.K71R"}
            for cid in expr.index
        ])
        weights = scorer.weights_from_carriers(
            carriers=carriers, cells=expr.index, genes=expr.columns,
        )
        analyzer = CascadeAnalyzer(mapk_kg)
        r_naive = analyzer.analyze_pathway(expr, pathway="MAPK")
        r_weighted = analyzer.analyze_pathway(
            expr, pathway="MAPK", gene_weights=weights
        )
        # ERK is an intermediate gene in the MAPK cascade. Down-weighting
        # its contribution should decrease the intermediate-layer mean
        # (the elevated ERK signal no longer fully counts).
        assert r_weighted.intermediate_score < r_naive.intermediate_score

    def test_empty_carriers_no_op(
        self, mapk_kg, healthy_expression
    ):
        scorer = AlphaMissenseScorer.empty()
        weights = scorer.weights_from_carriers(
            carriers=pd.DataFrame(columns=["cell_id", "gene", "variant_id"]),
            cells=healthy_expression.index,
            genes=healthy_expression.columns,
        )
        analyzer = CascadeAnalyzer(mapk_kg)
        r_baseline = analyzer.analyze(healthy_expression, pathways=["MAPK"])
        r_weighted = analyzer.analyze(
            healthy_expression, pathways=["MAPK"], gene_weights=weights
        )
        for a, b in zip(r_baseline.per_pathway, r_weighted.per_pathway):
            assert a.upstream_score == pytest.approx(b.upstream_score)
            assert a.downstream_score == pytest.approx(b.downstream_score)
