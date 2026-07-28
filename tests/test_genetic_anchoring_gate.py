"""
Tests for Gate 7 — the Genetic Anchoring Gate (feature-level).

Gate 7 is the positive-evidence counterpart to the negative controls: it tests
whether a subtype's *defining genes* are over-represented for disease
genetic-risk genes against a background-matched null. A germline-variant
enrichment cannot be manufactured by any postmortem/technical confound, so a pass
is confound-immune evidence the subtype axis is genetically implicated.

The reference scenario is Voineagu's discrimination reproduced from public data
(see the v0.7.0 reproduction bundle, ``gate7_feature_level_anchoring.py``): under
a brain-expressed-matched null, a NEURONAL/SYNAPTIC gene set is strongly enriched
for autism genetic risk (fold ~16x, p ~1e-23) while a GLIAL/IMMUNE set is not
(fold ~2.6x, n.s.). These tests reconstruct that discrimination synthetically and
deterministically — no network, no GEO files.
"""

import numpy as np
import pytest

from pathway_subtyping.genetics import (
    feature_level_anchoring,
    hypergeometric_enrichment,
)
from pathway_subtyping.validation import ValidationGates


class TestHypergeometricEnrichment:
    def test_strong_overlap_enriched(self):
        universe = {f"G{i}" for i in range(1000)}
        risk = {f"G{i}" for i in range(100)}  # 10% base rate
        test = {f"G{i}" for i in range(50)}  # entirely inside the risk block
        r = hypergeometric_enrichment(test, risk, universe, label="t", null="n")
        assert r.risk_hits == 50
        assert r.fold == pytest.approx(10.0)  # 50/50 hits vs 10% expected
        assert r.p_value < 1e-9

    def test_no_overlap_not_enriched(self):
        universe = {f"G{i}" for i in range(1000)}
        risk = {f"G{i}" for i in range(100)}
        test = {f"G{i}" for i in range(500, 560)}  # disjoint from the risk block
        r = hypergeometric_enrichment(test, risk, universe)
        assert r.risk_hits == 0
        assert r.fold == pytest.approx(0.0)
        assert r.p_value == pytest.approx(1.0)

    def test_only_universe_members_count(self):
        # Genes outside the universe are ignored on both the test and risk side.
        universe = {"A", "B", "C", "D"}
        r = hypergeometric_enrichment({"A", "B", "OUTSIDE"}, {"A", "OUTSIDE2"}, universe)
        assert r.testset_n == 2  # OUTSIDE dropped
        assert r.risk_in_universe == 1  # OUTSIDE2 dropped
        assert r.risk_hits == 1

    def test_undefined_test_returns_nan(self):
        universe = {f"G{i}" for i in range(100)}
        r = hypergeometric_enrichment(set(), {"G0"}, universe)
        assert r.testset_n == 0
        assert np.isnan(r.p_value)
        assert np.isnan(r.fold)


class TestGeneticAnchoringVoineaguLike:
    """Synthetic reconstruction of the neuronal-enriched / glial-null result."""

    def _make_scenario(self):
        # Brain-expressed-matched universe of 2000 genes; 200 carry ASD risk
        # (10% base rate). The neuronal set is drawn mostly from the risk block;
        # the glial set is drawn entirely from outside it.
        universe = {f"G{i}" for i in range(2000)}
        risk = {f"G{i}" for i in range(200)}
        neuronal = {f"G{i}" for i in range(40)} | {f"G{i}" for i in range(1000, 1020)}
        glial = {f"G{i}" for i in range(1500, 1560)}
        # A genome-wide reference (larger, dilutes the base rate) for contrast.
        genome_wide = {f"G{i}" for i in range(6000)}
        return universe, risk, neuronal, glial, genome_wide

    def test_discrimination_neuronal_enriched_glial_not(self):
        universe, risk, neuronal, glial, _ = self._make_scenario()
        res = feature_level_anchoring({"NEURONAL": neuronal, "GLIAL": glial}, risk, universe)
        neur = res["NEURONAL|background-matched"]
        gli = res["GLIAL|background-matched"]

        assert neur.fold > 3.0 and neur.p_value < 1e-3  # neuronal enriched
        assert gli.fold < 2.0 and gli.p_value > 0.05  # glial not

    def test_null_matters_matched_vs_genomewide(self):
        # Reproduces the "null matters" point: with a background-matched null the
        # base rate is higher, so a truly enriched set shows a *higher* fold than
        # against a diluted genome-wide null.
        universe, risk, neuronal, _, genome_wide = self._make_scenario()
        res = feature_level_anchoring(
            {"NEURONAL": neuronal}, risk, universe, reference_universe=genome_wide
        )
        matched = res["NEURONAL|background-matched"]
        wide = res["NEURONAL|genome-wide"]
        # Same numerator (risk hits), lower background base rate genome-wide ->
        # matched fold is the smaller number, both still enriched.
        assert matched.risk_hits == wide.risk_hits
        assert matched.fold < wide.fold
        assert matched.p_value < 1e-3

    def test_gate_passes_with_neuronal_axis_anchored(self):
        universe, risk, neuronal, glial, genome_wide = self._make_scenario()
        gates = ValidationGates(show_progress=False)
        result = gates.genetic_anchoring_gate(
            subtype_gene_sets={"NEURONAL": neuronal, "GLIAL": glial},
            risk_genes=risk,
            background_universe=universe,
            reference_universe=genome_wide,
        )
        assert result.passed is True
        assert result.details["anchored_subtypes"] == ["NEURONAL"]
        assert result.details["best_subtype"] == "NEURONAL"
        assert result.metric_value > 3.0
        assert result.details["gate_polarity"] == "positive_evidence"

        per = result.details["per_subtype"]
        assert per["NEURONAL"]["anchored"] is True
        assert per["GLIAL"]["anchored"] is False
        # Reference null reported per subtype but does not decide the gate.
        assert "reference_null" in per["NEURONAL"]
        assert per["NEURONAL"]["reference_null"]["null"] == "genome-wide"

    def test_gate_fails_when_no_axis_anchored(self):
        # A partition whose defining genes carry no genetic risk fails Gate 7.
        universe = {f"G{i}" for i in range(2000)}
        risk = {f"G{i}" for i in range(200)}
        subtype_a = {f"G{i}" for i in range(500, 560)}  # disjoint from risk
        subtype_b = {f"G{i}" for i in range(600, 660)}
        gates = ValidationGates(show_progress=False)
        result = gates.genetic_anchoring_gate(
            subtype_gene_sets={"A": subtype_a, "B": subtype_b},
            risk_genes=risk,
            background_universe=universe,
        )
        assert result.passed is False
        assert result.details["anchored_subtypes"] == []
        assert "low-power" in result.details["interpretation"]


class TestGeneticAnchoringInRunAll:
    def test_run_all_includes_gate_when_anchoring_present(self):
        import pandas as pd

        rng = np.random.default_rng(20260716)
        n = 60
        pathway_scores = pd.DataFrame(rng.normal(size=(n, 8)), columns=[f"PW{i}" for i in range(8)])
        gene_burdens = pd.DataFrame(rng.normal(size=(n, 12)), columns=[f"G{i}" for i in range(12)])
        pathways = {f"PW{i}": [f"G{i}", f"G{(i + 1) % 12}"] for i in range(8)}
        clusters = rng.integers(0, 2, n)

        universe = {f"GENE{i}" for i in range(2000)}
        risk = {f"GENE{i}" for i in range(200)}
        anchored_set = {f"GENE{i}" for i in range(40)} | {f"GENE{i}" for i in range(1000, 1020)}
        null_set = {f"GENE{i}" for i in range(1500, 1560)}

        gates = ValidationGates(
            seed=20260716, n_permutations=10, n_bootstrap=10, show_progress=False
        )
        res = gates.run_all(
            pathway_scores=pathway_scores,
            cluster_labels=clusters,
            pathways=pathways,
            gene_burdens=gene_burdens,
            n_clusters=2,
            gmm_seed=20260716,
            genetic_anchoring={
                "subtype_gene_sets": {0: anchored_set, 1: null_set},
                "risk_genes": risk,
                "background_universe": universe,
            },
        )
        names = [r.name for r in res.results]
        assert "Gate 7: Genetic Anchoring" in names
        gate7 = next(r for r in res.results if r.name == "Gate 7: Genetic Anchoring")
        assert gate7.passed is True
        assert gate7.details["anchored_subtypes"] == ["0"]
