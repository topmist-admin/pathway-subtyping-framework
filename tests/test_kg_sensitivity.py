"""Tests for the time-sliced KG sensitivity gate (Gate K).

Each verdict is exercised on a graph with known ground truth, using
partition functions whose KG-dependence is constructed rather than emergent,
so the expected verdict follows from the setup rather than from a fit.
"""

from __future__ import annotations

import numpy as np
import pytest

from pathway_subtyping.kg_sensitivity import (
    REWIRING_DEGREE,
    REWIRING_MODULE,
    REWIRING_UNIFORM,
    VERDICT_FRAGILE,
    VERDICT_KG_SENSITIVE,
    VERDICT_ROBUST,
    kg_timeslice_sensitivity,
    module_map_from_pathways,
    rewire_kg,
    within_module_fractions,
)
from pathway_subtyping.knowledge_graph.builder import KnowledgeGraph
from pathway_subtyping.knowledge_graph.schema import EdgeType, NodeType

N_SAMPLES = 40

# Two labelings with near-zero mutual agreement, so "the partition changed"
# is unambiguous rather than a matter of degree.
LABELS_A = np.array([i % 2 for i in range(N_SAMPLES)])
LABELS_B = np.array([0] * (N_SAMPLES // 2) + [1] * (N_SAMPLES // 2))

CRITICAL = ("G0", "G1")


def _chain_kg(n_genes: int = 101) -> KnowledgeGraph:
    """A gene chain with ``n_genes - 1`` GENE_INTERACTS edges.

    Sized so a single random removal has ~1% chance of hitting the critical
    edge; that is what separates 'kg-sensitive' from 'generically-fragile'
    with a 50-draw null.
    """
    kg = KnowledgeGraph()
    for i in range(n_genes):
        kg.add_node(f"G{i}", NodeType.GENE)
    for i in range(n_genes - 1):
        kg.add_edge(f"G{i}", f"G{i + 1}", EdgeType.GENE_INTERACTS)
    return kg


def _hub_kg(n_genes: int = 40) -> KnowledgeGraph:
    """A star-plus-chain graph with a pronounced degree skew.

    G0 is a hub wired to a third of the graph; the rest form a chain. Uniform
    rewiring flattens that skew, degree-preserving rewiring must not.
    """
    kg = KnowledgeGraph()
    for i in range(n_genes):
        kg.add_node(f"G{i}", NodeType.GENE)
    for i in range(1, n_genes // 3):
        kg.add_edge("G0", f"G{i}", EdgeType.GENE_INTERACTS)
    for i in range(n_genes // 3, n_genes - 1):
        kg.add_edge(f"G{i}", f"G{i + 1}", EdgeType.GENE_INTERACTS)
    return kg


def _pathway_kg(n_per: int = 10) -> KnowledgeGraph:
    """Two disjoint pathway modules, each an internally-wired gene clique-ish chain.

    Genes G0..G9 belong to P0, G10..G19 to P1, with GENE_IN_PATHWAY edges so
    ``module_map_from_pathways`` has something to read.
    """
    kg = KnowledgeGraph()
    kg.add_node("P0", NodeType.PATHWAY)
    kg.add_node("P1", NodeType.PATHWAY)
    for i in range(2 * n_per):
        kg.add_node(f"G{i}", NodeType.GENE)
        kg.add_edge(f"G{i}", "P0" if i < n_per else "P1", EdgeType.GENE_IN_PATHWAY)
    for block in (range(n_per - 1), range(n_per, 2 * n_per - 1)):
        for i in block:
            kg.add_edge(f"G{i}", f"G{i + 1}", EdgeType.GENE_INTERACTS)
    return kg


def _drop_edge(kg: KnowledgeGraph, src: str, tgt: str) -> KnowledgeGraph:
    """Return a copy of ``kg`` with every edge between src and tgt removed."""
    out = KnowledgeGraph(schema=kg.schema)
    for node_id, node_type in kg._node_types.items():
        out.add_node(node_id, node_type)
    for (s, t, key), edge_type in kg._edge_types.items():
        if (s, t) == (src, tgt):
            continue
        out.add_edge(s, t, edge_type, weight=float(kg.graph.edges[s, t, key].get("weight", 1.0)))
    return out


# --------------------------------------------------------------------------- #
# rewire_kg
# --------------------------------------------------------------------------- #


def test_rewire_preserves_nodes_and_edge_count_and_is_deterministic():
    kg = _chain_kg(20)
    budget_rm = {"gene_interacts_gene": 3}
    budget_add = {"gene_interacts_gene": 3}

    a = rewire_kg(kg, budget_rm, budget_add, np.random.default_rng(7))
    b = rewire_kg(kg, budget_rm, budget_add, np.random.default_rng(7))

    assert a.n_nodes == kg.n_nodes
    assert a.n_edges == kg.n_edges  # 3 out, 3 in
    assert kg.n_edges == 19, "input graph must not be mutated"

    edges_a = sorted((s, t) for (s, t, _) in a._edge_types)
    edges_b = sorted((s, t) for (s, t, _) in b._edge_types)
    assert edges_a == edges_b, "same seed must reproduce the same rewiring"


def test_rewire_actually_changes_something():
    kg = _chain_kg(20)
    p = rewire_kg(
        kg, {"gene_interacts_gene": 4}, {"gene_interacts_gene": 4}, np.random.default_rng(1)
    )
    assert sorted((s, t) for (s, t, _) in p._edge_types) != sorted(
        (s, t) for (s, t, _) in kg._edge_types
    )


def test_rewire_rejects_unknown_mode():
    kg = _chain_kg(20)
    with pytest.raises(ValueError, match="rewiring must be one of"):
        rewire_kg(kg, {}, {}, np.random.default_rng(0), rewiring="shuffle")


def test_module_mode_requires_a_module_map():
    kg = _chain_kg(20)
    with pytest.raises(ValueError, match="requires a non-empty module_map"):
        rewire_kg(
            kg,
            {"gene_interacts_gene": 1},
            {"gene_interacts_gene": 1},
            np.random.default_rng(0),
            rewiring=REWIRING_MODULE,
        )


# --------------------------------------------------------------------------- #
# Degree-preserving null
# --------------------------------------------------------------------------- #


def _degree_sequence(kg):
    return sorted((n, kg.graph.in_degree(n), kg.graph.out_degree(n)) for n in kg.graph.nodes())


def test_degree_mode_preserves_the_degree_sequence_exactly():
    """The whole point of the refinement: hubs must stay hubs."""
    kg = _hub_kg()
    before = _degree_sequence(kg)

    p = rewire_kg(
        kg,
        {"gene_interacts_gene": 5},
        {"gene_interacts_gene": 5},  # balanced -> fully swap-able
        np.random.default_rng(3),
        rewiring=REWIRING_DEGREE,
    )

    assert _degree_sequence(p) == before
    assert p.n_edges == kg.n_edges
    # ...and it did actually rewire, so preservation is not vacuous.
    assert sorted((s, t) for (s, t, _) in p._edge_types) != sorted(
        (s, t) for (s, t, _) in kg._edge_types
    )


def test_uniform_mode_does_not_preserve_degrees():
    """Contrast case — this is why 'degree' is the stronger null."""
    kg = _hub_kg()
    before = _degree_sequence(kg)
    p = rewire_kg(
        kg,
        {"gene_interacts_gene": 5},
        {"gene_interacts_gene": 5},
        np.random.default_rng(3),
        rewiring=REWIRING_UNIFORM,
    )
    assert _degree_sequence(p) != before


# --------------------------------------------------------------------------- #
# Module-aware null
# --------------------------------------------------------------------------- #


def test_module_map_derives_from_pathway_membership():
    kg = _pathway_kg()
    mm = module_map_from_pathways(kg)
    assert mm["G0"] == frozenset({"P0"})
    assert mm["G10"] == frozenset({"P1"})
    assert mm["G0"].isdisjoint(mm["G10"]), "different pathways are different modules"


def test_within_module_fractions_reads_the_observed_diff():
    from pathway_subtyping.knowledge_graph.diff import diff_kgs

    v1 = _pathway_kg()
    v2 = _drop_edge(v1, "G0", "G1")  # both in P0 -> a within-module removal
    mm = module_map_from_pathways(v1)
    frac = within_module_fractions(diff_kgs(v1, v2), mm)
    assert frac["gene_interacts_gene"] == pytest.approx(1.0)


def test_module_mode_runs_and_preserves_degrees():
    kg = _pathway_kg()
    before = _degree_sequence(kg)
    mm = module_map_from_pathways(kg)
    p = rewire_kg(
        kg,
        {"gene_interacts_gene": 2},
        {"gene_interacts_gene": 2},
        np.random.default_rng(5),
        rewiring=REWIRING_MODULE,
        module_map=mm,
        within_module_frac={"gene_interacts_gene": 1.0},
    )
    assert _degree_sequence(p) == before


def test_null_choice_can_change_the_verdict():
    """The refinement is load-bearing, not cosmetic.

    Identical data, identical observed ARI, opposite conclusions -- separated
    only by whether the null respects the degree sequence. The real change
    removes a *peripheral* edge; uniform rewiring hits peripheral edges often
    (so the disruption looks ordinary), while a degree-aware null hits them
    rarely (so the disruption looks specific).

    This is why ``degree`` is the default: ``uniform`` misreads changes at both
    ends of the degree distribution.
    """
    v1 = _hub_kg(60)
    peripheral = ("G45", "G46")
    v2 = _drop_edge(v1, *peripheral)

    def partition_fn(kg, cohort):
        return LABELS_A if kg.graph.has_edge(*peripheral) else LABELS_B

    uniform = kg_timeslice_sensitivity(
        v1, v2, partition_fn, None, n_null=60, seed=42, rewiring=REWIRING_UNIFORM
    )
    degree = kg_timeslice_sensitivity(
        v1, v2, partition_fn, None, n_null=60, seed=42, rewiring=REWIRING_DEGREE
    )

    assert uniform.observed_ari == pytest.approx(
        degree.observed_ari
    ), "the observed statistic is identical; only the reference differs"
    assert uniform.verdict == VERDICT_FRAGILE
    assert degree.verdict == VERDICT_KG_SENSITIVE
    assert degree.empirical_p < uniform.empirical_p
    assert degree.rewiring == REWIRING_DEGREE  # recorded on the result


def test_module_mode_abstains_without_pathway_edges():
    """A chain graph has no GENE_IN_PATHWAY edges, so modules are undefined."""
    v1 = _chain_kg(30)
    v2 = _drop_edge(v1, *CRITICAL)
    res = kg_timeslice_sensitivity(
        v1,
        v2,
        lambda kg, c: LABELS_A,
        None,
        n_null=5,
        seed=42,
        rewiring=REWIRING_MODULE,
    )
    assert not res.testable
    assert "module-aware" in res.verdict


# --------------------------------------------------------------------------- #
# Verdicts
# --------------------------------------------------------------------------- #


def test_robust_when_partition_ignores_the_graph():
    """A finding that does not depend on the KG must come back robust."""
    v1 = _chain_kg()
    v2 = _drop_edge(v1, *CRITICAL)

    res = kg_timeslice_sensitivity(v1, v2, lambda kg, cohort: LABELS_A, None, n_null=20, seed=42)

    assert res.testable
    assert res.verdict == VERDICT_ROBUST
    assert res.passed
    assert res.observed_ari == pytest.approx(1.0)


def test_kg_sensitive_when_the_specific_curated_edge_drives_the_partition():
    """Observed disruption >> size-matched random disruption."""
    v1 = _chain_kg()
    v2 = _drop_edge(v1, *CRITICAL)

    def partition_fn(kg, cohort):
        return LABELS_A if kg.has_edge(*CRITICAL) else LABELS_B

    res = kg_timeslice_sensitivity(v1, v2, partition_fn, None, n_null=50, seed=42)

    assert res.testable
    assert res.verdict == VERDICT_KG_SENSITIVE
    assert not res.passed
    # The real swap broke it; random swaps of equal size mostly did not.
    assert res.observed_ari < res.null_median
    assert res.empirical_p < res.alpha


def test_generically_fragile_when_any_perturbation_breaks_it():
    """Observed disruption ~ random disruption -> the partition is the problem."""
    v1 = _chain_kg()
    v2 = _drop_edge(v1, *CRITICAL)

    def partition_fn(kg, cohort):
        # Sensitive to the edge *count*, so every size-matched perturbation
        # disrupts it exactly as much as the real change does.
        return LABELS_A if kg.n_edges % 2 == 0 else LABELS_B

    res = kg_timeslice_sensitivity(v1, v2, partition_fn, None, n_null=30, seed=42)

    assert res.testable
    assert res.verdict == VERDICT_FRAGILE
    assert not res.passed
    assert res.empirical_p >= res.alpha


# --------------------------------------------------------------------------- #
# Abstentions -- must be distinguishable from rejections
# --------------------------------------------------------------------------- #


def test_abstains_when_graphs_are_identical():
    v1 = _chain_kg(20)
    res = kg_timeslice_sensitivity(
        v1, _chain_kg(20), lambda kg, c: LABELS_A, None, n_null=5, seed=42
    )
    assert not res.testable
    assert res.verdict.startswith("not-testable")
    assert not res.passed


def test_abstains_on_degenerate_partition():
    v1 = _chain_kg(20)
    v2 = _drop_edge(v1, *CRITICAL)
    single = np.zeros(N_SAMPLES, dtype=int)

    res = kg_timeslice_sensitivity(v1, v2, lambda kg, c: single, None, n_null=5, seed=42)
    assert not res.testable
    assert "degenerate" in res.verdict


def test_abstention_is_not_a_rejection():
    """The distinction Gate A got wrong: abstentions must not read as failures."""
    v1 = _chain_kg(20)
    res = kg_timeslice_sensitivity(
        v1, _chain_kg(20), lambda kg, c: LABELS_A, None, n_null=5, seed=42
    )
    # passed is False, but testable is ALSO False -- so any failure rate
    # computed over many runs can exclude this one from its denominator.
    assert res.passed is False
    assert res.testable is False
    assert np.isnan(res.observed_ari)


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #


def test_result_serializes_and_carries_the_diff():
    v1 = _chain_kg()
    v2 = _drop_edge(v1, *CRITICAL)
    res = kg_timeslice_sensitivity(v1, v2, lambda kg, c: LABELS_A, None, n_null=10, seed=42)
    d = res.to_dict()
    assert d["verdict"] == VERDICT_ROBUST
    assert d["diff_summary"]["edges_removed_total"] == 1
    assert d["n_null_draws"] == 10
    assert d["regression"] is None


def test_optional_scalar_regression_is_reported_but_does_not_decide():
    """A flagged scalar regression must not change the verdict."""
    v1 = _chain_kg()
    v2 = _drop_edge(v1, *CRITICAL)

    def score_fn(kg, cohort):
        return {"n_edges": float(kg.n_edges)}  # differs between v1 and v2

    res = kg_timeslice_sensitivity(
        v1,
        v2,
        lambda kg, c: LABELS_A,
        None,
        n_null=10,
        seed=42,
        score_fns=[score_fn],
        tolerance=0.0,
    )
    assert res.regression is not None
    assert not res.regression.passed, "scalar score did move"
    assert res.verdict == VERDICT_ROBUST, "but the partition verdict is unchanged"
    assert res.passed
