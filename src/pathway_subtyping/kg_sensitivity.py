"""
Time-sliced knowledge-graph sensitivity (Gate K).

Asks one question about a finding that was derived with help from a knowledge
graph: **does it survive being recomputed against a different KG version?**

A subtype that appears under Reactome 2026 and vanishes under Reactome 2025 is a
knowledge artifact, not a property of the cohort. Nothing in the existing gate
battery checks for that, because every other gate holds the KG fixed and varies
the data.

WHY A RAW AGREEMENT NUMBER IS NOT ENOUGH
----------------------------------------
The obvious implementation is: partition the cohort under KG v1, partition it
under KG v2, report the ARI between the two labelings. That number is
uninterpretable on its own, for exactly the reason the pre-v0.8.0 stability gate
was uninterpretable -- there is no reference telling you what value to expect.

An ARI of 0.6 between two KG versions could mean either:

  (a) the curated change was substantive and specifically disrupted this finding; or
  (b) this partition is fragile to *any* perturbation of comparable size, and the
      KG version is incidental.

These have opposite implications. (a) is a statement about the knowledge base;
(b) is a statement about the partition, and a much worse one. Distinguishing them
requires a null, so this module builds one: rewire the baseline KG at random,
matching the observed diff's edge-addition and edge-removal counts, and measure
how much the partition moves under that size-matched random perturbation.

  observed agreement >> null  ->  robust
  observed agreement ~  null  ->  generically fragile (any perturbation breaks it)
  observed agreement << null  ->  KG-sensitive (this specific curated change matters)

Reuses ``diff_kgs`` to size the perturbation and, optionally, ``run_kg_regression``
for scalar score deltas alongside the partition-level verdict.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .knowledge_graph.builder import KnowledgeGraph
from .knowledge_graph.diff import KGDiff, diff_kgs
from .knowledge_graph.regression import KGRegressionReport, ScoreFn, run_kg_regression
from .knowledge_graph.schema import EdgeType
from .utils.metrics import safe_adjusted_rand_score

logger = logging.getLogger(__name__)


#: A caller-supplied function that turns (knowledge graph, cohort payload) into
#: cluster labels. Mirrors ``regression.ScoreFn``, but returns a labeling rather
#: than a dict of scalars.
PartitionFn = Callable[[KnowledgeGraph, Any], np.ndarray]

VERDICT_ROBUST = "robust"
VERDICT_KG_SENSITIVE = "kg-sensitive"
VERDICT_FRAGILE = "generically-fragile"

#: Null models, weakest to strongest. See :func:`rewire_kg`.
REWIRING_UNIFORM = "uniform"
REWIRING_DEGREE = "degree"
REWIRING_MODULE = "module"


@dataclass
class KGSensitivityResult:
    """Outcome of one time-sliced sensitivity run.

    Attributes:
        verdict: One of ``robust``, ``kg-sensitive``, ``generically-fragile``, or
            a ``not-testable (...)`` string. Four outcomes, not two -- an
            abstention is not a failure and must never be counted as one.
        passed: True only for ``robust``. Abstentions are False but are NOT
            rejections; see ``testable``.
        testable: False when the run abstained. Any false-positive or
            failure rate computed over many runs MUST use the testable subset as
            its denominator.
        observed_ari: Agreement between the v1 and v2 partitions.
        null_aris: Agreement between the v1 partition and each size-matched
            randomly-rewired partition.
        null_p05: 5th percentile of ``null_aris`` (the low tail is the one that
            matters -- we ask whether the real change disrupted *more* than chance).
        empirical_p: Add-one empirical p that the observed agreement is at or
            below the null distribution.
        diff: The structural ``KGDiff`` between the two graphs.
        regression: Optional scalar-score report from ``run_kg_regression``.
    """

    verdict: str
    testable: bool
    observed_ari: float = float("nan")
    null_aris: List[float] = field(default_factory=list)
    null_median: float = float("nan")
    null_p05: float = float("nan")
    empirical_p: float = float("nan")
    n_clusters_v1: int = 0
    n_clusters_v2: int = 0
    ari_min: float = 0.0
    alpha: float = 0.05
    rewiring: str = REWIRING_DEGREE
    diff: Optional[KGDiff] = None
    regression: Optional[KGRegressionReport] = None

    @property
    def passed(self) -> bool:
        return self.verdict == VERDICT_ROBUST

    def summary(self) -> str:
        if not self.testable:
            return f"KGSensitivity({self.verdict})"
        return (
            f"KGSensitivity({self.verdict}, observed_ari={self.observed_ari:.3f}, "
            f"null_median={self.null_median:.3f}, p={self.empirical_p:.3f})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "verdict": self.verdict,
            "passed": self.passed,
            "testable": self.testable,
            "observed_ari": self.observed_ari,
            "null_median": self.null_median,
            "null_p05": self.null_p05,
            "empirical_p": self.empirical_p,
            "n_null_draws": len(self.null_aris),
            "n_clusters_v1": self.n_clusters_v1,
            "n_clusters_v2": self.n_clusters_v2,
            "ari_min": self.ari_min,
            "alpha": self.alpha,
            "rewiring": self.rewiring,
            "diff_summary": dict(self.diff.summary) if self.diff else None,
            "regression": self.regression.to_dict() if self.regression else None,
        }


# --------------------------------------------------------------------------- #
# Size-matched random rewiring (the null)
# --------------------------------------------------------------------------- #

def _edges_by_type(kg: KnowledgeGraph) -> Dict[str, List[Tuple[str, str]]]:
    out: Dict[str, List[Tuple[str, str]]] = {}
    for (src, tgt, _key), edge_type in kg._edge_types.items():
        out.setdefault(edge_type.value, []).append((src, tgt))
    return out


def _edge_key_index(kg: KnowledgeGraph) -> Dict[Tuple[str, str, str], List[int]]:
    """{(src, tgt, edge_type_value): [multigraph keys]} for O(1) typed removal."""
    idx: Dict[Tuple[str, str, str], List[int]] = {}
    for (src, tgt, key), edge_type in kg._edge_types.items():
        idx.setdefault((src, tgt, edge_type.value), []).append(key)
    return idx


def _drop_typed_edge(kg: KnowledgeGraph, src: str, tgt: str, et_value: str) -> None:
    """Remove every edge of ``et_value`` between src and tgt, in place."""
    keys = [
        k for (s, t, k), et in list(kg._edge_types.items())
        if s == src and t == tgt and et.value == et_value
    ]
    for k in keys:
        if kg.graph.has_edge(src, tgt, k):
            kg.graph.remove_edge(src, tgt, k)
        kg._edge_types.pop((src, tgt, k), None)


def module_map_from_pathways(kg: KnowledgeGraph) -> Dict[str, frozenset]:
    """Derive a module assignment from ``GENE_IN_PATHWAY`` membership.

    Two genes count as same-module when they share at least one pathway. This is
    the KG's own native notion of a module, so module-aware rewiring needs no
    external community detection and no caller-supplied partition.

    Returns:
        {gene_id: frozenset(pathway_ids)}. Genes in no pathway map to an empty
        frozenset and are treated as sharing a module with nobody.
    """
    membership: Dict[str, set] = {}
    for (src, tgt, _key), edge_type in kg._edge_types.items():
        if edge_type is EdgeType.GENE_IN_PATHWAY:
            membership.setdefault(src, set()).add(tgt)
    return {gene: frozenset(paths) for gene, paths in membership.items()}


def _same_module(module_map: Dict[str, frozenset], a: str, b: str) -> bool:
    ma = module_map.get(a)
    mb = module_map.get(b)
    if not ma or not mb:
        return False
    return not ma.isdisjoint(mb)


def _degree_preserving_swaps(
    kg: KnowledgeGraph,
    et_value: str,
    n_swaps: int,
    rng: np.random.Generator,
    module_map: Optional[Dict[str, frozenset]] = None,
    within_module: Optional[bool] = None,
) -> int:
    """Maslov–Sneppen double-edge swaps on one edge type, in place.

    Takes two edges ``a->b`` and ``c->d`` and rewrites them as ``a->d`` and
    ``c->b``. Every node keeps its exact in-degree and out-degree, so the null
    perturbs topology without touching the degree sequence — which is what stops
    hub genes from being stripped or peripheral genes from being inflated.

    Each swap changes two edges (two removed, two added), so it consumes one unit
    of the removal budget and one of the addition budget.

    Args:
        module_map: When supplied with ``within_module``, only accept swaps whose
            *resulting* edges are within-module (True) or cross-module (False).

    Returns:
        Number of swaps actually performed. May be fewer than requested on dense
        or highly constrained graphs; the caller logs the shortfall.
    """
    performed = 0
    for _ in range(max(1, int(n_swaps)) * 40):
        if performed >= int(n_swaps):
            break
        present = _edges_by_type(kg).get(et_value, [])
        if len(present) < 2:
            break
        i, j = rng.choice(len(present), size=2, replace=False)
        a, b = present[int(i)]
        c, d = present[int(j)]
        if len({a, b, c, d}) < 4:
            continue
        if kg.graph.has_edge(a, d) or kg.graph.has_edge(c, b):
            continue
        if module_map is not None and within_module is not None:
            keep = (
                _same_module(module_map, a, d) == within_module
                and _same_module(module_map, c, b) == within_module
            )
            if not keep:
                continue
        try:
            edge_type = EdgeType(et_value)
        except ValueError:
            break
        try:
            kg.add_edge(a, d, edge_type, weight=1.0)
            kg.add_edge(c, b, edge_type, weight=1.0)
        except ValueError:
            _drop_typed_edge(kg, a, d, et_value)
            _drop_typed_edge(kg, c, b, et_value)
            continue
        _drop_typed_edge(kg, a, b, et_value)
        _drop_typed_edge(kg, c, d, et_value)
        performed += 1
    return performed


def _degree_weighted_choice(
    kg: KnowledgeGraph,
    candidates: Sequence[str],
    rng: np.random.Generator,
) -> str:
    """Sample a node with probability proportional to its total degree.

    Used for the residual when a diff is unbalanced (more additions than
    removals, or vice versa). Strict degree preservation is impossible there —
    adding an edge necessarily raises somebody's degree — so the residual instead
    preserves the *shape* of the degree distribution by attaching preferentially
    to already-connected nodes, which is how curated releases actually grow.
    """
    degrees = np.array([kg.graph.degree(n) for n in candidates], dtype=float)
    total = degrees.sum()
    if total <= 0:
        return candidates[int(rng.integers(len(candidates)))]
    return candidates[int(rng.choice(len(candidates), p=degrees / total))]


def _clone_kg(kg: KnowledgeGraph) -> KnowledgeGraph:
    """Structural copy preserving node types, edge types and weights."""
    clone = KnowledgeGraph(schema=kg.schema)
    for node_id, node_type in kg._node_types.items():
        node = kg.get_node(node_id)
        clone.add_node(
            node_id,
            node_type,
            attributes=dict(node.attributes) if node is not None else None,
        )
    for (src, tgt, key), edge_type in kg._edge_types.items():
        data = kg.graph.edges[src, tgt, key]
        clone.add_edge(src, tgt, edge_type, weight=float(data.get("weight", 1.0)))
    return clone


def rewire_kg(
    kg: KnowledgeGraph,
    n_remove_by_type: Dict[str, int],
    n_add_by_type: Dict[str, int],
    rng: np.random.Generator,
    *,
    rewiring: str = REWIRING_DEGREE,
    module_map: Optional[Dict[str, frozenset]] = None,
    within_module_frac: Optional[Dict[str, float]] = None,
) -> KnowledgeGraph:
    """Perturb ``kg`` by a per-edge-type budget of removals and additions.

    The perturbation is *size-matched* to an observed diff: it changes the same
    number of edges of the same types, but chooses which edges at random. That
    is what makes it a reference for "how much would a curated change of this
    magnitude move the partition if it carried no information?"

    Node set and edge-type composition are always preserved. What varies is how
    much *else* is preserved, and that choice sets how strict the null is:

    ``uniform``
        Remove random edges; attach additions to uniformly-chosen endpoints. The
        weakest null, and **anti-conservative for the ``kg-sensitive`` verdict**:
        it destroys the degree sequence, so it disrupts the partition more than a
        real release would, making the real change look comparatively benign.

    ``degree`` (default)
        Maslov–Sneppen double-edge swaps for the balanced part of the budget, so
        every node keeps its exact in- and out-degree. Hubs stay hubs. Any
        residual (when additions and removals are unequal, as in a release that
        mostly grows) is degree-weighted, preserving the shape of the degree
        distribution but not its exact values — strict preservation is impossible
        there, because adding an edge necessarily raises somebody's degree.

    ``module``
        Degree-preserving swaps additionally constrained to reproduce the
        observed diff's within-module / cross-module split. Modules come from
        ``GENE_IN_PATHWAY`` membership via :func:`module_map_from_pathways`, so a
        release that concentrated its edits inside a few actively-researched
        pathways is matched by a null that does the same.

    Args:
        kg: Baseline graph to perturb.
        n_remove_by_type: {edge_type_value: count} edges to delete.
        n_add_by_type: {edge_type_value: count} edges to introduce.
        rng: Seeded generator; the same seed reproduces the same rewiring.
        rewiring: One of ``uniform``, ``degree``, ``module``.
        module_map: Required for ``module``; see :func:`module_map_from_pathways`.
        within_module_frac: {edge_type_value: fraction in [0, 1]} of changes to
            place within-module. Required for ``module``.

    Returns:
        A new ``KnowledgeGraph``. The input is not modified.

    Raises:
        ValueError: on an unknown ``rewiring`` mode, or ``module`` without a
            ``module_map``.
    """
    if rewiring not in (REWIRING_UNIFORM, REWIRING_DEGREE, REWIRING_MODULE):
        raise ValueError(
            f"rewiring must be one of 'uniform', 'degree', 'module'; got {rewiring!r}"
        )
    if rewiring == REWIRING_MODULE and not module_map:
        raise ValueError("rewiring='module' requires a non-empty module_map")

    clone = _clone_kg(kg)

    for et_value in sorted(set(n_remove_by_type) | set(n_add_by_type)):
        n_remove = int(n_remove_by_type.get(et_value, 0))
        n_add = int(n_add_by_type.get(et_value, 0))
        if n_remove <= 0 and n_add <= 0:
            continue

        if rewiring in (REWIRING_DEGREE, REWIRING_MODULE):
            # Each swap consumes one removal AND one addition, so the balanced
            # part of the budget is what can be done degree-preservingly.
            n_paired = min(n_remove, n_add)
            if n_paired > 0:
                if rewiring == REWIRING_MODULE:
                    frac = float((within_module_frac or {}).get(et_value, 0.0))
                    n_within = int(round(n_paired * max(0.0, min(1.0, frac))))
                    done_w = _degree_preserving_swaps(
                        clone, et_value, n_within, rng, module_map, True
                    )
                    done_c = _degree_preserving_swaps(
                        clone, et_value, n_paired - n_within, rng, module_map, False
                    )
                    done = done_w + done_c
                else:
                    done = _degree_preserving_swaps(clone, et_value, n_paired, rng)
                if done < n_paired:
                    logger.debug(
                        "[GateK] %s: %d/%d degree-preserving swaps on %s "
                        "(graph too constrained); residual falls back to weighted",
                        rewiring, done, n_paired, et_value,
                    )
                n_remove -= done
                n_add -= done

        # --- residual removals -------------------------------------------- #
        if n_remove > 0:
            present = _edges_by_type(clone).get(et_value, [])
            if present:
                n = min(n_remove, len(present))
                if rewiring == REWIRING_UNIFORM:
                    idx = rng.choice(len(present), size=n, replace=False)
                    chosen = [present[int(i)] for i in idx]
                else:
                    # Weight by endpoint degree so hubs shed edges in proportion
                    # to how many they have, rather than uniformly.
                    w = np.array(
                        [clone.graph.degree(s) + clone.graph.degree(t) for s, t in present],
                        dtype=float,
                    )
                    p = w / w.sum() if w.sum() > 0 else None
                    idx = rng.choice(len(present), size=n, replace=False, p=p)
                    chosen = [present[int(i)] for i in idx]
                for src, tgt in chosen:
                    _drop_typed_edge(clone, src, tgt, et_value)

        # --- residual additions -------------------------------------------- #
        if n_add > 0:
            present = _edges_by_type(clone).get(et_value, [])
            if not present:
                # No exemplar of this edge type -> no way to infer a
                # schema-valid endpoint pair. Skip rather than guess.
                logger.debug("[GateK] no exemplar for %s; skipping additions", et_value)
                continue
            try:
                edge_type = EdgeType(et_value)
            except ValueError:
                continue
            sources = sorted({s for s, _ in present})
            targets = sorted({t for _, t in present})
            added = 0
            # Bounded attempts: a dense graph may have few free slots, and we
            # would rather under-perturb (and log it) than spin.
            for _ in range(n_add * 20):
                if added >= n_add:
                    break
                if rewiring == REWIRING_UNIFORM:
                    src = sources[int(rng.integers(len(sources)))]
                    tgt = targets[int(rng.integers(len(targets)))]
                else:
                    src = _degree_weighted_choice(clone, sources, rng)
                    tgt = _degree_weighted_choice(clone, targets, rng)
                if src == tgt or clone.graph.has_edge(src, tgt):
                    continue
                if rewiring == REWIRING_MODULE and module_map is not None:
                    want = rng.random() < float(
                        (within_module_frac or {}).get(et_value, 0.0)
                    )
                    if _same_module(module_map, src, tgt) != want:
                        continue
                try:
                    clone.add_edge(src, tgt, edge_type, weight=1.0)
                except ValueError:
                    continue  # schema rejected this endpoint pair
                added += 1
            if added < n_add:
                logger.debug(
                    "[GateK] added %d/%d edges of type %s (graph too dense)",
                    added, n_add, et_value,
                )

    return clone


def within_module_fractions(
    diff: KGDiff, module_map: Dict[str, frozenset]
) -> Dict[str, float]:
    """Per-edge-type fraction of an observed diff that is within-module.

    This is what ``rewiring='module'`` matches. A release that concentrated 80%
    of its new edges inside shared pathways yields ``{edge_type: 0.8}``, and the
    null then places 80% of its own changes within-module — so a `kg-sensitive`
    verdict cannot be explained away by "the real change was concentrated."
    """
    out: Dict[str, float] = {}
    for et_value in set(diff.edges_added) | set(diff.edges_removed):
        pairs = list(diff.edges_added.get(et_value, [])) + list(
            diff.edges_removed.get(et_value, [])
        )
        if not pairs:
            continue
        within = sum(1 for s, t in pairs if _same_module(module_map, s, t))
        out[et_value] = within / len(pairs)
    return out


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def kg_timeslice_sensitivity(
    v1: KnowledgeGraph,
    v2: KnowledgeGraph,
    partition_fn: PartitionFn,
    cohort: Any,
    *,
    n_null: int = 50,
    ari_min: float = 0.80,
    alpha: float = 0.05,
    seed: int = 42,
    rewiring: str = REWIRING_DEGREE,
    score_fns: Optional[Iterable[ScoreFn]] = None,
    tolerance: float = 0.05,
) -> KGSensitivityResult:
    """Test whether a partition survives swapping the knowledge graph.

    THE DECISION RULE, in full::

        if observed_ari >= ari_min:              verdict = "robust"
        elif empirical_p < alpha:                verdict = "kg-sensitive"
        else:                                    verdict = "generically-fragile"

    Both terms are live: ``ari_min`` decides robustness and ``empirical_p``
    decides how a non-robust result is *explained*. Nothing else in this module
    enters the rule -- the ``KGDiff`` and the optional ``run_kg_regression``
    report are context, not criteria. (Gate A shipped documenting three criteria
    of which only one decided; that gap was found by hostile review. If a future
    revision wants the regression report to decide, it must change this rule in
    the same commit as the docstring.)

    Abstains, rather than ruling, when: the two graphs are identical (nothing to
    test), either partition has fewer than two clusters, or the two partitions
    cover different numbers of samples.

    Args:
        v1: Baseline knowledge graph (e.g. Reactome 2025).
        v2: Comparison knowledge graph (e.g. Reactome 2026).
        partition_fn: ``f(kg, cohort) -> labels``. Must be deterministic given
            its inputs; seed any internal clustering yourself.
        cohort: Opaque payload handed to ``partition_fn`` (expression matrix,
            pathway-score frame, whatever the caller's pipeline consumes).
        n_null: Size-matched random rewirings to draw.
        ari_min: Agreement at or above which the finding is called robust.
        alpha: Significance level for the empirical p.
        seed: Seeds the rewiring generator; the run is reproducible.
        score_fns: Optional scalar scoring callables, forwarded verbatim to
            ``run_kg_regression`` for a side-by-side scalar view.
        tolerance: Relative-delta tolerance for that regression report.

    Returns:
        A ``KGSensitivityResult``.
    """
    diff = diff_kgs(v1, v2)
    if diff.is_identical():
        return KGSensitivityResult(
            verdict="not-testable (knowledge graphs are identical)",
            testable=False,
            diff=diff,
            ari_min=ari_min,
            alpha=alpha,
            rewiring=rewiring,
        )

    labels_v1 = np.asarray(partition_fn(v1, cohort))
    labels_v2 = np.asarray(partition_fn(v2, cohort))

    if labels_v1.shape[0] != labels_v2.shape[0]:
        return KGSensitivityResult(
            verdict="not-testable (partitions cover different sample counts)",
            testable=False,
            diff=diff,
            ari_min=ari_min,
            alpha=alpha,
            rewiring=rewiring,
        )

    k1, k2 = len(np.unique(labels_v1)), len(np.unique(labels_v2))
    if k1 < 2 or k2 < 2:
        return KGSensitivityResult(
            verdict="not-testable (degenerate partition under at least one KG)",
            testable=False,
            n_clusters_v1=k1,
            n_clusters_v2=k2,
            diff=diff,
            ari_min=ari_min,
            alpha=alpha,
            rewiring=rewiring,
        )

    observed_ari = float(safe_adjusted_rand_score(labels_v1, labels_v2))
    if not np.isfinite(observed_ari):
        return KGSensitivityResult(
            verdict="not-testable (ARI undefined on these partitions)",
            testable=False,
            n_clusters_v1=k1,
            n_clusters_v2=k2,
            diff=diff,
            ari_min=ari_min,
            alpha=alpha,
            rewiring=rewiring,
        )

    # Size-match the perturbation to the real diff, per edge type.
    n_remove_by_type = {k: len(v) for k, v in diff.edges_removed.items()}
    n_add_by_type = {k: len(v) for k, v in diff.edges_added.items()}

    module_map: Optional[Dict[str, frozenset]] = None
    frac: Optional[Dict[str, float]] = None
    if rewiring == REWIRING_MODULE:
        module_map = module_map_from_pathways(v1)
        if not module_map:
            return KGSensitivityResult(
                verdict="not-testable (no GENE_IN_PATHWAY edges for module-aware null)",
                testable=False,
                observed_ari=observed_ari,
                n_clusters_v1=k1,
                n_clusters_v2=k2,
                diff=diff,
                ari_min=ari_min,
                alpha=alpha,
                rewiring=rewiring,
            )
        frac = within_module_fractions(diff, module_map)

    rng = np.random.default_rng(seed)
    null_aris: List[float] = []
    for _ in range(int(n_null)):
        perturbed = rewire_kg(
            v1, n_remove_by_type, n_add_by_type, rng,
            rewiring=rewiring, module_map=module_map, within_module_frac=frac,
        )
        labels_p = np.asarray(partition_fn(perturbed, cohort))
        if labels_p.shape[0] != labels_v1.shape[0]:
            continue
        value = float(safe_adjusted_rand_score(labels_v1, labels_p))
        if np.isfinite(value):
            null_aris.append(value)

    if len(null_aris) < 2:
        return KGSensitivityResult(
            verdict="not-testable (null distribution could not be constructed)",
            testable=False,
            observed_ari=observed_ari,
            n_clusters_v1=k1,
            n_clusters_v2=k2,
            diff=diff,
            ari_min=ari_min,
            alpha=alpha,
            rewiring=rewiring,
        )

    null_arr = np.asarray(null_aris, dtype=float)
    # Add-one empirical p in the LOW tail: how often does a size-matched random
    # rewiring disrupt the partition at least as much as the real KG change?
    empirical_p = float((np.sum(null_arr <= observed_ari) + 1) / (null_arr.size + 1))

    if observed_ari >= ari_min:
        verdict = VERDICT_ROBUST
    elif empirical_p < alpha:
        verdict = VERDICT_KG_SENSITIVE
    else:
        verdict = VERDICT_FRAGILE

    regression: Optional[KGRegressionReport] = None
    if score_fns is not None:
        score_fns = list(score_fns)
        if score_fns:
            regression = run_kg_regression(
                v1, v2, score_fns, [cohort], tolerance=tolerance
            )

    result = KGSensitivityResult(
        verdict=verdict,
        testable=True,
        observed_ari=observed_ari,
        null_aris=null_aris,
        null_median=float(np.median(null_arr)),
        null_p05=float(np.percentile(null_arr, 5)),
        empirical_p=empirical_p,
        n_clusters_v1=k1,
        n_clusters_v2=k2,
        ari_min=ari_min,
        alpha=alpha,
        rewiring=rewiring,
        diff=diff,
        regression=regression,
    )
    logger.info("[GateK] %s", result.summary())
    return result


__all__ = [
    "KGSensitivityResult",
    "PartitionFn",
    "kg_timeslice_sensitivity",
    "rewire_kg",
    "module_map_from_pathways",
    "within_module_fractions",
    "REWIRING_UNIFORM",
    "REWIRING_DEGREE",
    "REWIRING_MODULE",
    "VERDICT_ROBUST",
    "VERDICT_KG_SENSITIVE",
    "VERDICT_FRAGILE",
]
