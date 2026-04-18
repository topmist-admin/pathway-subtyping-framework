"""
Knowledge-graph diff utility for v0.6 F3.

Compares two ``KnowledgeGraph`` instances and reports what changed:
node additions/removals by type, edge additions/removals by type, and
edges whose direction flipped between the two snapshots. Used to audit
KG upgrades between releases (e.g. OmniPath 2024 vs 2025) before any
code depending on KG topology is re-run.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Set, Tuple

from .builder import KnowledgeGraph
from .schema import EdgeType, NodeType

logger = logging.getLogger(__name__)


EdgeTuple = Tuple[str, str, str]  # (source, target, edge_type_value)


@dataclass
class KGDiff:
    """Structural diff between two ``KnowledgeGraph`` instances.

    All lists are sorted for deterministic output.

    Attributes:
        nodes_added: {node_type_value: [node_id, ...]} present in v2, absent in v1.
        nodes_removed: {node_type_value: [node_id, ...]} present in v1, absent in v2.
        edges_added: {edge_type_value: [(src, tgt), ...]} present in v2, absent in v1.
        edges_removed: {edge_type_value: [(src, tgt), ...]} present in v1, absent in v2.
        edges_flipped: [(src, tgt, edge_type)] whose direction reversed.
        summary: counts for quick logging.
    """

    nodes_added: Dict[str, List[str]] = field(default_factory=dict)
    nodes_removed: Dict[str, List[str]] = field(default_factory=dict)
    edges_added: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)
    edges_removed: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)
    edges_flipped: List[Tuple[str, str, str]] = field(default_factory=list)
    summary: Dict[str, int] = field(default_factory=dict)

    def is_identical(self) -> bool:
        return (
            not self.nodes_added
            and not self.nodes_removed
            and not self.edges_added
            and not self.edges_removed
            and not self.edges_flipped
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "nodes_added": {k: sorted(v) for k, v in self.nodes_added.items()},
            "nodes_removed": {k: sorted(v) for k, v in self.nodes_removed.items()},
            "edges_added": {
                k: sorted(list(v)) for k, v in self.edges_added.items()
            },
            "edges_removed": {
                k: sorted(list(v)) for k, v in self.edges_removed.items()
            },
            "edges_flipped": sorted(self.edges_flipped),
            "summary": dict(self.summary),
        }


# --------------------------------------------------------------------------- #
# Internals
# --------------------------------------------------------------------------- #

def _nodes_by_type(kg: KnowledgeGraph) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    for node_id, node_type in kg._node_types.items():
        out.setdefault(node_type.value, set()).add(node_id)
    return out


def _edges_by_type(kg: KnowledgeGraph) -> Dict[str, Set[Tuple[str, str]]]:
    out: Dict[str, Set[Tuple[str, str]]] = {}
    for (src, tgt, key), edge_type in kg._edge_types.items():
        out.setdefault(edge_type.value, set()).add((src, tgt))
    return out


def _directed_edges(kg: KnowledgeGraph) -> Set[Tuple[str, str, str]]:
    out: Set[Tuple[str, str, str]] = set()
    for (src, tgt, _), edge_type in kg._edge_types.items():
        out.add((src, tgt, edge_type.value))
    return out


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def diff_kgs(v1: KnowledgeGraph, v2: KnowledgeGraph) -> KGDiff:
    """Compute a structural diff between two knowledge graphs.

    Args:
        v1: Baseline knowledge graph (e.g. v0.5 KG).
        v2: New knowledge graph (e.g. v0.6 KG).

    Returns:
        A KGDiff capturing node/edge additions, removals, and flipped edges.
    """
    v1_nodes = _nodes_by_type(v1)
    v2_nodes = _nodes_by_type(v2)
    v1_edges = _edges_by_type(v1)
    v2_edges = _edges_by_type(v2)

    diff = KGDiff()

    # Node diffs per type
    for node_type in set(v1_nodes) | set(v2_nodes):
        a = v1_nodes.get(node_type, set())
        b = v2_nodes.get(node_type, set())
        added = sorted(b - a)
        removed = sorted(a - b)
        if added:
            diff.nodes_added[node_type] = added
        if removed:
            diff.nodes_removed[node_type] = removed

    # Edge diffs per type
    for edge_type in set(v1_edges) | set(v2_edges):
        a = v1_edges.get(edge_type, set())
        b = v2_edges.get(edge_type, set())
        added = sorted(b - a)
        removed = sorted(a - b)
        if added:
            diff.edges_added[edge_type] = added
        if removed:
            diff.edges_removed[edge_type] = removed

    # Direction flips: edges present in both graphs with same edge_type where
    # v1 has (src, tgt) and v2 has (tgt, src). Only considered for directed
    # edge types.
    directed_types = {
        EdgeType.GENE_REGULATES,
        EdgeType.DRUG_TARGETS,
        EdgeType.GO_IS_A,
        EdgeType.GO_PART_OF,
        EdgeType.PATHWAY_CONTAINS,
        EdgeType.GENE_IN_PATHWAY,
        EdgeType.GENE_ENCODES,
        EdgeType.VARIANT_IN_GENE,
        EdgeType.VARIANT_AFFECTS,
        EdgeType.DRUG_IN_PATHWAY,
        EdgeType.GENE_ASSOCIATED_PHENOTYPE,
        EdgeType.PATHWAY_ASSOCIATED_PHENOTYPE,
        EdgeType.GENE_HAS_GO,
        EdgeType.GENE_EXPRESSED_IN,
    }
    for edge_type in directed_types:
        et_value = edge_type.value
        a = v1_edges.get(et_value, set())
        b = v2_edges.get(et_value, set())
        if not a or not b:
            continue
        for (src, tgt) in a:
            if (src, tgt) in b:
                continue
            if (tgt, src) in b:
                diff.edges_flipped.append((src, tgt, et_value))
    diff.edges_flipped.sort()

    diff.summary = {
        "nodes_added_total": sum(len(v) for v in diff.nodes_added.values()),
        "nodes_removed_total": sum(len(v) for v in diff.nodes_removed.values()),
        "edges_added_total": sum(len(v) for v in diff.edges_added.values()),
        "edges_removed_total": sum(len(v) for v in diff.edges_removed.values()),
        "edges_flipped_total": len(diff.edges_flipped),
    }

    logger.info(
        "[KGDiff] +%d/-%d nodes, +%d/-%d edges, %d direction flips",
        diff.summary["nodes_added_total"],
        diff.summary["nodes_removed_total"],
        diff.summary["edges_added_total"],
        diff.summary["edges_removed_total"],
        diff.summary["edges_flipped_total"],
    )
    return diff
