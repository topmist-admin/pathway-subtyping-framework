"""
Knowledge Graph Builder

Constructs a heterogeneous biological knowledge graph from various data sources
including pathways, GO terms, protein-protein interactions, and gene annotations.

Requires: pip install pathway-subtyping[graph]
"""

import json
import logging
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    import networkx as nx
except ImportError:
    raise ImportError(
        "networkx is required for the knowledge_graph module. "
        "Install with: pip install pathway-subtyping[graph]"
    )

from .schema import EdgeType, GraphSchema, Node, NodeType

logger = logging.getLogger(__name__)


@dataclass
class PPINetwork:
    """Protein-protein interaction network data."""

    interactions: List[Tuple[str, str, float]]  # (protein1, protein2, score)
    source: str = "STRING"
    score_threshold: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.interactions)

    def filter_by_score(self, min_score: float) -> "PPINetwork":
        """Filter interactions by minimum combined score."""
        filtered = [(p1, p2, score) for p1, p2, score in self.interactions if score >= min_score]
        return PPINetwork(
            interactions=filtered,
            source=self.source,
            score_threshold=min_score,
            metadata={**self.metadata, "filtered_from": len(self.interactions)},
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable format."""
        return {
            "n_interactions": len(self.interactions),
            "source": self.source,
            "score_threshold": self.score_threshold,
            "metadata": self.metadata,
        }


@dataclass
class KnowledgeGraphStats:
    """Statistics about the knowledge graph."""

    n_nodes: int
    n_edges: int
    node_type_counts: Dict[str, int]
    edge_type_counts: Dict[str, int]
    avg_degree: float
    density: float
    n_connected_components: int

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable format."""
        return {
            "n_nodes": self.n_nodes,
            "n_edges": self.n_edges,
            "node_type_counts": self.node_type_counts,
            "edge_type_counts": self.edge_type_counts,
            "avg_degree": round(self.avg_degree, 4),
            "density": round(self.density, 6),
            "n_connected_components": self.n_connected_components,
        }

    def format_report(self) -> str:
        """Generate human-readable report."""
        lines = [
            "Knowledge Graph Statistics",
            "=" * 40,
            f"Nodes: {self.n_nodes}",
            f"Edges: {self.n_edges}",
            f"Connected components: {self.n_connected_components}",
            f"Average degree: {self.avg_degree:.2f}",
            f"Density: {self.density:.6f}",
            "",
            "Node types:",
        ]
        for nt, count in sorted(self.node_type_counts.items()):
            lines.append(f"  {nt}: {count}")
        lines.append("")
        lines.append("Edge types:")
        for et, count in sorted(self.edge_type_counts.items()):
            lines.append(f"  {et}: {count}")
        return "\n".join(lines)


class KnowledgeGraph:
    """
    Heterogeneous knowledge graph for biological data.

    Uses NetworkX internally for graph operations with support for
    multiple node types and edge types.
    """

    def __init__(self, schema: Optional[GraphSchema] = None):
        """
        Initialize knowledge graph.

        Args:
            schema: Optional schema defining valid node/edge types
        """
        self._graph = nx.MultiDiGraph()
        self._schema = schema or GraphSchema.default_schema()
        self._node_types: Dict[str, NodeType] = {}
        self._edge_types: Dict[Tuple[str, str, int], EdgeType] = {}

    @property
    def graph(self) -> nx.MultiDiGraph:
        """Access underlying NetworkX graph."""
        return self._graph

    @property
    def schema(self) -> GraphSchema:
        """Get graph schema."""
        return self._schema

    def add_node(
        self,
        node_id: str,
        node_type: NodeType,
        attributes: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Add a node to the graph.

        Args:
            node_id: Unique node identifier
            node_type: Type of the node
            attributes: Optional node attributes
        """
        if node_type not in self._schema.node_types:
            raise ValueError(f"Invalid node type: {node_type}")

        attrs = attributes or {}
        self._graph.add_node(node_id, node_type=node_type.value, **attrs)
        self._node_types[node_id] = node_type

    def add_nodes(
        self,
        node_ids: List[str],
        node_type: NodeType,
        attributes: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        """
        Add multiple nodes of the same type.

        Args:
            node_ids: List of node identifiers
            node_type: Type of all nodes
            attributes: Optional dict mapping node_id to attributes
        """
        attrs = attributes or {}
        for node_id in node_ids:
            node_attrs = attrs.get(node_id, {})
            self.add_node(node_id, node_type, node_attrs)

    def add_edge(
        self,
        source: str,
        target: str,
        edge_type: EdgeType,
        weight: float = 1.0,
        attributes: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Add an edge to the graph.

        Args:
            source: Source node ID
            target: Target node ID
            edge_type: Type of the edge
            weight: Edge weight
            attributes: Optional edge attributes
        """
        if source not in self._node_types:
            raise ValueError(f"Source node not found: {source}")
        if target not in self._node_types:
            raise ValueError(f"Target node not found: {target}")

        source_type = self._node_types[source]
        target_type = self._node_types[target]

        if not self._schema.is_valid_edge(source_type, target_type, edge_type):
            raise ValueError(
                f"Invalid edge type {edge_type} for " f"{source_type} -> {target_type}"
            )

        attrs = attributes or {}
        key = self._graph.add_edge(
            source, target, edge_type=edge_type.value, weight=weight, **attrs
        )
        self._edge_types[(source, target, key)] = edge_type

    def add_edges(
        self,
        edges: List[Tuple[str, str]],
        edge_type: EdgeType,
        weights: Optional[List[float]] = None,
        skip_missing_nodes: bool = True,
    ) -> int:
        """
        Add multiple edges of the same type.

        Args:
            edges: List of (source, target) tuples
            edge_type: Type of all edges
            weights: Optional weights for each edge
            skip_missing_nodes: If True, skip edges with missing nodes

        Returns:
            Number of edges added
        """
        if weights is None:
            weights = [1.0] * len(edges)

        added = 0
        for (source, target), weight in zip(edges, weights):
            try:
                self.add_edge(source, target, edge_type, weight)
                added += 1
            except ValueError as e:
                if not skip_missing_nodes:
                    raise
                logger.debug(
                    "[KnowledgeGraph] Skipping edge %s->%s: %s",
                    source,
                    target,
                    e,
                )

        return added

    def get_node(self, node_id: str) -> Optional[Node]:
        """Get a node by ID."""
        if node_id not in self._graph:
            return None
        attrs = dict(self._graph.nodes[node_id])
        node_type = NodeType(attrs.pop("node_type"))
        return Node(id=node_id, node_type=node_type, attributes=attrs)

    def get_node_type(self, node_id: str) -> Optional[NodeType]:
        """Get the type of a node."""
        return self._node_types.get(node_id)

    def get_nodes_by_type(self, node_type: NodeType) -> List[str]:
        """Get all nodes of a specific type."""
        return [node_id for node_id, nt in self._node_types.items() if nt == node_type]

    def get_neighbors(
        self,
        node_id: str,
        edge_type: Optional[EdgeType] = None,
        direction: str = "out",
    ) -> List[str]:
        """
        Get neighbors of a node.

        Args:
            node_id: Node to get neighbors for
            edge_type: Optional edge type filter
            direction: "out" for successors, "in" for predecessors,
                       "both" for both

        Returns:
            List of neighbor node IDs
        """
        if node_id not in self._graph:
            return []

        neighbors: Set[str] = set()

        if direction in ("out", "both"):
            for _, target, data in self._graph.out_edges(node_id, data=True):
                if edge_type is None or data.get("edge_type") == edge_type.value:
                    neighbors.add(target)

        if direction in ("in", "both"):
            for source, _, data in self._graph.in_edges(node_id, data=True):
                if edge_type is None or data.get("edge_type") == edge_type.value:
                    neighbors.add(source)

        return list(neighbors)

    def get_edges_by_type(self, edge_type: EdgeType) -> List[Tuple[str, str, float]]:
        """Get all edges of a specific type."""
        edges = []
        for source, target, data in self._graph.edges(data=True):
            if data.get("edge_type") == edge_type.value:
                edges.append((source, target, data.get("weight", 1.0)))
        return edges

    def has_node(self, node_id: str) -> bool:
        """Check if node exists."""
        return node_id in self._graph

    def has_edge(
        self,
        source: str,
        target: str,
        edge_type: Optional[EdgeType] = None,
    ) -> bool:
        """Check if edge exists."""
        if not self._graph.has_edge(source, target):
            return False
        if edge_type is None:
            return True
        for _, _, data in self._graph.edges(source, data=True):
            if data.get("edge_type") == edge_type.value:
                return True
        return False

    @property
    def n_nodes(self) -> int:
        """Number of nodes in the graph."""
        return self._graph.number_of_nodes()

    @property
    def n_edges(self) -> int:
        """Number of edges in the graph."""
        return self._graph.number_of_edges()

    def get_stats(self) -> KnowledgeGraphStats:
        """Get statistics about the graph."""
        node_type_counts: Dict[str, int] = {}
        for node_type in NodeType:
            count = len(self.get_nodes_by_type(node_type))
            if count > 0:
                node_type_counts[node_type.value] = count

        edge_type_counts: Dict[str, int] = {}
        for source, target, data in self._graph.edges(data=True):
            et = data.get("edge_type", "unknown")
            edge_type_counts[et] = edge_type_counts.get(et, 0) + 1

        # Calculate average degree
        if self.n_nodes > 0:
            avg_degree = self.n_edges / self.n_nodes
        else:
            avg_degree = 0.0

        # Calculate density
        if self.n_nodes > 1:
            max_edges = self.n_nodes * (self.n_nodes - 1)
            density = self.n_edges / max_edges
        else:
            density = 0.0

        # Count connected components (treating as undirected)
        undirected = self._graph.to_undirected()
        n_components = nx.number_connected_components(undirected)

        return KnowledgeGraphStats(
            n_nodes=self.n_nodes,
            n_edges=self.n_edges,
            node_type_counts=node_type_counts,
            edge_type_counts=edge_type_counts,
            avg_degree=avg_degree,
            density=density,
            n_connected_components=n_components,
        )

    def subgraph(
        self,
        node_ids: List[str],
        include_edges: bool = True,
    ) -> "KnowledgeGraph":
        """
        Create a subgraph with specified nodes.

        Args:
            node_ids: Nodes to include
            include_edges: Whether to include edges between nodes

        Returns:
            New KnowledgeGraph with subset of nodes
        """
        sub = KnowledgeGraph(schema=self._schema)

        for node_id in node_ids:
            if node_id in self._node_types:
                node = self.get_node(node_id)
                if node:
                    sub.add_node(node_id, node.node_type, node.attributes)

        if include_edges:
            node_set = set(node_ids)
            for source, target, data in self._graph.edges(data=True):
                if source in node_set and target in node_set:
                    edge_type = EdgeType(data.get("edge_type"))
                    weight = data.get("weight", 1.0)
                    attrs = {k: v for k, v in data.items() if k not in ("edge_type", "weight")}
                    sub.add_edge(source, target, edge_type, weight, attrs)

        return sub

    # ── Topology-aware methods (v0.5 QC layer support) ────────────────

    def get_pathway_genes(self, pathway_id: str) -> List[str]:
        """Get all gene nodes connected to a pathway via GENE_IN_PATHWAY edges.

        Args:
            pathway_id: Pathway node ID.

        Returns:
            List of gene node IDs in the pathway.
        """
        return self.get_neighbors(pathway_id, EdgeType.GENE_IN_PATHWAY, direction="in")

    def get_pathway_subgraph(self, pathway_id: str) -> "KnowledgeGraph":
        """Extract a subgraph containing only genes in a pathway and their edges.

        Args:
            pathway_id: Pathway node ID.

        Returns:
            KnowledgeGraph containing only pathway genes and inter-gene edges.
        """
        genes = self.get_pathway_genes(pathway_id)
        return self.subgraph(genes, include_edges=True)

    def get_directed_edges_in_pathway(
        self,
        pathway_id: str,
        exclude_symmetric: bool = True,
    ) -> List[Tuple[str, str, float]]:
        """Get directed (non-symmetric) edges between genes in a pathway.

        Args:
            pathway_id: Pathway node ID.
            exclude_symmetric: If True, exclude PPI and co-expression edges
                (which are symmetric and don't indicate signal direction).

        Returns:
            List of (source, target, weight) tuples for directed edges.
        """
        from .schema import EDGE_TYPE_METADATA

        genes = set(self.get_pathway_genes(pathway_id))
        directed = []

        for source, target, data in self._graph.edges(data=True):
            if source not in genes or target not in genes:
                continue
            edge_type_val = data.get("edge_type")
            if edge_type_val is None:
                continue
            try:
                et = EdgeType(edge_type_val)
            except ValueError:
                continue

            if exclude_symmetric:
                meta = EDGE_TYPE_METADATA.get(et, {})
                if meta.get("symmetric", False):
                    continue

            weight = data.get("weight", 1.0)
            directed.append((source, target, weight))

        return directed

    def get_in_degree(
        self,
        node_id: str,
        edge_type: Optional[EdgeType] = None,
        within_nodes: Optional[Set[str]] = None,
    ) -> int:
        """Count incoming edges to a node.

        Args:
            node_id: Target node ID.
            edge_type: Optional filter by edge type.
            within_nodes: Optional set of nodes to restrict to (pathway scope).

        Returns:
            Number of incoming edges.
        """
        if node_id not in self._graph:
            return 0
        count = 0
        for source, _, data in self._graph.in_edges(node_id, data=True):
            if within_nodes is not None and source not in within_nodes:
                continue
            if edge_type is not None and data.get("edge_type") != edge_type.value:
                continue
            count += 1
        return count

    def get_out_degree(
        self,
        node_id: str,
        edge_type: Optional[EdgeType] = None,
        within_nodes: Optional[Set[str]] = None,
    ) -> int:
        """Count outgoing edges from a node.

        Args:
            node_id: Source node ID.
            edge_type: Optional filter by edge type.
            within_nodes: Optional set of nodes to restrict to (pathway scope).

        Returns:
            Number of outgoing edges.
        """
        if node_id not in self._graph:
            return 0
        count = 0
        for _, target, data in self._graph.out_edges(node_id, data=True):
            if within_nodes is not None and target not in within_nodes:
                continue
            if edge_type is not None and data.get("edge_type") != edge_type.value:
                continue
            count += 1
        return count

    def partition_pathway_genes(
        self,
        pathway_id: str,
        edge_type: Optional[EdgeType] = None,
    ) -> Dict[str, List[str]]:
        """Partition pathway genes into upstream, intermediate, and downstream layers.

        Uses in-degree/out-degree within the pathway subgraph:
        - **upstream** (initiators): genes with in-degree = 0 within the pathway
        - **downstream** (effectors): genes with out-degree = 0 within the pathway
        - **intermediate** (transducers): genes with both incoming and outgoing edges

        If the pathway has no directed edges (PPI-only), falls back to
        equal-thirds partitioning by node degree centrality.

        Args:
            pathway_id: Pathway node ID.
            edge_type: Optional edge type filter for degree computation.
                If None, uses all non-symmetric edge types.

        Returns:
            Dict with keys 'upstream', 'intermediate', 'downstream',
            each mapping to a list of gene node IDs.
        """
        genes = self.get_pathway_genes(pathway_id)
        if not genes:
            return {"upstream": [], "intermediate": [], "downstream": []}

        gene_set = set(genes)

        # Check if we have directed (non-symmetric) edges in this pathway
        directed_edges = self.get_directed_edges_in_pathway(pathway_id, exclude_symmetric=True)

        if directed_edges:
            # Use directed edges for topology-aware partitioning
            upstream = []
            downstream = []
            intermediate = []

            for gene in genes:
                in_deg = self.get_in_degree(gene, edge_type, within_nodes=gene_set)
                out_deg = self.get_out_degree(gene, edge_type, within_nodes=gene_set)

                if in_deg == 0 and out_deg > 0:
                    upstream.append(gene)
                elif out_deg == 0 and in_deg > 0:
                    downstream.append(gene)
                elif in_deg > 0 and out_deg > 0:
                    intermediate.append(gene)
                else:
                    # Isolated node (no directed edges) — classify by PPI degree
                    intermediate.append(gene)
        else:
            # Fallback: no directed edges, partition by PPI degree centrality
            sub = self.get_pathway_subgraph(pathway_id)
            if sub.n_edges > 0:
                centrality = nx.degree_centrality(sub.graph.to_undirected())
                sorted_genes = sorted(genes, key=lambda g: centrality.get(g, 0))
            else:
                sorted_genes = list(genes)

            n = len(sorted_genes)
            n_up = max(1, n // 3)
            n_down = max(1, n // 3)
            upstream = sorted_genes[:n_up]
            downstream = sorted_genes[-n_down:]
            intermediate = sorted_genes[n_up:-n_down] if n > n_up + n_down else []

        return {
            "upstream": upstream,
            "intermediate": intermediate,
            "downstream": downstream,
        }

    def topological_sort_pathway(self, pathway_id: str) -> List[str]:
        """Topologically sort genes within a pathway by directed edges.

        If the pathway subgraph has cycles, returns genes sorted by
        out-degree (most connections first) as a best-effort ordering.

        Args:
            pathway_id: Pathway node ID.

        Returns:
            List of gene node IDs in topological order (upstream first).
        """
        sub = self.get_pathway_subgraph(pathway_id)
        if sub.n_nodes == 0:
            return []

        try:
            return list(nx.topological_sort(sub.graph))
        except nx.NetworkXUnfeasible:
            # Graph has cycles — fallback to degree-based ordering
            logger.debug(
                "[KnowledgeGraph] Pathway %s has cycles; using degree-based sort",
                pathway_id,
            )
            genes = list(sub.graph.nodes())
            gene_set = set(genes)
            genes.sort(
                key=lambda g: self.get_out_degree(g, within_nodes=gene_set),
                reverse=True,
            )
            return genes

    def find_cascade_paths(
        self,
        source: str,
        target: str,
        pathway_id: Optional[str] = None,
        max_depth: int = 5,
    ) -> List[List[str]]:
        """Find all directed paths from source to target.

        Args:
            source: Source gene node ID.
            target: Target gene node ID.
            pathway_id: Optional pathway to restrict search to.
            max_depth: Maximum path length.

        Returns:
            List of paths, each a list of gene node IDs.
        """
        if pathway_id:
            sub = self.get_pathway_subgraph(pathway_id)
            graph = sub.graph
        else:
            graph = self._graph

        if source not in graph or target not in graph:
            return []

        try:
            paths = list(nx.all_simple_paths(graph, source, target, cutoff=max_depth))
            return paths
        except nx.NodeNotFound:
            return []

    def get_shared_genes(
        self,
        pathway_a: str,
        pathway_b: str,
    ) -> List[str]:
        """Find genes shared between two pathways.

        Args:
            pathway_a: First pathway node ID.
            pathway_b: Second pathway node ID.

        Returns:
            List of gene node IDs present in both pathways.
        """
        genes_a = set(self.get_pathway_genes(pathway_a))
        genes_b = set(self.get_pathway_genes(pathway_b))
        return list(genes_a & genes_b)

    def get_pathway_crosstalk(
        self,
        pathway_a: str,
        pathway_b: str,
    ) -> Dict[str, Any]:
        """Quantify crosstalk between two pathways.

        Measures shared genes and inter-pathway edge density.

        Args:
            pathway_a: First pathway node ID.
            pathway_b: Second pathway node ID.

        Returns:
            Dict with crosstalk metrics:
                - shared_genes: list of shared gene IDs
                - n_shared: count of shared genes
                - inter_pathway_edges: edges connecting genes across pathways
                - edge_density: inter-pathway edges / possible inter-pathway edges
                - jaccard_index: |A ∩ B| / |A ∪ B|
        """
        genes_a = set(self.get_pathway_genes(pathway_a))
        genes_b = set(self.get_pathway_genes(pathway_b))
        shared = genes_a & genes_b
        union = genes_a | genes_b

        # Count edges between pathway A genes and pathway B genes
        inter_edges = []
        for source, target, data in self._graph.edges(data=True):
            if (source in genes_a and target in genes_b and source not in shared) or (
                source in genes_b and target in genes_a and source not in shared
            ):
                inter_edges.append((source, target, data.get("weight", 1.0)))

        max_inter = len(genes_a - shared) * len(genes_b - shared)
        edge_density = len(inter_edges) / max(max_inter, 1)
        jaccard = len(shared) / max(len(union), 1)

        return {
            "shared_genes": list(shared),
            "n_shared": len(shared),
            "inter_pathway_edges": inter_edges,
            "n_inter_edges": len(inter_edges),
            "edge_density": edge_density,
            "jaccard_index": jaccard,
            "pathway_a_size": len(genes_a),
            "pathway_b_size": len(genes_b),
        }

    # ── Topology-weighted scoring ────────────────────────────────────────

    def compute_centrality(
        self,
        pathway_id: Optional[str] = None,
        method: str = "degree",
    ) -> Dict[str, float]:
        """Compute node centrality scores, optionally within a pathway.

        Args:
            pathway_id: If given, compute centrality only within the pathway
                subgraph. Otherwise uses the full gene-gene graph.
            method: Centrality method. One of "degree", "betweenness",
                "closeness", "pagerank".

        Returns:
            Dict mapping gene node ID to centrality score.
        """
        if pathway_id:
            sub = self.get_pathway_subgraph(pathway_id)
            g = sub.graph
        else:
            # Use only gene nodes
            gene_nodes = self.get_nodes_by_type(NodeType.GENE)
            g = self._graph.subgraph(gene_nodes)

        if g.number_of_nodes() == 0:
            return {}

        undirected = g.to_undirected()

        if method == "degree":
            return dict(nx.degree_centrality(undirected))
        elif method == "betweenness":
            return dict(nx.betweenness_centrality(undirected))
        elif method == "closeness":
            return dict(nx.closeness_centrality(undirected))
        elif method == "pagerank":
            return dict(nx.pagerank(g, alpha=0.85))
        else:
            raise ValueError(f"Unknown centrality method: {method}")

    def topology_weighted_pathway_score(
        self,
        expression: Any,  # pd.DataFrame
        pathway_id: str,
        centrality_method: str = "degree",
    ) -> Any:  # np.ndarray
        """Compute pathway scores weighted by gene centrality in the network.

        Instead of equal-weight mean of all pathway genes, each gene's
        contribution is weighted by its network centrality (degree,
        betweenness, or PageRank). Hub genes contribute more.

        Args:
            expression: Expression matrix (cells x genes), pd.DataFrame.
            pathway_id: Pathway node ID.
            centrality_method: Centrality method for weighting.

        Returns:
            1D numpy array of per-cell weighted pathway scores.
        """
        import numpy as np_local

        genes = self.get_pathway_genes(pathway_id)
        present = [g for g in genes if g in expression.columns]
        if not present:
            return np_local.zeros(len(expression))

        centrality = self.compute_centrality(pathway_id, method=centrality_method)

        weights = np_local.array([centrality.get(g, 0.0) for g in present])
        total_weight = weights.sum()
        if total_weight == 0:
            weights = np_local.ones(len(present))
            total_weight = len(present)
        weights = weights / total_weight

        vals = expression[present].values
        # Z-score per gene
        means = vals.mean(axis=0, keepdims=True)
        stds = vals.std(axis=0, keepdims=True)
        stds[stds == 0] = 1.0
        z = (vals - means) / stds

        return z @ weights

    # ── Hierarchical pathway queries ──────────────────────────────────

    def get_child_pathways(self, pathway_id: str) -> List[str]:
        """Get sub-pathways contained within a parent pathway.

        Uses PATHWAY_CONTAINS edges (parent -> child).

        Args:
            pathway_id: Parent pathway node ID.

        Returns:
            List of child pathway node IDs.
        """
        return self.get_neighbors(pathway_id, EdgeType.PATHWAY_CONTAINS, direction="out")

    def get_parent_pathways(self, pathway_id: str) -> List[str]:
        """Get parent pathways that contain this pathway.

        Args:
            pathway_id: Child pathway node ID.

        Returns:
            List of parent pathway node IDs.
        """
        return self.get_neighbors(pathway_id, EdgeType.PATHWAY_CONTAINS, direction="in")

    def get_pathway_hierarchy(self, root_pathway_id: str) -> Dict[str, Any]:
        """Get the full hierarchy tree rooted at a pathway.

        Recursively traverses PATHWAY_CONTAINS edges to build a tree.

        Args:
            root_pathway_id: Root pathway node ID.

        Returns:
            Nested dict: {"id": str, "name": str, "children": [...],
            "genes": [...]}
        """
        node = self.get_node(root_pathway_id)
        name = node.attributes.get("name", root_pathway_id) if node else root_pathway_id

        children = self.get_child_pathways(root_pathway_id)
        genes = self.get_pathway_genes(root_pathway_id)

        return {
            "id": root_pathway_id,
            "name": name,
            "n_genes": len(genes),
            "genes": genes,
            "children": [self.get_pathway_hierarchy(child) for child in children],
        }

    def get_all_descendant_genes(self, pathway_id: str) -> List[str]:
        """Get all genes in a pathway and its sub-pathways (recursively).

        Args:
            pathway_id: Root pathway node ID.

        Returns:
            Deduplicated list of all gene node IDs across the hierarchy.
        """
        all_genes: Set[str] = set()
        all_genes.update(self.get_pathway_genes(pathway_id))

        for child in self.get_child_pathways(pathway_id):
            all_genes.update(self.get_all_descendant_genes(child))

        return list(all_genes)

    # ── Cross-omics entity resolution ─────────────────────────────────

    def resolve_gene_to_protein(self, gene_id: str) -> List[str]:
        """Resolve a gene to its encoded protein(s).

        Uses GENE_ENCODES edges.

        Args:
            gene_id: Gene node ID.

        Returns:
            List of protein node IDs.
        """
        return self.get_neighbors(gene_id, EdgeType.GENE_ENCODES, direction="out")

    def resolve_protein_to_gene(self, protein_id: str) -> List[str]:
        """Resolve a protein back to its encoding gene(s).

        Args:
            protein_id: Protein node ID.

        Returns:
            List of gene node IDs.
        """
        return self.get_neighbors(protein_id, EdgeType.GENE_ENCODES, direction="in")

    def resolve_entity_chain(
        self,
        start_id: str,
        target_type: "NodeType",
        max_hops: int = 3,
    ) -> List[List[str]]:
        """Resolve an entity across omics layers by traversing edges.

        Finds all paths from start_id to any node of target_type within
        max_hops, following directed edges.

        Args:
            start_id: Starting node ID (e.g., a gene).
            target_type: Target node type (e.g., NodeType.DRUG).
            max_hops: Maximum path length.

        Returns:
            List of paths, each a list of node IDs.
        """
        if start_id not in self._graph:
            return []

        paths: List[List[str]] = []
        self._dfs_resolve(start_id, target_type, max_hops, [start_id], paths)
        return paths

    def _dfs_resolve(
        self,
        current: str,
        target_type: "NodeType",
        remaining_hops: int,
        path: List[str],
        results: List[List[str]],
    ) -> None:
        """DFS helper for entity resolution."""
        if remaining_hops <= 0:
            return

        for _, neighbor, _ in self._graph.out_edges(current, data=True):
            if neighbor in path:
                continue
            new_path = path + [neighbor]
            nt = self._node_types.get(neighbor)
            if nt == target_type:
                results.append(new_path)
            self._dfs_resolve(neighbor, target_type, remaining_hops - 1, new_path, results)

    def get_drug_targets_in_pathway(self, pathway_id: str) -> Dict[str, List[str]]:
        """Find all drugs targeting genes in a pathway.

        Args:
            pathway_id: Pathway node ID.

        Returns:
            Dict mapping gene_id to list of drug_ids that target it.
        """
        genes = self.get_pathway_genes(pathway_id)
        gene_drugs: Dict[str, List[str]] = {}

        for gene in genes:
            drugs = self.get_neighbors(gene, EdgeType.DRUG_TARGETS, direction="in")
            if drugs:
                gene_drugs[gene] = drugs

        return gene_drugs

    def save(self, path: str) -> None:
        """
        Save graph to file.

        Args:
            path: Output path (supports .gpickle, .json)
        """
        filepath = Path(path)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        if filepath.suffix == ".gpickle":
            with open(filepath, "wb") as f:
                pickle.dump(
                    {
                        "graph": self._graph,
                        "node_types": self._node_types,
                        "edge_types": self._edge_types,
                    },
                    f,
                )
        elif filepath.suffix == ".json":
            data = nx.node_link_data(self._graph)
            data["node_types"] = {k: v.value for k, v in self._node_types.items()}
            with open(filepath, "w") as f:
                json.dump(data, f, indent=2)
        else:
            raise ValueError(f"Unsupported format: {filepath.suffix}")

        logger.info("[KnowledgeGraph] Saved graph to %s", filepath)

    @classmethod
    def load(cls, path: str) -> "KnowledgeGraph":
        """
        Load graph from file.

        Args:
            path: Input path

        Returns:
            KnowledgeGraph
        """
        filepath = Path(path)
        kg = cls()

        if filepath.suffix == ".gpickle":
            with open(filepath, "rb") as f:
                data = pickle.load(f)  # noqa: S301
            kg._graph = data["graph"]
            kg._node_types = data["node_types"]
            kg._edge_types = data.get("edge_types", {})
        elif filepath.suffix == ".json":
            with open(filepath, "r") as f:
                data = json.load(f)
            node_types_data = data.pop("node_types", {})
            kg._graph = nx.node_link_graph(data)
            kg._node_types = {k: NodeType(v) for k, v in node_types_data.items()}
        else:
            raise ValueError(f"Unsupported format: {filepath.suffix}")

        logger.info(
            "[KnowledgeGraph] Loaded graph from %s: %d nodes, %d edges",
            filepath,
            kg.n_nodes,
            kg.n_edges,
        )
        return kg

    def get_citations(self) -> List[str]:
        """Return literature citations for methods used."""
        return [
            "Szklarczyk D, et al. The STRING database in 2023: "
            "protein-protein association networks and functional "
            "enrichment analyses for any sequenced genome of interest. "
            "Nucleic Acids Res. 2023;51(D1):D599-D606.",
            "Ashburner M, et al. Gene Ontology: tool for the "
            "unification of biology. Nat Genet. 2000;25(1):25-9.",
        ]


class KnowledgeGraphBuilder:
    """
    Builder for constructing knowledge graphs from various data sources.

    Supports:
    - Pathway databases (GO, Reactome)
    - Protein-protein interaction networks (STRING)
    - Gene annotations
    - Directed signaling edges (Reactome, KEGG, custom)
    """

    def __init__(self, schema: Optional[GraphSchema] = None):
        """
        Initialize builder.

        Args:
            schema: Optional schema for validation
        """
        self._schema = schema or GraphSchema.default_schema()
        self._genes: Set[str] = set()
        self._pathways: Dict[str, Set[str]] = {}
        self._pathway_names: Dict[str, str] = {}
        self._go_terms: Dict[str, Dict[str, str]] = {}
        self._gene_go_annotations: Dict[str, Set[str]] = {}
        self._ppi_edges: List[Tuple[str, str, float]] = []
        self._go_hierarchy: List[Tuple[str, str, str]] = []
        self._signaling_edges: List[Tuple[str, str, float]] = []

    def add_genes(self, gene_list: List[str]) -> "KnowledgeGraphBuilder":
        """
        Add genes to the graph.

        Args:
            gene_list: List of gene symbols/IDs

        Returns:
            Self for chaining
        """
        self._genes.update(gene_list)
        logger.info(
            "[KnowledgeGraph] Added %d genes (total: %d)",
            len(gene_list),
            len(self._genes),
        )
        return self

    def add_pathways(self, pathway_db: Any) -> "KnowledgeGraphBuilder":
        """
        Add pathways from a PathwayDatabase-like object.

        Args:
            pathway_db: Object with .pathways (Dict[str, Set[str]]) and
                        .pathway_names (Dict[str, str]) attributes

        Returns:
            Self for chaining
        """
        for pathway_id, genes in pathway_db.pathways.items():
            self._pathways[pathway_id] = genes.copy() if isinstance(genes, set) else set(genes)
            self._pathway_names[pathway_id] = pathway_db.pathway_names.get(pathway_id, pathway_id)
            self._genes.update(genes)

        logger.info(
            "[KnowledgeGraph] Added %d pathways from database",
            len(pathway_db.pathways),
        )
        return self

    def add_pathways_from_dict(
        self,
        pathways: Dict[str, List[str]],
        pathway_names: Optional[Dict[str, str]] = None,
    ) -> "KnowledgeGraphBuilder":
        """
        Add pathways from a plain dictionary (framework standard format).

        Args:
            pathways: Dict mapping pathway ID to list of gene symbols
            pathway_names: Optional dict mapping pathway ID to display name

        Returns:
            Self for chaining
        """
        names = pathway_names or {}
        for pathway_id, genes in pathways.items():
            self._pathways[pathway_id] = set(genes)
            self._pathway_names[pathway_id] = names.get(pathway_id, pathway_id)
            self._genes.update(genes)

        logger.info(
            "[KnowledgeGraph] Added %d pathways from dict",
            len(pathways),
        )
        return self

    def add_go_terms(
        self,
        go_terms: Dict[str, Dict[str, str]],
        annotations: Dict[str, Set[str]],
    ) -> "KnowledgeGraphBuilder":
        """
        Add GO terms and gene annotations.

        Args:
            go_terms: Dict mapping GO ID to term info (name, namespace)
            annotations: Dict mapping gene symbol to set of GO IDs

        Returns:
            Self for chaining
        """
        self._go_terms.update(go_terms)
        for gene, go_ids in annotations.items():
            if gene not in self._gene_go_annotations:
                self._gene_go_annotations[gene] = set()
            self._gene_go_annotations[gene].update(go_ids)
            self._genes.add(gene)

        logger.info(
            "[KnowledgeGraph] Added %d GO terms with %d gene annotations",
            len(go_terms),
            len(annotations),
        )
        return self

    def add_go_hierarchy(
        self,
        hierarchy: List[Tuple[str, str, str]],
    ) -> "KnowledgeGraphBuilder":
        """
        Add GO term hierarchy relationships.

        Args:
            hierarchy: List of (child_id, parent_id, relationship_type)
                      where relationship_type is 'is_a' or 'part_of'

        Returns:
            Self for chaining
        """
        self._go_hierarchy.extend(hierarchy)
        logger.info(
            "[KnowledgeGraph] Added %d GO hierarchy relationships",
            len(hierarchy),
        )
        return self

    def add_ppi(
        self,
        ppi_network: PPINetwork,
        gene_to_protein: Optional[Dict[str, str]] = None,
    ) -> "KnowledgeGraphBuilder":
        """
        Add protein-protein interactions.

        Args:
            ppi_network: PPINetwork with interaction data
            gene_to_protein: Optional mapping from gene symbols to protein IDs

        Returns:
            Self for chaining
        """
        for protein1, protein2, score in ppi_network.interactions:
            if gene_to_protein:
                protein_to_gene = {v: k for k, v in gene_to_protein.items()}
                gene1 = protein_to_gene.get(protein1, protein1)
                gene2 = protein_to_gene.get(protein2, protein2)
            else:
                gene1 = self._protein_id_to_gene(protein1)
                gene2 = self._protein_id_to_gene(protein2)

            self._ppi_edges.append((gene1, gene2, score))
            self._genes.add(gene1)
            self._genes.add(gene2)

        logger.info(
            "[KnowledgeGraph] Added %d PPI edges (threshold: %s)",
            len(ppi_network.interactions),
            ppi_network.score_threshold,
        )
        return self

    def add_ppi_from_dataframe(
        self,
        df: Any,
        protein1_col: str = "protein1",
        protein2_col: str = "protein2",
        score_col: str = "combined_score",
        min_score: float = 0.0,
    ) -> "KnowledgeGraphBuilder":
        """
        Add PPI from a pandas DataFrame.

        Args:
            df: DataFrame with PPI data
            protein1_col: Column name for first protein
            protein2_col: Column name for second protein
            score_col: Column name for interaction score
            min_score: Minimum score threshold

        Returns:
            Self for chaining
        """
        count = 0
        for _, row in df.iterrows():
            score = row[score_col]
            if score < min_score:
                continue

            protein1 = str(row[protein1_col])
            protein2 = str(row[protein2_col])

            gene1 = self._protein_id_to_gene(protein1)
            gene2 = self._protein_id_to_gene(protein2)

            self._ppi_edges.append((gene1, gene2, score))
            self._genes.add(gene1)
            self._genes.add(gene2)
            count += 1

        logger.info("[KnowledgeGraph] Added %d PPI edges from DataFrame", count)
        return self

    def add_signaling_edges(
        self,
        edges: List[Tuple[str, str]],
        weights: Optional[List[float]] = None,
    ) -> "KnowledgeGraphBuilder":
        """Add directed signaling/regulatory edges between genes.

        These represent directed signal flow (e.g., RAS activates RAF,
        kinase phosphorylates substrate). Unlike PPI edges, these are
        directional and used for cascade layer partitioning.

        Args:
            edges: List of (source_gene, target_gene) tuples indicating
                signal direction (source regulates/activates target).
            weights: Optional confidence weights per edge. Defaults to 1.0.

        Returns:
            Self for chaining.
        """
        if weights is None:
            weights = [1.0] * len(edges)

        for (source, target), weight in zip(edges, weights):
            self._signaling_edges.append((source, target, weight))
            self._genes.add(source)
            self._genes.add(target)

        logger.info(
            "[KnowledgeGraph] Added %d directed signaling edges (total: %d)",
            len(edges),
            len(self._signaling_edges),
        )
        return self

    def add_signaling_edges_from_dict(
        self,
        pathway_edges: Dict[str, List[Tuple[str, str]]],
        default_weight: float = 1.0,
    ) -> "KnowledgeGraphBuilder":
        """Add directed signaling edges organized by pathway.

        Convenience method for adding edges from pathway-keyed dictionaries,
        as commonly extracted from Reactome or KEGG.

        Args:
            pathway_edges: Dict mapping pathway_id to list of (source, target)
                directed edges within that pathway.
            default_weight: Weight for all edges.

        Returns:
            Self for chaining.
        """
        total = 0
        for pathway_id, edges in pathway_edges.items():
            for source, target in edges:
                self._signaling_edges.append((source, target, default_weight))
                self._genes.add(source)
                self._genes.add(target)
                total += 1

        logger.info(
            "[KnowledgeGraph] Added %d signaling edges across %d pathways",
            total,
            len(pathway_edges),
        )
        return self

    def _protein_id_to_gene(self, protein_id: str) -> str:
        """
        Convert protein ID to gene symbol.

        STRING uses format like '9606.ENSP00000269305'.
        """
        if "." in protein_id:
            protein_id = protein_id.split(".")[-1]
        return protein_id

    def build(self) -> KnowledgeGraph:
        """
        Build the knowledge graph from added data.

        Returns:
            Constructed KnowledgeGraph
        """
        kg = KnowledgeGraph(schema=self._schema)

        # Add gene nodes
        logger.info("[KnowledgeGraph] Adding %d gene nodes...", len(self._genes))
        for gene in self._genes:
            kg.add_node(gene, NodeType.GENE)

        # Add pathway nodes and gene-pathway edges
        logger.info(
            "[KnowledgeGraph] Adding %d pathway nodes...",
            len(self._pathways),
        )
        for pathway_id, genes in self._pathways.items():
            kg.add_node(
                pathway_id,
                NodeType.PATHWAY,
                {"name": self._pathway_names.get(pathway_id, pathway_id)},
            )
            for gene in genes:
                if kg.has_node(gene):
                    kg.add_edge(gene, pathway_id, EdgeType.GENE_IN_PATHWAY)

        # Add GO term nodes and gene-GO edges
        logger.info(
            "[KnowledgeGraph] Adding %d GO term nodes...",
            len(self._go_terms),
        )
        for go_id, term_info in self._go_terms.items():
            kg.add_node(
                go_id,
                NodeType.GO_TERM,
                {
                    "name": term_info.get("name", go_id),
                    "namespace": term_info.get("namespace", "unknown"),
                },
            )

        for gene, go_ids in self._gene_go_annotations.items():
            if not kg.has_node(gene):
                continue
            for go_id in go_ids:
                if kg.has_node(go_id):
                    kg.add_edge(gene, go_id, EdgeType.GENE_HAS_GO)

        # Add GO hierarchy edges
        logger.info(
            "[KnowledgeGraph] Adding %d GO hierarchy edges...",
            len(self._go_hierarchy),
        )
        for child_id, parent_id, rel_type in self._go_hierarchy:
            if kg.has_node(child_id) and kg.has_node(parent_id):
                if rel_type == "is_a":
                    kg.add_edge(child_id, parent_id, EdgeType.GO_IS_A)
                elif rel_type == "part_of":
                    kg.add_edge(child_id, parent_id, EdgeType.GO_PART_OF)

        # Add PPI edges
        logger.info("[KnowledgeGraph] Adding %d PPI edges...", len(self._ppi_edges))
        ppi_added = 0
        for gene1, gene2, score in self._ppi_edges:
            if kg.has_node(gene1) and kg.has_node(gene2) and gene1 != gene2:
                normalized_score = score / 1000.0 if score > 1 else score
                kg.add_edge(gene1, gene2, EdgeType.GENE_INTERACTS, normalized_score)
                ppi_added += 1

        # Add directed signaling edges
        if self._signaling_edges:
            logger.info(
                "[KnowledgeGraph] Adding %d signaling edges...",
                len(self._signaling_edges),
            )
            sig_added = 0
            for source, target, weight in self._signaling_edges:
                if kg.has_node(source) and kg.has_node(target) and source != target:
                    kg.add_edge(source, target, EdgeType.GENE_REGULATES, weight)
                    sig_added += 1

        stats = kg.get_stats()
        logger.info(
            "[KnowledgeGraph] Built knowledge graph: %d nodes, %d edges, " "%d components",
            stats.n_nodes,
            stats.n_edges,
            stats.n_connected_components,
        )

        return kg

    def clear(self) -> "KnowledgeGraphBuilder":
        """Clear all added data."""
        self._genes.clear()
        self._pathways.clear()
        self._pathway_names.clear()
        self._go_terms.clear()
        self._gene_go_annotations.clear()
        self._ppi_edges.clear()
        self._go_hierarchy.clear()
        self._signaling_edges.clear()
        return self


def load_ppi_from_file(
    file_path: str,
    min_score: float = 400.0,
    max_edges: Optional[int] = None,
) -> PPINetwork:
    """
    Load PPI network from a STRING-format file.

    Args:
        file_path: Path to STRING PPI file
        min_score: Minimum combined score (STRING uses 0-1000)
        max_edges: Optional maximum number of edges to load

    Returns:
        PPINetwork
    """
    interactions: List[Tuple[str, str, float]] = []
    filepath = Path(file_path)

    with open(filepath, "r") as f:
        f.readline()  # Skip header

        for i, line in enumerate(f):
            if max_edges and i >= max_edges:
                break

            parts = line.strip().split()
            if len(parts) >= 3:
                protein1 = parts[0]
                protein2 = parts[1]
                score = float(parts[-1])

                if score >= min_score:
                    interactions.append((protein1, protein2, score))

    logger.info(
        "[KnowledgeGraph] Loaded %d PPI interactions from %s " "(min_score=%s)",
        len(interactions),
        file_path,
        min_score,
    )

    return PPINetwork(
        interactions=interactions,
        source="STRING",
        score_threshold=min_score,
        metadata={"file": str(filepath)},
    )
