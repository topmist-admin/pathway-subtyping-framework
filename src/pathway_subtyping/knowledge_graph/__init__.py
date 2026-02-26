"""
Knowledge Graph Module

Provides tools for building and exporting heterogeneous biological
knowledge graphs from pathway databases, GO terms, and protein-protein
interaction networks.

Requires: pip install pathway-subtyping[graph]
"""

try:
    import networkx  # noqa: F401
except ImportError:
    raise ImportError(
        "The knowledge_graph module requires networkx. "
        "Install with: pip install pathway-subtyping[graph]"
    )

from .builder import (
    KnowledgeGraph,
    KnowledgeGraphBuilder,
    KnowledgeGraphStats,
    PPINetwork,
    load_ppi_from_file,
)
from .exporters import (
    NodeMapping,
    create_node_mapping,
    to_adjacency_matrix,
    to_csv,
    to_dgl,
    to_edge_list,
    to_neo4j_cypher,
    to_pyg,
)
from .schema import (
    EDGE_TYPE_METADATA,
    Edge,
    EdgeType,
    GraphSchema,
    Node,
    NodeFeatures,
    NodeType,
)

__all__ = [
    # Builder
    "KnowledgeGraph",
    "KnowledgeGraphBuilder",
    "KnowledgeGraphStats",
    "PPINetwork",
    "load_ppi_from_file",
    # Exporters
    "NodeMapping",
    "create_node_mapping",
    "to_adjacency_matrix",
    "to_csv",
    "to_dgl",
    "to_edge_list",
    "to_neo4j_cypher",
    "to_pyg",
    # Schema
    "EDGE_TYPE_METADATA",
    "Edge",
    "EdgeType",
    "GraphSchema",
    "Node",
    "NodeFeatures",
    "NodeType",
]
