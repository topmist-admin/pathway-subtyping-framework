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
from .diff import KGDiff, diff_kgs
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
from .regression import KGRegressionReport, ScoreDelta, run_kg_regression
from .schema import (
    EDGE_TYPE_METADATA,
    Edge,
    EdgeType,
    GraphSchema,
    Node,
    NodeFeatures,
    NodeType,
)
from .sources import (
    KG_SOURCES,
    KG_SOURCES_V05,
    KG_SOURCES_V06,
    KGSource,
    get_source,
    list_sources,
    manifest_digest,
)

__all__ = [
    # Builder
    "KnowledgeGraph",
    "KnowledgeGraphBuilder",
    "KnowledgeGraphStats",
    "PPINetwork",
    "load_ppi_from_file",
    # Diff
    "KGDiff",
    "diff_kgs",
    # Exporters
    "NodeMapping",
    "create_node_mapping",
    "to_adjacency_matrix",
    "to_csv",
    "to_dgl",
    "to_edge_list",
    "to_neo4j_cypher",
    "to_pyg",
    # Regression
    "KGRegressionReport",
    "ScoreDelta",
    "run_kg_regression",
    # Schema
    "EDGE_TYPE_METADATA",
    "Edge",
    "EdgeType",
    "GraphSchema",
    "Node",
    "NodeFeatures",
    "NodeType",
    # Sources
    "KG_SOURCES",
    "KG_SOURCES_V05",
    "KG_SOURCES_V06",
    "KGSource",
    "get_source",
    "list_sources",
    "manifest_digest",
]
