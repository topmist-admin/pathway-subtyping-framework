"""
Knowledge Graph Schema

Defines node and edge types for the heterogeneous biological knowledge graph.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Set


class NodeType(Enum):
    """Types of nodes in the knowledge graph."""

    GENE = "gene"
    PATHWAY = "pathway"
    GO_TERM = "go_term"
    CELL_TYPE = "cell_type"
    DRUG = "drug"
    PROTEIN = "protein"
    VARIANT = "variant"
    PHENOTYPE = "phenotype"


class EdgeType(Enum):
    """Types of edges in the knowledge graph."""

    # Gene-gene relationships
    GENE_INTERACTS = "gene_interacts_gene"  # PPI network
    GENE_COEXPRESSED = "gene_coexpressed_gene"  # Co-expression

    # Gene-pathway/GO relationships
    GENE_IN_PATHWAY = "gene_in_pathway"
    GENE_HAS_GO = "gene_has_go"

    # Hierarchical relationships
    PATHWAY_CONTAINS = "pathway_contains_pathway"
    GO_IS_A = "go_is_a"  # GO term hierarchy
    GO_PART_OF = "go_part_of"  # GO term part_of relationship

    # Gene-protein relationships
    GENE_ENCODES = "gene_encodes_protein"
    PROTEIN_INTERACTS = "protein_interacts_protein"

    # Drug relationships
    DRUG_TARGETS = "drug_targets_gene"
    DRUG_IN_PATHWAY = "drug_in_pathway"

    # Variant relationships
    VARIANT_IN_GENE = "variant_in_gene"
    VARIANT_AFFECTS = "variant_affects_protein"

    # Cell type relationships
    GENE_EXPRESSED_IN = "gene_expressed_in_cell_type"

    # Phenotype relationships
    GENE_ASSOCIATED_PHENOTYPE = "gene_associated_phenotype"
    PATHWAY_ASSOCIATED_PHENOTYPE = "pathway_associated_phenotype"


@dataclass
class Node:
    """A node in the knowledge graph."""

    id: str
    node_type: NodeType
    attributes: Dict[str, Any] = field(default_factory=dict)

    def __hash__(self) -> int:
        return hash((self.id, self.node_type))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Node):
            return False
        return self.id == other.id and self.node_type == other.node_type

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable format."""
        return {
            "id": self.id,
            "node_type": self.node_type.value,
            "attributes": self.attributes,
        }


@dataclass
class Edge:
    """An edge in the knowledge graph."""

    source: str
    target: str
    edge_type: EdgeType
    weight: float = 1.0
    attributes: Dict[str, Any] = field(default_factory=dict)

    def __hash__(self) -> int:
        return hash((self.source, self.target, self.edge_type))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Edge):
            return False
        return (
            self.source == other.source
            and self.target == other.target
            and self.edge_type == other.edge_type
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable format."""
        return {
            "source": self.source,
            "target": self.target,
            "edge_type": self.edge_type.value,
            "weight": self.weight,
            "attributes": self.attributes,
        }


@dataclass
class NodeFeatures:
    """Feature vector for a node."""

    node_id: str
    node_type: NodeType
    features: List[float] = field(default_factory=list)
    feature_names: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable format."""
        return {
            "node_id": self.node_id,
            "node_type": self.node_type.value,
            "n_features": len(self.features),
            "feature_names": self.feature_names,
        }


@dataclass
class GraphSchema:
    """
    Schema definition for the knowledge graph.

    Defines allowed node types, edge types, and valid source-target combinations.
    """

    node_types: Set[NodeType] = field(default_factory=set)
    edge_types: Set[EdgeType] = field(default_factory=set)
    valid_edges: Dict[EdgeType, tuple] = field(default_factory=dict)

    @classmethod
    def default_schema(cls) -> "GraphSchema":
        """Create the default schema for biological knowledge graph."""
        schema = cls()

        # Add all node types
        schema.node_types = {nt for nt in NodeType}

        # Add all edge types
        schema.edge_types = {et for et in EdgeType}

        # Define valid source-target combinations for each edge type
        schema.valid_edges = {
            EdgeType.GENE_INTERACTS: (NodeType.GENE, NodeType.GENE),
            EdgeType.GENE_COEXPRESSED: (NodeType.GENE, NodeType.GENE),
            EdgeType.GENE_IN_PATHWAY: (NodeType.GENE, NodeType.PATHWAY),
            EdgeType.GENE_HAS_GO: (NodeType.GENE, NodeType.GO_TERM),
            EdgeType.PATHWAY_CONTAINS: (NodeType.PATHWAY, NodeType.PATHWAY),
            EdgeType.GO_IS_A: (NodeType.GO_TERM, NodeType.GO_TERM),
            EdgeType.GO_PART_OF: (NodeType.GO_TERM, NodeType.GO_TERM),
            EdgeType.GENE_ENCODES: (NodeType.GENE, NodeType.PROTEIN),
            EdgeType.PROTEIN_INTERACTS: (NodeType.PROTEIN, NodeType.PROTEIN),
            EdgeType.DRUG_TARGETS: (NodeType.DRUG, NodeType.GENE),
            EdgeType.DRUG_IN_PATHWAY: (NodeType.DRUG, NodeType.PATHWAY),
            EdgeType.VARIANT_IN_GENE: (NodeType.VARIANT, NodeType.GENE),
            EdgeType.VARIANT_AFFECTS: (NodeType.VARIANT, NodeType.PROTEIN),
            EdgeType.GENE_EXPRESSED_IN: (NodeType.GENE, NodeType.CELL_TYPE),
            EdgeType.GENE_ASSOCIATED_PHENOTYPE: (
                NodeType.GENE,
                NodeType.PHENOTYPE,
            ),
            EdgeType.PATHWAY_ASSOCIATED_PHENOTYPE: (
                NodeType.PATHWAY,
                NodeType.PHENOTYPE,
            ),
        }

        return schema

    def is_valid_edge(
        self,
        source_type: NodeType,
        target_type: NodeType,
        edge_type: EdgeType,
    ) -> bool:
        """Check if an edge type is valid for given source and target node types."""
        if edge_type not in self.valid_edges:
            return False
        expected_source, expected_target = self.valid_edges[edge_type]
        return source_type == expected_source and target_type == expected_target


# Edge type metadata for analysis
EDGE_TYPE_METADATA = {
    EdgeType.GENE_INTERACTS: {
        "symmetric": True,
        "source": "STRING",
        "description": "Protein-protein interaction",
    },
    EdgeType.GENE_COEXPRESSED: {
        "symmetric": True,
        "source": "Co-expression atlas",
        "description": "Co-expression in tissue",
    },
    EdgeType.GENE_IN_PATHWAY: {
        "symmetric": False,
        "source": "Reactome/GO",
        "description": "Gene membership in pathway",
    },
    EdgeType.GENE_HAS_GO: {
        "symmetric": False,
        "source": "GO",
        "description": "Gene annotation with GO term",
    },
    EdgeType.GO_IS_A: {
        "symmetric": False,
        "source": "GO",
        "description": "GO term hierarchy (is_a relationship)",
    },
    EdgeType.GO_PART_OF: {
        "symmetric": False,
        "source": "GO",
        "description": "GO term part_of relationship",
    },
    EdgeType.DRUG_TARGETS: {
        "symmetric": False,
        "source": "DrugBank",
        "description": "Drug-gene targeting relationship",
    },
    EdgeType.GENE_EXPRESSED_IN: {
        "symmetric": False,
        "source": "Single-cell atlas",
        "description": "Gene expression in cell type",
    },
    EdgeType.GENE_ASSOCIATED_PHENOTYPE: {
        "symmetric": False,
        "source": "Disease gene database",
        "description": "Gene-phenotype association",
    },
}
