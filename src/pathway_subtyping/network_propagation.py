"""
Network Propagation

Implements diffusion-based signal refinement on biological networks.
Spreads gene-level signals through protein-protein interaction and
pathway networks to improve signal detection.

Supports four propagation methods:
- Random Walk with Restart (RWR)
- Heat Diffusion (Laplacian kernel)
- Insulated Diffusion (normalized Laplacian, boundary-preserving)
- Personalized PageRank
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

if TYPE_CHECKING:
    from pathway_subtyping.knowledge_graph import KnowledgeGraph

logger = logging.getLogger(__name__)


class PropagationMethod(Enum):
    """Network propagation methods."""

    RANDOM_WALK = "random_walk"
    HEAT_DIFFUSION = "heat_diffusion"
    INSULATED = "insulated"
    PAGERANK = "pagerank"


@dataclass
class PropagationConfig:
    """Configuration for network propagation."""

    method: PropagationMethod = PropagationMethod.RANDOM_WALK

    # Random walk parameters
    restart_prob: float = 0.5
    n_iterations: int = 100
    convergence_threshold: float = 1e-6

    # Heat diffusion parameters
    diffusion_time: float = 0.1

    # Edge weight handling
    normalize_edges: bool = True
    use_edge_weights: bool = True

    # Node filtering
    min_degree: int = 1
    max_degree: int = 1000

    # Output processing
    normalize_output: bool = True
    preserve_zeros: bool = False

    # Reproducibility
    seed: Optional[int] = None


@dataclass
class PropagationResult:
    """Result of network propagation."""

    gene_scores: Dict[str, float]
    n_iterations: int
    converged: bool
    method: str
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable format."""
        top_genes = dict(sorted(self.gene_scores.items(), key=lambda x: -x[1])[:10])
        return {
            "method": self.method,
            "n_iterations": self.n_iterations,
            "converged": self.converged,
            "n_genes_scored": len(self.gene_scores),
            "top_genes": {k: round(v, 6) for k, v in top_genes.items()},
            "metadata": self.metadata,
        }

    def format_report(self) -> str:
        """Generate human-readable report."""
        lines = [
            "Network Propagation Result",
            "=" * 40,
            f"Method: {self.method}",
            f"Converged: {self.converged} ({self.n_iterations} iterations)",
            f"Genes scored: {len(self.gene_scores)}",
            "",
            "Top 10 genes by propagated score:",
        ]
        for gene, score in sorted(self.gene_scores.items(), key=lambda x: -x[1])[:10]:
            lines.append(f"  {gene}: {score:.4f}")
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return literature citations for methods used."""
        citations: List[str] = []
        if self.method in ("random_walk", "pagerank"):
            citations.append(
                "Kohler S, et al. Walking the interactome for "
                "prioritization of candidate disease genes. "
                "Am J Hum Genet. 2008;82(4):949-58."
            )
        if self.method in ("heat_diffusion", "insulated"):
            citations.append(
                "Vandin F, Upfal E, Raphael BJ. Algorithms for "
                "detecting significantly mutated pathways in cancer. "
                "J Comput Biol. 2011;18(3):507-22."
            )
        return citations


class NetworkPropagator:
    """
    Propagates gene-level signals through biological networks.

    Uses random walk with restart (RWR) or heat diffusion to spread
    signals from seed genes to their network neighbors.
    """

    def __init__(self, config: Optional[PropagationConfig] = None):
        """
        Initialize network propagator.

        Args:
            config: Propagation configuration
        """
        self.config = config or PropagationConfig()
        self._adj_matrix: Optional[sparse.csr_matrix] = None
        self._node_to_idx: Dict[str, int] = {}
        self._idx_to_node: Dict[int, str] = {}

    def build_network(
        self,
        knowledge_graph: "KnowledgeGraph",
        edge_types: Optional[List[str]] = None,
    ) -> None:
        """
        Build propagation network from knowledge graph.

        Args:
            knowledge_graph: KnowledgeGraph instance
            edge_types: Edge types to include (default: gene_interacts_gene)
        """
        if edge_types is None:
            edge_types = ["gene_interacts_gene"]

        # Get gene nodes via duck-typing
        gene_nodes = [
            n
            for n in knowledge_graph.graph.nodes()
            if knowledge_graph.graph.nodes[n].get("node_type") == "gene"
        ]

        if not gene_nodes:
            logger.warning("[NetworkPropagation] No gene nodes found in knowledge graph")
            return

        # Build node index
        self._node_to_idx = {node: idx for idx, node in enumerate(gene_nodes)}
        self._idx_to_node = {idx: node for node, idx in self._node_to_idx.items()}
        n_nodes = len(gene_nodes)

        logger.info(
            "[NetworkPropagation] Building network with %d gene nodes",
            n_nodes,
        )

        # Build adjacency matrix
        row_indices: List[int] = []
        col_indices: List[int] = []
        weights: List[float] = []

        for source, target, data in knowledge_graph.graph.edges(data=True):
            edge_type = data.get("edge_type", "")
            if edge_type not in edge_types:
                continue

            if source not in self._node_to_idx or target not in self._node_to_idx:
                continue

            src_idx = self._node_to_idx[source]
            tgt_idx = self._node_to_idx[target]
            weight = data.get("weight", 1.0) if self.config.use_edge_weights else 1.0

            # Add edge in both directions (symmetric)
            row_indices.extend([src_idx, tgt_idx])
            col_indices.extend([tgt_idx, src_idx])
            weights.extend([weight, weight])

        self._adj_matrix = sparse.csr_matrix(
            (weights, (row_indices, col_indices)),
            shape=(n_nodes, n_nodes),
        )

        # Remove self-loops
        self._adj_matrix.setdiag(0)
        self._adj_matrix.eliminate_zeros()

        # Apply degree filtering
        degrees = np.array(self._adj_matrix.sum(axis=1)).flatten()
        valid_nodes = (degrees >= self.config.min_degree) & (degrees <= self.config.max_degree)

        if not np.all(valid_nodes):
            n_filtered = int(np.sum(~valid_nodes))
            logger.info(
                "[NetworkPropagation] Filtered %d nodes by degree",
                n_filtered,
            )

        logger.info(
            "[NetworkPropagation] Built network: %d nodes, %d edges",
            n_nodes,
            self._adj_matrix.nnz // 2,
        )

    def build_network_from_edges(
        self,
        edges: List[Tuple[str, str, float]],
        nodes: Optional[List[str]] = None,
    ) -> None:
        """
        Build propagation network from edge list.

        Args:
            edges: List of (source, target, weight) tuples
            nodes: Optional list of nodes (inferred from edges if None)
        """
        if nodes is None:
            node_set: Set[str] = set()
            for src, tgt, _ in edges:
                node_set.add(src)
                node_set.add(tgt)
            nodes = sorted(node_set)

        self._node_to_idx = {node: idx for idx, node in enumerate(nodes)}
        self._idx_to_node = {idx: node for node, idx in self._node_to_idx.items()}
        n_nodes = len(nodes)

        row_indices: List[int] = []
        col_indices: List[int] = []
        weights: List[float] = []

        for src, tgt, weight in edges:
            if src not in self._node_to_idx or tgt not in self._node_to_idx:
                continue

            src_idx = self._node_to_idx[src]
            tgt_idx = self._node_to_idx[tgt]

            row_indices.extend([src_idx, tgt_idx])
            col_indices.extend([tgt_idx, src_idx])
            weights.extend([weight, weight])

        self._adj_matrix = sparse.csr_matrix(
            (weights, (row_indices, col_indices)),
            shape=(n_nodes, n_nodes),
        )

        self._adj_matrix.setdiag(0)
        self._adj_matrix.eliminate_zeros()

        logger.info(
            "[NetworkPropagation] Built network: %d nodes, %d edges",
            n_nodes,
            len(edges),
        )

    def propagate(
        self,
        gene_scores: Dict[str, float],
    ) -> PropagationResult:
        """
        Propagate gene scores through the network.

        Args:
            gene_scores: Dictionary mapping gene ID to score

        Returns:
            PropagationResult with propagated scores
        """
        if self._adj_matrix is None:
            raise RuntimeError("Network not built. Call build_network() first.")

        method = self.config.method

        if method == PropagationMethod.RANDOM_WALK:
            return self._random_walk_with_restart(gene_scores)
        elif method == PropagationMethod.HEAT_DIFFUSION:
            return self._heat_diffusion(gene_scores)
        elif method == PropagationMethod.INSULATED:
            return self._insulated_diffusion(gene_scores)
        elif method == PropagationMethod.PAGERANK:
            return self._personalized_pagerank(gene_scores)
        else:
            raise ValueError(f"Unknown propagation method: {method}")

    def propagate_score_matrix(
        self,
        score_matrix: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Propagate all samples' gene scores through the network.

        Args:
            score_matrix: DataFrame with samples as rows, genes as columns

        Returns:
            DataFrame with propagated scores (same index, expanded columns)
        """
        if self._adj_matrix is None:
            raise RuntimeError("Network not built. Call build_network() first.")

        logger.info(
            "[NetworkPropagation] Propagating scores for %d samples",
            len(score_matrix),
        )

        network_genes = list(self._node_to_idx.keys())
        result_data = np.zeros((len(score_matrix), len(network_genes)))

        for s_idx, (sample_id, row) in enumerate(score_matrix.iterrows()):
            sample_scores = {
                gene: score
                for gene, score in row.items()
                if score != 0 and gene in self._node_to_idx
            }

            if not sample_scores:
                continue

            prop_result = self.propagate(sample_scores)

            for gene, score in prop_result.gene_scores.items():
                if gene in self._node_to_idx:
                    g_idx = network_genes.index(gene)
                    result_data[s_idx, g_idx] = score

        return pd.DataFrame(
            result_data,
            index=score_matrix.index,
            columns=network_genes,
        )

    def _random_walk_with_restart(
        self,
        gene_scores: Dict[str, float],
    ) -> PropagationResult:
        """
        Random walk with restart propagation.

        p_{t+1} = (1 - alpha) * W * p_t + alpha * p_0
        """
        n_nodes = self._adj_matrix.shape[0]
        alpha = self.config.restart_prob

        W = self._normalize_adjacency(self._adj_matrix)

        p0 = np.zeros(n_nodes)
        for gene, score in gene_scores.items():
            if gene in self._node_to_idx:
                p0[self._node_to_idx[gene]] = score

        if np.sum(p0) > 0:
            p0 = p0 / np.sum(p0)

        p = p0.copy()
        converged = False
        iteration = 0

        # W = D^-1 A is ROW-stochastic, so the walk on a column probability
        # vector must apply its TRANSPOSE: mass leaves node j split by j's own
        # out-degree. Using W.dot(p) instead makes each node divide its INCOMING
        # mass by its own degree, which does not conserve probability (a 5-node
        # test summed to 0.719 instead of 1.0) and systematically SUPPRESSES hubs
        # rather than boosting them -- the opposite of what propagation is for.
        # Verified against the closed form alpha*inv(I-(1-alpha)W^T)p0 and
        # networkx.pagerank; pinned by test_rwr_matches_closed_form_and_networkx.
        W_t = W.T.tocsr() if sparse.issparse(W) else W.T

        for iteration in range(self.config.n_iterations):
            p_new = (1 - alpha) * W_t.dot(p) + alpha * p0

            diff = np.max(np.abs(p_new - p))
            if diff < self.config.convergence_threshold:
                converged = True
                p = p_new
                break

            p = p_new

        if not converged:
            logger.warning(
                "[NetworkPropagation] random-walk did NOT converge after %d "
                "iterations (last delta %.3g > tol %.3g). The returned scores are "
                "un-converged; raise n_iterations or relax convergence_threshold. "
                "Check `result.converged` before using these values.",
                self.config.n_iterations,
                float(diff),
                self.config.convergence_threshold,
            )

        propagated = self._vector_to_scores(p, gene_scores)

        return PropagationResult(
            gene_scores=propagated,
            n_iterations=iteration + 1,
            converged=converged,
            method="random_walk",
            metadata={"alpha": alpha},
        )

    def _heat_diffusion(
        self,
        gene_scores: Dict[str, float],
    ) -> PropagationResult:
        """
        Heat diffusion kernel propagation.

        p = exp(-t * L) * p_0
        where L is the Laplacian matrix.
        """
        n_nodes = self._adj_matrix.shape[0]
        t = self.config.diffusion_time

        degrees = np.array(self._adj_matrix.sum(axis=1)).flatten()
        D = sparse.diags(degrees)
        L = D - self._adj_matrix

        p0 = np.zeros(n_nodes)
        for gene, score in gene_scores.items():
            if gene in self._node_to_idx:
                p0[self._node_to_idx[gene]] = score

        p = expm_multiply(-t * L, p0)

        propagated = self._vector_to_scores(p, gene_scores)

        return PropagationResult(
            gene_scores=propagated,
            n_iterations=1,
            converged=True,
            method="heat_diffusion",
            metadata={"diffusion_time": t},
        )

    def _insulated_diffusion(
        self,
        gene_scores: Dict[str, float],
    ) -> PropagationResult:
        """
        Insulated heat diffusion (preserves signal at boundaries).

        Uses normalized Laplacian: L_norm = I - D^(-1/2) * A * D^(-1/2)
        """
        n_nodes = self._adj_matrix.shape[0]
        t = self.config.diffusion_time

        degrees = np.array(self._adj_matrix.sum(axis=1)).flatten()
        degrees[degrees == 0] = 1  # Avoid division by zero

        D_inv_sqrt = sparse.diags(1.0 / np.sqrt(degrees))
        L_norm = sparse.eye(n_nodes) - D_inv_sqrt @ self._adj_matrix @ D_inv_sqrt

        p0 = np.zeros(n_nodes)
        for gene, score in gene_scores.items():
            if gene in self._node_to_idx:
                p0[self._node_to_idx[gene]] = score

        p = expm_multiply(-t * L_norm, p0)

        propagated = self._vector_to_scores(p, gene_scores)

        return PropagationResult(
            gene_scores=propagated,
            n_iterations=1,
            converged=True,
            method="insulated",
            metadata={"diffusion_time": t},
        )

    def _personalized_pagerank(
        self,
        gene_scores: Dict[str, float],
    ) -> PropagationResult:
        """
        Personalized PageRank propagation.

        Equivalent to RWR with alpha as damping factor.
        """
        result = self._random_walk_with_restart(gene_scores)
        return PropagationResult(
            gene_scores=result.gene_scores,
            n_iterations=result.n_iterations,
            converged=result.converged,
            method="pagerank",
            metadata=result.metadata,
        )

    def _normalize_adjacency(
        self,
        adj: sparse.csr_matrix,
    ) -> sparse.csr_matrix:
        """Row-normalize adjacency matrix to create transition matrix."""
        row_sums = np.array(adj.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1

        if self.config.normalize_edges:
            D_inv = sparse.diags(1.0 / row_sums)
            return D_inv @ adj
        else:
            return adj

    def _vector_to_scores(
        self,
        p: np.ndarray,
        original_scores: Dict[str, float],
    ) -> Dict[str, float]:
        """Convert propagated vector back to gene score dictionary."""
        if self.config.normalize_output and np.max(p) > 0:
            p = p / np.max(p)

        result: Dict[str, float] = {}
        for idx, score in enumerate(p):
            gene = self._idx_to_node.get(idx)
            if gene is None:
                continue

            if self.config.preserve_zeros:
                if gene not in original_scores or original_scores[gene] == 0:
                    if score < 1e-10:
                        continue

            if score > 1e-10:
                result[gene] = float(score)

        return result

    def get_network_stats(self) -> Dict[str, Any]:
        """Get statistics about the propagation network."""
        if self._adj_matrix is None:
            return {"error": "Network not built"}

        degrees = np.array(self._adj_matrix.sum(axis=1)).flatten()

        return {
            "n_nodes": self._adj_matrix.shape[0],
            "n_edges": self._adj_matrix.nnz // 2,
            "avg_degree": float(np.mean(degrees)),
            "max_degree": int(np.max(degrees)),
            "min_degree": (int(np.min(degrees[degrees > 0])) if np.any(degrees > 0) else 0),
            "density": float(self._adj_matrix.nnz / (self._adj_matrix.shape[0] ** 2)),
        }


def propagate_network_scores(
    gene_scores: Dict[str, float],
    edges: List[Tuple[str, str, float]],
    method: PropagationMethod = PropagationMethod.RANDOM_WALK,
    restart_prob: float = 0.5,
    diffusion_time: float = 0.1,
    seed: Optional[int] = None,
) -> PropagationResult:
    """
    One-shot network propagation (convenience wrapper).

    Args:
        gene_scores: Dictionary mapping gene ID to score
        edges: List of (source, target, weight) tuples
        method: Propagation method
        restart_prob: Restart probability for RWR (alpha)
        diffusion_time: Diffusion time for heat kernel methods
        seed: Random seed for reproducibility

    Returns:
        PropagationResult with propagated scores
    """
    config = PropagationConfig(
        method=method,
        restart_prob=restart_prob,
        diffusion_time=diffusion_time,
        seed=seed,
    )
    propagator = NetworkPropagator(config)
    propagator.build_network_from_edges(edges)
    return propagator.propagate(gene_scores)
