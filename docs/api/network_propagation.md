# Network Propagation API

> **Module**: `pathway_subtyping.network_propagation`

Diffusion-based signal refinement on biological networks. Spreads gene-level signals through protein-protein interaction and pathway networks using random walk with restart, heat diffusion, insulated diffusion, or personalized PageRank.

---

## Quick Example

```python
from pathway_subtyping import (
    NetworkPropagator,
    PropagationConfig,
    PropagationMethod,
    propagate_network_scores,
)

edges = [("A", "B", 1.0), ("B", "C", 1.0), ("C", "D", 0.5)]

# One-shot convenience wrapper
result = propagate_network_scores(
    {"A": 1.0},
    edges,
    method=PropagationMethod.RANDOM_WALK,
    restart_prob=0.5,
)
print(result.format_report())

# Or drive the propagator directly for repeated use
propagator = NetworkPropagator(
    PropagationConfig(method=PropagationMethod.HEAT_DIFFUSION, diffusion_time=0.1)
)
propagator.build_network_from_edges(edges)
print(propagator.get_network_stats())
print(propagator.propagate({"A": 1.0}).gene_scores)
```

`NetworkPropagator`, `PropagationConfig`, `PropagationMethod`, `PropagationResult`, and `propagate_network_scores` are re-exported from the top-level `pathway_subtyping` package.

---

## Enums

### `PropagationMethod`

| Value | String | Description |
|-------|--------|-------------|
| `RANDOM_WALK` | `"random_walk"` | Random walk with restart: `p_{t+1} = (1 - alpha) * W * p_t + alpha * p_0` |
| `HEAT_DIFFUSION` | `"heat_diffusion"` | Heat kernel `exp(-t * L)` on the unnormalized Laplacian `L = D - A` |
| `INSULATED` | `"insulated"` | Heat kernel on the normalized Laplacian `I - D^(-1/2) A D^(-1/2)` |
| `PAGERANK` | `"pagerank"` | Personalized PageRank; delegates to the random-walk solver, reporting `method="pagerank"` |

---

## Functions

### `propagate_network_scores()`

```python
def propagate_network_scores(
    gene_scores: Dict[str, float],
    edges: List[Tuple[str, str, float]],
    method: PropagationMethod = PropagationMethod.RANDOM_WALK,
    restart_prob: float = 0.5,
    diffusion_time: float = 0.1,
    seed: Optional[int] = None,
) -> PropagationResult
```

One-shot network propagation (convenience wrapper). Builds a `PropagationConfig`, constructs a `NetworkPropagator`, calls `build_network_from_edges(edges)`, and returns `propagate(gene_scores)`.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_scores` | Dict[str, float] | required | Gene ID to score mapping |
| `edges` | List[Tuple[str, str, float]] | required | `(source, target, weight)` tuples |
| `method` | PropagationMethod | `RANDOM_WALK` | Propagation method |
| `restart_prob` | float | 0.5 | Restart probability for RWR (alpha) |
| `diffusion_time` | float | 0.1 | Diffusion time for heat-kernel methods |
| `seed` | int or None | None | Stored on the config; see note below |

**Returns:** `PropagationResult`.

> Note: the wrapper only forwards `method`, `restart_prob`, `diffusion_time`, and `seed` to `PropagationConfig`. All other config fields keep their defaults. The four propagation solvers are deterministic and no code in this module reads `config.seed`.

---

## Classes

### `NetworkPropagator`

```python
NetworkPropagator(config: Optional[PropagationConfig] = None)
```

Propagates gene-level signals through a biological network. When `config` is `None`, a default `PropagationConfig()` is used.

A network must be built before propagating; calling `propagate()` or `propagate_score_matrix()` first raises `RuntimeError("Network not built. Call build_network() first.")`.

#### `build_network()`

```python
def build_network(
    self,
    knowledge_graph: "KnowledgeGraph",
    edge_types: Optional[List[str]] = None,
) -> None
```

Build the propagation network from a `KnowledgeGraph` (see `pathway_subtyping.knowledge_graph`). Nodes are those whose `node_type` attribute equals `"gene"`; edges are kept when their `edge_type` attribute is in `edge_types`, which defaults to `["gene_interacts_gene"]`.

Edge weights are read from each edge's `weight` attribute (default `1.0`) when `config.use_edge_weights` is `True`, otherwise every kept edge gets weight `1.0`. Edges are inserted symmetrically, self-loops are zeroed out, and the result is stored as a SciPy CSR matrix. If no gene nodes are found, a warning is logged and the method returns without building a matrix.

> Degree filtering: `config.min_degree` / `config.max_degree` are used to count out-of-range nodes and log that count. The current implementation does not remove those nodes from the adjacency matrix.

#### `build_network_from_edges()`

```python
def build_network_from_edges(
    self,
    edges: List[Tuple[str, str, float]],
    nodes: Optional[List[str]] = None,
) -> None
```

Build the network directly from an edge list. When `nodes` is `None`, the node set is inferred from the edges and sorted. Edges referencing unknown nodes are skipped. Edges are inserted symmetrically and self-loops removed. `config.use_edge_weights` is not consulted here — the weight in each tuple is used as given.

#### `propagate()`

```python
def propagate(self, gene_scores: Dict[str, float]) -> PropagationResult
```

Propagate a single score vector through the network, dispatching on `config.method`. Raises `ValueError` for an unknown method and `RuntimeError` if no network has been built.

#### `propagate_score_matrix()`

```python
def propagate_score_matrix(self, score_matrix: pd.DataFrame) -> pd.DataFrame
```

Propagate every row of a samples x genes DataFrame. Per sample, only non-zero entries whose gene is present in the network are used as seeds; samples with no usable seeds yield an all-zero row. Returns a DataFrame with the same index as `score_matrix` and one column per network node.

#### `get_network_stats()`

```python
def get_network_stats(self) -> Dict[str, Any]
```

Returns `{"n_nodes", "n_edges", "avg_degree", "max_degree", "min_degree", "density"}`. If no network has been built, returns `{"error": "Network not built"}`.

---

## Dataclasses

### `PropagationConfig`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `method` | PropagationMethod | `RANDOM_WALK` | Propagation method |
| `restart_prob` | float | 0.5 | RWR restart probability (alpha) |
| `n_iterations` | int | 100 | Maximum RWR iterations |
| `convergence_threshold` | float | 1e-6 | Max absolute change between iterations for convergence |
| `diffusion_time` | float | 0.1 | Time parameter `t` for heat / insulated diffusion |
| `normalize_edges` | bool | True | Row-normalize the adjacency into a transition matrix (RWR / PageRank) |
| `use_edge_weights` | bool | True | Read edge `weight` attributes in `build_network()` |
| `min_degree` | int | 1 | Lower degree bound (counted and logged only — see above) |
| `max_degree` | int | 1000 | Upper degree bound (counted and logged only — see above) |
| `normalize_output` | bool | True | Divide propagated scores by their maximum |
| `preserve_zeros` | bool | False | Drop near-zero scores for genes that were not seeds |
| `seed` | int or None | None | Stored for reproducibility; not read by any solver in this module |

Scores below `1e-10` are omitted from the returned dictionary regardless of `preserve_zeros`.

### `PropagationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `gene_scores` | Dict[str, float] | Propagated score per gene |
| `n_iterations` | int | Iterations run (always 1 for the heat-kernel methods) |
| `converged` | bool | Whether the solver converged (always `True` for the heat-kernel methods) |
| `method` | str | Method string, e.g. `"random_walk"` |
| `metadata` | Dict[str, Any] | Solver parameters — `{"alpha": ...}` for RWR/PageRank, `{"diffusion_time": ...}` for the heat-kernel methods |

**Methods:**

- `to_dict()` — JSON-serializable summary with `method`, `n_iterations`, `converged`, `n_genes_scored`, `top_genes` (top 10, rounded to 6 places), and `metadata`.
- `format_report()` — human-readable report listing the top 10 genes by propagated score.
- `get_citations()` — literature citations for the method used: Kohler et al. 2008 for `"random_walk"` / `"pagerank"`, Vandin et al. 2011 for `"heat_diffusion"` / `"insulated"`.

---

## Notes

- Internal solver helpers (`_random_walk_with_restart`, `_heat_diffusion`, `_insulated_diffusion`, `_personalized_pagerank`, `_normalize_adjacency`, `_vector_to_scores`) are implementation details and are not documented here.
- For RWR and PageRank the seed vector `p_0` is L1-normalized when its sum is positive. The heat-kernel methods use the raw seed values.
- Heat and insulated diffusion are single-shot matrix-exponential applications via `scipy.sparse.linalg.expm_multiply`; `n_iterations` and `convergence_threshold` do not apply to them.

---

## See Also

- [GNN & Graph Embeddings](gnn.md) — graph neural network scoring over a knowledge graph
- [Signaling Databases](signaling_databases.md) — interaction databases usable as edge sources
- [Framework Overview](../framework_overview.md) — architecture and pipeline flow
