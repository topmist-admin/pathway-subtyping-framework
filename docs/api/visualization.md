# Visualization API

> **Module**: `pathway_subtyping.visualization`

Generates interactive and publication-quality visualizations for pathway subtyping results. Interactive features require the `[viz]` extra (`pip install pathway-subtyping[viz]`). Static matplotlib fallbacks work with the base install.

---

## Quick Example

```python
from pathway_subtyping import (
    create_interactive_report,
    ReportConfig,
    DimReductionMethod,
)

# Generate a self-contained HTML report
result = create_interactive_report(
    pathway_scores=pathway_scores_df,
    labels=cluster_labels,
    output_path="outputs/report.html",
    config=ReportConfig(
        title="My Cohort Analysis",
        dim_reduction=DimReductionMethod.UMAP,
    ),
    seed=42,
)
# Open outputs/report.html in any browser — no server needed
```

---

## Enums

### `DimReductionMethod`

| Value | Description | Notes |
|-------|-------------|-------|
| `PCA` | Principal Component Analysis | Always available, fast |
| `TSNE` | t-distributed Stochastic Neighbor Embedding | Good for local structure, slower |
| `UMAP` | Uniform Manifold Approximation and Projection | Requires `[viz]` extra, best for large datasets |

### `FigureFormat`

| Value | Description | Interactive |
|-------|-------------|-------------|
| `PNG` | Raster image | No |
| `SVG` | Vector graphic | No |
| `PDF` | Portable document | No |
| `HTML` | Interactive Plotly | Yes |

---

## Functions

### `create_interactive_report()`

```python
def create_interactive_report(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    output_path: str,
    config: Optional[ReportConfig] = None,
    pathways: Optional[Dict[str, List[str]]] = None,
    seed: Optional[int] = None,
) -> VisualizationResult
```

Generate a self-contained interactive HTML report with scatter plot, heatmap, distribution chart, and radar chart.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_scores` | pd.DataFrame | required | Samples x pathways matrix |
| `labels` | np.ndarray | required | Cluster label per sample |
| `output_path` | str | required | Path for HTML file |
| `config` | ReportConfig or None | None | Report configuration |
| `pathways` | Dict or None | None | Pathway definitions (for enrichment display) |
| `seed` | int or None | None | Random seed |

---

### `compute_dim_reduction()`

```python
def compute_dim_reduction(
    pathway_scores: pd.DataFrame,
    method: DimReductionMethod = DimReductionMethod.PCA,
    n_components: int = 2,
    seed: Optional[int] = None,
) -> np.ndarray
```

Compute dimensionality reduction for visualization. Returns (n_samples, n_components) array.

---

### `plot_interactive_scatter()`

```python
def plot_interactive_scatter(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    method: DimReductionMethod = DimReductionMethod.PCA,
    title: str = "Subtype Scatter Plot",
    seed: Optional[int] = None,
) -> Any  # Plotly Figure
```

Interactive Plotly scatter plot with hover information showing sample ID and subtype.

---

### `plot_interactive_heatmap()`

```python
def plot_interactive_heatmap(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    title: str = "Pathway Heatmap",
    top_n: Optional[int] = None,
) -> Any  # Plotly Figure
```

Interactive heatmap of mean pathway Z-scores per subtype. Use `top_n` to limit to the most variable pathways.

---

### `plot_cluster_distribution()`

```python
def plot_cluster_distribution(
    labels: np.ndarray,
    title: str = "Cluster Distribution",
) -> Any  # Plotly Figure
```

Bar chart of cluster sizes.

---

### `plot_subtype_trajectories()`

```python
def plot_subtype_trajectories(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    top_n: int = 8,
    title: str = "Subtype Trajectories",
) -> Any  # Plotly Figure
```

Radar chart of subtype pathway profiles showing the `top_n` most variable pathways.

---

### `plot_static_scatter()`

```python
def plot_static_scatter(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    method: DimReductionMethod = DimReductionMethod.PCA,
    output_path: Optional[str] = None,
    seed: Optional[int] = None,
) -> Any  # matplotlib Figure
```

Static matplotlib scatter plot. Works with base install (no Plotly required).

---

### `export_figure()`

```python
def export_figure(
    figure: Any,
    output_path: str,
    formats: Optional[List[FigureFormat]] = None,
) -> Dict[str, str]
```

Export a Plotly or matplotlib figure in multiple formats.

**Returns:** Dict mapping format to file path.

---

### `generate_all_figures()`

```python
def generate_all_figures(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    output_dir: str,
    config: Optional[ReportConfig] = None,
    seed: Optional[int] = None,
) -> List[VisualizationResult]
```

Generate all visualizations (scatter, heatmap, distribution, radar) and save to `output_dir`.

---

## Dataclasses

### `ReportConfig`

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `title` | str | `"Pathway Subtyping Report"` | Report title |
| `dim_reduction` | DimReductionMethod | `PCA` | Dimensionality reduction method |
| `include_heatmap` | bool | True | Include pathway heatmap |
| `include_distribution` | bool | True | Include cluster distribution |
| `include_trajectories` | bool | True | Include radar chart |
| `top_n_pathways` | int or None | None | Limit heatmap/radar to top N pathways |
| `disclaimer` | str | research disclaimer | Footer disclaimer text |

### `VisualizationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `figure_type` | str | Type of figure (e.g., `"scatter"`, `"heatmap"`) |
| `output_paths` | Dict[str, str] | Format to file path mapping |
| `interactive` | bool | Whether Plotly figure was generated |
| `static_fallback` | bool | Whether matplotlib fallback was used |
| `metadata` | Dict[str, Any] | Additional metadata |

**Methods:** `to_dict()`

---

## Installation

```bash
# Base install (static matplotlib only)
pip install pathway-subtyping

# With interactive visualizations
pip install pathway-subtyping[viz]
# Installs: plotly, umap-learn, kaleido (for PNG/SVG/PDF export)
```

---

## See Also

- [Framework Overview](../framework_overview.md) — pipeline architecture
- [Performance & Hardware](../guides/performance-and-hardware.md) — memory considerations for large datasets
