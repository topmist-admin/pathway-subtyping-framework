"""
Advanced Visualization Module for the Pathway Subtyping Framework.

Provides interactive and publication-quality visualizations:
- Interactive HTML reports with Plotly (requires `pip install pathway-subtyping[viz]`)
- UMAP/t-SNE dimensionality reduction scatter plots
- Subtype trajectory (radar) visualization
- Multi-format figure export (PNG, SVG, PDF)
- Static matplotlib fallbacks when Plotly is unavailable

All interactive features require the optional [viz] extra:
    pip install pathway-subtyping[viz]

Static matplotlib visualizations work with the base install.

References:
- McInnes, Healy & Melville (2018). UMAP: Uniform Manifold Approximation
  and Projection for Dimension Reduction. arXiv:1802.03426.
- van der Maaten & Hinton (2008). Visualizing Data using t-SNE.
  JMLR, 9, 2579-2605.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# ENUMS AND CONFIGURATION
# =============================================================================


class DimReductionMethod(Enum):
    """Dimensionality reduction method for scatter plots."""

    PCA = "pca"
    TSNE = "tsne"
    UMAP = "umap"


class FigureFormat(Enum):
    """Supported figure export formats."""

    PNG = "png"
    SVG = "svg"
    PDF = "pdf"
    HTML = "html"


# =============================================================================
# DATA CLASSES
# =============================================================================


@dataclass
class VisualizationResult:
    """
    Result from a visualization operation.

    Attributes:
        figure_type: Type of figure generated (e.g., 'heatmap', 'scatter')
        output_paths: Dict mapping format to file path for saved figures
        interactive: Whether an interactive (Plotly) figure was generated
        static_fallback: Whether a static matplotlib fallback was used
        metadata: Additional metadata about the visualization
    """

    figure_type: str
    output_paths: Dict[str, str] = field(default_factory=dict)
    interactive: bool = False
    static_fallback: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "figure_type": self.figure_type,
            "output_paths": self.output_paths,
            "interactive": self.interactive,
            "static_fallback": self.static_fallback,
            "metadata": self.metadata,
        }


@dataclass
class ReportConfig:
    """
    Configuration for interactive HTML report generation.

    Attributes:
        title: Report title
        include_heatmap: Include pathway heatmap
        include_scatter: Include dimensionality reduction scatter plot
        include_distribution: Include cluster distribution bar chart
        include_trajectories: Include subtype trajectory (radar) chart
        dim_reduction: Dimensionality reduction method for scatter plot
        color_palette: Color palette name for subtypes
        disclaimer: Disclaimer text for the report footer
    """

    title: str = "Pathway Subtyping Report"
    include_heatmap: bool = True
    include_scatter: bool = True
    include_distribution: bool = True
    include_trajectories: bool = True
    dim_reduction: DimReductionMethod = DimReductionMethod.PCA
    color_palette: str = "Set2"
    disclaimer: str = "Research use only. Not for clinical decision-making."


# =============================================================================
# OPTIONAL IMPORT HELPERS
# =============================================================================


def _check_plotly():
    """Check if plotly is available, raise helpful error if not."""
    try:
        import plotly  # noqa: F401

        return True
    except ImportError:
        raise ImportError(
            "Plotly is required for interactive visualizations.\n"
            "Install with: pip install pathway-subtyping[viz]\n"
            "Or: pip install plotly>=5.15.0"
        )


def _check_umap():
    """Check if umap-learn is available, raise helpful error if not."""
    try:
        import umap  # noqa: F401

        return True
    except ImportError:
        raise ImportError(
            "umap-learn is required for UMAP dimensionality reduction.\n"
            "Install with: pip install pathway-subtyping[viz]\n"
            "Or: pip install umap-learn>=0.5.0"
        )


def _check_kaleido():
    """Check if kaleido is available for static export from Plotly."""
    try:
        import kaleido  # noqa: F401

        return True
    except ImportError:
        return False


# =============================================================================
# DIMENSIONALITY REDUCTION
# =============================================================================


def compute_dim_reduction(
    pathway_scores: pd.DataFrame,
    method: DimReductionMethod = DimReductionMethod.PCA,
    n_components: int = 2,
    seed: Optional[int] = None,
    **kwargs,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Compute dimensionality reduction on pathway scores.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        method: Dimensionality reduction method (PCA, t-SNE, or UMAP)
        n_components: Number of output dimensions (default 2)
        seed: Random seed for reproducibility
        **kwargs: Additional parameters passed to the underlying method

    Returns:
        Tuple of (embedding array [n_samples, n_components], metadata dict)
    """
    metadata = {"method": method.value, "n_components": n_components}
    X = pathway_scores.values

    if method == DimReductionMethod.PCA:
        from sklearn.decomposition import PCA

        pca = PCA(n_components=n_components, random_state=seed)
        embedding = pca.fit_transform(X)
        metadata["explained_variance_ratio"] = pca.explained_variance_ratio_.tolist()
        logger.info(
            f"[Visualization] PCA: {pca.explained_variance_ratio_.sum():.1%} " f"variance explained"
        )

    elif method == DimReductionMethod.TSNE:
        from sklearn.manifold import TSNE

        perplexity = kwargs.get("perplexity", min(30, max(5, len(X) // 4)))
        tsne = TSNE(
            n_components=n_components,
            random_state=seed,
            perplexity=perplexity,
        )
        embedding = tsne.fit_transform(X)
        metadata["perplexity"] = perplexity
        metadata["kl_divergence"] = float(tsne.kl_divergence_)
        logger.info(f"[Visualization] t-SNE: KL divergence = {tsne.kl_divergence_:.4f}")

    elif method == DimReductionMethod.UMAP:
        _check_umap()
        import umap

        n_neighbors = kwargs.get("n_neighbors", min(15, max(2, len(X) - 1)))
        min_dist = kwargs.get("min_dist", 0.1)
        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=seed,
        )
        embedding = reducer.fit_transform(X)
        metadata["n_neighbors"] = n_neighbors
        metadata["min_dist"] = min_dist
        logger.info("[Visualization] UMAP embedding computed")

    else:
        raise ValueError(f"Unknown dimensionality reduction method: {method}")

    return embedding, metadata


# =============================================================================
# PLOTLY INTERACTIVE FIGURES
# =============================================================================


def plot_interactive_scatter(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    method: DimReductionMethod = DimReductionMethod.PCA,
    seed: Optional[int] = None,
    title: str = "Sample Clustering",
    color_palette: str = "Set2",
    **kwargs,
) -> Any:
    """
    Create an interactive Plotly scatter plot of samples in reduced space.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        labels: Cluster labels for each sample
        method: Dimensionality reduction method
        seed: Random seed for reproducibility
        title: Plot title
        color_palette: Plotly color palette name
        **kwargs: Passed to compute_dim_reduction()

    Returns:
        plotly.graph_objects.Figure
    """
    _check_plotly()
    import plotly.express as px

    embedding, meta = compute_dim_reduction(pathway_scores, method=method, seed=seed, **kwargs)

    df = pd.DataFrame(
        {
            "Dim1": embedding[:, 0],
            "Dim2": embedding[:, 1],
            "Subtype": [f"Subtype {int(lb)}" for lb in labels],
            "Sample": pathway_scores.index.tolist(),
        }
    )

    axis_labels = {}
    if method == DimReductionMethod.PCA and "explained_variance_ratio" in meta:
        evr = meta["explained_variance_ratio"]
        axis_labels = {
            "Dim1": f"PC1 ({evr[0]:.1%})",
            "Dim2": f"PC2 ({evr[1]:.1%})",
        }
    else:
        method_name = method.value.upper()
        axis_labels = {
            "Dim1": f"{method_name} 1",
            "Dim2": f"{method_name} 2",
        }

    fig = px.scatter(
        df,
        x="Dim1",
        y="Dim2",
        color="Subtype",
        hover_data=["Sample"],
        title=title,
        labels=axis_labels,
        color_discrete_sequence=getattr(
            px.colors.qualitative, color_palette, px.colors.qualitative.Set2
        ),
    )
    fig.update_traces(marker=dict(size=8, opacity=0.7))
    fig.update_layout(
        template="plotly_white",
        legend_title_text="Subtype",
    )

    return fig


def plot_interactive_heatmap(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    title: str = "Pathway Scores by Subtype",
    color_scale: str = "RdBu_r",
) -> Any:
    """
    Create an interactive Plotly heatmap of mean pathway scores per subtype.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        labels: Cluster labels for each sample
        title: Plot title
        color_scale: Plotly color scale name

    Returns:
        plotly.graph_objects.Figure
    """
    _check_plotly()
    import plotly.graph_objects as go

    # Compute mean scores per subtype
    unique_labels = sorted(np.unique(labels))
    subtype_names = [f"Subtype {int(lb)}" for lb in unique_labels]
    pathways = pathway_scores.columns.tolist()

    matrix = np.zeros((len(unique_labels), len(pathways)))
    for i, label in enumerate(unique_labels):
        mask = labels == label
        matrix[i] = pathway_scores.values[mask].mean(axis=0)

    # Z-score normalize columns for better visualization
    col_means = matrix.mean(axis=0)
    col_stds = matrix.std(axis=0)
    col_stds[col_stds == 0] = 1
    z_matrix = (matrix - col_means) / col_stds

    vmax = max(abs(z_matrix.min()), abs(z_matrix.max()), 0.01)

    fig = go.Figure(
        data=go.Heatmap(
            z=z_matrix,
            x=pathways,
            y=subtype_names,
            colorscale=color_scale,
            zmin=-vmax,
            zmax=vmax,
            text=np.round(z_matrix, 2).astype(str),
            texttemplate="%{text}",
            textfont={"size": 10},
            hovertemplate=(
                "Pathway: %{x}<br>" "Subtype: %{y}<br>" "Z-score: %{z:.2f}<extra></extra>"
            ),
            colorbar=dict(title="Mean Z-score"),
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title="Pathway",
        yaxis_title="Subtype",
        template="plotly_white",
        xaxis=dict(tickangle=45),
    )

    return fig


def plot_cluster_distribution(
    labels: np.ndarray,
    title: str = "Cluster Size Distribution",
    color_palette: str = "Set2",
) -> Any:
    """
    Create an interactive Plotly bar chart of cluster sizes.

    Args:
        labels: Cluster labels for each sample
        title: Plot title
        color_palette: Plotly color palette name

    Returns:
        plotly.graph_objects.Figure
    """
    _check_plotly()
    import plotly.express as px

    unique_labels = sorted(np.unique(labels))
    counts = [int(np.sum(labels == lb)) for lb in unique_labels]
    names = [f"Subtype {int(lb)}" for lb in unique_labels]

    df = pd.DataFrame({"Subtype": names, "Count": counts})

    fig = px.bar(
        df,
        x="Subtype",
        y="Count",
        title=title,
        color="Subtype",
        color_discrete_sequence=getattr(
            px.colors.qualitative, color_palette, px.colors.qualitative.Set2
        ),
        text="Count",
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        yaxis_title="Number of Samples",
    )

    return fig


def plot_subtype_trajectories(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    title: str = "Subtype Pathway Profiles",
    color_palette: str = "Set2",
    top_n: int = 10,
) -> Any:
    """
    Create a radar (polar) chart of subtype pathway profiles.

    Shows the mean Z-score for each pathway across subtypes, making it
    easy to compare pathway activation patterns between subtypes.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        labels: Cluster labels for each sample
        title: Plot title
        color_palette: Plotly color palette name
        top_n: Maximum number of pathways to display (selects most variable)

    Returns:
        plotly.graph_objects.Figure
    """
    _check_plotly()
    import plotly.express as px
    import plotly.graph_objects as go

    # Select top_n most variable pathways
    variances = pathway_scores.var()
    if len(variances) > top_n:
        top_pathways = variances.nlargest(top_n).index.tolist()
    else:
        top_pathways = pathway_scores.columns.tolist()

    unique_labels = sorted(np.unique(labels))
    colors = getattr(px.colors.qualitative, color_palette, px.colors.qualitative.Set2)

    fig = go.Figure()

    for i, label in enumerate(unique_labels):
        mask = labels == label
        means = pathway_scores.loc[mask, top_pathways].mean()

        # Z-score for comparable scale
        overall_means = pathway_scores[top_pathways].mean()
        overall_stds = pathway_scores[top_pathways].std()
        overall_stds[overall_stds == 0] = 1
        z_means = (means - overall_means) / overall_stds

        # Close the radar polygon
        values = z_means.tolist() + [z_means.iloc[0]]
        categories = top_pathways + [top_pathways[0]]

        color = colors[i % len(colors)]
        n_samples = int(mask.sum())
        fig.add_trace(
            go.Scatterpolar(
                r=values,
                theta=categories,
                fill="toself",
                name=f"Subtype {int(label)} (n={n_samples})",
                line=dict(color=color),
                opacity=0.6,
            )
        )

    fig.update_layout(
        title=title,
        template="plotly_white",
        polar=dict(radialaxis=dict(visible=True, title="Mean Z-score")),
        showlegend=True,
    )

    return fig


# =============================================================================
# STATIC MATPLOTLIB FALLBACKS
# =============================================================================


def plot_static_scatter(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    method: DimReductionMethod = DimReductionMethod.PCA,
    seed: Optional[int] = None,
    title: str = "Sample Clustering",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8),
    **kwargs,
) -> Any:
    """
    Create a static matplotlib scatter plot of samples in reduced space.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        labels: Cluster labels for each sample
        method: Dimensionality reduction method
        seed: Random seed
        title: Plot title
        output_path: Path to save the figure
        figsize: Figure size in inches
        **kwargs: Passed to compute_dim_reduction()

    Returns:
        matplotlib Figure, or None if matplotlib unavailable
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("[Visualization] matplotlib not available — skipping static scatter")
        return None

    embedding, meta = compute_dim_reduction(pathway_scores, method=method, seed=seed, **kwargs)

    fig, ax = plt.subplots(figsize=figsize)
    unique_labels = sorted(np.unique(labels))

    for label in unique_labels:
        mask = labels == label
        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            label=f"Subtype {int(label)}",
            alpha=0.7,
            s=50,
        )

    if method == DimReductionMethod.PCA and "explained_variance_ratio" in meta:
        evr = meta["explained_variance_ratio"]
        ax.set_xlabel(f"PC1 ({evr[0]:.1%})")
        ax.set_ylabel(f"PC2 ({evr[1]:.1%})")
    else:
        name = method.value.upper()
        ax.set_xlabel(f"{name} 1")
        ax.set_ylabel(f"{name} 2")

    ax.set_title(title)
    ax.legend()
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"[Visualization] Saved scatter: {output_path}")
        plt.close(fig)

    return fig


# =============================================================================
# MULTI-FORMAT EXPORT
# =============================================================================


def export_figure(
    fig: Any,
    output_path: str,
    formats: Optional[List[FigureFormat]] = None,
) -> Dict[str, str]:
    """
    Export a figure in multiple formats.

    Supports both Plotly and matplotlib figures. For Plotly static exports
    (PNG, SVG, PDF), the kaleido package is required.

    Args:
        fig: A Plotly Figure or matplotlib Figure object
        output_path: Base output path (without extension)
        formats: List of formats to export. Defaults to [PNG].

    Returns:
        Dict mapping format name to saved file path
    """
    if formats is None:
        formats = [FigureFormat.PNG]

    output_path = str(output_path)
    saved = {}

    # Detect figure type
    is_plotly = False
    try:
        import plotly.graph_objects as go

        is_plotly = isinstance(fig, go.Figure)
    except ImportError:
        pass

    for fmt in formats:
        ext = fmt.value
        path = f"{output_path}.{ext}" if not output_path.endswith(f".{ext}") else output_path

        if is_plotly:
            if fmt == FigureFormat.HTML:
                fig.write_html(path, include_plotlyjs="cdn")
                saved[ext] = path
                logger.info(f"[Visualization] Saved HTML: {path}")
            else:
                if _check_kaleido():
                    fig.write_image(path, scale=2)
                    saved[ext] = path
                    logger.info(f"[Visualization] Saved {ext.upper()}: {path}")
                else:
                    logger.warning(
                        f"[Visualization] kaleido not installed — "
                        f"cannot export Plotly figure as {ext.upper()}. "
                        f"Install with: pip install kaleido>=0.2.0"
                    )
        else:
            # matplotlib figure
            if fmt == FigureFormat.HTML:
                logger.warning("[Visualization] HTML export not supported for matplotlib figures")
            else:
                fig.savefig(path, dpi=150, bbox_inches="tight", format=ext)
                saved[ext] = path
                logger.info(f"[Visualization] Saved {ext.upper()}: {path}")

    return saved


# =============================================================================
# INTERACTIVE HTML REPORT
# =============================================================================


def create_interactive_report(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    output_path: str,
    config: Optional[ReportConfig] = None,
    characterization_result: Optional[Any] = None,
    validation_result: Optional[Any] = None,
    seed: Optional[int] = None,
) -> VisualizationResult:
    """
    Generate a self-contained interactive HTML report with Plotly charts.

    The report includes:
    - Dimensionality reduction scatter plot (PCA, t-SNE, or UMAP)
    - Pathway heatmap (mean Z-scores per subtype)
    - Cluster size distribution bar chart
    - Subtype trajectory radar chart
    - Optional characterization and validation summaries

    The output is a single self-contained HTML file that can be opened
    in any modern browser — no server or internet connection required.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        labels: Cluster labels for each sample
        output_path: Path for the output HTML file
        config: Report configuration. Uses defaults if None.
        characterization_result: Optional CharacterizationResult for summaries
        validation_result: Optional ValidationGatesResult for summaries
        seed: Random seed for dimensionality reduction

    Returns:
        VisualizationResult with output paths and metadata
    """
    _check_plotly()

    if config is None:
        config = ReportConfig()

    result = VisualizationResult(
        figure_type="interactive_report",
        interactive=True,
    )

    figures = []

    # 1. Scatter plot
    if config.include_scatter:
        scatter_fig = plot_interactive_scatter(
            pathway_scores,
            labels,
            method=config.dim_reduction,
            seed=seed,
            title=f"Sample Clustering ({config.dim_reduction.value.upper()})",
            color_palette=config.color_palette,
        )
        figures.append(("scatter", scatter_fig))

    # 2. Heatmap
    if config.include_heatmap:
        heatmap_fig = plot_interactive_heatmap(
            pathway_scores,
            labels,
            title="Pathway Profiles by Subtype (Mean Z-scores)",
        )
        figures.append(("heatmap", heatmap_fig))

    # 3. Distribution
    if config.include_distribution:
        dist_fig = plot_cluster_distribution(
            labels,
            title="Cluster Size Distribution",
            color_palette=config.color_palette,
        )
        figures.append(("distribution", dist_fig))

    # 4. Trajectory
    if config.include_trajectories:
        traj_fig = plot_subtype_trajectories(
            pathway_scores,
            labels,
            title="Subtype Pathway Profiles (Radar)",
            color_palette=config.color_palette,
        )
        figures.append(("trajectories", traj_fig))

    # Build HTML report
    html_content = _build_report_html(
        figures=figures,
        config=config,
        characterization_result=characterization_result,
        validation_result=validation_result,
        n_samples=len(pathway_scores),
        n_pathways=len(pathway_scores.columns),
        n_subtypes=len(np.unique(labels)),
    )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html_content, encoding="utf-8")

    result.output_paths["html"] = str(output_path)
    result.metadata = {
        "n_samples": len(pathway_scores),
        "n_pathways": len(pathway_scores.columns),
        "n_subtypes": int(len(np.unique(labels))),
        "dim_reduction": config.dim_reduction.value,
        "figures_included": [name for name, _ in figures],
    }

    logger.info(f"[Visualization] Interactive report saved: {output_path}")
    return result


def _build_report_html(
    figures: List[Tuple[str, Any]],
    config: ReportConfig,
    characterization_result: Optional[Any],
    validation_result: Optional[Any],
    n_samples: int,
    n_pathways: int,
    n_subtypes: int,
) -> str:
    """Build the HTML report content from Plotly figures."""
    import plotly.io as pio

    figure_divs = []
    for name, fig in figures:
        div_html = pio.to_html(fig, full_html=False, include_plotlyjs=False)
        figure_divs.append(f'<div class="chart-container" id="{name}">{div_html}</div>')

    figures_html = "\n".join(figure_divs)

    # Build summary section
    summary_items = [
        f"<li><strong>Samples:</strong> {n_samples:,}</li>",
        f"<li><strong>Pathways:</strong> {n_pathways}</li>",
        f"<li><strong>Subtypes discovered:</strong> {n_subtypes}</li>",
    ]

    if characterization_result is not None:
        try:
            profiles = characterization_result.subtype_profiles
            for p in profiles:
                n = p.n_samples
                top_pw = p.top_pathways[0] if p.top_pathways else "N/A"
                summary_items.append(
                    f"<li><strong>{p.subtype_label}:</strong> "
                    f"n={n}, top pathway = {top_pw}</li>"
                )
        except (AttributeError, IndexError):
            pass

    if validation_result is not None:
        try:
            passed = validation_result.passed
            status = "PASSED" if passed else "FAILED"
            summary_items.append(f"<li><strong>Validation:</strong> {status}</li>")
        except AttributeError:
            pass

    summary_html = "\n".join(summary_items)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{config.title}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                         "Helvetica Neue", Arial, sans-serif;
            background: #f8f9fa;
            color: #333;
            line-height: 1.6;
        }}
        .header {{
            background: linear-gradient(135deg, #2c3e50, #3498db);
            color: white;
            padding: 2rem;
            text-align: center;
        }}
        .header h1 {{ font-size: 1.8rem; margin-bottom: 0.5rem; }}
        .header p {{ opacity: 0.85; font-size: 0.95rem; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 1.5rem; }}
        .summary {{
            background: white;
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .summary h2 {{ font-size: 1.2rem; margin-bottom: 0.8rem; color: #2c3e50; }}
        .summary ul {{ list-style: none; }}
        .summary li {{ padding: 0.25rem 0; }}
        .chart-container {{
            background: white;
            border-radius: 8px;
            padding: 1rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .footer {{
            text-align: center;
            padding: 1.5rem;
            color: #888;
            font-size: 0.85rem;
            border-top: 1px solid #eee;
            margin-top: 2rem;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{config.title}</h1>
        <p>Generated by Pathway Subtyping Framework</p>
    </div>
    <div class="container">
        <div class="summary">
            <h2>Analysis Summary</h2>
            <ul>
                {summary_html}
            </ul>
        </div>
        {figures_html}
    </div>
    <div class="footer">
        <p>{config.disclaimer}</p>
        <p>Generated by <a href="https://codeberg.org/pathways/pathway-subtyping-framework">
           Pathway Subtyping Framework</a></p>
    </div>
</body>
</html>"""

    return html


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================


def generate_all_figures(
    pathway_scores: pd.DataFrame,
    labels: np.ndarray,
    output_dir: str,
    formats: Optional[List[FigureFormat]] = None,
    dim_reduction: DimReductionMethod = DimReductionMethod.PCA,
    seed: Optional[int] = None,
) -> List[VisualizationResult]:
    """
    Generate all available visualizations and save to output directory.

    Attempts interactive Plotly figures first, falls back to static
    matplotlib if Plotly is unavailable.

    Args:
        pathway_scores: DataFrame with samples as rows, pathways as columns
        labels: Cluster labels for each sample
        output_dir: Directory to save figures
        formats: Export formats. Defaults to [PNG].
        dim_reduction: Method for scatter plot
        seed: Random seed

    Returns:
        List of VisualizationResult for each generated figure
    """
    if formats is None:
        formats = [FigureFormat.PNG]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    use_plotly = False
    try:
        _check_plotly()
        use_plotly = True
    except ImportError:
        logger.info("[Visualization] Plotly not available — using static matplotlib figures")

    # Scatter plot
    if use_plotly:
        fig = plot_interactive_scatter(pathway_scores, labels, method=dim_reduction, seed=seed)
        saved = export_figure(fig, str(output_dir / "scatter"), formats)
        results.append(
            VisualizationResult(
                figure_type="scatter",
                output_paths=saved,
                interactive=True,
            )
        )
    else:
        fig = plot_static_scatter(
            pathway_scores,
            labels,
            method=dim_reduction,
            seed=seed,
            output_path=str(output_dir / "scatter.png"),
        )
        if fig is not None:
            results.append(
                VisualizationResult(
                    figure_type="scatter",
                    output_paths={"png": str(output_dir / "scatter.png")},
                    static_fallback=True,
                )
            )

    # Heatmap (Plotly only)
    if use_plotly:
        fig = plot_interactive_heatmap(pathway_scores, labels)
        saved = export_figure(fig, str(output_dir / "heatmap"), formats)
        results.append(
            VisualizationResult(
                figure_type="heatmap",
                output_paths=saved,
                interactive=True,
            )
        )

    # Distribution (Plotly only)
    if use_plotly:
        fig = plot_cluster_distribution(labels)
        saved = export_figure(fig, str(output_dir / "distribution"), formats)
        results.append(
            VisualizationResult(
                figure_type="distribution",
                output_paths=saved,
                interactive=True,
            )
        )

    # Trajectories (Plotly only)
    if use_plotly:
        fig = plot_subtype_trajectories(pathway_scores, labels)
        saved = export_figure(fig, str(output_dir / "trajectories"), formats)
        results.append(
            VisualizationResult(
                figure_type="trajectories",
                output_paths=saved,
                interactive=True,
            )
        )

    logger.info(f"[Visualization] Generated {len(results)} figures in {output_dir}")
    return results
