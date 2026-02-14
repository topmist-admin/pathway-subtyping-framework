"""
Tests for the advanced visualization module.
"""

import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.visualization import (
    DimReductionMethod,
    FigureFormat,
    ReportConfig,
    VisualizationResult,
    compute_dim_reduction,
    create_interactive_report,
    export_figure,
    generate_all_figures,
    plot_cluster_distribution,
    plot_interactive_heatmap,
    plot_interactive_scatter,
    plot_static_scatter,
    plot_subtype_trajectories,
)

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def sample_data():
    """Generate sample pathway scores and labels for testing."""
    np.random.seed(42)
    n_samples = 60
    n_pathways = 8

    pathway_names = [f"pathway_{i}" for i in range(n_pathways)]
    sample_names = [f"sample_{i}" for i in range(n_samples)]

    # Create 3 clusters with distinct pathway patterns
    scores = np.random.randn(n_samples, n_pathways)
    labels = np.array([0] * 20 + [1] * 20 + [2] * 20)

    # Make clusters separable
    scores[:20, :3] += 2.0  # Cluster 0: high on first 3 pathways
    scores[20:40, 3:6] += 2.0  # Cluster 1: high on middle pathways
    scores[40:, 6:] += 2.0  # Cluster 2: high on last 2 pathways

    pathway_scores = pd.DataFrame(scores, index=sample_names, columns=pathway_names)
    return pathway_scores, labels


@pytest.fixture
def output_dir(tmp_path):
    """Create a temporary output directory."""
    out = tmp_path / "viz_output"
    out.mkdir()
    return str(out)


# =============================================================================
# DATA CLASSES
# =============================================================================


class TestVisualizationResult:
    def test_to_dict(self):
        result = VisualizationResult(
            figure_type="scatter",
            output_paths={"png": "/tmp/scatter.png"},
            interactive=True,
            metadata={"n_samples": 100},
        )
        d = result.to_dict()
        assert d["figure_type"] == "scatter"
        assert d["interactive"]
        assert d["output_paths"]["png"] == "/tmp/scatter.png"
        assert d["metadata"]["n_samples"] == 100

    def test_defaults(self):
        result = VisualizationResult(figure_type="heatmap")
        assert result.output_paths == {}
        assert not result.interactive
        assert not result.static_fallback
        assert result.metadata == {}


class TestReportConfig:
    def test_defaults(self):
        config = ReportConfig()
        assert config.title == "Pathway Subtyping Report"
        assert config.include_heatmap
        assert config.include_scatter
        assert config.include_distribution
        assert config.include_trajectories
        assert config.dim_reduction == DimReductionMethod.PCA

    def test_custom_config(self):
        config = ReportConfig(
            title="My Report",
            dim_reduction=DimReductionMethod.TSNE,
            include_trajectories=False,
        )
        assert config.title == "My Report"
        assert config.dim_reduction == DimReductionMethod.TSNE
        assert not config.include_trajectories


# =============================================================================
# ENUMS
# =============================================================================


class TestEnums:
    def test_dim_reduction_values(self):
        assert DimReductionMethod.PCA.value == "pca"
        assert DimReductionMethod.TSNE.value == "tsne"
        assert DimReductionMethod.UMAP.value == "umap"

    def test_figure_format_values(self):
        assert FigureFormat.PNG.value == "png"
        assert FigureFormat.SVG.value == "svg"
        assert FigureFormat.PDF.value == "pdf"
        assert FigureFormat.HTML.value == "html"


# =============================================================================
# DIMENSIONALITY REDUCTION
# =============================================================================


class TestDimReduction:
    def test_pca(self, sample_data):
        scores, labels = sample_data
        embedding, meta = compute_dim_reduction(scores, DimReductionMethod.PCA, seed=42)
        assert embedding.shape == (60, 2)
        assert "explained_variance_ratio" in meta
        assert len(meta["explained_variance_ratio"]) == 2
        assert meta["method"] == "pca"

    def test_tsne(self, sample_data):
        scores, labels = sample_data
        embedding, meta = compute_dim_reduction(scores, DimReductionMethod.TSNE, seed=42)
        assert embedding.shape == (60, 2)
        assert "kl_divergence" in meta
        assert meta["method"] == "tsne"

    def test_umap(self, sample_data):
        scores, labels = sample_data
        embedding, meta = compute_dim_reduction(scores, DimReductionMethod.UMAP, seed=42)
        assert embedding.shape == (60, 2)
        assert "n_neighbors" in meta
        assert meta["method"] == "umap"

    def test_pca_3d(self, sample_data):
        scores, labels = sample_data
        embedding, meta = compute_dim_reduction(
            scores, DimReductionMethod.PCA, n_components=3, seed=42
        )
        assert embedding.shape == (60, 3)
        assert len(meta["explained_variance_ratio"]) == 3

    def test_reproducibility(self, sample_data):
        scores, labels = sample_data
        e1, _ = compute_dim_reduction(scores, DimReductionMethod.PCA, seed=42)
        e2, _ = compute_dim_reduction(scores, DimReductionMethod.PCA, seed=42)
        np.testing.assert_array_equal(e1, e2)

    def test_small_dataset(self):
        """Test with very small dataset (edge case for perplexity/neighbors)."""
        scores = pd.DataFrame(
            np.random.randn(10, 3),
            index=[f"s{i}" for i in range(10)],
            columns=["p1", "p2", "p3"],
        )
        embedding, meta = compute_dim_reduction(scores, DimReductionMethod.TSNE, seed=42)
        assert embedding.shape == (10, 2)

        embedding, meta = compute_dim_reduction(scores, DimReductionMethod.UMAP, seed=42)
        assert embedding.shape == (10, 2)


# =============================================================================
# PLOTLY INTERACTIVE FIGURES
# =============================================================================


class TestInteractiveScatter:
    def test_basic(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        assert fig is not None
        # Plotly figure should have data
        assert len(fig.data) > 0

    def test_pca_axes_labels(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, method=DimReductionMethod.PCA, seed=42)
        x_title = fig.layout.xaxis.title.text
        assert "PC1" in x_title

    def test_tsne_method(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, method=DimReductionMethod.TSNE, seed=42)
        assert fig is not None

    def test_umap_method(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, method=DimReductionMethod.UMAP, seed=42)
        assert fig is not None

    def test_custom_title(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, title="Custom Title", seed=42)
        assert fig.layout.title.text == "Custom Title"

    def test_traces_per_subtype(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        # Should have one trace per unique label
        assert len(fig.data) == 3


class TestInteractiveHeatmap:
    def test_basic(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_heatmap(scores, labels)
        assert fig is not None
        assert len(fig.data) > 0

    def test_dimensions(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_heatmap(scores, labels)
        heatmap_data = fig.data[0]
        # 3 subtypes x 8 pathways
        assert heatmap_data.z.shape == (3, 8)

    def test_custom_title(self, sample_data):
        scores, labels = sample_data
        fig = plot_interactive_heatmap(scores, labels, title="My Heatmap")
        assert fig.layout.title.text == "My Heatmap"


class TestClusterDistribution:
    def test_basic(self, sample_data):
        _, labels = sample_data
        fig = plot_cluster_distribution(labels)
        assert fig is not None
        assert len(fig.data) > 0

    def test_correct_counts(self, sample_data):
        _, labels = sample_data
        fig = plot_cluster_distribution(labels)
        # All bars should sum to total samples
        total = sum(trace.y[0] if hasattr(trace.y, "__getitem__") else 0 for trace in fig.data)
        # Each trace is one bar
        assert total == 60

    def test_custom_title(self, sample_data):
        _, labels = sample_data
        fig = plot_cluster_distribution(labels, title="My Distribution")
        assert fig.layout.title.text == "My Distribution"


class TestSubtypeTrajectories:
    def test_basic(self, sample_data):
        scores, labels = sample_data
        fig = plot_subtype_trajectories(scores, labels)
        assert fig is not None
        # One trace per subtype
        assert len(fig.data) == 3

    def test_top_n_pathways(self, sample_data):
        scores, labels = sample_data
        fig = plot_subtype_trajectories(scores, labels, top_n=5)
        assert fig is not None
        # Radar should use 5 pathways (+1 to close polygon)
        assert len(fig.data[0].theta) == 6

    def test_all_pathways(self, sample_data):
        scores, labels = sample_data
        fig = plot_subtype_trajectories(scores, labels, top_n=100)
        # Should use all 8 pathways (+1 to close)
        assert len(fig.data[0].theta) == 9


# =============================================================================
# STATIC MATPLOTLIB FALLBACK
# =============================================================================


class TestStaticScatter:
    def test_basic(self, sample_data):
        scores, labels = sample_data
        fig = plot_static_scatter(scores, labels, seed=42)
        import matplotlib.pyplot as plt

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_save_file(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "scatter.png")
        fig = plot_static_scatter(scores, labels, seed=42, output_path=path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_tsne_method(self, sample_data):
        scores, labels = sample_data
        fig = plot_static_scatter(scores, labels, method=DimReductionMethod.TSNE, seed=42)
        import matplotlib.pyplot as plt

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_umap_method(self, sample_data):
        scores, labels = sample_data
        fig = plot_static_scatter(scores, labels, method=DimReductionMethod.UMAP, seed=42)
        import matplotlib.pyplot as plt

        assert isinstance(fig, plt.Figure)
        plt.close(fig)


# =============================================================================
# MULTI-FORMAT EXPORT
# =============================================================================


class TestExportFigure:
    def test_export_plotly_html(self, sample_data, output_dir):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        saved = export_figure(fig, os.path.join(output_dir, "scatter"), [FigureFormat.HTML])
        assert "html" in saved
        assert os.path.exists(saved["html"])
        # HTML should be a valid file with plotly content
        content = Path(saved["html"]).read_text()
        assert "plotly" in content.lower()

    def test_export_plotly_png(self, sample_data, output_dir):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        saved = export_figure(fig, os.path.join(output_dir, "scatter"), [FigureFormat.PNG])
        assert "png" in saved
        assert os.path.exists(saved["png"])
        assert os.path.getsize(saved["png"]) > 0

    def test_export_plotly_svg(self, sample_data, output_dir):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        saved = export_figure(fig, os.path.join(output_dir, "scatter"), [FigureFormat.SVG])
        assert "svg" in saved
        content = Path(saved["svg"]).read_text()
        assert "<svg" in content

    def test_export_matplotlib_png(self, sample_data, output_dir):
        scores, labels = sample_data
        fig = plot_static_scatter(scores, labels, seed=42)
        saved = export_figure(fig, os.path.join(output_dir, "scatter"), [FigureFormat.PNG])
        assert "png" in saved
        assert os.path.exists(saved["png"])
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_export_multiple_formats(self, sample_data, output_dir):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        saved = export_figure(
            fig,
            os.path.join(output_dir, "scatter"),
            [FigureFormat.HTML, FigureFormat.PNG, FigureFormat.SVG],
        )
        assert len(saved) == 3
        for fmt in ["html", "png", "svg"]:
            assert fmt in saved
            assert os.path.exists(saved[fmt])

    def test_default_format_is_png(self, sample_data, output_dir):
        scores, labels = sample_data
        fig = plot_interactive_scatter(scores, labels, seed=42)
        saved = export_figure(fig, os.path.join(output_dir, "scatter"))
        assert "png" in saved


# =============================================================================
# INTERACTIVE HTML REPORT
# =============================================================================


class TestInteractiveReport:
    def test_basic_report(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        result = create_interactive_report(scores, labels, path, seed=42)
        assert os.path.exists(path)
        assert result.interactive
        assert result.figure_type == "interactive_report"
        assert "html" in result.output_paths

    def test_report_content(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        create_interactive_report(scores, labels, path, seed=42)
        content = Path(path).read_text()
        assert "Pathway Subtyping Report" in content
        assert "plotly" in content.lower()
        assert "Analysis Summary" in content

    def test_report_metadata(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        result = create_interactive_report(scores, labels, path, seed=42)
        assert result.metadata["n_samples"] == 60
        assert result.metadata["n_pathways"] == 8
        assert result.metadata["n_subtypes"] == 3
        assert "scatter" in result.metadata["figures_included"]
        assert "heatmap" in result.metadata["figures_included"]

    def test_custom_config(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        config = ReportConfig(
            title="My Custom Report",
            include_trajectories=False,
            dim_reduction=DimReductionMethod.TSNE,
        )
        result = create_interactive_report(scores, labels, path, config=config, seed=42)
        content = Path(path).read_text()
        assert "My Custom Report" in content
        assert "trajectories" not in result.metadata["figures_included"]

    def test_report_self_contained(self, sample_data, output_dir):
        """Report should include plotly.js CDN link."""
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        create_interactive_report(scores, labels, path, seed=42)
        content = Path(path).read_text()
        assert "cdn.plot.ly" in content

    def test_report_disclaimer(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        config = ReportConfig(disclaimer="Custom disclaimer text")
        create_interactive_report(scores, labels, path, config=config, seed=42)
        content = Path(path).read_text()
        assert "Custom disclaimer text" in content

    def test_creates_parent_dirs(self, sample_data, tmp_path):
        scores, labels = sample_data
        path = str(tmp_path / "nested" / "dir" / "report.html")
        result = create_interactive_report(scores, labels, path, seed=42)
        assert os.path.exists(path)

    def test_report_with_umap(self, sample_data, output_dir):
        scores, labels = sample_data
        path = os.path.join(output_dir, "report.html")
        config = ReportConfig(dim_reduction=DimReductionMethod.UMAP)
        result = create_interactive_report(scores, labels, path, config=config, seed=42)
        assert result.metadata["dim_reduction"] == "umap"


# =============================================================================
# GENERATE ALL FIGURES
# =============================================================================


class TestGenerateAllFigures:
    def test_basic(self, sample_data, output_dir):
        scores, labels = sample_data
        results = generate_all_figures(
            scores, labels, output_dir, formats=[FigureFormat.PNG], seed=42
        )
        # Should generate scatter + heatmap + distribution + trajectories
        assert len(results) == 4
        types = [r.figure_type for r in results]
        assert "scatter" in types
        assert "heatmap" in types
        assert "distribution" in types
        assert "trajectories" in types

    def test_html_export(self, sample_data, output_dir):
        scores, labels = sample_data
        results = generate_all_figures(
            scores, labels, output_dir, formats=[FigureFormat.HTML], seed=42
        )
        for r in results:
            assert "html" in r.output_paths
            assert os.path.exists(r.output_paths["html"])

    def test_creates_output_dir(self, sample_data, tmp_path):
        scores, labels = sample_data
        new_dir = str(tmp_path / "new_output")
        results = generate_all_figures(scores, labels, new_dir, seed=42)
        assert os.path.isdir(new_dir)
        assert len(results) > 0


# =============================================================================
# EDGE CASES
# =============================================================================


class TestEdgeCases:
    def test_single_cluster(self):
        """All samples in one cluster."""
        scores = pd.DataFrame(
            np.random.randn(20, 5),
            index=[f"s{i}" for i in range(20)],
            columns=[f"p{i}" for i in range(5)],
        )
        labels = np.zeros(20, dtype=int)
        fig = plot_interactive_scatter(scores, labels, seed=42)
        assert fig is not None
        assert len(fig.data) == 1  # One trace for one cluster

    def test_two_samples(self):
        """Minimum viable dataset."""
        scores = pd.DataFrame(
            np.random.randn(2, 3),
            index=["s0", "s1"],
            columns=["p0", "p1", "p2"],
        )
        labels = np.array([0, 1])
        fig = plot_interactive_scatter(scores, labels, seed=42)
        assert fig is not None

    def test_many_pathways(self):
        """Test with many pathways (trajectory top_n selection)."""
        scores = pd.DataFrame(
            np.random.randn(30, 50),
            index=[f"s{i}" for i in range(30)],
            columns=[f"pathway_{i}" for i in range(50)],
        )
        labels = np.array([0] * 10 + [1] * 10 + [2] * 10)
        fig = plot_subtype_trajectories(scores, labels, top_n=8)
        # Should only use 8 pathways (+1 to close)
        assert len(fig.data[0].theta) == 9

    def test_single_pathway(self):
        """Test with only one pathway."""
        scores = pd.DataFrame(
            np.random.randn(20, 1),
            index=[f"s{i}" for i in range(20)],
            columns=["only_pathway"],
        )
        labels = np.array([0] * 10 + [1] * 10)
        fig = plot_interactive_heatmap(scores, labels)
        assert fig is not None


class TestImportExports:
    """Test that all public functions are importable from the package."""

    def test_import_from_package(self):
        from pathway_subtyping import (
            DimReductionMethod,
            FigureFormat,
            ReportConfig,
            VisualizationResult,
            compute_dim_reduction,
            create_interactive_report,
            export_figure,
            generate_all_figures,
            plot_cluster_distribution,
            plot_interactive_heatmap,
            plot_interactive_scatter,
            plot_static_scatter,
            plot_subtype_trajectories,
        )

        assert DimReductionMethod.PCA.value == "pca"
        assert FigureFormat.PNG.value == "png"
