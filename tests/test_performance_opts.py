"""
Tests for performance optimization parameters (Issue #8).

Verifies that show_progress=False suppresses tqdm output without errors,
and that PipelineConfig accepts chunked processing options.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.expression import (
    ExpressionScoringMethod,
    score_pathways_from_expression,
)
from pathway_subtyping.pipeline import PipelineConfig
from pathway_subtyping.sensitivity import vary_feature_subset
from pathway_subtyping.simulation import (
    SimulationConfig,
    estimate_type_i_error,
    generate_synthetic_data,
    run_power_analysis,
)
from pathway_subtyping.validation import ValidationGates

# -- Shared fixtures ----------------------------------------------------------


@pytest.fixture
def small_synthetic():
    """Generate a small synthetic dataset for quick tests."""
    config = SimulationConfig(
        n_samples=30,
        n_pathways=5,
        n_subtypes=2,
        effect_size=2.0,
        seed=42,
    )
    data = generate_synthetic_data(config)

    from sklearn.mixture import GaussianMixture

    gmm = GaussianMixture(n_components=2, random_state=42)
    labels = gmm.fit_predict(data.pathway_scores.values)
    return data, labels


# -- ValidationGates ----------------------------------------------------------


class TestValidationGatesProgress:
    """Tests for show_progress parameter on ValidationGates."""

    def test_show_progress_false_runs(self, small_synthetic):
        """ValidationGates(show_progress=False) completes without error."""
        data, labels = small_synthetic
        vg = ValidationGates(
            seed=42,
            n_permutations=3,
            n_bootstrap=3,
            show_progress=False,
        )
        result = vg.run_all(
            pathway_scores=data.pathway_scores,
            cluster_labels=labels,
            pathways=data.pathways,
            gene_burdens=data.gene_burdens,
            n_clusters=2,
            gmm_seed=42,
        )
        assert result.summary is not None
        assert len(result.results) > 0

    def test_show_progress_default_is_true(self):
        """show_progress defaults to True."""
        vg = ValidationGates(seed=1)
        assert vg.show_progress is True

    def test_show_progress_stored(self):
        """show_progress value is stored on the instance."""
        vg = ValidationGates(seed=1, show_progress=False)
        assert vg.show_progress is False


# -- Simulation functions -----------------------------------------------------


class TestSimulationProgress:
    """Tests for show_progress on simulation analysis functions."""

    def test_type_i_error_no_progress(self):
        """estimate_type_i_error runs with show_progress=False."""
        result = estimate_type_i_error(
            n_simulations=2,
            n_samples=20,
            n_pathways=3,
            seed=42,
            show_progress=False,
        )
        assert 0.0 <= result.type_i_rate <= 1.0

    def test_power_analysis_no_progress(self):
        """run_power_analysis runs with show_progress=False."""
        result = run_power_analysis(
            n_samples=20,
            n_pathways=3,
            n_subtypes=2,
            effect_sizes=[0.5, 1.0],
            n_simulations_per_effect=2,
            seed=42,
            show_progress=False,
        )
        assert len(result.effect_sizes) == 2
        assert len(result.power_at_threshold) == 2


# -- Expression scoring -------------------------------------------------------


class TestExpressionProgress:
    """Tests for show_progress on expression pathway scoring."""

    def test_ssgsea_no_progress(self):
        """score_pathways_from_expression runs with show_progress=False."""
        rng = np.random.default_rng(42)
        genes = [f"GENE{i}" for i in range(20)]
        samples = [f"S{i}" for i in range(10)]
        expr = pd.DataFrame(
            rng.standard_normal((10, 20)),
            index=samples,
            columns=genes,
        )
        pathways = {
            "pathway_A": genes[:8],
            "pathway_B": genes[8:16],
        }
        result = score_pathways_from_expression(
            gene_expression=expr,
            pathways=pathways,
            method=ExpressionScoringMethod.MEAN_Z,
            show_progress=False,
        )
        assert result.pathway_scores.shape[0] == 10
        assert result.pathway_scores.shape[1] == 2


# -- Sensitivity analysis ------------------------------------------------------


class TestSensitivityProgress:
    """Tests for show_progress on sensitivity analysis."""

    def test_vary_feature_subset_no_progress(self, small_synthetic):
        """vary_feature_subset runs with show_progress=False."""
        data, _ = small_synthetic
        result = vary_feature_subset(
            pathway_scores=data.pathway_scores,
            n_clusters=2,
            seed=42,
            show_progress=False,
        )
        assert result.parameter is not None
        # configurations = all_features baseline + one per pathway LOO
        assert len(result.configurations) == data.pathway_scores.shape[1] + 1


# -- PipelineConfig chunked processing ----------------------------------------


class TestPipelineConfigChunked:
    """Tests for chunked processing config options."""

    def test_default_chunked_off(self):
        """use_chunked_processing defaults to False."""
        config = PipelineConfig(name="test")
        assert config.use_chunked_processing is False
        assert config.chunk_size == 1000

    def test_chunked_enabled(self):
        """PipelineConfig accepts use_chunked_processing=True."""
        config = PipelineConfig(
            name="test",
            use_chunked_processing=True,
            chunk_size=500,
        )
        assert config.use_chunked_processing is True
        assert config.chunk_size == 500
