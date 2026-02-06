"""
Tests for the simulation framework module.

Tests synthetic data generation, ground truth recovery,
type I error estimation, and power analysis.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.simulation import (
    PowerAnalysisResult,
    RecoveryResult,
    SampleSizeAnalysisResult,
    SimulatedData,
    SimulationConfig,
    TypeIErrorResult,
    estimate_type_i_error,
    evaluate_recovery,
    generate_synthetic_data,
    run_gmm_clustering,
    run_power_analysis,
    run_sample_size_analysis,
    validate_framework,
)


class TestSimulationConfig:
    """Tests for SimulationConfig."""

    def test_default_config(self):
        """Test default configuration values."""
        config = SimulationConfig()

        assert config.n_samples == 200
        assert config.n_subtypes == 3
        assert config.effect_size == 1.0
        assert len(config.subtype_proportions) == 3

    def test_equal_proportions(self):
        """Test that proportions sum to 1."""
        config = SimulationConfig(n_subtypes=4)

        assert len(config.subtype_proportions) == 4
        assert abs(sum(config.subtype_proportions) - 1.0) < 0.001

    def test_custom_proportions(self):
        """Test custom subtype proportions."""
        config = SimulationConfig(
            n_subtypes=3,
            subtype_proportions=[0.5, 0.3, 0.2],
        )

        assert abs(sum(config.subtype_proportions) - 1.0) < 0.001
        assert config.subtype_proportions[0] == 0.5


class TestGenerateSyntheticData:
    """Tests for synthetic data generation."""

    def test_basic_generation(self):
        """Test basic data generation."""
        config = SimulationConfig(
            n_samples=100,
            n_pathways=10,
            n_subtypes=3,
            seed=42,
        )

        data = generate_synthetic_data(config)

        assert isinstance(data, SimulatedData)
        assert len(data.true_labels) == 100
        assert len(data.pathways) == 10
        assert data.pathway_scores.shape[0] == 100

    def test_correct_subtype_counts(self):
        """Test that subtype counts match proportions."""
        config = SimulationConfig(
            n_samples=90,
            n_subtypes=3,
            seed=42,
        )

        data = generate_synthetic_data(config)

        unique, counts = np.unique(data.true_labels, return_counts=True)
        assert len(unique) == 3
        # Should be approximately equal (90/3 = 30 each)
        assert all(c > 20 and c < 40 for c in counts)

    def test_effect_size_impact(self):
        """Test that higher effect size creates more separation between subtypes."""
        from sklearn.metrics import silhouette_score

        # Low effect
        low_config = SimulationConfig(
            n_samples=100,
            n_subtypes=2,
            effect_size=0.1,
            seed=42,
        )
        low_data = generate_synthetic_data(low_config)

        # High effect
        high_config = SimulationConfig(
            n_samples=100,
            n_subtypes=2,
            effect_size=2.0,
            seed=42,
        )
        high_data = generate_synthetic_data(high_config)

        # High effect should have better cluster separation (higher silhouette)
        # Note: Z-score normalization sets variance=1, so we measure
        # true label separation instead of raw variance
        low_silhouette = silhouette_score(low_data.pathway_scores.values, low_data.true_labels)
        high_silhouette = silhouette_score(high_data.pathway_scores.values, high_data.true_labels)

        assert high_silhouette > low_silhouette

    def test_reproducibility(self):
        """Test that same seed produces same data."""
        config = SimulationConfig(n_samples=50, seed=123)

        data1 = generate_synthetic_data(config)
        data2 = generate_synthetic_data(config)

        np.testing.assert_array_equal(data1.true_labels, data2.true_labels)
        np.testing.assert_array_almost_equal(
            data1.pathway_scores.values,
            data2.pathway_scores.values,
        )

    def test_pathway_effects_structure(self):
        """Test that subtype_pathway_effects is populated."""
        config = SimulationConfig(n_samples=60, n_subtypes=3)
        data = generate_synthetic_data(config)

        assert len(data.subtype_pathway_effects) == 3
        for subtype, pathways in data.subtype_pathway_effects.items():
            assert len(pathways) > 0


class TestEvaluateRecovery:
    """Tests for ground truth recovery evaluation."""

    def test_perfect_recovery(self):
        """Test perfect label recovery."""
        true_labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        pred_labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])

        result = evaluate_recovery(pred_labels, true_labels)

        assert isinstance(result, RecoveryResult)
        assert result.ari == 1.0
        assert result.nmi == 1.0
        assert result.correct_k

    def test_random_labels(self):
        """Test with random predicted labels."""
        np.random.seed(42)
        true_labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        pred_labels = np.random.randint(0, 3, 9)

        result = evaluate_recovery(pred_labels, true_labels)

        # ARI should be low for random labels
        assert result.ari < 0.5

    def test_label_permutation(self):
        """Test that label permutation still gives high ARI."""
        true_labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        # Same structure, different label values
        pred_labels = np.array([2, 2, 2, 0, 0, 0, 1, 1, 1])

        result = evaluate_recovery(pred_labels, true_labels)

        assert result.ari == 1.0
        assert result.correct_k


class TestGMMClustering:
    """Tests for GMM clustering helper."""

    def test_basic_clustering(self):
        """Test basic clustering functionality."""
        np.random.seed(42)
        data = pd.DataFrame(np.random.randn(50, 5))

        labels = run_gmm_clustering(data, n_clusters=3, seed=42)

        assert len(labels) == 50
        assert len(np.unique(labels)) == 3

    def test_reproducibility(self):
        """Test clustering reproducibility with seed."""
        np.random.seed(42)
        data = pd.DataFrame(np.random.randn(50, 5))

        labels1 = run_gmm_clustering(data, n_clusters=3, seed=123)
        labels2 = run_gmm_clustering(data, n_clusters=3, seed=123)

        np.testing.assert_array_equal(labels1, labels2)


class TestTypeIError:
    """Tests for Type I error estimation."""

    def test_basic_estimation(self):
        """Test Type I error estimation runs."""
        result = estimate_type_i_error(
            n_samples=50,
            n_pathways=5,
            n_clusters=3,
            n_simulations=10,
            seed=42,
        )

        assert isinstance(result, TypeIErrorResult)
        assert result.n_simulations == 10
        assert 0 <= result.type_i_rate <= 1

    def test_null_ari_near_zero(self):
        """Test that null ARI is close to zero."""
        result = estimate_type_i_error(
            n_samples=50,
            n_pathways=5,
            n_clusters=3,
            n_simulations=20,
            seed=42,
        )

        # Mean ARI under null should be close to 0
        assert abs(result.null_ari_mean) < 0.3

    def test_result_to_dict(self):
        """Test serialization to dict."""
        result = estimate_type_i_error(
            n_samples=30,
            n_simulations=5,
            seed=42,
        )

        d = result.to_dict()

        assert "null_ari_mean" in d
        assert "type_i_error_rate" in d
        assert "n_simulations" in d


class TestPowerAnalysis:
    """Tests for power analysis."""

    def test_basic_power_analysis(self):
        """Test power analysis runs."""
        result = run_power_analysis(
            n_samples=50,
            n_pathways=5,
            n_subtypes=2,
            effect_sizes=[0.5, 1.0, 2.0],
            n_simulations_per_effect=5,
            seed=42,
        )

        assert isinstance(result, PowerAnalysisResult)
        assert len(result.effect_sizes) == 3
        assert len(result.power_at_threshold) == 3

    def test_power_increases_with_effect(self):
        """Test that power increases with effect size."""
        result = run_power_analysis(
            n_samples=100,
            n_pathways=10,
            n_subtypes=2,
            effect_sizes=[0.1, 1.0, 3.0],
            n_simulations_per_effect=10,
            seed=42,
        )

        # Power should generally increase with effect size
        # (allowing for simulation variance)
        assert result.power_at_threshold[3.0] >= result.power_at_threshold[0.1]

    def test_result_to_dict(self):
        """Test serialization."""
        result = run_power_analysis(
            n_samples=30,
            effect_sizes=[0.5, 1.0],
            n_simulations_per_effect=3,
            seed=42,
        )

        d = result.to_dict()

        assert "effect_sizes" in d
        assert "power_at_threshold" in d


class TestSampleSizeAnalysis:
    """Tests for sample size analysis."""

    def test_basic_analysis(self):
        """Test sample size analysis runs."""
        result = run_sample_size_analysis(
            effect_size=1.0,
            sample_sizes=[30, 50, 100],
            n_simulations=5,
            seed=42,
        )

        assert isinstance(result, SampleSizeAnalysisResult)
        assert len(result.sample_sizes) == 3
        assert len(result.power_by_size) == 3

    def test_power_increases_with_n(self):
        """Test that power increases with sample size."""
        result = run_sample_size_analysis(
            effect_size=1.0,
            sample_sizes=[30, 100, 300],
            n_simulations=10,
            seed=42,
        )

        # Power should generally increase with sample size
        assert result.power_by_size[300] >= result.power_by_size[30]


class TestValidateFramework:
    """Tests for comprehensive framework validation."""

    def test_basic_validation(self):
        """Test framework validation runs."""
        result = validate_framework(
            n_samples=50,
            n_pathways=5,
            n_subtypes=2,
            effect_size=1.5,
            n_runs=3,
            seed=42,
        )

        assert "recovery" in result
        assert "type_i_error" in result
        assert "mean_ari" in result["recovery"]

    def test_high_effect_recovery(self):
        """Test recovery with high effect size."""
        result = validate_framework(
            n_samples=100,
            n_pathways=10,
            n_subtypes=2,
            effect_size=2.0,
            n_runs=5,
            seed=42,
        )

        # Should achieve reasonable recovery with high effect
        assert result["recovery"]["mean_ari"] > 0.3

    def test_config_in_result(self):
        """Test that config is included in result."""
        result = validate_framework(
            n_samples=50,
            n_subtypes=3,
            effect_size=1.0,
            n_runs=2,
            seed=42,
        )

        assert "config" in result
        assert result["config"]["n_samples"] == 50
        assert result["config"]["n_subtypes"] == 3
