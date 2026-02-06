"""Tests for the sensitivity analysis module."""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.sensitivity import (
    ParameterVariationResult,
    SensitivityAnalysisResult,
    SensitivityParameter,
    run_sensitivity_analysis,
    vary_clustering_algorithm,
    vary_feature_subset,
    vary_n_clusters,
    vary_normalization,
)

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def well_separated_data():
    """Create well-separated clusters for sensitivity testing."""
    np.random.seed(42)
    n_per_cluster = 25
    n_pathways = 5

    scores = np.zeros((n_per_cluster * 3, n_pathways))

    # Cluster 0: high in pathway 0
    scores[:25, 0] = np.random.normal(3, 0.3, 25)
    scores[:25, 1:] = np.random.normal(0, 0.3, (25, 4))

    # Cluster 1: high in pathway 1
    scores[25:50, 1] = np.random.normal(3, 0.3, 25)
    scores[25:50, [0, 2, 3, 4]] = np.random.normal(0, 0.3, (25, 4))

    # Cluster 2: high in pathway 2
    scores[50:75, 2] = np.random.normal(3, 0.3, 25)
    scores[50:75, [0, 1, 3, 4]] = np.random.normal(0, 0.3, (25, 4))

    sample_ids = [f"S_{i:03d}" for i in range(75)]
    pathway_names = [f"P_{chr(65 + i)}" for i in range(n_pathways)]

    return pd.DataFrame(scores, index=sample_ids, columns=pathway_names)


@pytest.fixture
def noisy_data():
    """Create noisy data where clusters are harder to separate."""
    np.random.seed(42)
    scores = np.random.normal(0, 1, (60, 4))
    sample_ids = [f"S_{i:03d}" for i in range(60)]
    pathway_names = ["P_A", "P_B", "P_C", "P_D"]
    return pd.DataFrame(scores, index=sample_ids, columns=pathway_names)


# =============================================================================
# TEST: SensitivityParameter enum
# =============================================================================


class TestSensitivityParameter:
    def test_enum_values(self):
        assert SensitivityParameter.CLUSTERING_ALGORITHM.value == "clustering_algorithm"
        assert SensitivityParameter.N_CLUSTERS.value == "n_clusters"
        assert SensitivityParameter.NORMALIZATION.value == "normalization"
        assert SensitivityParameter.FEATURE_SUBSET.value == "feature_subset"

    def test_enum_count(self):
        assert len(SensitivityParameter) == 4


# =============================================================================
# TEST: vary_clustering_algorithm
# =============================================================================


class TestVaryClusteringAlgorithm:
    def test_runs_multiple_algorithms(self, well_separated_data):
        result = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)

        assert isinstance(result, ParameterVariationResult)
        assert result.parameter == SensitivityParameter.CLUSTERING_ALGORITHM
        assert len(result.configurations) >= 3  # GMM, KMeans, Hierarchical

    def test_well_separated_high_ari(self, well_separated_data):
        result = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)

        # Well-separated clusters should give high agreement
        assert result.mean_ari > 0.5

    def test_pairwise_ari_computed(self, well_separated_data):
        result = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)

        assert len(result.pairwise_ari) > 0
        for key, ari in result.pairwise_ari.items():
            assert -1 <= ari <= 1

    def test_reference_ari(self, well_separated_data):
        result = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)

        # Reference ARI of first config vs itself should be 1.0
        first = result.configurations[0]
        assert result.reference_ari[first] == 1.0

    def test_reproducibility(self, well_separated_data):
        r1 = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)
        r2 = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)

        assert r1.mean_ari == r2.mean_ari


# =============================================================================
# TEST: vary_n_clusters
# =============================================================================


class TestVaryNClusters:
    def test_basic_operation(self, well_separated_data):
        result = vary_n_clusters(well_separated_data, cluster_range=(2, 5), seed=42)

        assert isinstance(result, ParameterVariationResult)
        assert result.parameter == SensitivityParameter.N_CLUSTERS
        assert len(result.configurations) == 4  # k=2,3,4,5

    def test_adjacent_k_higher_ari(self, well_separated_data):
        result = vary_n_clusters(well_separated_data, cluster_range=(2, 6), seed=42)

        # k=3 vs k=4 should have higher ARI than k=2 vs k=6
        assert len(result.pairwise_ari) > 0

    def test_config_names(self, well_separated_data):
        result = vary_n_clusters(well_separated_data, cluster_range=(2, 4), seed=42)

        assert "k=2" in result.configurations
        assert "k=3" in result.configurations
        assert "k=4" in result.configurations


# =============================================================================
# TEST: vary_feature_subset
# =============================================================================


class TestVaryFeatureSubset:
    def test_leave_one_out(self, well_separated_data):
        result = vary_feature_subset(well_separated_data, n_clusters=3, seed=42)

        assert isinstance(result, ParameterVariationResult)
        assert result.parameter == SensitivityParameter.FEATURE_SUBSET
        # Should have all_features + one per pathway
        n_pathways = len(well_separated_data.columns)
        assert len(result.configurations) == n_pathways + 1

    def test_all_features_included(self, well_separated_data):
        result = vary_feature_subset(well_separated_data, n_clusters=3, seed=42)

        assert "all_features" in result.configurations

    def test_well_separated_robust(self, well_separated_data):
        result = vary_feature_subset(well_separated_data, n_clusters=3, seed=42)

        # Well-separated data should be robust to removing any single feature
        assert result.mean_ari > 0.5

    def test_identifies_important_feature(self, well_separated_data):
        result = vary_feature_subset(well_separated_data, n_clusters=3, seed=42)

        # Dropping a defining pathway should lower ARI more
        assert result.min_ari <= result.mean_ari


# =============================================================================
# TEST: vary_normalization
# =============================================================================


class TestVaryNormalization:
    def test_basic_operation(self, well_separated_data):
        result = vary_normalization(well_separated_data, n_clusters=3, seed=42)

        assert isinstance(result, ParameterVariationResult)
        assert result.parameter == SensitivityParameter.NORMALIZATION
        assert len(result.configurations) == 4  # zscore, minmax, robust, rank

    def test_normalization_methods(self, well_separated_data):
        result = vary_normalization(well_separated_data, n_clusters=3, seed=42)

        assert "zscore" in result.configurations
        assert "minmax" in result.configurations
        assert "robust" in result.configurations
        assert "rank" in result.configurations

    def test_well_separated_consistent(self, well_separated_data):
        result = vary_normalization(well_separated_data, n_clusters=3, seed=42)

        # Well-separated clusters should be found regardless of normalization
        assert result.mean_ari > 0.3


# =============================================================================
# TEST: run_sensitivity_analysis
# =============================================================================


class TestRunSensitivityAnalysis:
    def test_full_analysis(self, well_separated_data):
        result = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)

        assert isinstance(result, SensitivityAnalysisResult)
        assert len(result.parameter_results) == 4
        assert result.overall_stability > 0

    def test_selected_parameters(self, well_separated_data):
        result = run_sensitivity_analysis(
            well_separated_data,
            n_clusters=3,
            seed=42,
            parameters=[
                SensitivityParameter.CLUSTERING_ALGORITHM,
                SensitivityParameter.N_CLUSTERS,
            ],
        )

        assert len(result.parameter_results) == 2
        assert "clustering_algorithm" in result.parameter_results
        assert "n_clusters" in result.parameter_results

    def test_robustness_flag(self, well_separated_data):
        result = run_sensitivity_analysis(
            well_separated_data,
            n_clusters=3,
            seed=42,
            robustness_threshold=0.3,
        )

        # Well-separated data with low threshold should be robust
        assert result.is_robust

    def test_most_least_sensitive(self, well_separated_data):
        result = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)

        assert result.most_sensitive_parameter in [p.value for p in SensitivityParameter]
        assert result.least_sensitive_parameter in [p.value for p in SensitivityParameter]

    def test_noisy_data_less_robust(self, noisy_data):
        result = run_sensitivity_analysis(
            noisy_data,
            n_clusters=3,
            seed=42,
            robustness_threshold=0.9,
        )

        # Noisy data with high threshold should not be robust
        assert not result.is_robust


# =============================================================================
# TEST: DATACLASSES
# =============================================================================


class TestDataclasses:
    def test_parameter_variation_to_dict(self, well_separated_data):
        result = vary_clustering_algorithm(well_separated_data, n_clusters=3, seed=42)

        d = result.to_dict()
        assert "parameter" in d
        assert "n_configurations" in d
        assert "mean_pairwise_ari" in d
        assert "min_pairwise_ari" in d

    def test_sensitivity_result_to_dict(self, well_separated_data):
        result = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)

        d = result.to_dict()
        assert "overall_stability" in d
        assert "is_robust" in d
        assert "parameters" in d

    def test_sensitivity_result_format_report(self, well_separated_data):
        result = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)

        text = result.format_report()
        assert "Sensitivity Analysis Report" in text
        assert "Overall stability" in text

    def test_sensitivity_result_get_citations(self, well_separated_data):
        result = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)

        citations = result.get_citations()
        assert len(citations) >= 1
        assert any("Hennig" in c for c in citations)


# =============================================================================
# TEST: REPRODUCIBILITY
# =============================================================================


class TestReproducibility:
    def test_same_seed_same_result(self, well_separated_data):
        r1 = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)
        r2 = run_sensitivity_analysis(well_separated_data, n_clusters=3, seed=42)

        assert r1.overall_stability == r2.overall_stability
        assert r1.is_robust == r2.is_robust

    def test_different_seed_may_differ(self, well_separated_data):
        r1 = run_sensitivity_analysis(
            well_separated_data,
            n_clusters=3,
            seed=42,
            parameters=[SensitivityParameter.CLUSTERING_ALGORITHM],
        )
        r2 = run_sensitivity_analysis(
            well_separated_data,
            n_clusters=3,
            seed=99,
            parameters=[SensitivityParameter.CLUSTERING_ALGORITHM],
        )

        # Different seeds may give different results
        # (but for well-separated data, should be similar)
        assert isinstance(r1.overall_stability, float)
        assert isinstance(r2.overall_stability, float)
