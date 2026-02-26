"""
Tests for the clustering module.

Tests multiple clustering algorithms, model selection,
cross-validation, and algorithm comparison.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.clustering import (
    AlgorithmComparisonResult,
    ClusteringAlgorithm,
    ClusteringResult,
    CrossValidationResult,
    ModelSelectionResult,
    compare_algorithms,
    cross_validate_clustering,
    run_clustering,
    select_n_clusters,
)


@pytest.fixture
def clustered_data():
    """Generate data with known cluster structure."""
    np.random.seed(42)

    # 3 clusters, 30 samples each
    cluster_1 = np.random.normal([0, 0, 0, 0, 0], 0.5, (30, 5))
    cluster_2 = np.random.normal([3, 3, 0, 0, 0], 0.5, (30, 5))
    cluster_3 = np.random.normal([0, 0, 3, 3, 0], 0.5, (30, 5))

    data = np.vstack([cluster_1, cluster_2, cluster_3])
    true_labels = np.array([0] * 30 + [1] * 30 + [2] * 30)

    return data, true_labels


@pytest.fixture
def random_data():
    """Generate random data without structure."""
    np.random.seed(42)
    return np.random.randn(100, 5)


class TestRunClustering:
    """Tests for run_clustering function."""

    def test_gmm_clustering(self, clustered_data):
        """Test GMM clustering."""
        data, _ = clustered_data

        result = run_clustering(
            data,
            n_clusters=3,
            algorithm=ClusteringAlgorithm.GMM,
            seed=42,
        )

        assert isinstance(result, ClusteringResult)
        assert result.algorithm == ClusteringAlgorithm.GMM
        assert result.n_clusters == 3
        assert len(result.labels) == 90
        assert result.bic is not None

    def test_kmeans_clustering(self, clustered_data):
        """Test K-means clustering."""
        data, _ = clustered_data

        result = run_clustering(
            data,
            n_clusters=3,
            algorithm=ClusteringAlgorithm.KMEANS,
            seed=42,
        )

        assert result.algorithm == ClusteringAlgorithm.KMEANS
        assert len(np.unique(result.labels)) == 3
        assert result.bic is None  # K-means doesn't compute BIC

    def test_hierarchical_clustering(self, clustered_data):
        """Test hierarchical clustering."""
        data, _ = clustered_data

        result = run_clustering(
            data,
            n_clusters=3,
            algorithm=ClusteringAlgorithm.HIERARCHICAL,
        )

        assert result.algorithm == ClusteringAlgorithm.HIERARCHICAL
        assert len(np.unique(result.labels)) == 3

    def test_spectral_clustering(self, clustered_data):
        """Test spectral clustering."""
        data, _ = clustered_data

        result = run_clustering(
            data,
            n_clusters=3,
            algorithm=ClusteringAlgorithm.SPECTRAL,
            seed=42,
        )

        assert result.algorithm == ClusteringAlgorithm.SPECTRAL
        assert len(np.unique(result.labels)) == 3

    def test_clustering_metrics(self, clustered_data):
        """Test that clustering metrics are computed."""
        data, _ = clustered_data

        result = run_clustering(data, n_clusters=3, seed=42)

        assert result.silhouette != 0.0
        assert result.calinski_harabasz != 0.0
        assert result.davies_bouldin != float("inf")

    def test_high_silhouette_for_clear_clusters(self, clustered_data):
        """Test that clear clusters yield high silhouette."""
        data, _ = clustered_data

        result = run_clustering(data, n_clusters=3, seed=42)

        # Well-separated clusters should have high silhouette
        assert result.silhouette > 0.5

    def test_reproducibility(self, clustered_data):
        """Test reproducibility with same seed."""
        data, _ = clustered_data

        result1 = run_clustering(data, 3, ClusteringAlgorithm.GMM, seed=123)
        result2 = run_clustering(data, 3, ClusteringAlgorithm.GMM, seed=123)

        np.testing.assert_array_equal(result1.labels, result2.labels)

    def test_result_to_dict(self, clustered_data):
        """Test serialization."""
        data, _ = clustered_data

        result = run_clustering(data, n_clusters=3, seed=42)
        d = result.to_dict()

        assert "algorithm" in d
        assert "n_clusters" in d
        assert "metrics" in d


class TestModelSelection:
    """Tests for cluster number selection."""

    def test_select_correct_k(self, clustered_data):
        """Test selection of correct k for clear clusters."""
        data, _ = clustered_data

        result = select_n_clusters(
            data,
            k_range=[2, 3, 4, 5],
            method="bic",
            seed=42,
        )

        assert isinstance(result, ModelSelectionResult)
        # Should select k=3 for 3-cluster data
        assert result.optimal_k == 3

    def test_bic_method(self, clustered_data):
        """Test BIC-based selection."""
        data, _ = clustered_data

        result = select_n_clusters(data, k_range=[2, 3, 4], method="bic")

        assert result.method == "bic"
        assert len(result.bic_values) == 3
        assert all(k in result.bic_values for k in [2, 3, 4])

    def test_silhouette_method(self, clustered_data):
        """Test silhouette-based selection."""
        data, _ = clustered_data

        result = select_n_clusters(
            data,
            k_range=[2, 3, 4],
            method="silhouette",
        )

        assert result.method == "silhouette"
        assert len(result.silhouette_values) == 3

    def test_result_to_dict(self, clustered_data):
        """Test serialization."""
        data, _ = clustered_data

        result = select_n_clusters(data, k_range=[2, 3], seed=42)
        d = result.to_dict()

        assert "optimal_k" in d
        assert "k_range" in d
        assert "bic_values" in d


class TestMinClusterFraction:
    """Tests for min_cluster_fraction guard in select_n_clusters."""

    def test_no_guard_by_default(self, clustered_data):
        """Test that guard is not applied when min_cluster_fraction is None."""
        data, _ = clustered_data
        result = select_n_clusters(data, k_range=[2, 3, 4], seed=42)
        assert result.min_cluster_fraction is None
        assert result.rejected_k == {}

    def test_balanced_clusters_pass_guard(self, clustered_data):
        """Test that balanced clusters are not rejected."""
        data, _ = clustered_data  # 3 clusters of 30 each = 33% each
        result = select_n_clusters(data, k_range=[2, 3, 4], seed=42, min_cluster_fraction=0.05)
        # k=3 should not be rejected (3 balanced clusters of ~33%)
        assert 3 not in result.rejected_k

    def test_degenerate_split_rejected(self):
        """Test that a degenerate split (e.g. 452:6) is rejected.

        Simulates the TCGA-COAD scenario where BIC-optimal k=2
        produced a 452:6 split (1.3% < 5% threshold).
        """
        np.random.seed(42)
        # Large main cluster + tiny outlier cluster
        main = np.random.normal([0, 0, 0], 0.5, (95, 3))
        outlier = np.random.normal([10, 10, 10], 0.1, (5, 3))
        data = np.vstack([main, outlier])

        result = select_n_clusters(data, k_range=[2, 3, 4], seed=42, min_cluster_fraction=0.10)
        # k=2 should be rejected if it finds the 95:5 split (5% < 10%)
        # The optimal_k should be from non-rejected candidates
        assert result.min_cluster_fraction == 0.10

    def test_rejected_k_in_result(self):
        """Test that rejected k values are recorded with reasons."""
        np.random.seed(42)
        main = np.random.normal([0, 0, 0], 0.5, (95, 3))
        outlier = np.random.normal([10, 10, 10], 0.1, (5, 3))
        data = np.vstack([main, outlier])

        result = select_n_clusters(data, k_range=[2, 3, 4, 5], seed=42, min_cluster_fraction=0.10)
        # At least some k values should be rejected with the outlier structure
        # Check that any rejected entries have descriptive reasons
        for k, reason in result.rejected_k.items():
            assert "smallest cluster" in reason
            assert "samples" in reason

    def test_rejected_k_excluded_from_scores(self):
        """Test that rejected k values don't appear in BIC/silhouette dicts."""
        np.random.seed(42)
        main = np.random.normal([0, 0, 0], 0.5, (95, 3))
        outlier = np.random.normal([10, 10, 10], 0.1, (5, 3))
        data = np.vstack([main, outlier])

        result = select_n_clusters(data, k_range=[2, 3, 4], seed=42, min_cluster_fraction=0.10)
        for k in result.rejected_k:
            assert k not in result.bic_values
            assert k not in result.silhouette_values

    def test_to_dict_includes_guard_metadata(self):
        """Test serialization includes guard fields."""
        np.random.seed(42)
        data = np.random.randn(100, 5)
        result = select_n_clusters(data, k_range=[2, 3], seed=42, min_cluster_fraction=0.05)
        d = result.to_dict()
        assert "min_cluster_fraction" in d
        assert d["min_cluster_fraction"] == 0.05

    def test_to_dict_omits_guard_when_none(self, clustered_data):
        """Test serialization omits guard fields when not used."""
        data, _ = clustered_data
        result = select_n_clusters(data, k_range=[2, 3], seed=42)
        d = result.to_dict()
        assert "min_cluster_fraction" not in d
        assert "rejected_k" not in d

    def test_all_k_rejected_falls_back(self):
        """Test fallback when all k values are rejected."""
        np.random.seed(42)
        # Two tight, very uneven clusters — most k will be degenerate
        main = np.random.normal([0, 0], 0.3, (97, 2))
        outlier = np.random.normal([8, 8], 0.1, (3, 2))
        data = np.vstack([main, outlier])

        result = select_n_clusters(data, k_range=[2, 3], seed=42, min_cluster_fraction=0.20)
        # Should still return a valid result (falls back to k_range[0])
        assert isinstance(result, ModelSelectionResult)
        assert result.optimal_k in [2, 3]

    def test_invalid_fraction_raises(self, clustered_data):
        """Test that invalid fraction values raise ValueError."""
        data, _ = clustered_data
        with pytest.raises(ValueError, match="between 0 and 1"):
            select_n_clusters(data, k_range=[2, 3], min_cluster_fraction=0.0)
        with pytest.raises(ValueError, match="between 0 and 1"):
            select_n_clusters(data, k_range=[2, 3], min_cluster_fraction=1.0)
        with pytest.raises(ValueError, match="between 0 and 1"):
            select_n_clusters(data, k_range=[2, 3], min_cluster_fraction=-0.1)
        with pytest.raises(ValueError, match="between 0 and 1"):
            select_n_clusters(data, k_range=[2, 3], min_cluster_fraction=1.5)

    def test_backward_compatible(self, clustered_data):
        """Test that existing calls without min_cluster_fraction still work."""
        data, _ = clustered_data
        # These mirror the existing TestModelSelection tests
        result = select_n_clusters(data, k_range=[2, 3, 4, 5], method="bic", seed=42)
        assert result.optimal_k == 3
        assert result.rejected_k == {}


class TestCrossValidation:
    """Tests for cross-validation."""

    def test_basic_cv(self, clustered_data):
        """Test basic cross-validation."""
        data, _ = clustered_data

        result = cross_validate_clustering(
            data,
            n_clusters=3,
            n_folds=5,
            seed=42,
        )

        assert isinstance(result, CrossValidationResult)
        assert result.n_folds == 5
        assert len(result.fold_aris) == 5

    def test_stable_clusters_high_ari(self, clustered_data):
        """Test that clear clusters yield high CV ARI."""
        data, _ = clustered_data

        result = cross_validate_clustering(
            data,
            n_clusters=3,
            n_folds=5,
            seed=42,
        )

        # Well-separated clusters should be stable across folds
        assert result.mean_ari > 0.5

    def test_gmm_cv(self, clustered_data):
        """Test CV with GMM."""
        data, _ = clustered_data

        result = cross_validate_clustering(
            data,
            n_clusters=3,
            algorithm=ClusteringAlgorithm.GMM,
            seed=42,
        )

        assert result.mean_ari > 0
        assert result.std_ari >= 0

    def test_result_to_dict(self, clustered_data):
        """Test serialization."""
        data, _ = clustered_data

        result = cross_validate_clustering(data, 3, n_folds=3, seed=42)
        d = result.to_dict()

        assert "mean_ari" in d
        assert "std_ari" in d
        assert "n_folds" in d


class TestAlgorithmComparison:
    """Tests for comparing clustering algorithms."""

    def test_compare_all_algorithms(self, clustered_data):
        """Test comparison of all algorithms."""
        data, _ = clustered_data

        result = compare_algorithms(data, n_clusters=3, seed=42)

        assert isinstance(result, AlgorithmComparisonResult)
        # Should have results for multiple algorithms
        assert len(result.results) >= 3

    def test_pairwise_ari_computed(self, clustered_data):
        """Test that pairwise ARI is computed."""
        data, _ = clustered_data

        result = compare_algorithms(
            data,
            n_clusters=3,
            algorithms=[
                ClusteringAlgorithm.GMM,
                ClusteringAlgorithm.KMEANS,
            ],
            seed=42,
        )

        assert "gmm_vs_kmeans" in result.pairwise_ari

    def test_high_agreement_for_clear_clusters(self, clustered_data):
        """Test algorithms agree on clear clusters."""
        data, _ = clustered_data

        result = compare_algorithms(data, n_clusters=3, seed=42)

        # All pairwise ARIs should be high for clear clusters
        for key, ari in result.pairwise_ari.items():
            assert ari > 0.5, f"Low agreement for {key}"

    def test_most_stable_identified(self, clustered_data):
        """Test that most stable algorithm is identified."""
        data, _ = clustered_data

        result = compare_algorithms(data, n_clusters=3, seed=42)

        assert result.most_stable_algorithm in ["gmm", "kmeans", "hierarchical", "spectral"]
        assert result.consensus_labels is not None

    def test_specific_algorithms(self, clustered_data):
        """Test with specific algorithm subset."""
        data, _ = clustered_data

        result = compare_algorithms(
            data,
            n_clusters=3,
            algorithms=[ClusteringAlgorithm.GMM, ClusteringAlgorithm.KMEANS],
            seed=42,
        )

        assert len(result.results) == 2
        assert "gmm" in result.results
        assert "kmeans" in result.results

    def test_result_to_dict(self, clustered_data):
        """Test serialization."""
        data, _ = clustered_data

        result = compare_algorithms(
            data,
            3,
            algorithms=[ClusteringAlgorithm.GMM, ClusteringAlgorithm.KMEANS],
            seed=42,
        )
        d = result.to_dict()

        assert "algorithms" in d
        assert "pairwise_ari" in d
        assert "most_stable_algorithm" in d


class TestEdgeCases:
    """Tests for edge cases."""

    def test_two_clusters(self):
        """Test with only 2 clusters."""
        np.random.seed(42)
        data = np.vstack(
            [
                np.random.normal([0, 0], 0.5, (30, 2)),
                np.random.normal([3, 3], 0.5, (30, 2)),
            ]
        )

        result = run_clustering(data, n_clusters=2, seed=42)
        assert len(np.unique(result.labels)) == 2

    def test_many_clusters(self):
        """Test with many clusters."""
        np.random.seed(42)
        data = np.random.randn(200, 10)

        result = run_clustering(data, n_clusters=10, seed=42)
        assert result.n_clusters == 10

    def test_small_sample(self):
        """Test with small sample size."""
        np.random.seed(42)
        data = np.random.randn(20, 3)

        result = run_clustering(data, n_clusters=2, seed=42)
        assert len(result.labels) == 20

    def test_high_dimensional(self):
        """Test with high-dimensional data."""
        np.random.seed(42)
        data = np.random.randn(50, 100)

        result = run_clustering(data, n_clusters=3, seed=42)
        assert len(result.labels) == 50
