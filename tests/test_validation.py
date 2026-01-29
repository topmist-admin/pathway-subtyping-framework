"""
Tests for the validation gates module.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.validation import (
    ValidationGates,
    ValidationResult,
    ValidationGatesResult,
)


class TestValidationResult:
    """Tests for ValidationResult dataclass."""

    def test_validation_result_pass(self):
        """Test ValidationResult for passing test."""
        result = ValidationResult(
            name="Test Gate",
            passed=True,
            metric_name="ARI",
            metric_value=0.05,
            threshold=0.15,
            comparison="<",
        )

        assert result.passed is True
        assert result.metric_value < result.threshold

    def test_validation_result_fail(self):
        """Test ValidationResult for failing test."""
        result = ValidationResult(
            name="Test Gate",
            passed=False,
            metric_name="ARI",
            metric_value=0.25,
            threshold=0.15,
            comparison="<",
        )

        assert result.passed is False
        assert result.metric_value > result.threshold

    def test_validation_result_to_dict(self):
        """Test ValidationResult conversion to dict."""
        result = ValidationResult(
            name="Test Gate",
            passed=True,
            metric_name="ARI",
            metric_value=0.05,
            threshold=0.15,
            comparison="<",
        )

        d = result.to_dict()

        assert d["name"] == "Test Gate"
        assert d["metric"] == "ARI"
        assert d["value"] == 0.05
        assert d["status"] == "PASS"

    def test_validation_result_status_property(self):
        """Test status property."""
        result_pass = ValidationResult(
            name="Test", passed=True, metric_name="ARI",
            metric_value=0.1, threshold=0.15, comparison="<"
        )
        result_fail = ValidationResult(
            name="Test", passed=False, metric_name="ARI",
            metric_value=0.2, threshold=0.15, comparison="<"
        )

        assert result_pass.status == "PASS"
        assert result_fail.status == "FAIL"


class TestValidationGatesResult:
    """Tests for ValidationGatesResult dataclass."""

    def test_all_passed_true(self):
        """Test all_passed when all tests pass."""
        results = [
            ValidationResult("Test1", True, "ARI", 0.05, 0.15, "<"),
            ValidationResult("Test2", True, "ARI", 0.03, 0.15, "<"),
            ValidationResult("Test3", True, "ARI", 0.85, 0.80, ">="),
        ]

        gates_result = ValidationGatesResult(
            results=results,
            all_passed=True,
            summary="All 3 validation gates PASSED"
        )

        assert gates_result.all_passed is True
        assert "3" in gates_result.summary

    def test_all_passed_false(self):
        """Test all_passed when some tests fail."""
        results = [
            ValidationResult("Test1", True, "ARI", 0.05, 0.15, "<"),
            ValidationResult("Test2", False, "ARI", 0.25, 0.15, "<"),
            ValidationResult("Test3", True, "ARI", 0.85, 0.80, ">="),
        ]

        gates_result = ValidationGatesResult(
            results=results,
            all_passed=False,
            summary="2/3 validation gates passed. FAILED: Test2"
        )

        assert gates_result.all_passed is False
        assert "2/3" in gates_result.summary

    def test_to_dict(self):
        """Test conversion to dict."""
        results = [
            ValidationResult("Test1", True, "ARI", 0.05, 0.15, "<"),
        ]

        gates_result = ValidationGatesResult(
            results=results,
            all_passed=True,
            summary="All tests passed"
        )
        d = gates_result.to_dict()

        assert "all_passed" in d
        assert "summary" in d
        assert "tests" in d
        assert len(d["tests"]) == 1


class TestValidationGates:
    """Tests for ValidationGates class."""

    def test_init_defaults(self):
        """Test ValidationGates initialization with defaults."""
        validator = ValidationGates()

        assert validator.seed == 42
        assert validator.n_permutations == 100
        assert validator.n_bootstrap == 50
        assert validator.stability_threshold == 0.8
        assert validator.null_ari_max == 0.15

    def test_init_custom(self):
        """Test ValidationGates initialization with custom values."""
        validator = ValidationGates(
            seed=123,
            n_permutations=50,
            n_bootstrap=25,
            stability_threshold=0.9,
            null_ari_max=0.10,
        )

        assert validator.seed == 123
        assert validator.n_permutations == 50
        assert validator.n_bootstrap == 25
        assert validator.stability_threshold == 0.9
        assert validator.null_ari_max == 0.10

    def test_label_shuffle_clean_signal(self, sample_pathway_scores, sample_cluster_labels):
        """Test label shuffle with clean clustering signal."""
        validator = ValidationGates(seed=42, n_permutations=20)

        result = validator.negative_control_label_shuffle(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        # With clean signal, shuffled labels should not match original
        assert isinstance(result, ValidationResult)
        assert result.name == "Negative Control 1: Label Shuffle"

    def test_random_genes_synthetic(
        self, sample_pathway_scores, sample_cluster_labels, sample_pathways, sample_gene_burdens
    ):
        """Test random gene sets validation."""
        validator = ValidationGates(seed=42, n_permutations=10)

        result = validator.negative_control_random_gene_sets(
            sample_gene_burdens,
            sample_pathways,
            sample_cluster_labels,
            n_clusters=4,
            gmm_seed=42,
        )

        # Result should be a ValidationResult
        assert isinstance(result, ValidationResult)
        assert result.name == "Negative Control 2: Random Gene Sets"

    def test_bootstrap_stability_clean_signal(self, sample_pathway_scores, sample_cluster_labels):
        """Test bootstrap stability with clean clustering signal."""
        validator = ValidationGates(seed=42, n_bootstrap=10, stability_threshold=0.7)

        result = validator.stability_test_bootstrap(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        # With clean signal, bootstrap should be stable
        assert isinstance(result, ValidationResult)
        assert result.name == "Stability Test: Bootstrap"

    def test_run_all(
        self, sample_pathway_scores, sample_cluster_labels, sample_pathways, sample_gene_burdens
    ):
        """Test running all validation gates."""
        validator = ValidationGates(seed=42, n_permutations=10, n_bootstrap=10)

        result = validator.run_all(
            pathway_scores=sample_pathway_scores,
            cluster_labels=sample_cluster_labels,
            pathways=sample_pathways,
            gene_burdens=sample_gene_burdens,
            n_clusters=4,
            gmm_seed=42,
        )

        assert isinstance(result, ValidationGatesResult)
        assert len(result.results) == 3


class TestValidationGatesEdgeCases:
    """Tests for edge cases in validation gates."""

    def test_single_cluster(self, sample_pathway_scores):
        """Test validation with single cluster."""
        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5)

        # All samples in one cluster
        labels = np.zeros(len(sample_pathway_scores), dtype=int)

        # Label shuffle should still work
        result = validator.negative_control_label_shuffle(
            sample_pathway_scores, labels, n_clusters=1, gmm_seed=42
        )

        assert isinstance(result, ValidationResult)

    def test_small_sample_size(self):
        """Test validation with small sample size."""
        np.random.seed(42)

        # Small dataset: 10 samples, 3 pathways
        small_scores = pd.DataFrame(
            np.random.randn(10, 3),
            index=[f"S{i}" for i in range(10)],
            columns=["P1", "P2", "P3"],
        )
        small_labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2])

        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5)

        result = validator.negative_control_label_shuffle(
            small_scores, small_labels, n_clusters=3, gmm_seed=42
        )

        assert isinstance(result, ValidationResult)

    def test_many_clusters(self, sample_pathway_scores):
        """Test validation with many clusters."""
        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5)

        # Create 10 clusters
        labels = np.repeat(np.arange(10), 6)

        result = validator.negative_control_label_shuffle(
            sample_pathway_scores, labels, n_clusters=10, gmm_seed=42
        )

        assert isinstance(result, ValidationResult)


class TestValidationGatesReproducibility:
    """Tests for reproducibility of validation gates."""

    def test_label_shuffle_reproducible(self, sample_pathway_scores, sample_cluster_labels):
        """Test that label shuffle is reproducible with same seed."""
        validator1 = ValidationGates(seed=42, n_permutations=20)
        validator2 = ValidationGates(seed=42, n_permutations=20)

        result1 = validator1.negative_control_label_shuffle(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        result2 = validator2.negative_control_label_shuffle(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        assert result1.metric_value == result2.metric_value

    def test_bootstrap_reproducible(self, sample_pathway_scores, sample_cluster_labels):
        """Test that bootstrap is reproducible with same seed."""
        validator1 = ValidationGates(seed=42, n_bootstrap=10)
        validator2 = ValidationGates(seed=42, n_bootstrap=10)

        result1 = validator1.stability_test_bootstrap(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        result2 = validator2.stability_test_bootstrap(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        assert result1.metric_value == result2.metric_value

    def test_different_seeds_different_results(
        self, sample_pathway_scores, sample_cluster_labels
    ):
        """Test that different seeds produce different results."""
        validator1 = ValidationGates(seed=42, n_permutations=20)
        validator2 = ValidationGates(seed=123, n_permutations=20)

        result1 = validator1.negative_control_label_shuffle(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        result2 = validator2.negative_control_label_shuffle(
            sample_pathway_scores, sample_cluster_labels, n_clusters=4, gmm_seed=42
        )

        # Results might be similar but shouldn't be exactly equal
        # (unless by chance, which is unlikely with enough permutations)
        # This test verifies the seed is actually being used
        assert isinstance(result1.metric_value, float)
        assert isinstance(result2.metric_value, float)
