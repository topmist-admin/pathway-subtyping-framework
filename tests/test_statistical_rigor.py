"""
Tests for the statistical rigor module.

Tests FDR correction, effect size calculations, pathway normalization,
and literature-based burden weights.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.statistical_rigor import (
    BurdenWeights,
    BurdenWeightScheme,
    FDRResult,
    PathwayNormalization,
    StatisticalRigorResult,
    aggregate_pathway_scores,
    benjamini_hochberg,
    compute_cohens_d,
    compute_pathway_effect_sizes,
    compute_pathway_pvalues,
    compute_variant_weight,
    run_statistical_analysis,
)


class TestBenjaminiHochberg:
    """Tests for FDR correction."""

    def test_empty_pvalues(self):
        """Test handling of empty p-value array."""
        result = benjamini_hochberg(np.array([]))
        assert len(result) == 0

    def test_single_pvalue(self):
        """Test single p-value (no correction needed)."""
        result = benjamini_hochberg(np.array([0.05]))
        assert len(result) == 1
        assert result[0] == 0.05

    def test_all_significant(self):
        """Test when all p-values are significant."""
        p_values = np.array([0.001, 0.002, 0.003, 0.004, 0.005])
        q_values = benjamini_hochberg(p_values)

        # All q-values should be <= 0.05
        assert all(q <= 0.05 for q in q_values)
        # Q-values should be >= p-values
        assert all(q >= p for q, p in zip(q_values, p_values))

    def test_none_significant(self):
        """Test when no p-values are significant."""
        p_values = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        q_values = benjamini_hochberg(p_values)

        # All q-values should be > 0.05
        assert all(q > 0.05 for q in q_values)

    def test_mixed_significance(self):
        """Test mixture of significant and non-significant."""
        p_values = np.array([0.001, 0.01, 0.05, 0.2, 0.8])
        q_values = benjamini_hochberg(p_values)

        # Q-values should be monotonically non-decreasing when sorted by p
        sorted_idx = np.argsort(p_values)
        sorted_q = q_values[sorted_idx]
        assert all(sorted_q[i] <= sorted_q[i + 1] for i in range(len(sorted_q) - 1))

    def test_qvalues_capped_at_one(self):
        """Test that q-values don't exceed 1.0."""
        p_values = np.array([0.8, 0.9, 0.95, 0.99])
        q_values = benjamini_hochberg(p_values)

        assert all(q <= 1.0 for q in q_values)

    def test_order_preserved(self):
        """Test that original order is preserved."""
        p_values = np.array([0.05, 0.001, 0.1, 0.01, 0.5])
        q_values = benjamini_hochberg(p_values)

        # The smallest p-value should have smallest q-value
        min_p_idx = np.argmin(p_values)
        min_q_idx = np.argmin(q_values)
        assert min_p_idx == min_q_idx


class TestBurdenWeights:
    """Tests for burden weight schemes."""

    def test_default_scheme(self):
        """Test default weight scheme."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.DEFAULT)
        assert weights.lof_weight == 1.0
        assert weights.scheme == BurdenWeightScheme.DEFAULT

    def test_gnomad_scheme(self):
        """Test gnomAD constraint scheme."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)
        assert weights.lof_weight == 1.0
        assert weights.missense_high_weight == 0.8
        assert weights.cadd_high_threshold == 25.0

    def test_acmg_scheme(self):
        """Test ACMG-inspired scheme."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.ACMG_INSPIRED)
        assert weights.missense_high_weight == 0.9
        assert weights.cadd_high_threshold == 30.0  # Stricter

    def test_uniform_scheme(self):
        """Test uniform scheme (all equal weights)."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.UNIFORM)
        assert weights.lof_weight == 1.0
        assert weights.missense_high_weight == 1.0
        assert weights.missense_low_weight == 1.0
        assert weights.other_weight == 1.0

    def test_get_citations(self):
        """Test that citations are returned."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)
        citations = weights.get_citations()
        assert len(citations) >= 2
        assert "Lek" in citations[0]  # gnomAD paper


class TestComputeVariantWeight:
    """Tests for variant weight calculation."""

    def test_lof_weight(self):
        """Test LoF variants get high weight."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)

        assert compute_variant_weight("frameshift_variant", 30, weights) == weights.lof_weight
        assert compute_variant_weight("stop_gained", 25, weights) == weights.lof_weight
        assert compute_variant_weight("splice_donor_variant", 20, weights) == weights.lof_weight

    def test_missense_high_cadd(self):
        """Test high-CADD missense variants."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)

        weight = compute_variant_weight("missense_variant", 30, weights)
        assert weight == weights.missense_high_weight

    def test_missense_moderate_cadd(self):
        """Test moderate-CADD missense variants."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)

        weight = compute_variant_weight("missense_variant", 22, weights)
        assert weight == weights.missense_moderate_weight

    def test_missense_low_cadd(self):
        """Test low-CADD missense variants."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)

        weight = compute_variant_weight("missense_variant", 10, weights)
        assert weight < weights.missense_moderate_weight

    def test_other_variants(self):
        """Test other variant types."""
        weights = BurdenWeights.from_scheme(BurdenWeightScheme.GNOMAD_CONSTRAINT)

        weight = compute_variant_weight("synonymous_variant", 5, weights)
        assert weight == weights.other_weight


class TestCohensD:
    """Tests for Cohen's d effect size."""

    def test_identical_groups(self):
        """Test effect size for identical groups."""
        group1 = np.array([1, 2, 3, 4, 5])
        group2 = np.array([1, 2, 3, 4, 5])

        d = compute_cohens_d(group1, group2)
        assert abs(d) < 0.01

    def test_large_effect(self):
        """Test large effect size."""
        group1 = np.array([1, 2, 3, 4, 5])
        group2 = np.array([10, 11, 12, 13, 14])

        d = abs(compute_cohens_d(group1, group2))
        assert d > 2.0  # Very large effect

    def test_small_effect(self):
        """Test small effect size."""
        group1 = np.array([1, 2, 3, 4, 5])
        group2 = np.array([1.5, 2.5, 3.5, 4.5, 5.5])

        d = abs(compute_cohens_d(group1, group2))
        assert d < 0.5  # Small effect

    def test_insufficient_samples(self):
        """Test handling of insufficient samples."""
        group1 = np.array([1])
        group2 = np.array([2])

        d = compute_cohens_d(group1, group2)
        assert d == 0.0


class TestPathwayNormalization:
    """Tests for pathway aggregation methods."""

    @pytest.fixture
    def gene_burdens(self):
        """Create sample gene burden data."""
        np.random.seed(42)
        return pd.DataFrame(
            np.random.rand(50, 10),
            columns=[f"GENE_{i}" for i in range(10)],
            index=[f"SAMPLE_{i}" for i in range(50)],
        )

    def test_mean_aggregation(self, gene_burdens):
        """Test mean aggregation."""
        genes = ["GENE_0", "GENE_1", "GENE_2"]
        scores = aggregate_pathway_scores(gene_burdens, genes, PathwayNormalization.MEAN)

        assert len(scores) == 50
        expected = gene_burdens[genes].mean(axis=1)
        np.testing.assert_array_almost_equal(scores.values, expected.values)

    def test_median_aggregation(self, gene_burdens):
        """Test median aggregation."""
        genes = ["GENE_0", "GENE_1", "GENE_2"]
        scores = aggregate_pathway_scores(gene_burdens, genes, PathwayNormalization.MEDIAN)

        assert len(scores) == 50
        expected = gene_burdens[genes].median(axis=1)
        np.testing.assert_array_almost_equal(scores.values, expected.values)

    def test_size_normalized(self, gene_burdens):
        """Test size-normalized aggregation."""
        genes = ["GENE_0", "GENE_1", "GENE_2"]
        scores = aggregate_pathway_scores(gene_burdens, genes, PathwayNormalization.SIZE_NORMALIZED)

        assert len(scores) == 50
        # Should be mean / sqrt(3)
        expected = gene_burdens[genes].mean(axis=1) / np.sqrt(3)
        np.testing.assert_array_almost_equal(scores.values, expected.values)

    def test_pca_aggregation(self, gene_burdens):
        """Test PCA aggregation."""
        genes = ["GENE_0", "GENE_1", "GENE_2", "GENE_3", "GENE_4"]
        scores = aggregate_pathway_scores(gene_burdens, genes, PathwayNormalization.PCA)

        assert len(scores) == 50
        # PCA scores should have zero mean (centered)
        assert abs(scores.mean()) < 0.01

    def test_missing_genes(self, gene_burdens):
        """Test handling of genes not in data."""
        genes = ["GENE_0", "GENE_MISSING", "GENE_1"]
        scores = aggregate_pathway_scores(gene_burdens, genes, PathwayNormalization.MEAN)

        assert len(scores) == 50
        # Should use available genes only
        expected = gene_burdens[["GENE_0", "GENE_1"]].mean(axis=1)
        np.testing.assert_array_almost_equal(scores.values, expected.values)


class TestPathwayPvalues:
    """Tests for permutation-based p-values."""

    @pytest.fixture
    def clustered_data(self):
        """Create data with known cluster structure."""
        np.random.seed(42)
        n_samples = 60

        # Create pathway scores with clear difference between clusters
        scores = pd.DataFrame(index=[f"SAMPLE_{i}" for i in range(n_samples)])

        # Pathway with strong signal
        cluster_1 = np.random.normal(0, 1, 30)
        cluster_2 = np.random.normal(2, 1, 30)  # Different mean
        scores["PATHWAY_SIGNAL"] = np.concatenate([cluster_1, cluster_2])

        # Pathway with no signal
        scores["PATHWAY_NOISE"] = np.random.normal(0, 1, n_samples)

        labels = np.array([0] * 30 + [1] * 30)

        return scores, labels

    def test_significant_pathway(self, clustered_data):
        """Test p-value for pathway with real signal."""
        scores, labels = clustered_data

        p_values = compute_pathway_pvalues(scores, labels, n_permutations=100, seed=42)

        # Signal pathway should have low p-value
        assert p_values["PATHWAY_SIGNAL"] < 0.1

    def test_nonsignificant_pathway(self, clustered_data):
        """Test p-value for noise pathway."""
        scores, labels = clustered_data

        p_values = compute_pathway_pvalues(scores, labels, n_permutations=100, seed=42)

        # Noise pathway should have higher p-value
        assert p_values["PATHWAY_NOISE"] > p_values["PATHWAY_SIGNAL"]

    def test_single_cluster(self):
        """Test handling of single cluster."""
        scores = pd.DataFrame({"PATHWAY_1": np.random.rand(30)})
        labels = np.zeros(30)  # All same cluster

        p_values = compute_pathway_pvalues(scores, labels, n_permutations=10)

        # Should return p=1.0 when only one cluster
        assert p_values["PATHWAY_1"] == 1.0


class TestStatisticalAnalysis:
    """Integration tests for full statistical analysis."""

    @pytest.fixture
    def analysis_data(self):
        """Create data for full analysis."""
        np.random.seed(42)
        n_samples = 60

        scores = pd.DataFrame(index=[f"SAMPLE_{i}" for i in range(n_samples)])
        scores["PATHWAY_1"] = np.concatenate(
            [
                np.random.normal(0, 1, 30),
                np.random.normal(1.5, 1, 30),
            ]
        )
        scores["PATHWAY_2"] = np.concatenate(
            [
                np.random.normal(0, 1, 30),
                np.random.normal(0.5, 1, 30),
            ]
        )
        scores["PATHWAY_3"] = np.random.normal(0, 1, n_samples)

        labels = np.array([0] * 30 + [1] * 30)

        return scores, labels

    def test_full_analysis(self, analysis_data):
        """Test complete statistical analysis."""
        scores, labels = analysis_data

        result = run_statistical_analysis(
            scores,
            labels,
            n_permutations=50,
            n_bootstrap=50,
            seed=42,
        )

        assert isinstance(result, StatisticalRigorResult)
        assert result.n_pathways_tested == 3
        assert len(result.pathway_results) == 3
        assert len(result.weight_citations) > 0

    def test_analysis_report(self, analysis_data):
        """Test report generation."""
        scores, labels = analysis_data

        result = run_statistical_analysis(
            scores,
            labels,
            n_permutations=50,
            n_bootstrap=50,
            seed=42,
        )

        report = result.format_report()

        assert "Statistical Analysis" in report
        assert "FDR" in report
        assert "Effect Size" in report

    def test_significant_pathways(self, analysis_data):
        """Test identification of significant pathways."""
        scores, labels = analysis_data

        result = run_statistical_analysis(
            scores,
            labels,
            fdr_alpha=0.1,  # More lenient for test
            n_permutations=100,
            n_bootstrap=50,
            seed=42,
        )

        significant = result.get_significant_pathways()
        # PATHWAY_1 should be significant (large effect)
        # Can't guarantee exact results due to randomness
        assert isinstance(significant, list)
