"""
Tests for the ancestry/population stratification correction module.
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.ancestry import (
    AncestryAdjustmentResult,
    AncestryMethod,
    AncestryPCs,
    AncestryStratificationReport,
    adjust_pathway_scores,
    check_ancestry_independence,
    compute_ancestry_correlation,
    compute_ancestry_pcs,
    stratified_analysis,
)

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def genotype_matrix_3pop():
    """Genotype matrix with 3 population clusters."""
    rng = np.random.RandomState(42)
    n_samples = 60
    n_variants = 200

    genotypes = np.zeros((n_samples, n_variants))

    # Population 1 (samples 0-19): higher freq at variants 0-50
    genotypes[:20, :50] = rng.binomial(2, 0.3, (20, 50))
    genotypes[:20, 50:] = rng.binomial(2, 0.05, (20, 150))

    # Population 2 (samples 20-39): higher freq at variants 50-100
    genotypes[20:40, :50] = rng.binomial(2, 0.05, (20, 50))
    genotypes[20:40, 50:100] = rng.binomial(2, 0.3, (20, 50))
    genotypes[20:40, 100:] = rng.binomial(2, 0.05, (20, 100))

    # Population 3 (samples 40-59): higher freq at variants 100-150
    genotypes[40:60, :100] = rng.binomial(2, 0.05, (20, 100))
    genotypes[40:60, 100:150] = rng.binomial(2, 0.3, (20, 50))
    genotypes[40:60, 150:] = rng.binomial(2, 0.05, (20, 50))

    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    variant_ids = [f"VAR_{i}" for i in range(n_variants)]

    return pd.DataFrame(genotypes, index=sample_ids, columns=variant_ids)


@pytest.fixture
def ancestry_pcs_3pop(genotype_matrix_3pop):
    """Pre-computed ancestry PCs from 3-population genotype matrix."""
    return compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)


@pytest.fixture
def confounded_pathway_scores():
    """Pathway scores confounded with 3 ancestry groups."""
    rng = np.random.RandomState(42)
    n_samples = 60
    n_pathways = 5

    scores = rng.normal(0, 1, (n_samples, n_pathways))

    # Add ancestry-correlated signal: group 1 high in pathway 0, etc.
    scores[:20, 0] += 2.0
    scores[20:40, 1] += 2.0
    scores[40:60, 2] += 2.0

    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    pathway_names = [f"PATHWAY_{i}" for i in range(n_pathways)]

    return pd.DataFrame(scores, index=sample_ids, columns=pathway_names)


@pytest.fixture
def non_confounded_pathway_scores():
    """Pathway scores NOT confounded with ancestry."""
    rng = np.random.RandomState(99)
    n_samples = 60
    n_pathways = 5

    scores = rng.normal(0, 1, (n_samples, n_pathways))

    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    pathway_names = [f"PATHWAY_{i}" for i in range(n_pathways)]

    return pd.DataFrame(scores, index=sample_ids, columns=pathway_names)


# =============================================================================
# TEST: AncestryMethod Enum
# =============================================================================


class TestAncestryMethod:
    """Tests for AncestryMethod enum."""

    def test_regress_out_value(self):
        assert AncestryMethod.REGRESS_OUT.value == "regress_out"

    def test_covariate_aware_value(self):
        assert AncestryMethod.COVARIATE_AWARE.value == "covariate_aware"

    def test_stratified_value(self):
        assert AncestryMethod.STRATIFIED.value == "stratified"

    def test_from_string(self):
        method = AncestryMethod("regress_out")
        assert method == AncestryMethod.REGRESS_OUT


# =============================================================================
# TEST: compute_ancestry_pcs
# =============================================================================


class TestComputeAncestryPCs:
    """Tests for compute_ancestry_pcs function."""

    def test_basic_computation(self, genotype_matrix_3pop):
        result = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)
        assert isinstance(result, AncestryPCs)
        assert result.n_components == 5
        assert result.components.shape == (60, 5)
        assert len(result.sample_ids) == 60

    def test_variance_explained_decreasing(self, genotype_matrix_3pop):
        result = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)
        ratios = result.explained_variance_ratio
        for i in range(len(ratios) - 1):
            assert ratios[i] >= ratios[i + 1]

    def test_variance_explained_positive(self, genotype_matrix_3pop):
        result = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)
        assert all(v >= 0 for v in result.explained_variance_ratio)
        assert np.sum(result.explained_variance_ratio) <= 1.0

    def test_n_components_capping(self):
        """n_components should be capped at min(n_samples, n_variants)."""
        rng = np.random.RandomState(42)
        small_matrix = pd.DataFrame(
            rng.binomial(2, 0.3, (10, 5)),
            index=[f"S{i}" for i in range(10)],
            columns=[f"V{i}" for i in range(5)],
        )
        result = compute_ancestry_pcs(small_matrix, n_components=20, seed=42)
        assert result.n_components <= 5

    def test_monomorphic_variants(self):
        """Monomorphic (zero-variance) variants should not crash PCA."""
        rng = np.random.RandomState(42)
        n_samples = 30
        n_variants = 50
        genotypes = rng.binomial(2, 0.2, (n_samples, n_variants))
        # Make some columns all-zero (monomorphic)
        genotypes[:, :10] = 0

        df = pd.DataFrame(
            genotypes,
            index=[f"S{i}" for i in range(n_samples)],
            columns=[f"V{i}" for i in range(n_variants)],
        )
        result = compute_ancestry_pcs(df, n_components=5, seed=42)
        assert result.n_components == 5
        assert not np.any(np.isnan(result.components.values))

    def test_reproducibility(self, genotype_matrix_3pop):
        r1 = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)
        r2 = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)
        np.testing.assert_array_almost_equal(r1.components.values, r2.components.values)

    def test_different_seeds(self, genotype_matrix_3pop):
        r1 = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=42)
        r2 = compute_ancestry_pcs(genotype_matrix_3pop, n_components=5, seed=99)
        # PCA is deterministic for exact same data, but seeds may differ slightly
        # depending on SVD solver. Check that results are produced.
        assert r1.components.shape == r2.components.shape

    def test_nan_handling(self):
        """NaN genotype values should be handled gracefully."""
        rng = np.random.RandomState(42)
        genotypes = rng.binomial(2, 0.3, (20, 30)).astype(float)
        genotypes[0, 0] = np.nan
        genotypes[5, 10] = np.nan

        df = pd.DataFrame(
            genotypes,
            index=[f"S{i}" for i in range(20)],
            columns=[f"V{i}" for i in range(30)],
        )
        result = compute_ancestry_pcs(df, n_components=3, seed=42)
        assert not np.any(np.isnan(result.components.values))

    def test_to_dict(self, ancestry_pcs_3pop):
        d = ancestry_pcs_3pop.to_dict()
        assert d["n_components"] == 5
        assert d["n_samples"] == 60
        assert len(d["explained_variance_ratio"]) == 5
        assert 0.0 <= d["total_variance_explained"] <= 1.0


# =============================================================================
# TEST: adjust_pathway_scores
# =============================================================================


class TestAdjustPathwayScores:
    """Tests for adjust_pathway_scores function."""

    def test_regress_out_basic(self, confounded_pathway_scores, ancestry_pcs_3pop):
        result = adjust_pathway_scores(
            confounded_pathway_scores,
            ancestry_pcs_3pop,
            method=AncestryMethod.REGRESS_OUT,
        )
        assert isinstance(result, AncestryAdjustmentResult)
        assert result.method == AncestryMethod.REGRESS_OUT
        assert result.adjusted_scores.shape == confounded_pathway_scores.shape

    def test_r_squared_in_range(self, confounded_pathway_scores, ancestry_pcs_3pop):
        result = adjust_pathway_scores(confounded_pathway_scores, ancestry_pcs_3pop)
        for pathway, r2 in result.r_squared_per_pathway.items():
            assert 0.0 <= r2 <= 1.0, f"R^2 for {pathway} out of range: {r2}"

    def test_confounded_pathways_detected(self, confounded_pathway_scores, ancestry_pcs_3pop):
        result = adjust_pathway_scores(confounded_pathway_scores, ancestry_pcs_3pop)
        # Pathways 0, 1, 2 should have higher R^2 since they're confounded
        assert len(result.highly_confounded_pathways) > 0

    def test_covariate_aware_method(self, confounded_pathway_scores, ancestry_pcs_3pop):
        result = adjust_pathway_scores(
            confounded_pathway_scores,
            ancestry_pcs_3pop,
            method=AncestryMethod.COVARIATE_AWARE,
        )
        assert result.method == AncestryMethod.COVARIATE_AWARE
        assert result.adjusted_scores.shape == confounded_pathway_scores.shape

    def test_stratified_raises_error(self, confounded_pathway_scores, ancestry_pcs_3pop):
        with pytest.raises(ValueError, match="STRATIFIED"):
            adjust_pathway_scores(
                confounded_pathway_scores,
                ancestry_pcs_3pop,
                method=AncestryMethod.STRATIFIED,
            )

    def test_partial_sample_overlap(self, ancestry_pcs_3pop):
        """Should handle pathway scores with a subset of samples."""
        rng = np.random.RandomState(42)
        # Only use first 40 samples in pathway scores
        sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, 41)]
        scores = pd.DataFrame(
            rng.normal(0, 1, (40, 3)),
            index=sample_ids,
            columns=["P1", "P2", "P3"],
        )
        result = adjust_pathway_scores(scores, ancestry_pcs_3pop)
        assert result.adjusted_scores.shape[0] == 40

    def test_no_common_samples_raises(self, ancestry_pcs_3pop):
        """Should raise if no samples overlap."""
        scores = pd.DataFrame(
            np.random.normal(0, 1, (5, 3)),
            index=["X1", "X2", "X3", "X4", "X5"],
            columns=["P1", "P2", "P3"],
        )
        with pytest.raises(ValueError, match="No common samples"):
            adjust_pathway_scores(scores, ancestry_pcs_3pop)

    def test_n_pcs_parameter(self, confounded_pathway_scores, ancestry_pcs_3pop):
        result = adjust_pathway_scores(confounded_pathway_scores, ancestry_pcs_3pop, n_pcs=2)
        assert result.n_pcs_used == 2

    def test_adjusted_scores_normalized(self, confounded_pathway_scores, ancestry_pcs_3pop):
        """Adjusted scores should be approximately z-normalized."""
        result = adjust_pathway_scores(confounded_pathway_scores, ancestry_pcs_3pop)
        means = result.adjusted_scores.mean()
        stds = result.adjusted_scores.std()
        for pathway in means.index:
            assert abs(means[pathway]) < 0.1, f"Mean of {pathway} not near 0"
            assert abs(stds[pathway] - 1.0) < 0.2, f"Std of {pathway} not near 1"


# =============================================================================
# TEST: check_ancestry_independence
# =============================================================================


class TestAncestryIndependence:
    """Tests for check_ancestry_independence function."""

    def test_confounded_clusters_fail(self, ancestry_pcs_3pop):
        """Clusters aligned with ancestry groups should fail independence."""
        # Labels that match ancestry groups exactly
        labels = np.array([0] * 20 + [1] * 20 + [2] * 20)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop)
        assert isinstance(report, AncestryStratificationReport)
        assert not report.overall_independence_passed

    def test_random_clusters_pass(self, ancestry_pcs_3pop):
        """Random cluster labels should pass independence test."""
        rng = np.random.RandomState(42)
        labels = rng.randint(0, 3, size=60)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop)
        assert report.overall_independence_passed

    def test_single_cluster(self, ancestry_pcs_3pop):
        """Single cluster should pass (no comparison possible -> p=1)."""
        labels = np.zeros(60, dtype=int)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop)
        assert report.overall_independence_passed

    def test_pvalues_in_range(self, ancestry_pcs_3pop):
        labels = np.array([0] * 20 + [1] * 20 + [2] * 20)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop)
        for pval in report.independence_test_pvalues.values():
            assert 0.0 <= pval <= 1.0

    def test_n_pcs_to_test(self, ancestry_pcs_3pop):
        labels = np.array([0] * 20 + [1] * 20 + [2] * 20)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop, n_pcs_to_test=2)
        assert len(report.independence_test_pvalues) == 2

    def test_to_dict(self, ancestry_pcs_3pop):
        labels = np.array([0] * 20 + [1] * 20 + [2] * 20)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop)
        d = report.to_dict()
        assert "overall_independence_passed" in d
        assert "independence_test_pvalues" in d
        assert "n_pcs_tested" in d

    def test_format_report(self, ancestry_pcs_3pop):
        labels = np.array([0] * 20 + [1] * 20 + [2] * 20)
        report = check_ancestry_independence(labels, ancestry_pcs_3pop)
        text = report.format_report()
        assert "Ancestry Independence Test" in text
        assert "Kruskal-Wallis" in text


# =============================================================================
# TEST: stratified_analysis
# =============================================================================


class TestStratifiedAnalysis:
    """Tests for stratified_analysis function."""

    def test_basic_operation(self, non_confounded_pathway_scores):
        ancestry_groups = np.array([0] * 20 + [1] * 20 + [2] * 20)
        result = stratified_analysis(
            non_confounded_pathway_scores,
            ancestry_groups,
            n_clusters=2,
            seed=42,
        )
        assert result["n_groups"] == 3
        assert result["groups_analyzed"] == 3
        assert "cross_group_concordance" in result

    def test_small_group_skipped(self, non_confounded_pathway_scores):
        """Groups with too few samples should be skipped."""
        # Group 2 has only 5 samples (less than 2*3=6 needed for 2 clusters)
        ancestry_groups = np.array([0] * 27 + [1] * 28 + [2] * 5)
        result = stratified_analysis(
            non_confounded_pathway_scores,
            ancestry_groups,
            n_clusters=2,
            seed=42,
        )
        assert result["groups_analyzed"] == 2  # Group 2 skipped

    def test_concordance_present(self, non_confounded_pathway_scores):
        ancestry_groups = np.array([0] * 30 + [1] * 30)
        result = stratified_analysis(
            non_confounded_pathway_scores,
            ancestry_groups,
            n_clusters=2,
            seed=42,
        )
        concordance = result["cross_group_concordance"]
        assert "concordance_score" in concordance
        assert "sufficient_groups" in concordance


# =============================================================================
# TEST: compute_ancestry_correlation
# =============================================================================


class TestAncestryCorrelation:
    """Tests for compute_ancestry_correlation function."""

    def test_output_shape(self, confounded_pathway_scores, ancestry_pcs_3pop):
        corr = compute_ancestry_correlation(confounded_pathway_scores, ancestry_pcs_3pop, n_pcs=3)
        assert corr.shape == (5, 3)  # 5 pathways x 3 PCs

    def test_values_in_range(self, confounded_pathway_scores, ancestry_pcs_3pop):
        corr = compute_ancestry_correlation(confounded_pathway_scores, ancestry_pcs_3pop)
        assert (corr.abs() <= 1.0).all().all()

    def test_confounded_pathways_high_correlation(
        self, confounded_pathway_scores, ancestry_pcs_3pop
    ):
        """Confounded pathways should show higher correlation with ancestry PCs."""
        corr = compute_ancestry_correlation(confounded_pathway_scores, ancestry_pcs_3pop)
        # At least one pathway should have |r| > 0.3 with at least one PC
        max_abs_corr = corr.abs().max().max()
        assert max_abs_corr > 0.3

    def test_too_few_samples_raises(self, ancestry_pcs_3pop):
        """Should raise with fewer than 3 common samples."""
        scores = pd.DataFrame(
            np.random.normal(0, 1, (2, 3)),
            index=["SAMPLE_001", "SAMPLE_002"],
            columns=["P1", "P2", "P3"],
        )
        with pytest.raises(ValueError, match="at least 3"):
            compute_ancestry_correlation(scores, ancestry_pcs_3pop)


# =============================================================================
# TEST: Dataclass serialization
# =============================================================================


class TestDataclasses:
    """Tests for dataclass to_dict, format_report, and get_citations methods."""

    def test_adjustment_result_to_dict(self):
        result = AncestryAdjustmentResult(
            adjusted_scores=pd.DataFrame({"P1": [1, 2], "P2": [3, 4]}),
            method=AncestryMethod.REGRESS_OUT,
            r_squared_per_pathway={"P1": 0.15, "P2": 0.03},
            n_pcs_used=5,
            highly_confounded_pathways=["P1"],
        )
        d = result.to_dict()
        assert d["method"] == "regress_out"
        assert d["n_pcs_used"] == 5
        assert d["n_highly_confounded"] == 1
        assert d["mean_r_squared"] == round((0.15 + 0.03) / 2, 4)

    def test_adjustment_result_format_report(self):
        result = AncestryAdjustmentResult(
            adjusted_scores=pd.DataFrame({"P1": [1, 2]}),
            method=AncestryMethod.REGRESS_OUT,
            r_squared_per_pathway={"P1": 0.15},
            n_pcs_used=5,
            highly_confounded_pathways=["P1"],
        )
        report = result.format_report()
        assert "Ancestry Correction" in report
        assert "regress_out" in report
        assert "P1" in report

    def test_adjustment_result_citations(self):
        result = AncestryAdjustmentResult(
            adjusted_scores=pd.DataFrame(),
            method=AncestryMethod.REGRESS_OUT,
            r_squared_per_pathway={},
            n_pcs_used=5,
        )
        citations = result.get_citations()
        assert len(citations) == 2
        assert "Price" in citations[0]
        assert "Patterson" in citations[1]

    def test_stratification_report_to_dict(self):
        report = AncestryStratificationReport(
            ancestry_cluster_correlation={"PC1": {"0": 0.5, "1": -0.3}},
            independence_test_pvalues={"PC1": 0.01, "PC2": 0.5},
            overall_independence_passed=False,
        )
        d = report.to_dict()
        assert d["overall_independence_passed"] is False
        assert d["n_pcs_tested"] == 2
        assert d["n_significant_pcs"] == 1

    def test_stratification_report_format_report(self):
        report = AncestryStratificationReport(
            ancestry_cluster_correlation={},
            independence_test_pvalues={"PC1": 0.01},
            overall_independence_passed=False,
        )
        text = report.format_report()
        assert "FAIL" in text

    def test_empty_r_squared(self):
        """AncestryAdjustmentResult.to_dict with empty R^2 should not crash."""
        result = AncestryAdjustmentResult(
            adjusted_scores=pd.DataFrame(),
            method=AncestryMethod.REGRESS_OUT,
            r_squared_per_pathway={},
            n_pcs_used=5,
        )
        d = result.to_dict()
        assert d["mean_r_squared"] == 0.0
        assert d["max_r_squared"] == 0.0


# =============================================================================
# TEST: End-to-end correction
# =============================================================================


class TestEndToEnd:
    """End-to-end tests: ancestry correction improves subtype recovery."""

    def test_correction_reduces_ancestry_correlation(
        self, confounded_pathway_scores, ancestry_pcs_3pop
    ):
        """After correction, pathway-ancestry correlations should decrease."""
        corr_before = compute_ancestry_correlation(confounded_pathway_scores, ancestry_pcs_3pop)
        result = adjust_pathway_scores(confounded_pathway_scores, ancestry_pcs_3pop)
        corr_after = compute_ancestry_correlation(result.adjusted_scores, ancestry_pcs_3pop)
        # Mean absolute correlation should decrease
        mean_before = corr_before.abs().mean().mean()
        mean_after = corr_after.abs().mean().mean()
        assert mean_after < mean_before

    def test_non_confounded_not_harmed(self, non_confounded_pathway_scores, ancestry_pcs_3pop):
        """Correction on non-confounded data should not destroy signal."""
        result = adjust_pathway_scores(non_confounded_pathway_scores, ancestry_pcs_3pop)
        # Variance should still exist after correction
        stds = result.adjusted_scores.std()
        assert (stds > 0.5).all()
