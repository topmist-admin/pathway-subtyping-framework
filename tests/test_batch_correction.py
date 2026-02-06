"""Tests for the batch_correction module."""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.batch_correction import (
    BatchCorrectionMethod,
    BatchCorrectionResult,
    BatchEffectReport,
    correct_batch_effects,
    detect_batch_effects,
    validate_batch_correction,
)

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def pathway_scores_with_batch():
    """Create pathway scores with known batch effects."""
    np.random.seed(42)
    n_per_batch = 30
    n_pathways = 5

    scores = np.zeros((n_per_batch * 3, n_pathways))

    # Batch 1: shifted up in pathways 0, 1
    scores[:30, 0] = np.random.normal(2.0, 0.5, 30)
    scores[:30, 1] = np.random.normal(1.5, 0.5, 30)
    scores[:30, 2:] = np.random.normal(0, 0.5, (30, 3))

    # Batch 2: shifted down in pathways 0, 1
    scores[30:60, 0] = np.random.normal(-1.0, 0.5, 30)
    scores[30:60, 1] = np.random.normal(-0.5, 0.5, 30)
    scores[30:60, 2:] = np.random.normal(0, 0.5, (30, 3))

    # Batch 3: centered
    scores[60:90, :] = np.random.normal(0, 0.5, (30, 5))

    sample_ids = [f"SAMPLE_{i:03d}" for i in range(90)]
    pathway_names = [f"PATHWAY_{i}" for i in range(n_pathways)]

    df = pd.DataFrame(scores, index=sample_ids, columns=pathway_names)
    batch_labels = np.array([0] * 30 + [1] * 30 + [2] * 30)

    return df, batch_labels


@pytest.fixture
def pathway_scores_no_batch():
    """Create pathway scores without batch effects."""
    np.random.seed(42)
    scores = np.random.normal(0, 1, (60, 4))
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(60)]
    pathway_names = ["P_A", "P_B", "P_C", "P_D"]
    df = pd.DataFrame(scores, index=sample_ids, columns=pathway_names)
    batch_labels = np.array([0] * 30 + [1] * 30)
    return df, batch_labels


@pytest.fixture
def biological_labels():
    """Known biological group labels (3 groups across 90 samples)."""
    return np.array([0] * 30 + [1] * 30 + [2] * 30)


# =============================================================================
# TEST: BatchCorrectionMethod enum
# =============================================================================


class TestBatchCorrectionMethod:
    def test_enum_values(self):
        assert BatchCorrectionMethod.COMBAT.value == "combat"
        assert BatchCorrectionMethod.MEAN_CENTER.value == "mean_center"
        assert BatchCorrectionMethod.STANDARDIZE.value == "standardize"

    def test_enum_count(self):
        assert len(BatchCorrectionMethod) == 3


# =============================================================================
# TEST: detect_batch_effects
# =============================================================================


class TestDetectBatchEffects:
    def test_detects_known_batch_effects(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels)

        assert report.overall_batch_effect is True
        assert report.n_batches == 3
        assert len(report.significant_features) > 0
        assert "PATHWAY_0" in report.significant_features

    def test_no_batch_effects(self, pathway_scores_no_batch):
        scores, batch_labels = pathway_scores_no_batch
        report = detect_batch_effects(scores, batch_labels)

        # With random data, should have few or no significant features
        assert report.n_batches == 2

    def test_batch_sizes(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels)

        assert report.batch_sizes["0"] == 30
        assert report.batch_sizes["1"] == 30
        assert report.batch_sizes["2"] == 30

    def test_variance_explained_range(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels)

        for feat, eta2 in report.variance_explained.items():
            assert 0 <= eta2 <= 1

    def test_f_statistics_positive(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels)

        for feat, f_stat in report.f_statistics.items():
            assert f_stat >= 0

    def test_mismatched_lengths(self, pathway_scores_with_batch):
        scores, _ = pathway_scores_with_batch
        with pytest.raises(ValueError, match="must match"):
            detect_batch_effects(scores, np.array([0, 1, 2]))

    def test_single_batch_error(self, pathway_scores_with_batch):
        scores, _ = pathway_scores_with_batch
        with pytest.raises(ValueError, match="at least 2"):
            detect_batch_effects(scores, np.zeros(len(scores), dtype=int))

    def test_custom_threshold(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        strict = detect_batch_effects(scores, batch_labels, significance_threshold=0.001)
        lenient = detect_batch_effects(scores, batch_labels, significance_threshold=0.5)

        assert len(strict.significant_features) <= len(lenient.significant_features)

    def test_custom_batch_variable_name(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels, batch_variable="sequencing_site")
        assert report.batch_variable == "sequencing_site"


# =============================================================================
# TEST: correct_batch_effects
# =============================================================================


class TestCorrectBatchEffects:
    def test_combat_reduces_batch_variance(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels, method=BatchCorrectionMethod.COMBAT)

        assert isinstance(result, BatchCorrectionResult)
        mean_pre = np.mean(list(result.pre_correction_variance.values()))
        mean_post = np.mean(list(result.post_correction_variance.values()))
        assert mean_post < mean_pre

    def test_mean_center_reduces_batch_variance(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(
            scores, batch_labels, method=BatchCorrectionMethod.MEAN_CENTER
        )

        mean_pre = np.mean(list(result.pre_correction_variance.values()))
        mean_post = np.mean(list(result.post_correction_variance.values()))
        assert mean_post < mean_pre

    def test_standardize_reduces_batch_variance(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(
            scores, batch_labels, method=BatchCorrectionMethod.STANDARDIZE
        )

        mean_pre = np.mean(list(result.pre_correction_variance.values()))
        mean_post = np.mean(list(result.post_correction_variance.values()))
        assert mean_post < mean_pre

    def test_output_shape_preserved(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        assert result.corrected_scores.shape == scores.shape
        assert list(result.corrected_scores.columns) == list(scores.columns)
        assert list(result.corrected_scores.index) == list(scores.index)

    def test_no_nan_in_output(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        for method in BatchCorrectionMethod:
            result = correct_batch_effects(scores, batch_labels, method=method)
            assert not result.corrected_scores.isna().any().any(), f"NaN in {method.value} output"

    def test_mismatched_lengths(self, pathway_scores_with_batch):
        scores, _ = pathway_scores_with_batch
        with pytest.raises(ValueError, match="must match"):
            correct_batch_effects(scores, np.array([0, 1]))

    def test_single_batch_error(self, pathway_scores_with_batch):
        scores, _ = pathway_scores_with_batch
        with pytest.raises(ValueError, match="at least 2"):
            correct_batch_effects(scores, np.zeros(len(scores), dtype=int))

    def test_variance_reduction_dict(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        for col in scores.columns:
            assert col in result.variance_reduction
            # Variance reduction should be non-negative for affected features
            assert isinstance(result.variance_reduction[col], float)


# =============================================================================
# TEST: validate_batch_correction
# =============================================================================


class TestValidateBatchCorrection:
    def test_basic_validation(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        validation = validate_batch_correction(scores, result.corrected_scores, batch_labels)

        assert "batch_variance_before" in validation
        assert "batch_variance_after" in validation
        assert "batch_variance_reduced" in validation
        assert "mean_correlation_with_original" in validation

    def test_batch_variance_reduced(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        validation = validate_batch_correction(scores, result.corrected_scores, batch_labels)

        assert validation["batch_variance_reduced"]

    def test_with_biological_labels(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        bio_labels = np.array([0] * 15 + [1] * 15 + [0] * 15 + [1] * 15 + [0] * 15 + [1] * 15)
        result = correct_batch_effects(scores, batch_labels)

        validation = validate_batch_correction(
            scores, result.corrected_scores, batch_labels, bio_labels
        )

        assert "biological_variance_before" in validation
        assert "biological_variance_after" in validation
        assert "biological_signal_preserved" in validation

    def test_correlation_positive(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        validation = validate_batch_correction(scores, result.corrected_scores, batch_labels)

        assert validation["mean_correlation_with_original"] > 0


# =============================================================================
# TEST: DATACLASSES
# =============================================================================


class TestDataclasses:
    def test_batch_effect_report_to_dict(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels)

        d = report.to_dict()
        assert "batch_variable" in d
        assert "n_batches" in d
        assert "n_significant_features" in d
        assert "overall_batch_effect" in d

    def test_batch_effect_report_format_report(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        report = detect_batch_effects(scores, batch_labels)

        text = report.format_report()
        assert "Batch Effect Detection Report" in text
        assert "batch" in text.lower()

    def test_correction_result_to_dict(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        d = result.to_dict()
        assert d["method"] == "combat"
        assert "mean_variance_reduction" in d

    def test_correction_result_format_report(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        text = result.format_report()
        assert "Batch Correction Report" in text

    def test_correction_result_get_citations(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(scores, batch_labels)

        citations = result.get_citations()
        assert len(citations) >= 1
        assert any("Johnson" in c for c in citations)

    def test_mean_center_citations(self, pathway_scores_with_batch):
        scores, batch_labels = pathway_scores_with_batch
        result = correct_batch_effects(
            scores, batch_labels, method=BatchCorrectionMethod.MEAN_CENTER
        )

        citations = result.get_citations()
        assert any("Leek" in c for c in citations)
        # MEAN_CENTER should not cite Johnson (ComBat-specific)
        assert not any("Johnson" in c for c in citations)


# =============================================================================
# TEST: EDGE CASES
# =============================================================================


class TestEdgeCases:
    def test_two_batches(self):
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.normal(0, 1, (40, 3)),
            columns=["A", "B", "C"],
        )
        batch_labels = np.array([0] * 20 + [1] * 20)

        report = detect_batch_effects(scores, batch_labels)
        assert report.n_batches == 2

        result = correct_batch_effects(scores, batch_labels)
        assert result.corrected_scores.shape == (40, 3)

    def test_unequal_batch_sizes(self):
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.normal(0, 1, (50, 3)),
            columns=["A", "B", "C"],
        )
        batch_labels = np.array([0] * 10 + [1] * 40)

        result = correct_batch_effects(scores, batch_labels)
        assert result.corrected_scores.shape == (50, 3)
        assert not result.corrected_scores.isna().any().any()

    def test_string_batch_labels(self, pathway_scores_with_batch):
        scores, _ = pathway_scores_with_batch
        batch_labels = np.array(["site_A"] * 30 + ["site_B"] * 30 + ["site_C"] * 30)

        report = detect_batch_effects(scores, batch_labels)
        assert report.n_batches == 3

        result = correct_batch_effects(scores, batch_labels)
        assert result.corrected_scores.shape == scores.shape

    def test_zero_variance_column(self):
        np.random.seed(42)
        scores = pd.DataFrame(
            {
                "A": np.random.normal(0, 1, 40),
                "B": np.zeros(40),
                "C": np.random.normal(0, 1, 40),
            }
        )
        batch_labels = np.array([0] * 20 + [1] * 20)

        # Should not crash
        report = detect_batch_effects(scores, batch_labels)
        assert "B" in report.f_statistics

        result = correct_batch_effects(scores, batch_labels)
        assert not result.corrected_scores.isna().any().any()
