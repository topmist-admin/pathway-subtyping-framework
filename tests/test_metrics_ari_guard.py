"""
Regression tests for the guarded adjusted-Rand-index metrics.

These lock in the fix for the empty-input / single-true-cluster artifact that
produced 14 degenerate rows (13 spurious ARI=1.0) in the 47-dataset calibration
benchmark and drove the retracted R^2=0.889 adaptive-threshold model.
See CORRECTION_2026-07/ERRATUM_2026-07-08.md and KNOWN-ISSUES.md (NB17).
"""

import math

import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

from pathway_subtyping.utils.metrics import (
    ari_degenerate_reason,
    ari_with_validity,
    safe_adjusted_rand_score,
)


class TestDegenerateInputArtifact:
    """The raw sklearn behaviours these guards exist to neutralize."""

    def test_sklearn_returns_spurious_one_on_empty(self):
        # This is the bug: empty input -> "perfect" score.
        assert adjusted_rand_score([], []) == 1.0

    def test_sklearn_returns_spurious_one_on_single_true_cluster(self):
        assert adjusted_rand_score([0, 0, 0, 0], [0, 0, 0, 0]) == 1.0

    @pytest.mark.parametrize(
        "true,pred,reason_prefix",
        [
            ([], [], "empty_input"),
            ([0], [0], "too_few_samples"),
            ([0, 0, 0, 0], [0, 0, 0, 0], "degenerate_ground_truth"),
            ([0, 0, 0, 0], [0, 1, 2, 3], "degenerate_ground_truth"),
        ],
    )
    def test_guard_flags_degenerate(self, true, pred, reason_prefix):
        reason = ari_degenerate_reason(true, pred)
        assert reason is not None and reason.startswith(reason_prefix)
        value, valid, r = ari_with_validity(true, pred)
        assert valid is False
        assert math.isnan(value)
        assert safe_adjusted_rand_score(true, pred, default=-99.0) == -99.0


class TestGSE136196Fixture:
    """
    GEO GSE136196 was the 14th degenerate benchmark row: 1 true cluster, a
    multi-way prediction, and sklearn returned ARI = 0.0 (NOT 1.0). A guard
    keyed on the ARI *value* (== 1.0) would miss it; the structural guard
    (n_true_clusters < 2) must catch it. This test is the reason the fix keys
    on ground-truth structure, not the output value.
    """

    def test_gse136196_returns_zero_from_sklearn(self):
        true = [0] * 20  # single true cluster (the degeneracy)
        pred = [i % 4 for i in range(20)]  # multi-way prediction
        assert adjusted_rand_score(true, pred) == 0.0  # not 1.0 -> value guard fails

    def test_structural_guard_catches_gse136196(self):
        true = [0] * 20
        pred = [i % 4 for i in range(20)]
        reason = ari_degenerate_reason(true, pred)
        assert reason == "degenerate_ground_truth(n_true_clusters=1)"
        value, valid, _ = ari_with_validity(true, pred)
        assert valid is False and math.isnan(value)


class TestValidInputUnaffected:
    """A real ground truth must pass through unchanged."""

    @pytest.mark.parametrize(
        "true,pred",
        [
            ([0, 0, 1, 1], [1, 1, 0, 0]),  # perfect (label-swapped)
            ([0, 0, 1, 1], [0, 1, 0, 1]),  # poor agreement
            ([0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2]),
        ],
    )
    def test_matches_sklearn_when_valid(self, true, pred):
        expected = adjusted_rand_score(true, pred)
        value, valid, reason = ari_with_validity(true, pred)
        assert valid is True
        assert reason is None
        assert value == pytest.approx(expected)
        assert safe_adjusted_rand_score(true, pred) == pytest.approx(expected)


class TestLengthMismatchRaises:
    """A length mismatch is a caller bug, not degeneracy — it must raise."""

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="equal-length"):
            ari_with_validity([0, 0, 1], [0, 1])
        with pytest.raises(ValueError):
            safe_adjusted_rand_score([0, 0, 1], [0, 1])


class TestBenchmarkIntegration:
    """The benchmark scorer must surface ari=None (not 1.0) on degenerate truth."""

    def test_benchmark_single_true_cluster_gives_none(self):
        import pandas as pd

        from pathway_subtyping.benchmark import BenchmarkMethod, run_single_benchmark

        rng = np.random.default_rng(20260708)
        n = 40
        pathway_scores = pd.DataFrame(rng.normal(size=(n, 6)), columns=[f"PW{i}" for i in range(6)])
        gene_burdens = pd.DataFrame(rng.normal(size=(n, 10)), columns=[f"G{i}" for i in range(10)])
        true_labels = np.zeros(n, dtype=int)  # single true cluster = degenerate

        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            pathway_scores=pathway_scores,
            gene_burdens=gene_burdens,
            n_clusters=2,
            true_labels=true_labels,
            seed=20260708,
        )
        assert result.ari is None  # was a spurious 1.0 before the fix
        assert result.nmi is None
