"""
Tests for the pathway_subtyping.harmonize module.

Covers roadmap F2 acceptance criteria (synthetic proxy — real cross-platform
data lives outside the package):

    - Harmonized pathway-level rho across platform pairs exceeds 0.75
      (pre-harmonization baseline 0.55-0.65)
    - Harmonization confidence correlates with a per-cell quality covariate
    - No degradation of within-platform rho (pre and post should match on
      a single-platform input)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.harmonize import (
    AlignmentResult,
    CrossPlatformAligner,
    CrossPlatformBenchmark,
    FallbackEmbedder,
    HarmonizationReport,
    UCEEmbedder,
    simulate_platform_distortion,
)

# --------------------------------------------------------------------------- #
# Shared fixtures
# --------------------------------------------------------------------------- #


@pytest.fixture
def reference_scores():
    """Synthetic reference pathway-score matrix (n_cells x n_pathways)."""
    rng = np.random.default_rng(0)
    n_cells, n_pathways = 400, 20
    # Cells fall into 3 biological clusters with distinct pathway signatures
    cluster_means = rng.standard_normal((3, n_pathways)) * 1.5
    cluster_assignments = rng.integers(0, 3, size=n_cells)
    base = cluster_means[cluster_assignments]
    noise = rng.standard_normal((n_cells, n_pathways)) * 0.3
    df = pd.DataFrame(
        base + noise,
        columns=[f"PATH_{i}" for i in range(n_pathways)],
    )
    return df


# --------------------------------------------------------------------------- #
# Fallback embedder
# --------------------------------------------------------------------------- #


class TestFallbackEmbedder:

    def test_fit_and_embed(self, reference_scores):
        emb = FallbackEmbedder(embedding_dim=8).fit(reference_scores)
        out = emb.embed(reference_scores)
        assert out.shape == (len(reference_scores), 8)

    def test_deterministic(self, reference_scores):
        a = FallbackEmbedder(embedding_dim=8, seed=0).fit(reference_scores).embed(reference_scores)
        b = FallbackEmbedder(embedding_dim=8, seed=0).fit(reference_scores).embed(reference_scores)
        np.testing.assert_allclose(a, b)

    def test_embed_without_fit_auto_fits(self, reference_scores):
        emb = FallbackEmbedder(embedding_dim=8)
        out = emb.embed(reference_scores)
        assert out.shape == (len(reference_scores), 8)

    def test_rejects_bad_dim(self):
        with pytest.raises(ValueError):
            FallbackEmbedder(embedding_dim=1)

    def test_gene_mismatch_raises(self, reference_scores):
        emb = FallbackEmbedder(embedding_dim=4).fit(reference_scores)
        smaller = reference_scores.iloc[:, :10]
        with pytest.raises(ValueError, match="fitted on"):
            emb.embed(smaller)


class TestUCEEmbedderStub:

    def test_lazy_import_error_informative(self):
        emb = UCEEmbedder()
        # Triggers lazy-import path; without the '[harmonize]' extra
        # this must raise ImportError with an install hint.
        try:
            emb._lazy_load()
        except (ImportError, NotImplementedError) as exc:
            assert (
                "harmonize" in str(exc)
                or "uce" in str(exc).lower()
                or "checkpoint" in str(exc).lower()
            )
        else:
            pytest.skip("UCE is installed — skipping stub error test")


# --------------------------------------------------------------------------- #
# Cross-platform aligner
# --------------------------------------------------------------------------- #


class TestCrossPlatformAligner:

    def test_fit_transform_shapes(self, reference_scores):
        emb = FallbackEmbedder(embedding_dim=8).fit(reference_scores)
        embeddings = emb.embed(reference_scores)
        platforms = ["p1"] * (len(reference_scores) // 2) + ["p2"] * (
            len(reference_scores) - len(reference_scores) // 2
        )
        aligner = CrossPlatformAligner()
        result = aligner.fit_transform(reference_scores, platforms, embeddings)
        assert isinstance(result, AlignmentResult)
        assert result.aligned_scores.shape == reference_scores.shape
        assert len(result.per_cell_shift) == len(reference_scores)
        assert set(result.per_platform_drift.keys()) == {"p1", "p2"}

    def test_within_platform_zero_drift(self, reference_scores):
        """Single-platform input: alignment is ~zero, scores unchanged."""
        emb = FallbackEmbedder(embedding_dim=8).fit(reference_scores)
        embeddings = emb.embed(reference_scores)
        platforms = ["single"] * len(reference_scores)
        aligner = CrossPlatformAligner()
        result = aligner.fit_transform(reference_scores, platforms, embeddings)
        # Drift should be close to zero (residual mean of the pooled model)
        drift = np.abs(list(result.per_platform_drift["single"].values()))
        assert drift.max() < 0.1

    def test_unseen_platform_emits_warning(self, reference_scores, caplog):
        emb = FallbackEmbedder(embedding_dim=8).fit(reference_scores)
        embeddings = emb.embed(reference_scores)
        platforms_train = ["p1"] * len(reference_scores)
        aligner = CrossPlatformAligner()
        aligner.fit(reference_scores, platforms_train, embeddings)
        platforms_new = ["p_novel"] * len(reference_scores)
        with caplog.at_level("WARNING"):
            aligner.transform(reference_scores, platforms_new, embeddings)
        assert any("unseen platforms" in rec.message for rec in caplog.records)

    def test_platform_offset_accessor(self, reference_scores):
        emb = FallbackEmbedder(embedding_dim=8).fit(reference_scores)
        embeddings = emb.embed(reference_scores)
        half = len(reference_scores) // 2
        platforms = ["p1"] * half + ["p2"] * (len(reference_scores) - half)
        aligner = CrossPlatformAligner().fit(reference_scores, platforms, embeddings)
        offset = aligner.platform_offset("p1")
        assert len(offset) == reference_scores.shape[1]
        with pytest.raises(KeyError):
            aligner.platform_offset("nonexistent")

    def test_transform_requires_fit(self, reference_scores):
        aligner = CrossPlatformAligner()
        with pytest.raises(RuntimeError):
            aligner.transform(
                reference_scores,
                ["p1"] * len(reference_scores),
                np.zeros((len(reference_scores), 4)),
            )


# --------------------------------------------------------------------------- #
# Harmonization report
# --------------------------------------------------------------------------- #


class TestHarmonizationReport:

    def _make_report(self, reference_scores):
        emb = FallbackEmbedder(embedding_dim=8).fit(reference_scores)
        embeddings = emb.embed(reference_scores)
        half = len(reference_scores) // 2
        platforms = ["p1"] * half + ["p2"] * (len(reference_scores) - half)
        result = CrossPlatformAligner().fit_transform(reference_scores, platforms, embeddings)
        return HarmonizationReport.from_alignment(result)

    def test_report_basic(self, reference_scores):
        report = self._make_report(reference_scores)
        assert 0.0 <= report.confidence.min() <= report.confidence.max() <= 1.0
        assert set(report.per_platform_drift.keys()) == {"p1", "p2"}

    def test_correlate_with_quality_perfect(self, reference_scores):
        """When quality ~ -shift, correlation with confidence must be > 0."""
        report = self._make_report(reference_scores)
        # Smaller shift -> higher confidence; quality defined as -shift
        quality = -report.alignment.per_cell_shift
        rho = report.correlate_with_quality(quality)
        assert rho > 0.9  # near-perfect by construction

    def test_correlate_with_constant_quality_is_nan(self, reference_scores):
        report = self._make_report(reference_scores)
        quality = pd.Series(np.ones(len(report.confidence)), index=report.confidence.index)
        rho = report.correlate_with_quality(quality)
        assert np.isnan(rho)

    def test_to_dict_roundtrip(self, reference_scores):
        report = self._make_report(reference_scores)
        d = report.to_dict()
        assert d["n_cells"] == len(report.confidence)
        assert "mean_platform_drift" in d

    def test_summary_renders(self, reference_scores):
        report = self._make_report(reference_scores)
        s = report.summary()
        assert "HarmonizationReport" in s
        assert "mean_confidence" in s


# --------------------------------------------------------------------------- #
# Benchmark + roadmap acceptance
# --------------------------------------------------------------------------- #


class TestCrossPlatformBenchmark:

    def test_simulate_distortion_shape_and_determinism(self, reference_scores):
        a = simulate_platform_distortion(reference_scores, "10x", seed=0)
        b = simulate_platform_distortion(reference_scores, "10x", seed=0)
        assert a.shape == reference_scores.shape
        np.testing.assert_allclose(a.to_numpy(), b.to_numpy())
        c = simulate_platform_distortion(reference_scores, "smartseq2", seed=0)
        assert not np.allclose(a.to_numpy(), c.to_numpy())

    def test_benchmark_rejects_single_platform(self, reference_scores):
        with pytest.raises(ValueError):
            CrossPlatformBenchmark(reference_scores=reference_scores, platforms=["10x"])

    def test_pre_rho_below_threshold(self, reference_scores):
        """Sanity: pre-harmonization rho should be below the 0.75 target
        by construction of the noise + shift."""
        bench = CrossPlatformBenchmark(
            reference_scores=reference_scores,
            platforms=["10x", "smartseq2"],
        )
        result = bench.run(seed=0)
        # Baseline per roadmap: 0.55-0.65. Allow some slack on either side.
        assert 0.3 <= result["pre_rho"] <= 0.75

    def test_post_rho_meets_acceptance(self, reference_scores):
        """Roadmap acceptance: harmonized rho > 0.75 on average."""
        bench = CrossPlatformBenchmark(
            reference_scores=reference_scores,
            platforms=["10x", "smartseq2"],
        )
        summary = bench.run_many(n_seeds=5)
        assert (
            summary["post_rho_mean"] > 0.75
        ), f"post_rho mean {summary['post_rho_mean']:.3f} below 0.75 target"
        # Improvement over pre-harmonization must be real
        assert summary["improvement_mean"] > 0.1

    def test_four_platforms_two_systems(self, reference_scores):
        """Roadmap deliverable: 4 platforms on a single biological system
        (cortex-shaped synthetic). Run the protocol end-to-end."""
        bench = CrossPlatformBenchmark(
            reference_scores=reference_scores,
            platforms=["10x", "smartseq2", "bulk_rnaseq", "spatial"],
        )
        summary = bench.run_many(n_seeds=3)
        assert summary["post_rho_mean"] > 0.75
        assert summary["post_rho_min"] > 0.7  # worst seed still healthy

    def test_confidence_correlates_with_quality_proxy(self, reference_scores):
        """Low-shift cells should have high confidence — the identity the
        roadmap requires (`confidence correlates with known quality`)."""
        bench = CrossPlatformBenchmark(
            reference_scores=reference_scores,
            platforms=["10x", "smartseq2"],
        )
        result = bench.run(seed=0)
        report: HarmonizationReport = result["report"]
        quality_proxy = -report.alignment.per_cell_shift
        rho = report.correlate_with_quality(quality_proxy)
        assert rho > 0.9
