"""
Tests for the pathway_subtyping.perturb module.

Covers v0.6 Phase 2 F5 acceptance criteria on synthetic data:

    - Perturbing a synthetic "master regulator" produces the
      directionally expected MSV shift (roadmap proxy for MYC/MECP2).
    - Screen ranks the strongest perturbation first.
    - Conformal intervals from F1 wrap the perturbed MSV pipeline.
    - No breaking changes to existing scoring APIs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.perturb import (
    FallbackPerturber,
    GeneformerPerturber,
    MSVFromEmbedding,
    OfficialBackend,
    PerturbationMode,
    PerturbationReport,
    PerturbationResult,
    PerturbationScreen,
    ScreenResult,
)


# --------------------------------------------------------------------------- #
# Shared fixtures
# --------------------------------------------------------------------------- #

@pytest.fixture
def synthetic_cohort():
    """Synthetic cohort with an identifiable master-regulator gene.

    3 clusters x 100 cells, 50 genes. ``MARKER_A`` is high in cluster 0
    and drives pathway_0 in the reference score. Knocking out MARKER_A
    should push cluster-0 cells' pathway_0 score down.
    """
    rng = np.random.default_rng(0)
    n_per = 100
    n_genes = 50
    cluster_profiles = np.zeros((3, n_genes))
    # MARKER_A is index 0; cluster 0 overexpresses it
    cluster_profiles[0, 0] = 5.0
    cluster_profiles[1, 1] = 5.0  # MARKER_B drives cluster 1
    cluster_profiles[2, 2] = 5.0  # MARKER_C drives cluster 2

    cells: list = []
    cluster_ids: list = []
    for c in range(3):
        base = cluster_profiles[c]
        noise = rng.normal(0, 1.0, size=(n_per, n_genes))
        cells.append(base + noise)
        cluster_ids.extend([c] * n_per)

    X = np.vstack(cells)
    gene_names = [f"MARKER_A"] + [f"MARKER_B"] + [f"MARKER_C"] + [
        f"GENE_{i}" for i in range(3, n_genes)
    ]
    expression = pd.DataFrame(X, columns=gene_names)

    # Reference pathway scores: pathway_k is mean of cluster-k markers
    # plus small noise. This gives a clean signal that
    # MARKER_A controls pathway_0.
    pathway_scores = pd.DataFrame({
        "pathway_0": X[:, 0],  # driven by MARKER_A
        "pathway_1": X[:, 1],  # driven by MARKER_B
        "pathway_2": X[:, 2],  # driven by MARKER_C
        "pathway_generic": X[:, 3:10].mean(axis=1),
    })
    return expression, pathway_scores, cluster_ids


# --------------------------------------------------------------------------- #
# FallbackPerturber
# --------------------------------------------------------------------------- #

class TestFallbackPerturber:

    def test_fit_and_embed_shape(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        perturber = FallbackPerturber(embedding_dim=16)
        out = perturber.embed(expr)
        assert out.shape == (len(expr), 16)

    def test_deterministic(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        a = FallbackPerturber(embedding_dim=16, seed=0).fit(expr).embed(expr)
        b = FallbackPerturber(embedding_dim=16, seed=0).fit(expr).embed(expr)
        np.testing.assert_allclose(a, b)

    def test_rejects_non_dataframe(self):
        with pytest.raises(TypeError):
            FallbackPerturber().fit(np.zeros((10, 5)))  # type: ignore[arg-type]

    def test_rejects_bad_dim(self):
        with pytest.raises(ValueError):
            FallbackPerturber(embedding_dim=1)

    def test_knockout_zeros_the_gene_column(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        perturber = FallbackPerturber(embedding_dim=16).fit(expr)
        baseline = perturber.embed(expr)
        perturbed = perturber.perturb(expr, "MARKER_A", PerturbationMode.KNOCKOUT)
        # With MARKER_A zeroed, embeddings of cluster 0 cells should move
        # most (their baseline expression of MARKER_A was largest)
        delta = np.linalg.norm(perturbed - baseline, axis=1)
        cluster0_delta = delta[:100].mean()
        cluster_others_delta = delta[100:].mean()
        assert cluster0_delta > cluster_others_delta

    def test_overexpress_amplifies_the_gene(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        perturber = FallbackPerturber(embedding_dim=16).fit(expr)
        baseline = perturber.embed(expr)
        perturbed = perturber.perturb(expr, "MARKER_A", PerturbationMode.OVEREXPRESS)
        delta = np.linalg.norm(perturbed - baseline, axis=1)
        # Overexpression should produce non-trivial change on all cells
        assert delta.mean() > 0.0

    def test_unknown_gene_raises(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        perturber = FallbackPerturber(embedding_dim=16).fit(expr)
        with pytest.raises(KeyError, match="NONEXISTENT"):
            perturber.perturb(expr, "NONEXISTENT", PerturbationMode.KNOCKOUT)


class TestOfficialBackend:

    def test_lazy_load_surfaces_informative_error(self):
        """Either missing deps OR missing model_directory must be informative."""
        backend = OfficialBackend()
        try:
            backend._lazy_load()
        except ImportError as exc:
            msg = str(exc).lower()
            assert "perturb" in msg or "geneformer" in msg
        except ValueError as exc:
            # With deps installed, constructor + lazy_load with no model_directory
            # should surface a clear message pointing at model_directory.
            msg = str(exc).lower()
            assert "model_directory" in msg or "checkpoint" in msg
        else:
            pytest.skip("loaded without error — checkpoint must be present")

    def test_backend_id_is_stable_and_config_sensitive(self):
        """Backend ID changes with checkpoint path + emb_mode + max_input_len."""
        a = OfficialBackend(model_directory="/tmp/ckpt", emb_mode="cls", max_input_len=2048)
        b = OfficialBackend(model_directory="/tmp/ckpt", emb_mode="cls", max_input_len=2048)
        c = OfficialBackend(model_directory="/tmp/ckpt", emb_mode="mean", max_input_len=2048)
        d = OfficialBackend(model_directory="/tmp/ckpt2", emb_mode="cls", max_input_len=2048)
        e = OfficialBackend(model_directory="/tmp/ckpt", emb_mode="cls", max_input_len=4096)
        assert a._backend_id() == b._backend_id()
        assert a._backend_id() != c._backend_id()
        assert a._backend_id() != d._backend_id()
        assert a._backend_id() != e._backend_id()

    def test_cache_lookup_miss_returns_none(self, tmp_path):
        """Cache lookup on an empty cache directory returns None (no error)."""
        import pandas as pd
        backend = OfficialBackend(
            model_directory="/tmp/ckpt", cache_dir=str(tmp_path),
        )
        expr = pd.DataFrame(
            [[1.0, 2.0, 3.0]], columns=["A", "B", "C"], index=["cell_0"],
        )
        assert backend._cache_lookup(expr) is None

    def test_cache_roundtrip(self, tmp_path):
        """Storing an embedding then looking it up returns the same array."""
        import numpy as np
        import pandas as pd
        backend = OfficialBackend(
            model_directory="/tmp/ckpt", cache_dir=str(tmp_path),
        )
        expr = pd.DataFrame(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            columns=["A", "B", "C"], index=["c0", "c1"],
        )
        emb = np.arange(2 * 8, dtype=float).reshape(2, 8)
        backend._cache_store(expr, emb)
        out = backend._cache_lookup(expr)
        assert out is not None
        np.testing.assert_array_equal(out, emb)

    def test_cache_invalidates_on_input_change(self, tmp_path):
        """Changing one expression value must produce a different cache key."""
        import numpy as np
        import pandas as pd
        backend = OfficialBackend(
            model_directory="/tmp/ckpt", cache_dir=str(tmp_path),
        )
        expr_a = pd.DataFrame(
            [[1.0, 2.0]], columns=["A", "B"], index=["c0"],
        )
        expr_b = pd.DataFrame(
            [[1.0, 2.01]], columns=["A", "B"], index=["c0"],
        )
        backend._cache_store(expr_a, np.zeros((1, 4)))
        assert backend._cache_lookup(expr_a) is not None
        assert backend._cache_lookup(expr_b) is None

    def test_cache_disabled_when_cache_dir_none(self, tmp_path):
        """No cache_dir → store is a no-op; lookup returns None."""
        import numpy as np
        import pandas as pd
        backend = OfficialBackend(model_directory="/tmp/ckpt", cache_dir=None)
        expr = pd.DataFrame([[1.0]], columns=["A"], index=["c0"])
        backend._cache_store(expr, np.zeros((1, 4)))
        assert backend._cache_lookup(expr) is None
        # Ensure nothing was written anywhere adjacent
        assert not any(tmp_path.iterdir())


# --------------------------------------------------------------------------- #
# GeneformerPerturber (high-level wrapper)
# --------------------------------------------------------------------------- #

class TestGeneformerPerturber:

    def test_default_backend_is_fallback(self):
        p = GeneformerPerturber()
        assert isinstance(p.backend, FallbackPerturber)

    def test_perturb_returns_result_with_delta(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        p = GeneformerPerturber(FallbackPerturber(embedding_dim=16))
        result = p.perturb(expr, "MARKER_A", PerturbationMode.KNOCKOUT)
        assert isinstance(result, PerturbationResult)
        assert result.baseline_embedding.shape == result.perturbed_embedding.shape
        assert result.delta_embedding.shape == result.baseline_embedding.shape
        assert result.mean_l2_impact > 0.0

    def test_accepts_string_mode(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        p = GeneformerPerturber(FallbackPerturber(embedding_dim=16))
        result = p.perturb(expr, "MARKER_A", "knockout")
        assert result.mode == PerturbationMode.KNOCKOUT

    def test_to_dict_fields(self, synthetic_cohort):
        expr, _, _ = synthetic_cohort
        p = GeneformerPerturber(FallbackPerturber(embedding_dim=16))
        r = p.perturb(expr, "MARKER_A", PerturbationMode.KNOCKOUT).to_dict()
        assert r["gene"] == "MARKER_A"
        assert r["mode"] == "knockout"


# --------------------------------------------------------------------------- #
# MSVFromEmbedding
# --------------------------------------------------------------------------- #

class TestMSVFromEmbedding:

    def test_fit_transform_round_trip(self, synthetic_cohort):
        expr, scores, _ = synthetic_cohort
        p = FallbackPerturber(embedding_dim=16).fit(expr)
        emb = p.embed(expr)
        head = MSVFromEmbedding(ridge_alpha=1e-2).fit(emb, scores)
        pred = head.transform(emb)
        assert pred.shape == scores.shape
        assert list(pred.columns) == list(scores.columns)

    def test_fit_requires_matching_rows(self, synthetic_cohort):
        expr, scores, _ = synthetic_cohort
        emb = np.zeros((5, 16))
        with pytest.raises(ValueError, match="same number"):
            MSVFromEmbedding().fit(emb, scores)

    def test_transform_requires_fit(self):
        with pytest.raises(RuntimeError):
            MSVFromEmbedding().transform(np.zeros((3, 16)))

    def test_delta_helper(self, synthetic_cohort):
        expr, scores, _ = synthetic_cohort
        p = FallbackPerturber(embedding_dim=16).fit(expr)
        emb = p.embed(expr)
        head = MSVFromEmbedding().fit(emb, scores)
        perturbed = p.perturb(expr, "MARKER_A", PerturbationMode.KNOCKOUT)
        delta = head.delta(emb, perturbed)
        assert delta.shape == scores.shape

    def test_pathway_names_stored(self, synthetic_cohort):
        expr, scores, _ = synthetic_cohort
        p = FallbackPerturber(embedding_dim=16).fit(expr)
        emb = p.embed(expr)
        head = MSVFromEmbedding().fit(emb, scores)
        assert head.pathway_names == list(scores.columns)


# --------------------------------------------------------------------------- #
# PerturbationScreen
# --------------------------------------------------------------------------- #

class TestPerturbationScreen:

    def _make_screen(self, synthetic_cohort):
        expr, scores, _ = synthetic_cohort
        perturber = GeneformerPerturber(FallbackPerturber(embedding_dim=16))
        emb = perturber.embed(expr)
        head = MSVFromEmbedding().fit(emb, scores)
        return expr, scores, PerturbationScreen(perturber, head)

    def test_run_basic(self, synthetic_cohort):
        expr, scores, screen = self._make_screen(synthetic_cohort)
        result = screen.run(expr, ["MARKER_A", "MARKER_B", "GENE_10"])
        assert isinstance(result, ScreenResult)
        assert list(result.gene_panel) == ["MARKER_A", "MARKER_B", "GENE_10"]
        assert result.delta_msv_by_gene.shape == (3, scores.shape[1])

    def test_skip_missing(self, synthetic_cohort):
        expr, scores, screen = self._make_screen(synthetic_cohort)
        result = screen.run(expr, ["MARKER_A", "NONEXISTENT"], skip_missing=True)
        assert "NONEXISTENT" not in result.gene_panel

    def test_skip_missing_false_raises(self, synthetic_cohort):
        expr, scores, screen = self._make_screen(synthetic_cohort)
        with pytest.raises(KeyError):
            screen.run(expr, ["MARKER_A", "NONEXISTENT"], skip_missing=False)

    def test_empty_panel_raises(self, synthetic_cohort):
        expr, scores, screen = self._make_screen(synthetic_cohort)
        with pytest.raises(ValueError, match="no valid genes"):
            screen.run(expr, ["NOT_THERE_1", "NOT_THERE_2"])

    def test_rank_largest_first(self, synthetic_cohort):
        expr, _, screen = self._make_screen(synthetic_cohort)
        result = screen.run(
            expr, ["MARKER_A", "MARKER_B", "MARKER_C", "GENE_20"]
        )
        ranking = result.rank()
        assert ranking.iloc[0, 0] >= ranking.iloc[-1, 0]


# --------------------------------------------------------------------------- #
# PerturbationReport — directional signature (roadmap acceptance)
# --------------------------------------------------------------------------- #

class TestDirectionalSignature:

    def _screen_result(self, synthetic_cohort):
        expr, scores, _ = synthetic_cohort
        perturber = GeneformerPerturber(FallbackPerturber(embedding_dim=16))
        emb = perturber.embed(expr)
        head = MSVFromEmbedding().fit(emb, scores)
        return PerturbationScreen(perturber, head).run(
            expr, ["MARKER_A", "MARKER_B", "MARKER_C", "GENE_20"],
        )

    def test_master_regulator_directional_expectation(self, synthetic_cohort):
        """Knocking out MARKER_A (drives pathway_0) should DECREASE pathway_0.

        This is the roadmap-style acceptance criterion for F5:
        perturbing a known master regulator produces the expected sign
        of MSV shift.
        """
        result = self._screen_result(synthetic_cohort)
        report = PerturbationReport.from_screen(result)
        signature = {
            "MARKER_A": {"pathway_0": -1},
            "MARKER_B": {"pathway_1": -1},
            "MARKER_C": {"pathway_2": -1},
        }
        check = report.check_directional_signature(signature)
        # All three master regulators should pass (expected sign matches observed)
        assert check["passed"].all(), check

    def test_missing_gene_is_flagged(self, synthetic_cohort):
        result = self._screen_result(synthetic_cohort)
        report = PerturbationReport.from_screen(result)
        check = report.check_directional_signature({
            "NOT_IN_PANEL": {"pathway_0": -1},
        })
        assert not check["passed"].all()
        assert check["reason"].iloc[0] == "missing gene or pathway"

    def test_top_k(self, synthetic_cohort):
        result = self._screen_result(synthetic_cohort)
        report = PerturbationReport.from_screen(result)
        top2 = report.top_k(2)
        assert len(top2) == 2

    def test_dominant_pathway_per_gene(self, synthetic_cohort):
        result = self._screen_result(synthetic_cohort)
        report = PerturbationReport.from_screen(result)
        dom = report.dominant_pathway_per_gene()
        # Master regulator of cluster 0 should dominate pathway_0
        assert dom["MARKER_A"] == "pathway_0"

    def test_summary_renders(self, synthetic_cohort):
        result = self._screen_result(synthetic_cohort)
        report = PerturbationReport.from_screen(result)
        s = report.summary()
        assert "PerturbationReport" in s

    def test_to_dict_has_top5(self, synthetic_cohort):
        result = self._screen_result(synthetic_cohort)
        report = PerturbationReport.from_screen(result)
        d = report.to_dict()
        assert "top5" in d
        assert "dominant_pathway_per_gene" in d


# --------------------------------------------------------------------------- #
# F1 integration: conformal intervals on perturbed MSV
# --------------------------------------------------------------------------- #

class TestConformalIntegration:

    def test_conformal_wraps_msv_head(self, synthetic_cohort):
        """Acceptance: conformal intervals remain calibrated when wrapping
        the GeneformerPerturber -> MSVFromEmbedding chain."""
        from pathway_subtyping.uncertainty import ConformalPathwayPredictor

        expr, scores, _ = synthetic_cohort
        perturber = GeneformerPerturber(FallbackPerturber(embedding_dim=16))
        emb = perturber.embed(expr)
        head = MSVFromEmbedding().fit(emb, scores)

        # Calibrate a conformal predictor on pathway_0 using the head.
        def score_fn(X):
            # X is an embedding; return predicted pathway_0 score
            return head.transform(X)["pathway_0"].to_numpy()

        rng = np.random.default_rng(0)
        perm = rng.permutation(len(emb))
        cal_idx = perm[: len(emb) // 2]
        te_idx = perm[len(emb) // 2:]

        predictor = ConformalPathwayPredictor(score_fn=score_fn, coverage=0.9)
        predictor.calibrate(emb[cal_idx], scores["pathway_0"].to_numpy()[cal_idx])
        empirical = predictor.coverage_on(
            emb[te_idx], scores["pathway_0"].to_numpy()[te_idx]
        )
        # Should be at or above nominal (finite-sample oracle effect)
        assert empirical >= 0.85
