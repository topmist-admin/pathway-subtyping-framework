"""
Tests for the v0.6 F6 scGPT embedding layer.

Covers:
    - FallbackSCGPTEmbedder: deterministic output, shape checks, API conformance
    - OfficialSCGPTBackend: informative ImportError / NotImplementedError
      when optional deps are missing
    - scGPTEmbedder: high-level delegation
    - EmbeddingCache: content-hashed put/get/invalidation
    - cache_key_for: stable across reruns, invalidated by content change
    - Harmonize integration: scGPT embeddings substitute for UCE in
      CrossPlatformAligner without glue code
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.embed import (
    Embedder,
    EmbeddingCache,
    EmbeddingResult,
    FallbackSCGPTEmbedder,
    OfficialSCGPTBackend,
    cache_key_for,
    scGPTEmbedder,
)

# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #


@pytest.fixture
def expression():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 120, 40
    df = pd.DataFrame(
        rng.standard_normal((n_cells, n_genes)),
        columns=[f"GENE_{i}" for i in range(n_genes)],
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    df.index.name = "cell_id"
    return df


# --------------------------------------------------------------------------- #
# FallbackSCGPTEmbedder
# --------------------------------------------------------------------------- #


class TestFallback:

    def test_shape_and_api_surface(self, expression):
        emb = FallbackSCGPTEmbedder(embedding_dim=16).embed(expression)
        assert isinstance(emb, EmbeddingResult)
        assert emb.embeddings.shape == (len(expression), 16)
        assert emb.n_cells == len(expression)
        assert emb.embedding_dim == 16
        assert emb.cell_index.equals(expression.index)

    def test_backend_identifier_records_dim(self, expression):
        emb = FallbackSCGPTEmbedder(embedding_dim=16).embed(expression)
        assert emb.backend == "scgpt:fallback/dim=16"

    def test_deterministic_across_runs(self, expression):
        a = FallbackSCGPTEmbedder(embedding_dim=16, seed=0).embed(expression)
        b = FallbackSCGPTEmbedder(embedding_dim=16, seed=0).embed(expression)
        np.testing.assert_allclose(a.embeddings, b.embeddings)

    def test_rejects_non_dataframe(self):
        emb = FallbackSCGPTEmbedder(embedding_dim=8)
        with pytest.raises(TypeError):
            emb.fit(np.zeros((5, 5)))  # type: ignore[arg-type]

    def test_rejects_bad_dim(self):
        with pytest.raises(ValueError):
            FallbackSCGPTEmbedder(embedding_dim=0)

    def test_as_dataframe_alignment(self, expression):
        emb = FallbackSCGPTEmbedder(embedding_dim=8).embed(expression)
        df = emb.as_dataframe()
        assert df.shape == (len(expression), 8)
        assert df.index.equals(expression.index)


class TestOfficialBackendStub:

    def test_lazy_import_error_informative(self):
        be = OfficialSCGPTBackend()
        try:
            be._lazy_load()
        except (ImportError, NotImplementedError) as exc:
            msg = str(exc).lower()
            assert "embed" in msg or "scgpt" in msg or "checkpoint" in msg
        else:
            pytest.skip("scGPT installed — skipping stub error test")


# --------------------------------------------------------------------------- #
# scGPTEmbedder (high-level)
# --------------------------------------------------------------------------- #


class TestSCGPTEmbedder:

    def test_default_is_fallback(self):
        e = scGPTEmbedder()
        assert isinstance(e.backend, FallbackSCGPTEmbedder)

    def test_delegates_embed_call(self, expression):
        e = scGPTEmbedder(FallbackSCGPTEmbedder(embedding_dim=12))
        out = e.embed(expression)
        assert out.embedding_dim == 12

    def test_callable_protocol(self, expression):
        e = scGPTEmbedder(FallbackSCGPTEmbedder(embedding_dim=12))
        out = e(expression)
        assert isinstance(out, EmbeddingResult)

    def test_embedder_subclass(self):
        assert issubclass(scGPTEmbedder, Embedder)


# --------------------------------------------------------------------------- #
# EmbeddingCache
# --------------------------------------------------------------------------- #


class TestEmbeddingCache:

    def test_put_get_roundtrip(self, expression, tmp_path):
        cache = EmbeddingCache(tmp_path)
        emb = FallbackSCGPTEmbedder(embedding_dim=16).embed(expression)
        key = cache_key_for(backend=emb.backend, expression=expression)
        assert not cache.has(key)
        cache.put(key, emb)
        assert cache.has(key)
        reloaded = cache.get(key)
        np.testing.assert_allclose(reloaded.embeddings, emb.embeddings)
        assert reloaded.backend == emb.backend
        assert reloaded.cell_index.equals(emb.cell_index)

    def test_miss_raises(self, tmp_path):
        cache = EmbeddingCache(tmp_path)
        with pytest.raises(KeyError):
            cache.get("nope")

    def test_delete(self, expression, tmp_path):
        cache = EmbeddingCache(tmp_path)
        emb = FallbackSCGPTEmbedder(embedding_dim=8).embed(expression)
        key = cache_key_for(backend=emb.backend, expression=expression)
        cache.put(key, emb)
        assert cache.delete(key) is True
        assert not cache.has(key)
        assert cache.delete(key) is False

    def test_clear(self, expression, tmp_path):
        cache = EmbeddingCache(tmp_path)
        for dim in (4, 6, 8):
            emb = FallbackSCGPTEmbedder(embedding_dim=dim).embed(expression)
            cache.put(
                cache_key_for(backend=emb.backend, expression=expression),
                emb,
            )
        n_cleared = cache.clear()
        assert n_cleared >= 3  # one .npy + one .json per entry

    def test_key_invalidates_on_content_change(self, expression):
        k1 = cache_key_for(backend="scgpt:fallback/dim=16", expression=expression)
        shuffled = expression.sample(frac=1, random_state=1)
        k2 = cache_key_for(backend="scgpt:fallback/dim=16", expression=shuffled)
        assert k1 != k2

    def test_key_invalidates_on_backend_change(self, expression):
        k1 = cache_key_for(backend="scgpt:fallback/dim=16", expression=expression)
        k2 = cache_key_for(backend="scgpt:official/v1.2.0", expression=expression)
        assert k1 != k2

    def test_key_stable_across_reruns(self, expression):
        k1 = cache_key_for(backend="scgpt:fallback/dim=16", expression=expression)
        k2 = cache_key_for(backend="scgpt:fallback/dim=16", expression=expression)
        assert k1 == k2

    def test_key_includes_extra(self, expression):
        k1 = cache_key_for(
            backend="scgpt:fallback/dim=16",
            expression=expression,
            extra={"tokenizer": "v1"},
        )
        k2 = cache_key_for(
            backend="scgpt:fallback/dim=16",
            expression=expression,
            extra={"tokenizer": "v2"},
        )
        assert k1 != k2


# --------------------------------------------------------------------------- #
# Integration with F2 harmonize (scGPT as drop-in UCE substitute)
# --------------------------------------------------------------------------- #


class TestHarmonizeIntegration:

    def test_scgpt_embeddings_work_with_cross_platform_aligner(self, expression):
        """Acceptance: F2 aligner consumes scGPT embeddings without glue code."""
        from pathway_subtyping.harmonize import CrossPlatformAligner

        half = len(expression) // 2
        platforms = ["10x"] * half + ["smartseq2"] * (len(expression) - half)
        # Build pathway scores from reference for the aligner to work on.
        rng = np.random.default_rng(0)
        pathway_scores = pd.DataFrame(
            rng.standard_normal((len(expression), 5)),
            index=expression.index,
            columns=[f"PATH_{i}" for i in range(5)],
        )
        embeddings = FallbackSCGPTEmbedder(embedding_dim=16).embed(expression).embeddings
        aligner = CrossPlatformAligner()
        result = aligner.fit_transform(pathway_scores, platforms, embeddings)
        assert result.aligned_scores.shape == pathway_scores.shape
