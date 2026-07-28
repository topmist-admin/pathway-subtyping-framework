"""
Tests for v0.6 F8 Nicheformer spatial-joint embedding layer.

Covers:
    - FallbackNicheformerEmbedder shape, determinism, and cross-modality
      projection consistency
    - NicheformerEmbedder.embed_joint: dissociated + spatial inputs
      produce embeddings in a shared basis
    - OfficialNicheformerBackend stub raises informative ImportError
    - Roadmap acceptance: dissociated-vs-spatial pathway-score rho
      exceeds 0.7 on a paired cortex-shaped synthetic reference (paired
      cells share biological identity, differ only in modality-specific
      noise — which Nicheformer's shared basis should collapse).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from pathway_subtyping.embed import (
    Embedder,
    EmbeddingResult,
    FallbackNicheformerEmbedder,
    NicheformerEmbedder,
    OfficialNicheformerBackend,
)

# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #


@pytest.fixture
def paired_cortex():
    """Paired dissociated + spatial matrices on the same cells.

    Three simulated cell types. Each cell's true pathway score comes
    from cluster identity. Dissociated vs spatial differ only in
    modality-specific Gaussian noise — shared biology should collapse
    under a biology-invariant embedder.
    """
    rng = np.random.default_rng(0)
    n_cells = 150
    n_genes = 40
    cluster_profiles = rng.standard_normal((3, n_genes)) * 1.5
    cluster_ids = rng.integers(0, 3, size=n_cells)

    # Shared biological signal
    base = cluster_profiles[cluster_ids]
    diss_noise = rng.normal(0, 0.3, size=base.shape)
    spat_noise = rng.normal(0, 0.3, size=base.shape)

    gene_cols = [f"GENE_{i}" for i in range(n_genes)]
    cell_idx = [f"cell_{i}" for i in range(n_cells)]
    dissociated = pd.DataFrame(base + diss_noise, index=cell_idx, columns=gene_cols)
    spatial = pd.DataFrame(base + spat_noise, index=cell_idx, columns=gene_cols)

    # Pathway scores derived from dominant cluster-profile dimensions
    pathway_cols = [f"PATH_{i}" for i in range(5)]
    diss_path = pd.DataFrame(
        (base + diss_noise) @ rng.standard_normal((n_genes, 5)),
        index=cell_idx,
        columns=pathway_cols,
    )
    spat_path = pd.DataFrame(
        (base + spat_noise) @ rng.standard_normal((n_genes, 5)),
        index=cell_idx,
        columns=pathway_cols,
    )
    return dissociated, spatial, diss_path, spat_path


# --------------------------------------------------------------------------- #
# Fallback embedder
# --------------------------------------------------------------------------- #


class TestFallbackNicheformer:

    def test_shape_and_api(self, paired_cortex):
        diss, spat, _, _ = paired_cortex
        emb = FallbackNicheformerEmbedder(embedding_dim=16).embed(diss)
        assert isinstance(emb, EmbeddingResult)
        assert emb.embeddings.shape == (len(diss), 16)
        assert emb.backend.startswith("nicheformer:fallback")

    def test_deterministic(self, paired_cortex):
        diss, _, _, _ = paired_cortex
        a = FallbackNicheformerEmbedder(embedding_dim=16, seed=0).embed(diss)
        b = FallbackNicheformerEmbedder(embedding_dim=16, seed=0).embed(diss)
        np.testing.assert_allclose(a.embeddings, b.embeddings)

    def test_shared_basis_on_matching_genes(self, paired_cortex):
        """Fit on dissociated, project spatial through same basis."""
        diss, spat, _, _ = paired_cortex
        embedder = FallbackNicheformerEmbedder(embedding_dim=16).fit(diss)
        diss_emb = embedder.embed(diss)
        spat_emb = embedder.embed(spat)
        # Paired cells should be near each other in the embedding space
        per_cell_delta = np.linalg.norm(diss_emb.embeddings - spat_emb.embeddings, axis=1)
        # Compare against the magnitude of the dissociated embedding
        norm_diss = np.linalg.norm(diss_emb.embeddings, axis=1)
        # Inter-modality distance should be a fraction of within-modality magnitude
        assert (per_cell_delta < 1.2 * norm_diss).mean() > 0.8

    def test_rejects_bad_dim(self):
        with pytest.raises(ValueError):
            FallbackNicheformerEmbedder(embedding_dim=1)

    def test_rejects_non_dataframe(self):
        with pytest.raises(TypeError):
            FallbackNicheformerEmbedder().fit(np.zeros((5, 5)))  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# Official backend stub
# --------------------------------------------------------------------------- #


class TestOfficialStub:

    def test_informative_error(self):
        be = OfficialNicheformerBackend()
        try:
            be._lazy_load()
        except (ImportError, NotImplementedError) as exc:
            msg = str(exc).lower()
            assert "embed" in msg or "nicheformer" in msg or "checkpoint" in msg
        else:
            pytest.skip("Nicheformer installed — skipping stub error test")


# --------------------------------------------------------------------------- #
# NicheformerEmbedder (high-level)
# --------------------------------------------------------------------------- #


class TestNicheformerEmbedder:

    def test_default_backend_is_fallback(self):
        e = NicheformerEmbedder()
        assert isinstance(e.backend, FallbackNicheformerEmbedder)

    def test_is_embedder_subclass(self):
        assert issubclass(NicheformerEmbedder, Embedder)

    def test_embed_joint_requires_matching_columns(self, paired_cortex):
        diss, spat, _, _ = paired_cortex
        embedder = NicheformerEmbedder()
        # Drop a column from spat -> columns no longer match
        with pytest.raises(ValueError, match="same gene columns"):
            embedder.embed_joint(diss, spat.iloc[:, :-1])

    def test_embed_joint_shapes(self, paired_cortex):
        diss, spat, _, _ = paired_cortex
        embedder = NicheformerEmbedder(FallbackNicheformerEmbedder(embedding_dim=16))
        d, s = embedder.embed_joint(diss, spat)
        assert d.embeddings.shape == (len(diss), 16)
        assert s.embeddings.shape == (len(spat), 16)


# --------------------------------------------------------------------------- #
# Roadmap acceptance: dissociated vs spatial pathway-score rho > 0.7
# --------------------------------------------------------------------------- #


class TestJointSpaceAcceptance:

    def test_dissociated_spatial_rho_exceeds_threshold(self, paired_cortex):
        """With paired cells embedded in a shared basis, dissociated and
        spatial pathway scores (predicted from embeddings via a small
        linear head) should correlate strongly per pathway.

        This is the roadmap acceptance proxy: 'dissociated-vs-spatial
        pathway score rho > 0.7 on paired cortex reference'.
        """
        from pathway_subtyping.perturb import MSVFromEmbedding

        diss, spat, diss_path, spat_path = paired_cortex
        embedder = NicheformerEmbedder(FallbackNicheformerEmbedder(embedding_dim=16))
        d_emb, s_emb = embedder.embed_joint(diss, spat)

        # Train a shared head on the dissociated modality; evaluate on both.
        head = MSVFromEmbedding(ridge_alpha=1e-2).fit(
            d_emb.embeddings,
            diss_path,
        )
        diss_pred = head.transform(d_emb.embeddings)
        spat_pred = head.transform(s_emb.embeddings)

        # Per-pathway Spearman rho between the two modalities' predictions
        rhos = []
        for col in diss_pred.columns:
            rho, _ = spearmanr(diss_pred[col].to_numpy(), spat_pred[col].to_numpy())
            if not np.isnan(rho):
                rhos.append(float(rho))
        mean_rho = float(np.mean(rhos))
        assert mean_rho > 0.7, f"mean dissociated/spatial rho = {mean_rho:.3f} below 0.7"

    def test_integrates_with_cross_platform_aligner(self, paired_cortex):
        """Acceptance hook: dissociated/spatial joint alignment works via
        F2's CrossPlatformAligner with modality labels."""
        from pathway_subtyping.harmonize import CrossPlatformAligner

        diss, spat, diss_path, spat_path = paired_cortex
        embedder = NicheformerEmbedder(FallbackNicheformerEmbedder(embedding_dim=16))
        d_emb, s_emb = embedder.embed_joint(diss, spat)

        # Stack as if they were two "platforms"
        all_scores = pd.concat([diss_path, spat_path], axis=0).reset_index(drop=True)
        platforms = (["dissociated"] * len(diss_path)) + (["spatial"] * len(spat_path))
        embeddings = np.vstack([d_emb.embeddings, s_emb.embeddings])

        aligner = CrossPlatformAligner()
        aligned = aligner.fit_transform(all_scores, platforms, embeddings)
        assert aligned.aligned_scores.shape == all_scores.shape
        assert set(aligner.platforms) == {"dissociated", "spatial"}
