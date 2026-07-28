"""
Nicheformer embedding wrapper for joint dissociated + spatial analysis.

``NicheformerEmbedder`` follows the same pattern as the F6 scGPT
wrapper: a pluggable backend with an opt-in ``OfficialNicheformerBackend``
(requires ``pathway-subtyping[embed]``) and a deterministic
``FallbackNicheformerEmbedder`` for tests and CI.

F8 key capability (v0.6 roadmap): users who have both 10x dissociated
scRNA-seq and Visium spatial of the same tissue can embed both modalities
into a common biology-invariant space and then use the F2
``CrossPlatformAligner`` (with modality labels in place of platform
labels) to produce joint-space pathway scores.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from typing import Any, List, Optional

import numpy as np
import pandas as pd

from .base import Embedder, EmbeddingResult

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Official backend (opt-in)
# --------------------------------------------------------------------------- #


class OfficialNicheformerBackend(Embedder):
    """Production Nicheformer wrapper — requires ``pathway-subtyping[embed]``.

    Usage::

        backend = OfficialNicheformerBackend(checkpoint="nicheformer-v1.0")
        result = backend.embed(dissociated_df)   # or spatial_df

    The same backend is used for both modalities; Nicheformer is trained
    to produce a modality-invariant embedding so the aligner has a clean
    biology target.
    """

    embedding_dim = 512

    def __init__(
        self,
        checkpoint: str = "nicheformer-v1.0",
        device: str = "cpu",
    ) -> None:
        self.checkpoint = checkpoint
        self.device = device
        self._model: Any = None

    def _lazy_load(self) -> None:
        if self._model is not None:
            return
        try:
            import nicheformer  # type: ignore  # noqa: F401
            import torch  # noqa: F401
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "OfficialNicheformerBackend requires torch + nicheformer. "
                "Install with: pip install pathway-subtyping[embed]"
            ) from exc
        raise NotImplementedError(
            "Live Nicheformer inference is gated on the packaged checkpoint "
            "loader expected from the '[embed]' extra. Use "
            "FallbackNicheformerEmbedder for local testing."
        )

    @property
    def backend_id(self) -> str:
        return f"nicheformer:official/{self.checkpoint}"

    def embed(self, expression: pd.DataFrame) -> EmbeddingResult:  # pragma: no cover
        self._lazy_load()
        arr = np.zeros((len(expression), self.embedding_dim), dtype=float)
        return EmbeddingResult(
            embeddings=arr,
            cell_index=expression.index,
            backend=self.backend_id,
            meta={"checkpoint": self.checkpoint},
        )


# --------------------------------------------------------------------------- #
# Deterministic fallback
# --------------------------------------------------------------------------- #


class FallbackNicheformerEmbedder(Embedder):
    """Deterministic PCA substitute for Nicheformer.

    Treats both dissociated and spatial inputs as expression matrices
    (same gene columns) and projects through a shared PCA basis fitted
    once on the first modality seen. The shared basis is what makes
    dissociated/spatial scores comparable.
    """

    def __init__(
        self,
        embedding_dim: int = 64,
        seed: Optional[int] = 0,
    ) -> None:
        if embedding_dim < 2:
            raise ValueError("embedding_dim must be >= 2")
        self.embedding_dim = int(embedding_dim)
        self.seed = seed
        self._mean: Optional[np.ndarray] = None
        self._components: Optional[np.ndarray] = None
        self._columns: Optional[List[str]] = None
        self._fitted: bool = False

    @property
    def backend_id(self) -> str:
        return f"nicheformer:fallback/dim={self.embedding_dim}"

    # --------------------------------------------------------------- fit ---
    def fit(self, expression: pd.DataFrame) -> "FallbackNicheformerEmbedder":
        if not isinstance(expression, pd.DataFrame):
            raise TypeError("FallbackNicheformerEmbedder requires a DataFrame input")
        self._columns = list(expression.columns)
        X = expression.to_numpy(dtype=float)
        self._mean = X.mean(axis=0)
        centered = X - self._mean
        k = min(self.embedding_dim, centered.shape[1], centered.shape[0])
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        self._components = Vt[:k]
        if k < self.embedding_dim:
            rng = np.random.default_rng(self.seed)
            extra = rng.standard_normal((self.embedding_dim - k, centered.shape[1]))
            extra = np.linalg.qr(extra.T)[0].T
            self._components = np.vstack([self._components, extra])[: self.embedding_dim]
        self._fitted = True
        logger.info(
            "[FallbackNicheformerEmbedder] fit: dim=%d n_genes=%d",
            self.embedding_dim,
            centered.shape[1],
        )
        return self

    # ------------------------------------------------------------- embed ---
    def embed(self, expression: pd.DataFrame) -> EmbeddingResult:
        if not self._fitted:
            self.fit(expression)
        assert self._mean is not None and self._components is not None
        X = expression.to_numpy(dtype=float)
        if X.shape[1] != self._components.shape[1]:
            raise ValueError(
                f"expression has {X.shape[1]} genes; embedder was fitted on "
                f"{self._components.shape[1]}"
            )
        arr = (X - self._mean) @ self._components.T
        return EmbeddingResult(
            embeddings=arr,
            cell_index=expression.index,
            backend=self.backend_id,
            meta={"seed": self.seed},
        )


# --------------------------------------------------------------------------- #
# High-level wrapper
# --------------------------------------------------------------------------- #


class NicheformerEmbedder(Embedder):
    """Public Nicheformer interface — delegates to a backend.

    Usage (same backend used for both modalities)::

        embedder = NicheformerEmbedder()           # fallback by default
        diss_emb = embedder.embed(dissociated_df)
        spat_emb = embedder.embed(spatial_df)
        # diss_emb.embeddings and spat_emb.embeddings share the same basis
    """

    def __init__(self, backend: Optional[Embedder] = None) -> None:
        self.backend = backend if backend is not None else FallbackNicheformerEmbedder()

    @property
    def embedding_dim(self) -> int:  # type: ignore[override]
        return int(self.backend.embedding_dim)

    def embed(self, expression: pd.DataFrame) -> EmbeddingResult:
        return self.backend.embed(expression)

    # -------------------------------------------------------- joint frame ---
    def embed_joint(
        self,
        dissociated: pd.DataFrame,
        spatial: pd.DataFrame,
    ) -> tuple[EmbeddingResult, EmbeddingResult]:
        """Embed dissociated + spatial matrices in a shared basis.

        The backend is fitted/used on the dissociated matrix first so the
        shared basis is driven by the higher-coverage modality; the
        spatial matrix is then projected through the same basis.

        Gene columns must match between the two frames. Cells remain
        separate — use :class:`pathway_subtyping.harmonize.CrossPlatformAligner`
        with modality labels to produce a joint-frame score matrix.
        """
        if list(dissociated.columns) != list(spatial.columns):
            raise ValueError("dissociated and spatial inputs must share the same gene columns")
        diss = self.embed(dissociated)
        spat = self.embed(spatial)
        return diss, spat
