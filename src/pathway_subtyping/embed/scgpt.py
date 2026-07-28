"""
scGPT embedding wrapper.

``scGPTEmbedder`` is the stable public API. Like the UCE and Geneformer
wrappers before it, the embedder delegates inference to a pluggable
backend:

    - ``OfficialSCGPTBackend`` — wraps the published scGPT model
      (Cui et al. 2024). Requires ``pip install pathway-subtyping[embed]``
      which pulls in torch + transformers + the scgpt checkpoint.
    - ``FallbackSCGPTEmbedder`` — deterministic PCA-based substitute
      that implements the same API for tests and CI.

Downstream modules (``harmonize/align.py``, ``perturb/msv_from_embedding.py``)
consume :class:`EmbeddingResult` uniformly and do not depend on which
backend produced the embedding.

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


class OfficialSCGPTBackend(Embedder):
    """Production scGPT wrapper — requires ``pathway-subtyping[embed]``.

    Usage::

        embedder = OfficialSCGPTBackend(checkpoint="scgpt-whole-human-v1")
        result = embedder.embed(expression_df)
    """

    embedding_dim = 512

    def __init__(
        self,
        checkpoint: str = "scgpt-whole-human-v1",
        device: str = "cpu",
    ) -> None:
        self.checkpoint = checkpoint
        self.device = device
        self._model: Any = None

    def _lazy_load(self) -> None:
        if self._model is not None:
            return
        try:
            import scgpt  # type: ignore  # noqa: F401
            import torch  # noqa: F401
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "OfficialSCGPTBackend requires torch + scgpt. "
                "Install with: pip install pathway-subtyping[embed]"
            ) from exc
        raise NotImplementedError(
            "Live scGPT inference is gated on the packaged checkpoint "
            "loader expected from the '[embed]' extra. Use "
            "FallbackSCGPTEmbedder for local testing or CI that does "
            "not ship the checkpoint."
        )

    @property
    def backend_id(self) -> str:
        return f"scgpt:official/{self.checkpoint}"

    def embed(self, expression: pd.DataFrame) -> EmbeddingResult:  # pragma: no cover
        self._lazy_load()
        arr = np.zeros((len(expression), self.embedding_dim), dtype=float)
        return EmbeddingResult(
            embeddings=arr,
            cell_index=expression.index,
            backend=self.backend_id,
            meta={"checkpoint": self.checkpoint, "device": self.device},
        )


# --------------------------------------------------------------------------- #
# Deterministic fallback
# --------------------------------------------------------------------------- #


class FallbackSCGPTEmbedder(Embedder):
    """Deterministic PCA substitute for scGPT.

    Fits a PCA basis on the first expression matrix it sees and projects
    subsequent matrices through it. Preserves the :class:`Embedder` API,
    so tests and CI can run without any optional dependencies.
    """

    def __init__(
        self,
        embedding_dim: int = 128,
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
        return f"scgpt:fallback/dim={self.embedding_dim}"

    # --------------------------------------------------------------- fit ---
    def fit(self, expression: pd.DataFrame) -> "FallbackSCGPTEmbedder":
        if not isinstance(expression, pd.DataFrame):
            raise TypeError("FallbackSCGPTEmbedder requires a DataFrame input")
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
            "[FallbackSCGPTEmbedder] fit: embedding_dim=%d n_genes=%d",
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


class scGPTEmbedder(Embedder):
    """Public scGPT interface — delegates to an :class:`Embedder` backend.

    Usage::

        embedder = scGPTEmbedder()                 # fallback by default
        result = embedder.embed(expression_df)
        # result.embeddings is (n_cells, embedding_dim)
    """

    def __init__(self, backend: Optional[Embedder] = None) -> None:
        self.backend = backend if backend is not None else FallbackSCGPTEmbedder()

    @property
    def embedding_dim(self) -> int:  # type: ignore[override]
        return int(self.backend.embedding_dim)

    def embed(self, expression: pd.DataFrame) -> EmbeddingResult:
        return self.backend.embed(expression)
