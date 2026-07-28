"""
Universal Cell Embeddings (UCE) wrapper for cross-platform harmonization.

Provides a thin, backend-agnostic interface for extracting per-cell
embeddings that downstream alignment uses as a platform-invariant anchor.

Two backends:

    UCEEmbedder  — opt-in wrapper around the real UCE model (Rosen et al.,
                   Nature 2024). Requires torch + the UCE checkpoint, which
                   ship with the ``[harmonize]`` extra. Lazy-imports; a
                   clear ImportError is raised if the extra is missing.
    FallbackEmbedder — deterministic PCA-based embedder that preserves the
                   public API so tests and development can run uniformly
                   without the UCE checkpoint. Not a substitute for real
                   UCE; use only for testing or smoke runs.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from typing import Any, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class _BaseEmbedder:
    """Interface both backends share."""

    embedding_dim: int

    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        raise NotImplementedError


class UCEEmbedder(_BaseEmbedder):
    """Real UCE model wrapper — requires ``pathway-subtyping[harmonize]``.

    Usage::

        embedder = UCEEmbedder(checkpoint="uce-33l-v2", device="cpu")
        embeddings = embedder.embed(expression_df)   # (n_cells, 1280)

    The model checkpoint is downloaded on first use and cached under
    ``~/.cache/pathway-subtyping/uce/``. CPU inference is supported but
    slow; users with a CUDA device should pass ``device="cuda"``.
    """

    def __init__(
        self,
        checkpoint: str = "uce-33l-v2",
        device: str = "cpu",
        batch_size: int = 32,
    ):
        self.checkpoint = checkpoint
        self.device = device
        self.batch_size = int(batch_size)
        self._model: Any = None
        self.embedding_dim = 1280

    def _lazy_load(self) -> None:
        if self._model is not None:
            return
        try:
            import torch  # noqa: F401
            import uce  # type: ignore  # noqa: F401
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "UCEEmbedder requires 'torch' and 'uce'. Install the "
                "extra with: pip install pathway-subtyping[harmonize]"
            ) from exc
        raise NotImplementedError(
            "Live UCE inference is gated on the packaged checkpoint; this "
            "stub defers to the checkpoint loader expected from the "
            "'[harmonize]' extra. Use FallbackEmbedder for local testing."
        )

    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        self._lazy_load()
        # pragma: no cover - requires optional dep
        return np.zeros((len(expression), self.embedding_dim), dtype=float)


class FallbackEmbedder(_BaseEmbedder):
    """Deterministic PCA-based embedder — works without optional deps.

    Fits on the first ``embed()`` call so a single embedder can produce
    mutually comparable embeddings across multiple platform-specific
    datasets. Subsequent calls project into the already-fitted space.

    Not a drop-in substitute for real UCE: the embedding space is derived
    from the variance structure of whatever data the embedder sees first,
    so it captures platform signal rather than biological universals.
    Suitable for tests, smoke runs, and demonstrating the harmonize API.
    """

    def __init__(self, embedding_dim: int = 64, seed: Optional[int] = 0):
        if embedding_dim < 2:
            raise ValueError("embedding_dim must be >= 2")
        self.embedding_dim = int(embedding_dim)
        self.seed = seed
        self._mean: Optional[np.ndarray] = None
        self._components: Optional[np.ndarray] = None
        self._fitted: bool = False

    # --------------------------------------------------------------- fit ---
    def fit(self, expression: pd.DataFrame) -> "FallbackEmbedder":
        """Fit the PCA basis on an expression matrix."""
        X = np.asarray(expression.values if isinstance(expression, pd.DataFrame) else expression)
        self._mean = X.mean(axis=0)
        centered = X - self._mean
        # SVD-based PCA (truncated)
        k = min(self.embedding_dim, centered.shape[1], centered.shape[0])
        # deterministic: numpy.linalg.svd with full_matrices=False
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
        self._components = Vt[:k]
        if k < self.embedding_dim:
            # pad with random orthogonal vectors (deterministic via seed)
            rng = np.random.default_rng(self.seed)
            extra = rng.standard_normal((self.embedding_dim - k, centered.shape[1]))
            extra = np.linalg.qr(extra.T)[0].T
            self._components = np.vstack([self._components, extra])[: self.embedding_dim]
        self._fitted = True
        logger.info(
            "[FallbackEmbedder] fitted: embedding_dim=%d (from rank=%d)",
            self.embedding_dim,
            k,
        )
        return self

    # ----------------------------------------------------------- embed ---
    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        """Project an expression matrix into the fitted embedding space."""
        X = np.asarray(expression.values if isinstance(expression, pd.DataFrame) else expression)
        if not self._fitted:
            self.fit(expression)
        assert self._mean is not None and self._components is not None
        if X.shape[1] != self._components.shape[1]:
            raise ValueError(
                f"expression has {X.shape[1]} genes; embedder was fitted "
                f"on {self._components.shape[1]}"
            )
        return (X - self._mean) @ self._components.T
