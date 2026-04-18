"""
Learned head mapping Geneformer embeddings to PSF MSV scores.

``MSVFromEmbedding`` is the small trainable bridge between Geneformer's
per-cell embedding and PSF's pathway-level molecular state vector. The
production head is an MLP fit during package build against paired
(embedding, MSV) data on public benchmarks. This implementation ships a
closed-form ridge-regression head — zero extra dependencies, deterministic
across reruns, and sufficient for the wrapper + screen machinery to be
tested end-to-end. The MLP variant can drop in later behind the same API.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class MSVFromEmbedding:
    """Linear head: (n_cells, embedding_dim) -> (n_cells, n_pathways).

    Usage::

        head = MSVFromEmbedding(ridge_alpha=1e-2)
        head.fit(embeddings, pathway_scores)
        predicted = head.transform(new_embeddings)
        # predicted: DataFrame shape (n_cells, n_pathways)
    """

    def __init__(self, ridge_alpha: float = 1e-2) -> None:
        self.ridge_alpha = float(ridge_alpha)
        self._beta: Optional[np.ndarray] = None  # (embedding_dim + 1, n_pathways)
        self._pathway_names: Optional[List[str]] = None
        self._fitted: bool = False

    # -------------------------------------------------------------- fit ---
    def fit(
        self,
        embeddings: np.ndarray,
        pathway_scores: pd.DataFrame,
    ) -> "MSVFromEmbedding":
        if len(embeddings) != len(pathway_scores):
            raise ValueError(
                "embeddings and pathway_scores must have the same number of rows"
            )
        if len(embeddings) < 5:
            raise ValueError("need >= 5 cells to fit the MSV head")
        X = np.hstack([np.asarray(embeddings, dtype=float),
                       np.ones((len(embeddings), 1))])
        Y = pathway_scores.to_numpy(dtype=float)
        d = X.shape[1]
        reg = self.ridge_alpha * np.eye(d)
        reg[-1, -1] = 0.0  # do not regularize bias
        self._beta = np.linalg.solve(X.T @ X + reg, X.T @ Y)
        self._pathway_names = list(pathway_scores.columns)
        self._fitted = True
        logger.info(
            "[MSVFromEmbedding] fit: n_cells=%d embedding_dim=%d n_pathways=%d",
            len(embeddings), embeddings.shape[1], Y.shape[1],
        )
        return self

    # ---------------------------------------------------------- transform ---
    def transform(self, embeddings: np.ndarray) -> pd.DataFrame:
        if not self._fitted or self._beta is None:
            raise RuntimeError("MSVFromEmbedding.fit() must be called first")
        X = np.hstack([np.asarray(embeddings, dtype=float),
                       np.ones((len(embeddings), 1))])
        Y_hat = X @ self._beta
        return pd.DataFrame(Y_hat, columns=self._pathway_names)

    def fit_transform(
        self,
        embeddings: np.ndarray,
        pathway_scores: pd.DataFrame,
    ) -> pd.DataFrame:
        self.fit(embeddings, pathway_scores)
        return self.transform(embeddings)

    # ------------------------------------------------------------ delta ---
    def delta(
        self,
        baseline_embeddings: np.ndarray,
        perturbed_embeddings: np.ndarray,
    ) -> pd.DataFrame:
        """Return the per-cell per-pathway MSV delta.

        Convenience wrapper: ``transform(perturbed) - transform(baseline)``.
        """
        if baseline_embeddings.shape != perturbed_embeddings.shape:
            raise ValueError("baseline and perturbed embeddings must share shape")
        baseline = self.transform(baseline_embeddings)
        perturbed = self.transform(perturbed_embeddings)
        return perturbed - baseline

    # --------------------------------------------------------- accessors ---
    @property
    def pathway_names(self) -> List[str]:
        if not self._fitted:
            raise RuntimeError("MSVFromEmbedding.fit() must be called first")
        assert self._pathway_names is not None
        return list(self._pathway_names)
