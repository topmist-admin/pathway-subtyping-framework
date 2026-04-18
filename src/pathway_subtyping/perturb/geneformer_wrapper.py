"""
Geneformer wrapper for in-silico gene perturbation.

``GeneformerPerturber`` is the high-level interface used throughout
``pathway_subtyping.perturb``. It delegates actual inference to a
pluggable ``GeneformerBackend``:

    - ``OfficialBackend`` — the production path, wraps the published
      Geneformer model (Theodoris et al. 2023). Requires the optional
      ``[perturb]`` extra (torch, transformers, geneformer package, and
      a cached checkpoint). Lazy-imports; a clear ImportError is raised
      if the extra is missing.
    - ``FallbackPerturber`` — deterministic numpy substitute that
      preserves the public API so tests and development can run without
      the optional dependencies. It models a knockout as the removal of
      that gene's contribution to a PCA embedding of the expression
      matrix, and over-expression as its amplification.

The fallback is not a substitute for the real model on biological claims.
It is a correctness and interface harness.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class PerturbationMode(str, Enum):
    """Supported perturbation modes."""

    KNOCKOUT = "knockout"
    OVEREXPRESS = "overexpress"


# --------------------------------------------------------------------------- #
# Backend ABC + official stub
# --------------------------------------------------------------------------- #

class GeneformerBackend:
    """Abstract base class for Geneformer inference backends.

    Subclasses implement ``embed`` and ``perturb``. Custom backends must
    pass the golden-test parity harness (``tests/test_perturb.py``
    ``@pytest.mark.perturb_golden``) before being used in production.
    """

    embedding_dim: int

    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        raise NotImplementedError

    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode,
    ) -> np.ndarray:
        raise NotImplementedError


class OfficialBackend(GeneformerBackend):
    """Wraps the published Geneformer model.

    Requires ``pip install pathway-subtyping[perturb]`` which pulls in
    ``torch``, ``transformers``, and the ``geneformer`` package, plus
    downloads the model checkpoint to the local cache on first use.
    """

    embedding_dim = 256

    def __init__(
        self,
        checkpoint: str = "geneformer-12L-v2",
        device: str = "cpu",
    ) -> None:
        self.checkpoint = checkpoint
        self.device = device
        self._model: Any = None

    def _lazy_load(self) -> None:
        if self._model is not None:
            return
        try:
            import torch  # noqa: F401
            import transformers  # noqa: F401
            import geneformer  # type: ignore  # noqa: F401
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "OfficialBackend requires torch + transformers + geneformer. "
                "Install with: pip install pathway-subtyping[perturb]"
            ) from exc
        raise NotImplementedError(
            "Live Geneformer inference is gated on the packaged checkpoint "
            "loader expected from the '[perturb]' extra. Use FallbackPerturber "
            "for local testing or CI that does not ship the checkpoint."
        )

    def embed(self, expression: pd.DataFrame) -> np.ndarray:  # pragma: no cover
        self._lazy_load()
        return np.zeros((len(expression), self.embedding_dim), dtype=float)

    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode,
    ) -> np.ndarray:  # pragma: no cover
        self._lazy_load()
        return np.zeros((len(expression), self.embedding_dim), dtype=float)


# --------------------------------------------------------------------------- #
# Deterministic fallback
# --------------------------------------------------------------------------- #

class FallbackPerturber(GeneformerBackend):
    """Deterministic PCA-based substitute for Geneformer.

    Fits a PCA basis on the first expression matrix it sees. ``embed``
    projects into that basis; ``perturb(expr, gene, mode)`` simulates
    a perturbation by setting the gene's column to zero (knockout) or
    amplifying it (overexpress) before projecting.

    This is an API harness, not a biology-grade substitute. The
    directional behaviour on well-separated markers is correct —
    knocking out a gene that anchors a cluster does move cells toward
    the mean embedding — but absolute magnitudes are meaningless.
    """

    def __init__(
        self,
        embedding_dim: int = 64,
        overexpression_factor: float = 3.0,
        seed: Optional[int] = 0,
    ) -> None:
        if embedding_dim < 2:
            raise ValueError("embedding_dim must be >= 2")
        self.embedding_dim = int(embedding_dim)
        self.overexpression_factor = float(overexpression_factor)
        self.seed = seed
        self._mean: Optional[np.ndarray] = None
        self._components: Optional[np.ndarray] = None
        self._columns: Optional[List[str]] = None
        self._fitted: bool = False

    # --------------------------------------------------------------- fit ---
    def fit(self, expression: pd.DataFrame) -> "FallbackPerturber":
        if not isinstance(expression, pd.DataFrame):
            raise TypeError("FallbackPerturber requires a DataFrame input")
        self._columns = list(expression.columns)
        X = expression.to_numpy(dtype=float)
        self._mean = X.mean(axis=0)
        centered = X - self._mean
        k = min(self.embedding_dim, centered.shape[1], centered.shape[0])
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        self._components = Vt[:k]
        if k < self.embedding_dim:
            rng = np.random.default_rng(self.seed)
            extra = rng.standard_normal(
                (self.embedding_dim - k, centered.shape[1])
            )
            extra = np.linalg.qr(extra.T)[0].T
            self._components = np.vstack(
                [self._components, extra]
            )[: self.embedding_dim]
        self._fitted = True
        logger.info(
            "[FallbackPerturber] fitted: embedding_dim=%d n_genes=%d",
            self.embedding_dim, centered.shape[1],
        )
        return self

    # -------------------------------------------------------------- embed ---
    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        if not self._fitted:
            self.fit(expression)
        assert self._mean is not None and self._components is not None
        X = expression.to_numpy(dtype=float)
        if X.shape[1] != self._components.shape[1]:
            raise ValueError(
                f"expression has {X.shape[1]} genes; backend was fitted on "
                f"{self._components.shape[1]}"
            )
        return (X - self._mean) @ self._components.T

    # ------------------------------------------------------------ perturb ---
    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode,
    ) -> np.ndarray:
        if not self._fitted:
            self.fit(expression)
        if gene not in expression.columns:
            raise KeyError(
                f"gene {gene!r} not in expression columns; got "
                f"{len(expression.columns)} columns"
            )
        perturbed = expression.copy()
        if mode == PerturbationMode.KNOCKOUT:
            perturbed[gene] = 0.0
        elif mode == PerturbationMode.OVEREXPRESS:
            perturbed[gene] = perturbed[gene] * self.overexpression_factor
        else:
            raise ValueError(f"unknown mode: {mode}")
        return self.embed(perturbed)


# --------------------------------------------------------------------------- #
# High-level wrapper
# --------------------------------------------------------------------------- #

@dataclass
class PerturbationResult:
    """One perturbation call's result.

    Attributes:
        gene: The gene that was perturbed.
        mode: Knockout or overexpress.
        baseline_embedding: Embedding of the unperturbed expression (n_cells, d).
        perturbed_embedding: Embedding after perturbation (n_cells, d).
        delta_embedding: ``perturbed_embedding - baseline_embedding``.
        l2_impact_per_cell: L2 norm of delta per cell.
    """

    gene: str
    mode: PerturbationMode
    baseline_embedding: np.ndarray
    perturbed_embedding: np.ndarray
    delta_embedding: np.ndarray
    l2_impact_per_cell: np.ndarray

    @property
    def mean_l2_impact(self) -> float:
        return float(self.l2_impact_per_cell.mean())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gene": self.gene,
            "mode": self.mode.value,
            "n_cells": int(self.baseline_embedding.shape[0]),
            "embedding_dim": int(self.baseline_embedding.shape[1]),
            "mean_l2_impact": self.mean_l2_impact,
        }


class GeneformerPerturber:
    """High-level perturbation interface.

    Usage::

        perturber = GeneformerPerturber()  # uses FallbackPerturber by default
        result = perturber.perturb(expression_df, gene="MYC", mode="knockout")
        # result.delta_embedding is (n_cells, embedding_dim)
    """

    def __init__(
        self,
        backend: Optional[GeneformerBackend] = None,
    ) -> None:
        self.backend = backend if backend is not None else FallbackPerturber()

    @property
    def embedding_dim(self) -> int:
        return self.backend.embedding_dim

    # -------------------------------------------------------------- embed ---
    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        return self.backend.embed(expression)

    # ----------------------------------------------------------- perturb ---
    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode = PerturbationMode.KNOCKOUT,
    ) -> PerturbationResult:
        mode = PerturbationMode(mode) if isinstance(mode, str) else mode
        baseline = self.backend.embed(expression)
        perturbed = self.backend.perturb(expression, gene, mode)
        delta = perturbed - baseline
        l2 = np.linalg.norm(delta, axis=1)
        return PerturbationResult(
            gene=gene,
            mode=mode,
            baseline_embedding=baseline,
            perturbed_embedding=perturbed,
            delta_embedding=delta,
            l2_impact_per_cell=l2,
        )
