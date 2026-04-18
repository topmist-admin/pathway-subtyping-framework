"""
Abstract Embedder interface shared by all foundation-model wrappers.

Each concrete backend (scGPT, UCE already in harmonize, future Geneformer
direct embeddings, etc.) implements ``.embed(expression) -> EmbeddingResult``.
The uniform interface lets PSF's harmonize / screen / uncertainty layers
swap models without bespoke glue code.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd


@dataclass
class EmbeddingResult:
    """Output of any :class:`Embedder` backend.

    Attributes:
        embeddings: ``(n_cells, embedding_dim)`` float array.
        cell_index: DataFrame index from the input expression (for
            alignment with downstream score frames).
        backend: Short string identifying which embedder produced this
            (``'scgpt:official/v1.2.0'``, ``'scgpt:fallback'``, …).
        meta: Optional extra metadata the backend wants to surface.
    """

    embeddings: np.ndarray
    cell_index: pd.Index
    backend: str
    meta: Dict[str, Any] = field(default_factory=dict)

    @property
    def embedding_dim(self) -> int:
        return int(self.embeddings.shape[1])

    @property
    def n_cells(self) -> int:
        return int(self.embeddings.shape[0])

    def as_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(
            self.embeddings,
            index=self.cell_index,
            columns=[f"dim_{i}" for i in range(self.embedding_dim)],
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "backend": self.backend,
            "n_cells": self.n_cells,
            "embedding_dim": self.embedding_dim,
            "meta": dict(self.meta),
        }


class Embedder:
    """Abstract base class for cell-embedding backends.

    Concrete implementations override :meth:`embed`. The default
    :meth:`embedding_dim` property reads a class attribute of the same
    name — subclasses should set it at the class level (for models with
    a fixed embedding size) or override the property (if dim depends on
    checkpoint choice).
    """

    embedding_dim: int = 0

    def embed(self, expression: pd.DataFrame) -> EmbeddingResult:
        raise NotImplementedError

    def __call__(self, expression: pd.DataFrame) -> EmbeddingResult:
        return self.embed(expression)
