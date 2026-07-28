"""
Cross-platform alignment of pathway scores.

Given platform-labelled cells (or samples) with both pathway scores and a
platform-invariant embedding, ``CrossPlatformAligner`` learns a per-platform
linear map that removes platform-specific drift from pathway scores while
preserving biological variation visible in the embedding.

The alignment is ComBat-style but anchored on an embedding rather than on
hand-specified biological covariates:

    ref_score_i,p = platform_score_i,p - drift_p(embedding_i)

Per platform p, we fit ``drift_p`` as the residual pathway-score mean after
regressing pathway scores on the embedding pooled across platforms; this
captures the platform-specific offset that cannot be explained by shared
biology.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Output of CrossPlatformAligner.transform."""

    aligned_scores: pd.DataFrame
    platform: pd.Series
    per_cell_shift: pd.Series  # |aligned - raw| averaged over pathways
    per_platform_drift: Dict[str, Dict[str, float]]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_cells": int(len(self.aligned_scores)),
            "n_pathways": int(self.aligned_scores.shape[1]),
            "platforms": sorted(self.platform.unique().tolist()),
            "per_platform_drift": self.per_platform_drift,
            "per_cell_shift_mean": float(self.per_cell_shift.mean()),
            "per_cell_shift_max": float(self.per_cell_shift.max()),
        }


class CrossPlatformAligner:
    """Fit and apply an embedding-anchored per-platform alignment.

    Usage::

        aligner = CrossPlatformAligner()
        aligner.fit(scores, platforms, embeddings)
        result = aligner.transform(scores, platforms, embeddings)
        # result.aligned_scores is the drift-removed pathway-score DataFrame

    The fit pools cells from all platforms, regresses each pathway on the
    embedding, and records the mean residual by platform — that residual
    becomes the additive correction removed at transform time.
    """

    def __init__(
        self,
        ridge_alpha: float = 1e-3,
        reference_platform: Optional[str] = None,
    ):
        self.ridge_alpha = float(ridge_alpha)
        self.configured_reference_platform = reference_platform
        self._pathway_names: Optional[List[str]] = None
        self._platform_betas: Dict[str, np.ndarray] = {}  # plat -> (d_embed+1, n_pathways)
        self._platform_offsets: Dict[str, np.ndarray] = {}  # plat -> (n_pathways,) for reporting
        self._reference_platform: Optional[str] = None
        self._fitted: bool = False

    # --------------------------------------------------------- internals ---
    @staticmethod
    def _design(embeddings: np.ndarray) -> np.ndarray:
        return np.hstack([embeddings, np.ones((len(embeddings), 1))])

    def _ridge_fit(self, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        """Closed-form ridge regression: (X.T X + alpha I)^-1 X.T Y."""
        d = X.shape[1]
        reg = self.ridge_alpha * np.eye(d)
        reg[-1, -1] = 0.0  # don't regularize intercept
        return np.linalg.solve(X.T @ X + reg, X.T @ Y)

    # -------------------------------------------------------------- fit ---
    def fit(
        self,
        scores: pd.DataFrame,
        platforms: Sequence[str],
        embeddings: np.ndarray,
    ) -> "CrossPlatformAligner":
        if len(scores) != len(platforms) or len(scores) != len(embeddings):
            raise ValueError("scores, platforms, embeddings must all have same length")
        if len(scores) < 5:
            raise ValueError("need >= 5 cells to fit a cross-platform alignment")

        self._pathway_names = list(scores.columns)
        platforms_ = pd.Series([str(p) for p in platforms], index=scores.index)
        Y = scores.to_numpy(dtype=float)
        emb = np.asarray(embeddings, dtype=float)

        unique_platforms = sorted(platforms_.unique())
        if len(unique_platforms) == 0:
            raise ValueError("at least one platform required")

        # Choose reference platform: user-specified or the largest by count.
        if self.configured_reference_platform in unique_platforms:
            self._reference_platform = self.configured_reference_platform
        else:
            counts = platforms_.value_counts()
            self._reference_platform = str(counts.idxmax())

        # Fit a separate biology regressor per platform:
        #   scores_p = embedding @ beta_p + intercept_p + residual
        # The per-platform betas capture both platform-specific mean shift
        # (the intercept row) and platform-specific coupling between
        # embedding-space biology and pathway score (the embedding rows).
        # Transform step expresses any platform's data in the reference-
        # platform frame by swapping beta_p with beta_ref.
        self._platform_betas = {}
        self._platform_offsets = {}
        for plat in unique_platforms:
            mask = (platforms_ == plat).to_numpy()
            if mask.sum() < 3:
                logger.warning(
                    "[CrossPlatformAligner] platform %r has only %d cells; "
                    "skipping platform-specific fit (will fall back to reference)",
                    plat,
                    int(mask.sum()),
                )
                continue
            X_plat = self._design(emb[mask])
            theta = self._ridge_fit(X_plat, Y[mask])
            # theta shape: (d_emb + 1, n_pathways); last row is intercept
            self._platform_betas[plat] = theta
            # Reporting-only offset: mean deviation from the reference intercept
            self._platform_offsets[plat] = theta[-1]

        # Normalise reporting offsets relative to reference intercept.
        ref_intercept = self._platform_betas[self._reference_platform][-1]
        for plat in self._platform_offsets:
            self._platform_offsets[plat] = self._platform_offsets[plat] - ref_intercept

        self._fitted = True
        logger.info(
            "[CrossPlatformAligner] fit complete: ref=%s platforms=%s pathways=%d",
            self._reference_platform,
            unique_platforms,
            Y.shape[1],
        )
        return self

    # --------------------------------------------------------- transform ---
    def transform(
        self,
        scores: pd.DataFrame,
        platforms: Sequence[str],
        embeddings: np.ndarray,
    ) -> AlignmentResult:
        self._require_fitted()
        if list(scores.columns) != self._pathway_names:
            raise ValueError("scores columns must match those seen at fit time")

        assert self._reference_platform is not None
        platforms_ = pd.Series([str(p) for p in platforms], index=scores.index, name="platform")
        emb = np.asarray(embeddings, dtype=float)

        Y_raw = scores.to_numpy(dtype=float)
        Y_out = Y_raw.copy()

        beta_ref = self._platform_betas[self._reference_platform]

        # For every non-reference platform, shift scores by the embedding-
        # dependent delta (beta_ref - beta_p) applied to each cell's embedding.
        # This lands cells in the reference-platform frame:
        #   aligned = scores_p + (beta_ref - beta_p) @ emb_cell
        # Cells whose platform was never seen at fit time remain unchanged.
        unseen: List[str] = []
        for plat in platforms_.unique():
            mask = (platforms_ == plat).to_numpy()
            if mask.sum() == 0:
                continue
            if plat not in self._platform_betas:
                unseen.append(plat)
                continue
            if plat == self._reference_platform:
                continue
            beta_p = self._platform_betas[plat]
            design = self._design(emb[mask])
            delta = design @ (beta_ref - beta_p)
            Y_out[mask] = Y_raw[mask] + delta

        if unseen:
            logger.warning(
                "[CrossPlatformAligner] unseen platforms at transform time: "
                "%s (no shift applied)",
                sorted(set(unseen)),
            )

        aligned = pd.DataFrame(Y_out, index=scores.index, columns=scores.columns)

        per_cell_shift = pd.Series(
            np.abs(Y_out - Y_raw).mean(axis=1),
            index=scores.index,
            name="harmonization_shift",
        )

        per_platform_drift: Dict[str, Dict[str, float]] = {}
        for plat, offset in self._platform_offsets.items():
            per_platform_drift[plat] = {
                name: float(val) for name, val in zip(self._pathway_names, offset)
            }

        return AlignmentResult(
            aligned_scores=aligned,
            platform=platforms_,
            per_cell_shift=per_cell_shift,
            per_platform_drift=per_platform_drift,
        )

    def fit_transform(
        self,
        scores: pd.DataFrame,
        platforms: Sequence[str],
        embeddings: np.ndarray,
    ) -> AlignmentResult:
        self.fit(scores, platforms, embeddings)
        return self.transform(scores, platforms, embeddings)

    # --------------------------------------------------------- accessors ---
    @property
    def platforms(self) -> List[str]:
        return sorted(self._platform_offsets.keys())

    @property
    def pathway_names(self) -> List[str]:
        self._require_fitted()
        assert self._pathway_names is not None
        return list(self._pathway_names)

    def platform_offset(self, platform: str) -> pd.Series:
        self._require_fitted()
        if platform not in self._platform_offsets:
            raise KeyError(f"platform '{platform}' not seen at fit time")
        return pd.Series(
            self._platform_offsets[platform],
            index=self._pathway_names,
        )

    def _require_fitted(self) -> None:
        if not self._fitted:
            raise RuntimeError("CrossPlatformAligner.fit() must be called first")
