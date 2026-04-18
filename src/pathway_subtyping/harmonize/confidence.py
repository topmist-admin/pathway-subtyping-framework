"""
Harmonization confidence and per-platform drift reporting.

``HarmonizationReport`` wraps an ``AlignmentResult`` and adds:

    - per-cell confidence: 1 / (1 + z-scored shift magnitude)
    - per-platform absolute drift per pathway
    - correlation between confidence and a quality covariate (e.g., read depth)
    - summary + plotting helpers

A low-confidence cell is one whose raw pathway scores had to move a lot
to align with the cross-platform reference — a useful flag for triage
alongside downstream MSV scoring or cluster interpretation.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from .align import AlignmentResult

logger = logging.getLogger(__name__)


@dataclass
class HarmonizationReport:
    """Cross-platform harmonization diagnostics.

    Usage::

        report = HarmonizationReport.from_alignment(result)
        report.confidence            # pd.Series indexed by cell id
        report.per_platform_drift    # dict[platform -> dict[pathway -> |drift|]]
        fig = report.plot_drift()
    """

    alignment: AlignmentResult
    _confidence: pd.Series = field(init=False)
    _platform_abs_drift: Dict[str, Dict[str, float]] = field(init=False)

    def __post_init__(self) -> None:
        shift = self.alignment.per_cell_shift
        if len(shift) == 0:
            raise ValueError("alignment contains no cells")
        # z-score the shift; clamp negatives before inverse-sigmoid mapping
        mu = float(shift.mean())
        sigma = float(shift.std(ddof=1)) if len(shift) > 1 else 1.0
        sigma = max(sigma, 1e-9)
        z = (shift - mu) / sigma
        # confidence high (near 1) when shift is small (low z); low (near 0)
        # when shift is large (high z). Sigmoid with negative sign centers
        # the score near 0.5 on the mean shift.
        self._confidence = pd.Series(
            1.0 / (1.0 + np.exp(z)),
            index=shift.index,
            name="harmonization_confidence",
        )

        self._platform_abs_drift = {}
        for plat, per_pw in self.alignment.per_platform_drift.items():
            self._platform_abs_drift[plat] = {
                pw: float(abs(val)) for pw, val in per_pw.items()
            }

    # --------------------------------------------------------- factories ---
    @classmethod
    def from_alignment(cls, alignment: AlignmentResult) -> "HarmonizationReport":
        return cls(alignment=alignment)

    # --------------------------------------------------------- accessors ---
    @property
    def confidence(self) -> pd.Series:
        return self._confidence

    @property
    def per_platform_drift(self) -> Dict[str, Dict[str, float]]:
        return self._platform_abs_drift

    @property
    def mean_platform_drift(self) -> Dict[str, float]:
        return {
            plat: float(np.mean(list(per_pw.values())))
            for plat, per_pw in self._platform_abs_drift.items()
        }

    # --------------------------------------------------------- summaries ---
    def correlate_with_quality(self, quality: pd.Series) -> float:
        """Spearman rho between confidence and a per-cell quality covariate.

        Higher-quality cells are expected to harmonize more cleanly, so
        a positive correlation (quality increases -> confidence increases)
        is the healthy signal. Returns NaN if quality is constant.
        """
        aligned_quality = quality.reindex(self._confidence.index)
        valid = aligned_quality.notna() & self._confidence.notna()
        if valid.sum() < 3:
            return float("nan")
        q = aligned_quality[valid].to_numpy()
        c = self._confidence[valid].to_numpy()
        if np.std(q) == 0 or np.std(c) == 0:
            return float("nan")
        return float(
            pd.Series(q).corr(pd.Series(c), method="spearman")
        )

    def summary(self) -> str:
        platform_drifts = ", ".join(
            f"{p}={d:.3f}" for p, d in self.mean_platform_drift.items()
        )
        return (
            f"HarmonizationReport(n={len(self._confidence)}, "
            f"mean_confidence={self._confidence.mean():.3f}, "
            f"mean_platform_drift={{{platform_drifts}}})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_cells": int(len(self._confidence)),
            "mean_confidence": float(self._confidence.mean()),
            "min_confidence": float(self._confidence.min()),
            "max_confidence": float(self._confidence.max()),
            "per_platform_drift": self._platform_abs_drift,
            "mean_platform_drift": self.mean_platform_drift,
        }

    # ----------------------------------------------------------- plotting ---
    def plot_drift(self, ax: Optional[Any] = None):
        """Per-platform mean |drift| per pathway. Returns the matplotlib figure."""
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "matplotlib required for HarmonizationReport.plot_drift(); "
                "install pathway-subtyping[viz]"
            ) from exc

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))
        else:
            fig = ax.figure

        df = pd.DataFrame(self._platform_abs_drift)
        df.plot(kind="bar", ax=ax)
        ax.set_xlabel("pathway")
        ax.set_ylabel("|drift|")
        ax.set_title("Per-platform pathway drift")
        ax.grid(True, alpha=0.3, axis="y")
        return fig
