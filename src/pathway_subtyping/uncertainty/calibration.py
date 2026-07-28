"""
Calibration assessment for probabilistic pathway predictions.

Provides:
    - reliability_curve: binned observed frequency vs predicted probability
    - ece: Expected Calibration Error
    - brier_score: mean squared error between probabilities and outcomes
    - CalibrationReport: aggregator producing all three with plot support

A well-calibrated classifier that reports p=0.8 is correct 80% of the time
across all its p~=0.8 predictions. Scoring functions that look accurate at
a single threshold can still be badly miscalibrated across the probability
range — this module catches that.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Metric primitives
# --------------------------------------------------------------------------- #


def reliability_curve(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_bins: int = 10,
    strategy: str = "uniform",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Binned reliability curve.

    Args:
        y_true: Binary ground truth in {0, 1}.
        y_prob: Predicted probability of the positive class.
        n_bins: Number of bins for the probability axis.
        strategy: ``'uniform'`` — equal-width bins in [0, 1].
                  ``'quantile'`` — equal-population bins.

    Returns:
        (bin_centers, mean_predicted, observed_frequency). ``nan`` for empty
        bins so plotting code can skip them.
    """
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.clip(np.asarray(y_prob).astype(float), 0.0, 1.0)

    if strategy == "uniform":
        edges = np.linspace(0.0, 1.0, n_bins + 1)
    elif strategy == "quantile":
        edges = np.quantile(y_prob, np.linspace(0.0, 1.0, n_bins + 1))
        edges[0], edges[-1] = 0.0, 1.0
    else:
        raise ValueError("strategy must be 'uniform' or 'quantile'")

    centers = 0.5 * (edges[:-1] + edges[1:])
    mean_pred = np.full(n_bins, np.nan)
    observed = np.full(n_bins, np.nan)

    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        if i == n_bins - 1:
            mask = (y_prob >= lo) & (y_prob <= hi)
        else:
            mask = (y_prob >= lo) & (y_prob < hi)
        if mask.sum() == 0:
            continue
        mean_pred[i] = float(y_prob[mask].mean())
        observed[i] = float(y_true[mask].mean())

    return centers, mean_pred, observed


def ece(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_bins: int = 15,
    strategy: str = "uniform",
) -> float:
    """Expected Calibration Error (ECE).

    ECE = sum_b (n_b / N) * |mean_predicted_b - observed_freq_b|
    """
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.clip(np.asarray(y_prob).astype(float), 0.0, 1.0)

    if strategy == "uniform":
        edges = np.linspace(0.0, 1.0, n_bins + 1)
    elif strategy == "quantile":
        edges = np.quantile(y_prob, np.linspace(0.0, 1.0, n_bins + 1))
        edges[0], edges[-1] = 0.0, 1.0
    else:
        raise ValueError("strategy must be 'uniform' or 'quantile'")

    n = len(y_prob)
    if n == 0:
        return 0.0

    total = 0.0
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        if i == n_bins - 1:
            mask = (y_prob >= lo) & (y_prob <= hi)
        else:
            mask = (y_prob >= lo) & (y_prob < hi)
        m = int(mask.sum())
        if m == 0:
            continue
        bin_pred = float(y_prob[mask].mean())
        bin_obs = float(y_true[mask].mean())
        total += (m / n) * abs(bin_pred - bin_obs)

    return float(total)


def brier_score(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    """Brier score — mean squared error between probabilities and outcomes."""
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.clip(np.asarray(y_prob).astype(float), 0.0, 1.0)
    return float(np.mean((y_prob - y_true) ** 2))


# --------------------------------------------------------------------------- #
# Aggregator
# --------------------------------------------------------------------------- #


@dataclass
class CalibrationReport:
    """Full calibration assessment for a set of probabilistic predictions.

    Usage::

        report = CalibrationReport.from_predictions(y_true, y_prob)
        print(report.summary())
        fig = report.plot()
    """

    y_true: np.ndarray
    y_prob: np.ndarray
    n_bins: int = 10
    strategy: str = "uniform"

    _ece: float = field(init=False)
    _brier: float = field(init=False)
    _centers: np.ndarray = field(init=False)
    _mean_pred: np.ndarray = field(init=False)
    _observed: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.y_true = np.asarray(self.y_true).astype(float)
        self.y_prob = np.clip(np.asarray(self.y_prob).astype(float), 0.0, 1.0)
        self._ece = ece(self.y_true, self.y_prob, n_bins=self.n_bins, strategy=self.strategy)
        self._brier = brier_score(self.y_true, self.y_prob)
        self._centers, self._mean_pred, self._observed = reliability_curve(
            self.y_true, self.y_prob, n_bins=self.n_bins, strategy=self.strategy
        )

    @classmethod
    def from_predictions(
        cls,
        y_true: np.ndarray,
        y_prob: np.ndarray,
        n_bins: int = 10,
        strategy: str = "uniform",
    ) -> "CalibrationReport":
        return cls(
            y_true=np.asarray(y_true),
            y_prob=np.asarray(y_prob),
            n_bins=n_bins,
            strategy=strategy,
        )

    # --------------------------------------------------------- accessors ---
    @property
    def ece(self) -> float:
        return self._ece

    @property
    def brier(self) -> float:
        return self._brier

    @property
    def reliability(self) -> Dict[str, np.ndarray]:
        return {
            "bin_centers": self._centers,
            "mean_predicted": self._mean_pred,
            "observed_frequency": self._observed,
        }

    def summary(self) -> str:
        n = len(self.y_prob)
        return f"CalibrationReport(n={n}, ECE={self._ece:.4f}, Brier={self._brier:.4f})"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_samples": int(len(self.y_prob)),
            "ece": self._ece,
            "brier": self._brier,
            "n_bins": self.n_bins,
            "strategy": self.strategy,
            "reliability": {k: v.tolist() for k, v in self.reliability.items()},
        }

    # ------------------------------------------------------------ plotting ---
    def plot(
        self,
        ax: Optional[Any] = None,
        show_histogram: bool = True,
    ):
        """Reliability diagram. Returns the matplotlib figure.

        Raises ImportError if matplotlib is unavailable.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "matplotlib required for CalibrationReport.plot(); "
                "install pathway-subtyping[viz]"
            ) from exc

        created_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))
            created_fig = True
        else:
            fig = ax.figure

        ax.plot([0, 1], [0, 1], linestyle="--", color="gray", label="perfect")
        mask = ~np.isnan(self._mean_pred)
        ax.plot(
            self._mean_pred[mask],
            self._observed[mask],
            marker="o",
            label=f"observed (ECE={self._ece:.3f})",
        )
        ax.set_xlabel("Mean predicted probability")
        ax.set_ylabel("Observed frequency")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title("Reliability diagram")
        ax.legend(loc="upper left")
        ax.grid(True, alpha=0.3)

        if show_histogram and created_fig:
            # Overlay a small histogram of predicted probs on a twin axis
            ax2 = ax.twinx()
            ax2.hist(self.y_prob, bins=20, alpha=0.2, color="steelblue")
            ax2.set_ylabel("Predicted-prob count", color="steelblue")

        return fig
