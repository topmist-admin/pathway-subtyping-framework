"""
Split-conformal prediction intervals for pathway scores.

Wraps any point-prediction scoring function and returns distribution-free
prediction intervals with user-specified coverage (e.g., 90%). Works for
both regression targets (continuous pathway scores) and classification
outputs (subtype assignment probabilities).

Algorithmic reference:
    Angelopoulos & Bates (2021). A Gentle Introduction to Conformal
    Prediction and Distribution-Free Uncertainty Quantification. arXiv:2107.07511.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, List, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ConformalInterval:
    """A single conformal prediction interval.

    For regression: (lower, upper) bracket around the point prediction.
    For classification: the set of labels in the prediction set (as a list).
    """

    prediction: float
    lower: float
    upper: float
    coverage: float
    n_calibration: int
    label_set: Optional[List[Any]] = None

    @property
    def width(self) -> float:
        return float(self.upper - self.lower)

    def contains(self, y: float) -> bool:
        return bool(self.lower <= y <= self.upper)

    def to_dict(self) -> dict:
        return {
            "prediction": self.prediction,
            "lower": self.lower,
            "upper": self.upper,
            "width": self.width,
            "coverage": self.coverage,
            "n_calibration": self.n_calibration,
            "label_set": self.label_set,
        }


class ConformalPathwayPredictor:
    """Split-conformal wrapper for any pathway score function.

    The predictor takes a scoring function ``score_fn(X) -> np.ndarray`` and,
    given a held-out calibration set (X_cal, y_cal), computes prediction
    intervals that achieve the requested coverage in the marginal sense.

    The guarantee is distribution-free and finite-sample:
        P(y_test in interval) >= 1 - alpha,
    where alpha = 1 - coverage, under the exchangeability assumption.

    Usage::

        predictor = ConformalPathwayPredictor(score_fn=my_msv_scorer, coverage=0.9)
        predictor.calibrate(X_cal, y_cal)
        intervals = predictor.predict(X_test)
        empirical = predictor.coverage_on(X_test, y_test)
    """

    # ---------------------------------------------------------------- init ---
    def __init__(
        self,
        score_fn: Callable[[np.ndarray], np.ndarray],
        coverage: float = 0.9,
        mode: str = "regression",
    ):
        if not 0.0 < coverage < 1.0:
            raise ValueError("coverage must be in (0, 1)")
        if mode not in {"regression", "classification"}:
            raise ValueError("mode must be 'regression' or 'classification'")

        self.score_fn = score_fn
        self.coverage = float(coverage)
        self.alpha = 1.0 - self.coverage
        self.mode = mode

        self._quantile: Optional[float] = None
        self._nonconformity: Optional[np.ndarray] = None
        self._n_calibration: int = 0
        self._labels: Optional[np.ndarray] = None

    # ----------------------------------------------------------- calibration ---
    def calibrate(
        self,
        X_cal: np.ndarray,
        y_cal: np.ndarray,
        labels: Optional[Sequence[Any]] = None,
    ) -> "ConformalPathwayPredictor":
        """Fit the conformal quantile on held-out calibration data.

        Args:
            X_cal: Feature matrix (n_cal, d) or sequence consumable by score_fn.
            y_cal: Ground-truth targets. Regression: (n_cal,) floats.
                   Classification: (n_cal,) integer labels.
            labels: For classification mode, optional label ordering. If None,
                    taken from ``np.unique(y_cal)``.

        Returns:
            Self (for chaining).
        """
        X_cal = np.asarray(X_cal)
        y_cal = np.asarray(y_cal)
        if len(y_cal) < 10:
            raise ValueError(
                "Need at least 10 calibration samples for a meaningful quantile"
            )

        predictions = np.asarray(self.score_fn(X_cal))

        if self.mode == "regression":
            nonconformity = np.abs(y_cal - predictions)
        else:  # classification
            if predictions.ndim != 2:
                raise ValueError(
                    "Classification score_fn must return (n, n_classes) probabilities"
                )
            if labels is None:
                labels = list(np.unique(y_cal))
            self._labels = np.asarray(labels)
            label_index = {lbl: i for i, lbl in enumerate(self._labels)}
            # APS-style nonconformity: 1 - p(true class)
            idx = np.array([label_index[int(y)] for y in y_cal])
            true_probs = predictions[np.arange(len(y_cal)), idx]
            nonconformity = 1.0 - true_probs

        # Finite-sample correction: ceil((n+1)(1-alpha)) / n
        n = len(nonconformity)
        q_level = np.ceil((n + 1) * (1.0 - self.alpha)) / n
        q_level = min(q_level, 1.0)
        self._quantile = float(np.quantile(nonconformity, q_level))
        self._nonconformity = nonconformity
        self._n_calibration = n

        logger.info(
            "[Conformal] calibrated: n=%d alpha=%.3f q=%.4f mode=%s",
            n, self.alpha, self._quantile, self.mode,
        )
        return self

    # ------------------------------------------------------------- predict ---
    def predict(self, X: np.ndarray) -> List[ConformalInterval]:
        """Return one ConformalInterval per input row."""
        if self._quantile is None:
            raise RuntimeError("calibrate() must be called before predict()")

        X = np.asarray(X)
        predictions = np.asarray(self.score_fn(X))

        intervals: List[ConformalInterval] = []

        if self.mode == "regression":
            for p in predictions:
                intervals.append(
                    ConformalInterval(
                        prediction=float(p),
                        lower=float(p - self._quantile),
                        upper=float(p + self._quantile),
                        coverage=self.coverage,
                        n_calibration=self._n_calibration,
                    )
                )
        else:  # classification
            assert self._labels is not None
            for probs in predictions:
                label_set = [
                    self._labels[i].item() if hasattr(self._labels[i], "item") else self._labels[i]
                    for i, p in enumerate(probs)
                    if (1.0 - p) <= self._quantile
                ]
                top_idx = int(np.argmax(probs))
                intervals.append(
                    ConformalInterval(
                        prediction=float(probs[top_idx]),
                        lower=float(max(0.0, probs[top_idx] - self._quantile)),
                        upper=float(min(1.0, probs[top_idx] + self._quantile)),
                        coverage=self.coverage,
                        n_calibration=self._n_calibration,
                        label_set=label_set,
                    )
                )

        return intervals

    # ----------------------------------------------------- empirical coverage ---
    def coverage_on(self, X: np.ndarray, y: np.ndarray) -> float:
        """Empirical coverage on (X, y) — fraction of y inside intervals."""
        intervals = self.predict(X)
        y = np.asarray(y)
        if self.mode == "regression":
            hits = sum(1 for ival, yi in zip(intervals, y) if ival.contains(float(yi)))
        else:
            hits = 0
            for ival, yi in zip(intervals, y):
                y_cmp = int(yi)
                if ival.label_set is None:
                    continue
                normalized = [int(lbl) if isinstance(lbl, (np.integer,)) else lbl for lbl in ival.label_set]
                if y_cmp in normalized:
                    hits += 1
        return float(hits) / max(len(intervals), 1)

    # ------------------------------------------------------------ utilities ---
    @property
    def quantile(self) -> Optional[float]:
        return self._quantile

    @property
    def n_calibration(self) -> int:
        return self._n_calibration
