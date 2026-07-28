"""
Active-learning sample selection.

``ActiveSampleSelector`` picks the most informative samples to label
from a pool of unlabelled data given a fixed budget. Three strategies:

    - ``uncertainty`` — pick samples whose F1 Bayesian-GMM posterior
      is highest-entropy (most uncertain subtype assignment).
    - ``diversity`` — pick samples that maximise coverage of feature
      space via greedy k-center on pairwise distances.
    - ``hybrid`` — weighted combination of the two.

All three return deterministic results for a fixed random seed. The
uncertainty strategy plugs into the F1 ``BayesianPathwayGMM`` by
default; users can supply any per-sample score via the lower-level
constructor.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class SelectionStrategy(str, Enum):
    UNCERTAINTY = "uncertainty"
    DIVERSITY = "diversity"
    HYBRID = "hybrid"


# --------------------------------------------------------------------------- #
# Result type
# --------------------------------------------------------------------------- #


@dataclass
class SelectionResult:
    """Outcome of an :class:`ActiveSampleSelector` run.

    Attributes:
        selected_indices: Row indices of selected samples in the pool.
        scores: Per-sample score used by the chosen strategy; higher =
            more preferred. ``nan`` for samples where the score wasn't
            computed (e.g., diversity ignores uncertainty scores).
        strategy: Which strategy produced the selection.
        budget: Label budget honoured by the selection.
    """

    selected_indices: np.ndarray
    scores: pd.Series
    strategy: SelectionStrategy
    budget: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "strategy": self.strategy.value,
            "budget": int(self.budget),
            "n_selected": int(len(self.selected_indices)),
            "selected_indices": self.selected_indices.tolist(),
        }


# --------------------------------------------------------------------------- #
# Core algorithms
# --------------------------------------------------------------------------- #


def _categorical_entropy(probs: np.ndarray) -> np.ndarray:
    """Per-row Shannon entropy on a (n, k) probability matrix."""
    eps = 1e-12
    return -np.sum(probs * np.log(probs + eps), axis=1)


def _greedy_k_center(features: np.ndarray, k: int, seed_idx: int = 0) -> np.ndarray:
    """Greedy k-center: pick the point farthest from the current set.

    Runs in O(n*k) distance evaluations. Deterministic given ``seed_idx``.
    """
    n = len(features)
    if k >= n:
        return np.arange(n)
    # Initial min distance array (to a virtual anchor at seed_idx)
    picked = [int(seed_idx)]
    dists = np.linalg.norm(features - features[seed_idx], axis=1)
    for _ in range(k - 1):
        next_idx = int(np.argmax(dists))
        picked.append(next_idx)
        new_d = np.linalg.norm(features - features[next_idx], axis=1)
        dists = np.minimum(dists, new_d)
    return np.asarray(picked, dtype=int)


# --------------------------------------------------------------------------- #
# Selector
# --------------------------------------------------------------------------- #


class ActiveSampleSelector:
    """Select informative samples to label under a fixed budget.

    Usage::

        selector = ActiveSampleSelector(strategy="uncertainty")
        result = selector.select(features=pool_features, budget=50, probs=gmm_probs)

        # hybrid with a custom uncertainty source:
        selector = ActiveSampleSelector(strategy="hybrid", alpha=0.6)
        result = selector.select(
            features=pool_features,
            budget=50,
            uncertainty_scores=custom_scores,
        )
    """

    def __init__(
        self,
        strategy: SelectionStrategy | str = SelectionStrategy.UNCERTAINTY,
        alpha: float = 0.5,
        seed: Optional[int] = None,
    ) -> None:
        self.strategy = SelectionStrategy(strategy)
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("alpha must be in [0, 1]")
        self.alpha = float(alpha)
        self.seed = seed

    # -------------------------------------------------------------- select ---
    def select(
        self,
        features: np.ndarray | pd.DataFrame,
        budget: int,
        probs: Optional[np.ndarray] = None,
        uncertainty_scores: Optional[np.ndarray] = None,
    ) -> SelectionResult:
        """Pick at most ``budget`` samples from the pool.

        Args:
            features: Pool features (n_samples, d). Required for diversity
                and hybrid strategies.
            budget: Maximum number of samples to select.
            probs: (n_samples, n_components) Bayesian-GMM posterior.
                Converted to per-sample entropy. Used by uncertainty
                and hybrid strategies when ``uncertainty_scores`` is
                not supplied.
            uncertainty_scores: Optional custom per-sample uncertainty
                (higher = more uncertain). Overrides ``probs`` if given.

        Returns:
            SelectionResult with the selected indices, per-sample
            scores, strategy, and budget.
        """
        X = (
            features.to_numpy(dtype=float)
            if isinstance(features, pd.DataFrame)
            else np.asarray(features, dtype=float)
        )
        n = len(X)
        if budget < 1:
            raise ValueError("budget must be >= 1")
        budget = min(int(budget), n)

        # Score computation per strategy
        if self.strategy == SelectionStrategy.UNCERTAINTY:
            scores = self._uncertainty_scores(n, probs, uncertainty_scores)
            ranked = np.argsort(-scores)
            selected = ranked[:budget]
        elif self.strategy == SelectionStrategy.DIVERSITY:
            seed_idx = 0 if self.seed is None else int(self.seed) % n
            selected = _greedy_k_center(X, budget, seed_idx=seed_idx)
            scores = np.full(n, np.nan)
            scores[selected] = 1.0  # diversity has no continuous score
        elif self.strategy == SelectionStrategy.HYBRID:
            u = self._uncertainty_scores(n, probs, uncertainty_scores)
            # Normalise u to [0, 1]
            if u.max() - u.min() > 1e-12:
                u_norm = (u - u.min()) / (u.max() - u.min())
            else:
                u_norm = np.zeros_like(u)
            # Greedy hybrid: at each step, pick the sample that maximises
            # alpha * u_norm + (1 - alpha) * min_distance_to_selected.
            # Distances are renormalised each step to keep scales comparable.
            seed_idx = int(np.argmax(u_norm))
            selected_list: List[int] = [seed_idx]
            min_d = np.linalg.norm(X - X[seed_idx], axis=1)
            for _ in range(budget - 1):
                d_norm = (min_d - min_d.min()) / max(min_d.max() - min_d.min(), 1e-12)
                score = self.alpha * u_norm + (1.0 - self.alpha) * d_norm
                score[selected_list] = -np.inf
                next_idx = int(np.argmax(score))
                selected_list.append(next_idx)
                new_d = np.linalg.norm(X - X[next_idx], axis=1)
                min_d = np.minimum(min_d, new_d)
            selected = np.asarray(selected_list, dtype=int)
            scores = u
        else:  # pragma: no cover - enum exhausts above
            raise RuntimeError(f"unknown strategy: {self.strategy}")

        score_series = pd.Series(
            scores,
            index=(features.index if isinstance(features, pd.DataFrame) else pd.RangeIndex(n)),
            name="score",
        )
        logger.info(
            "[ActiveSampleSelector] strategy=%s budget=%d selected=%d",
            self.strategy.value,
            budget,
            len(selected),
        )
        return SelectionResult(
            selected_indices=np.asarray(selected),
            scores=score_series,
            strategy=self.strategy,
            budget=budget,
        )

    # ------------------------------------------------------ helpers ---
    @staticmethod
    def _uncertainty_scores(
        n: int,
        probs: Optional[np.ndarray],
        uncertainty_scores: Optional[np.ndarray],
    ) -> np.ndarray:
        if uncertainty_scores is not None:
            u = np.asarray(uncertainty_scores, dtype=float)
            if len(u) != n:
                raise ValueError(f"uncertainty_scores has length {len(u)}; expected {n}")
            return u
        if probs is None:
            raise ValueError(
                "strategy requires either probs (categorical posterior) "
                "or uncertainty_scores (per-sample array)"
            )
        probs_arr = np.asarray(probs, dtype=float)
        if probs_arr.shape[0] != n:
            raise ValueError(f"probs has {probs_arr.shape[0]} rows; expected {n}")
        return _categorical_entropy(probs_arr)
