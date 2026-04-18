"""
Tests for v0.6 F12 active sample selection.

Covers:
    - Validation: bad strategy, bad budget, bad alpha, probs/features
      length mismatches
    - Uncertainty strategy prefers high-entropy samples
    - Diversity strategy deterministically covers the feature space
    - Hybrid strategy combines both
    - Roadmap acceptance: on a simulated autism-style cohort, active
      learning using 40% of labels reaches >= 90% of the full-cohort
      1-NN classification accuracy
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.active import (
    ActiveSampleSelector,
    SelectionResult,
    SelectionStrategy,
)


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #

class TestValidation:

    def test_rejects_bad_strategy(self):
        with pytest.raises(ValueError):
            ActiveSampleSelector(strategy="nonexistent")

    def test_rejects_bad_alpha(self):
        with pytest.raises(ValueError, match="alpha"):
            ActiveSampleSelector(strategy="hybrid", alpha=1.5)

    def test_rejects_zero_budget(self):
        sel = ActiveSampleSelector(strategy="uncertainty")
        with pytest.raises(ValueError, match="budget"):
            sel.select(
                features=np.zeros((5, 3)), budget=0,
                probs=np.ones((5, 2)) / 2,
            )

    def test_uncertainty_requires_probs_or_scores(self):
        sel = ActiveSampleSelector(strategy="uncertainty")
        with pytest.raises(ValueError, match="probs"):
            sel.select(features=np.zeros((5, 3)), budget=2)

    def test_wrong_probs_length_raises(self):
        sel = ActiveSampleSelector(strategy="uncertainty")
        with pytest.raises(ValueError, match="probs has"):
            sel.select(
                features=np.zeros((5, 3)), budget=2,
                probs=np.ones((3, 2)) / 2,
            )

    def test_wrong_uncertainty_scores_length_raises(self):
        sel = ActiveSampleSelector(strategy="uncertainty")
        with pytest.raises(ValueError, match="uncertainty_scores has"):
            sel.select(
                features=np.zeros((5, 3)), budget=2,
                uncertainty_scores=np.zeros(3),
            )


# --------------------------------------------------------------------------- #
# Uncertainty strategy
# --------------------------------------------------------------------------- #

class TestUncertaintyStrategy:

    def test_prefers_high_entropy(self):
        """A sample with uniform posterior (max entropy) should rank
        above samples whose posterior collapses on one component."""
        probs = np.array([
            [0.9, 0.1],    # confident
            [0.5, 0.5],    # most uncertain
            [0.95, 0.05],  # most confident
            [0.6, 0.4],    # mid
        ])
        sel = ActiveSampleSelector(strategy="uncertainty")
        result = sel.select(
            features=np.zeros((4, 2)), budget=2, probs=probs,
        )
        # Top-2 should be the uniform posterior (idx 1) and the mid one (idx 3)
        top = set(result.selected_indices.tolist())
        assert 1 in top  # the most uncertain must be selected
        assert 2 not in top  # the most confident must NOT be selected

    def test_custom_uncertainty_scores_used(self):
        scores = np.array([0.1, 0.9, 0.2, 0.8])
        sel = ActiveSampleSelector(strategy="uncertainty")
        result = sel.select(
            features=np.zeros((4, 2)), budget=2, uncertainty_scores=scores,
        )
        top = set(result.selected_indices.tolist())
        assert top == {1, 3}  # top-2 uncertainty scores


# --------------------------------------------------------------------------- #
# Diversity strategy
# --------------------------------------------------------------------------- #

class TestDiversityStrategy:

    def test_returns_budget_unique_points(self):
        rng = np.random.default_rng(0)
        X = rng.standard_normal((100, 5))
        sel = ActiveSampleSelector(strategy="diversity", seed=0)
        result = sel.select(features=X, budget=10)
        assert len(result.selected_indices) == 10
        assert len(set(result.selected_indices.tolist())) == 10

    def test_deterministic_for_fixed_seed(self):
        rng = np.random.default_rng(0)
        X = rng.standard_normal((60, 4))
        a = ActiveSampleSelector(strategy="diversity", seed=42).select(
            features=X, budget=8,
        )
        b = ActiveSampleSelector(strategy="diversity", seed=42).select(
            features=X, budget=8,
        )
        np.testing.assert_array_equal(a.selected_indices, b.selected_indices)


# --------------------------------------------------------------------------- #
# Hybrid strategy
# --------------------------------------------------------------------------- #

class TestHybridStrategy:

    def test_hybrid_uses_both_signals(self):
        rng = np.random.default_rng(0)
        X = rng.standard_normal((50, 3))
        probs = rng.dirichlet(np.ones(2), size=50)
        a = ActiveSampleSelector(strategy="hybrid", alpha=1.0).select(
            features=X, budget=10, probs=probs,
        )
        b = ActiveSampleSelector(strategy="hybrid", alpha=0.0).select(
            features=X, budget=10, probs=probs,
        )
        # alpha=1.0 is pure uncertainty; alpha=0.0 is pure-distance
        # The selected sets should typically differ.
        assert set(a.selected_indices) != set(b.selected_indices)


# --------------------------------------------------------------------------- #
# Roadmap acceptance: 90% accuracy at 40% labels
# --------------------------------------------------------------------------- #

def _autism_style_cohort(n_per_class: int = 120, n_pathways: int = 10, seed: int = 0):
    """Three-subtype synthetic cohort modelled after autism-style PSF data."""
    rng = np.random.default_rng(seed)
    cluster_means = rng.standard_normal((3, n_pathways)) * 1.3
    cluster_ids = np.repeat(np.arange(3), n_per_class)
    base = cluster_means[cluster_ids]
    noise = rng.normal(0, 0.7, size=base.shape)
    X = base + noise
    features = pd.DataFrame(
        X, columns=[f"PATH_{i}" for i in range(n_pathways)],
    )
    labels = pd.Series(cluster_ids, name="cluster")
    return features, labels


def _nn_accuracy(
    train_X: np.ndarray, train_y: np.ndarray,
    test_X: np.ndarray, test_y: np.ndarray,
) -> float:
    """1-NN accuracy from train set onto test set (cosine similarity)."""
    if len(train_X) == 0:
        # Majority-class fallback
        return float(
            (np.full_like(test_y, np.bincount(test_y).argmax()) == test_y).mean()
        )
    tn = np.linalg.norm(train_X, axis=1, keepdims=True); tn[tn == 0] = 1.0
    en = np.linalg.norm(test_X, axis=1, keepdims=True); en[en == 0] = 1.0
    sim = (test_X / en) @ (train_X / tn).T
    predictions = train_y[sim.argmax(axis=1)]
    return float((predictions == test_y).mean())


class TestActiveLearningAcceptance:

    def test_40pct_active_reaches_90pct_full_accuracy(self):
        """Roadmap: active learning reaches 90% of full-cohort accuracy
        while using only 40% of the labels."""
        from sklearn.mixture import GaussianMixture  # already a PSF dep

        features, labels = _autism_style_cohort(n_per_class=120, seed=0)
        rng = np.random.default_rng(0)
        perm = rng.permutation(len(features))
        # 20% held-out test set, 80% pool
        test_n = len(features) // 5
        test_idx = perm[:test_n]
        pool_idx = perm[test_n:]

        X_test = features.iloc[test_idx].to_numpy(dtype=float)
        y_test = labels.iloc[test_idx].to_numpy(dtype=int)
        X_pool = features.iloc[pool_idx].to_numpy(dtype=float)
        y_pool = labels.iloc[pool_idx].to_numpy(dtype=int)

        # Full-cohort accuracy: train on ALL pool labels
        acc_full = _nn_accuracy(X_pool, y_pool, X_test, y_test)

        # Active learning budget = 40% of pool
        budget = int(0.40 * len(X_pool))

        # Uncertainty scores from a 3-component GMM posterior
        gmm = GaussianMixture(n_components=3, random_state=0).fit(X_pool)
        probs = gmm.predict_proba(X_pool)

        # Hybrid strategy — best of both worlds
        result = ActiveSampleSelector(strategy="hybrid", alpha=0.5).select(
            features=X_pool, budget=budget, probs=probs,
        )
        X_train_al = X_pool[result.selected_indices]
        y_train_al = y_pool[result.selected_indices]
        acc_al = _nn_accuracy(X_train_al, y_train_al, X_test, y_test)

        assert acc_al >= 0.90 * acc_full, (
            f"active learning acc {acc_al:.3f} < 0.90 * full acc "
            f"{acc_full:.3f} = {0.90 * acc_full:.3f}"
        )

    def test_selection_result_to_dict(self):
        probs = np.array([[0.7, 0.3], [0.4, 0.6], [0.9, 0.1]])
        sel = ActiveSampleSelector(strategy="uncertainty")
        result = sel.select(
            features=np.zeros((3, 2)), budget=2, probs=probs,
        )
        d = result.to_dict()
        assert d["strategy"] == "uncertainty"
        assert d["budget"] == 2
        assert d["n_selected"] == 2
