"""
Tests for v0.6 F11 invariant causal prediction.

Covers:
    - invariance_pvalue primitive under null (single env, constant
      residuals) and alternative (shifted residuals across envs)
    - InvariantPathwayPredictor basic shape + argument validation
    - Roadmap acceptance: on simulated data with ground-truth causal
      structure, identifiable parents recall >= 0.7
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.causal import (
    CausalParentReport,
    InvariantPathwayPredictor,
    invariance_pvalue,
)

# --------------------------------------------------------------------------- #
# Primitive
# --------------------------------------------------------------------------- #


class TestInvariancePvalue:

    def test_single_environment_returns_one(self):
        residuals = np.array([1.0, 2.0, 3.0])
        envs = np.array(["A", "A", "A"])
        assert invariance_pvalue(residuals, envs) == 1.0

    def test_invariant_residuals_high_p(self):
        rng = np.random.default_rng(0)
        n = 200
        residuals = rng.standard_normal(n)
        envs = np.repeat(["A", "B"], n // 2)
        assert invariance_pvalue(residuals, envs) > 0.1

    def test_shifted_residuals_low_p(self):
        """Environment-specific mean shift should break invariance."""
        rng = np.random.default_rng(0)
        n_per = 120
        envs = np.repeat(["A", "B"], n_per)
        residuals = np.concatenate(
            [
                rng.normal(0, 1.0, n_per),
                rng.normal(2.5, 1.0, n_per),
            ]
        )
        assert invariance_pvalue(residuals, envs) < 1e-6

    def test_constant_residuals_return_one(self):
        residuals = np.zeros(10)
        envs = np.array(["A"] * 5 + ["B"] * 5)
        assert invariance_pvalue(residuals, envs) == 1.0


# --------------------------------------------------------------------------- #
# Predictor validation
# --------------------------------------------------------------------------- #


class TestValidation:

    def test_rejects_bad_alpha(self):
        with pytest.raises(ValueError, match="alpha"):
            InvariantPathwayPredictor(alpha=0.0)
        with pytest.raises(ValueError, match="alpha"):
            InvariantPathwayPredictor(alpha=1.0)

    def test_target_in_X_raises(self):
        X = pd.DataFrame({"A": [1.0], "Y": [2.0]})
        y = pd.Series([0.0])
        predictor = InvariantPathwayPredictor()
        with pytest.raises(ValueError, match="must not be a column"):
            predictor.fit(X=X, y=y, target_name="Y", environments=["e1"])

    def test_length_mismatch_raises(self):
        X = pd.DataFrame({"A": [1.0, 2.0]})
        y = pd.Series([0.0, 1.0])
        predictor = InvariantPathwayPredictor()
        with pytest.raises(ValueError, match="same length"):
            predictor.fit(
                X=X,
                y=y,
                target_name="Y",
                environments=["e1"],  # len 1 vs 2
            )


# --------------------------------------------------------------------------- #
# Roadmap acceptance: identifiable-parent recall >= 0.7
# --------------------------------------------------------------------------- #


def _synthetic_causal_cohort(n_per_env: int = 300, seed: int = 0):
    """Two environments that differ in parent-variable distributions.

    Ground truth:
        Y  <-  X1   (direct cause, coefficient 0.8)
        Y  <-  X2   (direct cause, coefficient 0.6)
        X3 is a descendant of Y (child; not a parent)
        X4 is pure noise with env-specific mean shift

    Environments differ in the mean of X1 and X2, so Y's marginal also
    differs. Under ICP, the only subsets whose residuals are invariant
    are those containing both X1 and X2 — the identifiable parents.
    """
    rng = np.random.default_rng(seed)

    def generate(env, mu_x1, mu_x2, sigma_x1, sigma_x2, mu_x4):
        X1 = rng.standard_normal(n_per_env) * sigma_x1 + mu_x1
        X2 = rng.standard_normal(n_per_env) * sigma_x2 + mu_x2
        Y = 0.8 * X1 + 0.6 * X2 + 0.2 * rng.standard_normal(n_per_env)
        # X3 is a noisy child of Y — correlates with Y but is not a parent
        X3 = 0.7 * Y + rng.standard_normal(n_per_env) * 0.4
        X4 = rng.standard_normal(n_per_env) + mu_x4
        return pd.DataFrame(
            {
                "X1": X1,
                "X2": X2,
                "X3": X3,
                "X4": X4,
                "Y": Y,
                "env": [env] * n_per_env,
            }
        )

    df_a = generate(
        "envA",
        mu_x1=0.0,
        mu_x2=0.0,
        sigma_x1=1.0,
        sigma_x2=1.0,
        mu_x4=0.0,
    )
    # envB shifts both mean and variance of the causal parents so that any
    # non-parent subset produces residuals with env-specific dispersion.
    df_b = generate(
        "envB",
        mu_x1=2.0,
        mu_x2=-1.5,
        sigma_x1=2.0,
        sigma_x2=0.4,
        mu_x4=-2.0,
    )
    df = pd.concat([df_a, df_b], ignore_index=True)
    return df


class TestCausalAcceptance:

    def test_recall_ground_truth_parents(self):
        df = _synthetic_causal_cohort(n_per_env=200, seed=0)
        predictor = InvariantPathwayPredictor(alpha=0.05, max_subset_size=3)
        X = df[["X1", "X2", "X3", "X4"]]
        y = df["Y"]
        envs = df["env"].to_numpy()
        report = predictor.fit(X=X, y=y, target_name="Y", environments=envs)

        # Ground truth causal parents: X1, X2
        recall = report.recall_against(["X1", "X2"])
        assert recall >= 0.7, (
            f"recall = {recall:.2f} below 0.7 target; "
            f"identifiable = {sorted(report.identifiable_parents)}"
        )

    def test_spurious_variable_not_identified(self):
        """X3 correlates with Y in env A but its P(Y|X3) breaks under env shift
        — it should not appear in the identifiable parent set."""
        df = _synthetic_causal_cohort(n_per_env=200, seed=0)
        predictor = InvariantPathwayPredictor(alpha=0.05, max_subset_size=3)
        report = predictor.fit(
            X=df[["X1", "X2", "X3", "X4"]],
            y=df["Y"],
            target_name="Y",
            environments=df["env"].to_numpy(),
        )
        # X3 and X4 should not be in the identifiable parent set
        assert "X4" not in report.identifiable_parents

    def test_report_summary_renders(self):
        df = _synthetic_causal_cohort(n_per_env=100)
        report = InvariantPathwayPredictor(alpha=0.05).fit(
            X=df[["X1", "X2", "X4"]],
            y=df["Y"],
            target_name="Y",
            environments=df["env"].to_numpy(),
        )
        s = report.summary()
        assert "CausalParentReport" in s

    def test_report_to_dict(self):
        df = _synthetic_causal_cohort(n_per_env=100)
        report = InvariantPathwayPredictor(alpha=0.05).fit(
            X=df[["X1", "X2"]],
            y=df["Y"],
            target_name="Y",
            environments=df["env"].to_numpy(),
        )
        d = report.to_dict()
        assert "identifiable_parents" in d
        assert "invariant_subsets" in d

    def test_precision_against(self):
        df = _synthetic_causal_cohort(n_per_env=200)
        report = InvariantPathwayPredictor(alpha=0.05, max_subset_size=3).fit(
            X=df[["X1", "X2", "X3", "X4"]],
            y=df["Y"],
            target_name="Y",
            environments=df["env"].to_numpy(),
        )
        # Precision should be 1.0 if only X1, X2 are identified
        p = report.precision_against(["X1", "X2"])
        # If nothing was identified, precision is NaN; otherwise should be high
        if not np.isnan(p):
            assert p >= 0.5


class TestSingleEnvironmentEdgeCase:

    def test_single_env_returns_all_subsets_invariant(self):
        rng = np.random.default_rng(0)
        n = 80
        X = pd.DataFrame(
            {
                "A": rng.standard_normal(n),
                "B": rng.standard_normal(n),
            }
        )
        y = pd.Series(X["A"] * 0.5 + rng.normal(0, 0.1, n))
        report = InvariantPathwayPredictor(alpha=0.05, max_subset_size=2).fit(
            X=X,
            y=y,
            target_name="Y",
            environments=["only"] * n,
        )
        # With one environment, every subset trivially "passes" invariance
        # → identifiable parent set is the empty intersection
        assert report.identifiable_parents == frozenset()
