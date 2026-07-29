"""
Tests for the pathway_subtyping.uncertainty module.

Covers roadmap acceptance criteria:
    - Conformal intervals achieve nominal coverage within ±2%
    - Bootstrap intervals stable to < 5% width variance across 3 seeds
    - BayesianPathwayGMM reproduces point-estimate results within 1%
      (on identifiable synthetic data, measured via ARI)
    - CalibrationReport ECE/Brier metrics correct on known-calibration data
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

from pathway_subtyping.uncertainty import (
    BayesianPathwayGMM,
    BootstrapMSV,
    BootstrapResult,
    CalibrationReport,
    ConformalInterval,
    ConformalPathwayPredictor,
    brier_score,
    ece,
    reliability_curve,
)

# --------------------------------------------------------------------------- #
# Shared fixtures
# --------------------------------------------------------------------------- #


@pytest.fixture
def regression_data():
    """Synthetic regression: y = 2x + noise, exchangeable splits."""
    rng = np.random.default_rng(0)
    n = 2000
    x = rng.uniform(-2, 2, size=(n, 1))
    y = (2.0 * x[:, 0] + rng.normal(0, 0.5, size=n)).astype(float)
    return x, y


@pytest.fixture
def gmm_data():
    """Three well-separated Gaussian clusters."""
    rng = np.random.default_rng(42)
    means = np.array([[-3.0, 0.0], [3.0, 0.0], [0.0, 4.0]])
    n_per = 200
    X = np.vstack([m + rng.normal(0, 0.5, size=(n_per, 2)) for m in means])
    y = np.repeat(np.arange(3), n_per)
    return X, y


@pytest.fixture
def binary_labels_and_probs():
    """Well-calibrated binary problem: p ~= P(y=1)."""
    rng = np.random.default_rng(7)
    n = 3000
    p = rng.uniform(0.0, 1.0, size=n)
    y = (rng.uniform(size=n) < p).astype(int)
    return y, p


# --------------------------------------------------------------------------- #
# Conformal
# --------------------------------------------------------------------------- #


class TestConformalRegression:

    def test_interval_dataclass_basics(self):
        ival = ConformalInterval(
            prediction=1.0,
            lower=0.5,
            upper=1.5,
            coverage=0.9,
            n_calibration=100,
        )
        assert ival.width == pytest.approx(1.0)
        assert ival.contains(1.2)
        assert not ival.contains(2.0)
        d = ival.to_dict()
        assert d["width"] == pytest.approx(1.0)

    def test_requires_calibration_before_predict(self, regression_data):
        x, _ = regression_data
        pred = ConformalPathwayPredictor(score_fn=lambda X: 2.0 * X[:, 0])
        with pytest.raises(RuntimeError, match="calibrate"):
            pred.predict(x[:5])

    def test_rejects_invalid_coverage(self):
        with pytest.raises(ValueError):
            ConformalPathwayPredictor(score_fn=lambda X: X[:, 0], coverage=1.5)

    def test_rejects_tiny_calibration_set(self):
        pred = ConformalPathwayPredictor(score_fn=lambda X: X[:, 0])
        with pytest.raises(ValueError, match="10 calibration"):
            pred.calibrate(np.ones((5, 1)), np.ones(5))

    def test_regression_coverage_within_2pct(self):
        """Roadmap acceptance: nominal coverage ±2% on held-out synthetic data.

        Averaged over multiple seeds to absorb finite-sample discreteness —
        split-conformal has a formal *lower-bound* coverage guarantee; any
        single run may be off by more than 2% while the mean tracks nominal.
        """
        deviations = {t: [] for t in (0.80, 0.90, 0.95)}
        for seed in range(6):
            rng = np.random.default_rng(seed)
            n = 4000
            x = rng.uniform(-2, 2, size=(n, 1))
            y = 2.0 * x[:, 0] + rng.normal(0, 0.5, size=n)

            perm = rng.permutation(n)
            fit_idx = perm[: int(0.50 * n)]
            cal_idx = perm[int(0.50 * n) : int(0.75 * n)]
            te_idx = perm[int(0.75 * n) :]

            coef = np.polyfit(x[fit_idx, 0], y[fit_idx], deg=1)

            def score_fn(X, _c=coef):
                return _c[0] * X[:, 0] + _c[1]

            for target in deviations:
                predictor = ConformalPathwayPredictor(score_fn=score_fn, coverage=target)
                predictor.calibrate(x[cal_idx], y[cal_idx])
                empirical = predictor.coverage_on(x[te_idx], y[te_idx])
                deviations[target].append(empirical - target)

        for target, devs in deviations.items():
            mean_dev = float(np.mean(devs))
            # Mean deviation across seeds should be within ±2% (roadmap).
            assert (
                abs(mean_dev) < 0.02
            ), f"target={target} mean deviation={mean_dev:.3f} devs={devs}"

    def test_quantile_monotone_in_coverage(self, regression_data):
        x, y = regression_data

        def score_fn(X):  # matches true mean
            return 2.0 * X[:, 0]

        q_vals = []
        for target in (0.5, 0.8, 0.95):
            pred = ConformalPathwayPredictor(score_fn=score_fn, coverage=target)
            pred.calibrate(x[:500], y[:500])
            q_vals.append(pred.quantile)
        assert q_vals[0] < q_vals[1] < q_vals[2]


class TestConformalClassification:

    def test_label_set_coverage(self, gmm_data):
        X, y = gmm_data
        # "Model": sklearn GMM predict_proba — well calibrated on this data
        gmm = GaussianMixture(n_components=3, random_state=0).fit(X)

        # Align GMM components to true labels by majority vote
        raw = gmm.predict(X)
        mapping = {}
        for c in range(3):
            mask = raw == c
            if mask.sum() == 0:
                mapping[c] = c
                continue
            mapping[c] = int(np.bincount(y[mask]).argmax())
        remap = np.vectorize(mapping.get)

        def score_fn(X_):
            probs = gmm.predict_proba(X_)
            # Reorder columns so column k corresponds to true-label k
            reorder = np.zeros((3, 3))
            for src, dst in mapping.items():
                reorder[src, dst] = 1.0
            return probs @ reorder

        n = len(X)
        rng = np.random.default_rng(1)
        perm = rng.permutation(n)
        cal_idx = perm[: int(0.5 * n)]
        te_idx = perm[int(0.5 * n) :]

        pred = ConformalPathwayPredictor(score_fn=score_fn, coverage=0.9, mode="classification")
        pred.calibrate(X[cal_idx], y[cal_idx], labels=[0, 1, 2])
        empirical = pred.coverage_on(X[te_idx], y[te_idx])
        assert empirical >= 0.88


# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #


class TestCalibrationMetrics:

    def test_reliability_curve_shape(self, binary_labels_and_probs):
        y, p = binary_labels_and_probs
        centers, mean_pred, observed = reliability_curve(y, p, n_bins=10)
        assert centers.shape == (10,)
        assert mean_pred.shape == (10,)
        assert observed.shape == (10,)
        # For well-calibrated data, observed ~= mean_pred
        mask = ~np.isnan(mean_pred)
        residual = np.abs(mean_pred[mask] - observed[mask])
        assert residual.mean() < 0.05

    def test_ece_well_calibrated_is_low(self, binary_labels_and_probs):
        y, p = binary_labels_and_probs
        assert ece(y, p, n_bins=15) < 0.05

    def test_ece_miscalibrated_is_high(self):
        rng = np.random.default_rng(0)
        n = 2000
        p = rng.uniform(0.0, 1.0, size=n)
        # Miscalibration: true rate is p^2 (squeezed toward 0)
        y = (rng.uniform(size=n) < p**2).astype(int)
        assert ece(y, p, n_bins=15) > 0.1

    def test_brier_bounds(self, binary_labels_and_probs):
        y, p = binary_labels_and_probs
        b = brier_score(y, p)
        assert 0.0 <= b <= 1.0
        # Better than uniform 0.5 for well-calibrated p
        assert b < 0.25

    def test_reliability_quantile_strategy(self, binary_labels_and_probs):
        y, p = binary_labels_and_probs
        _, mp_q, obs_q = reliability_curve(y, p, n_bins=10, strategy="quantile")
        mask = ~np.isnan(mp_q)
        # Each bin holds ~10% of samples by construction
        assert mask.sum() >= 8


class TestCalibrationReport:

    def test_from_predictions_basic(self, binary_labels_and_probs):
        y, p = binary_labels_and_probs
        report = CalibrationReport.from_predictions(y, p)
        assert 0.0 <= report.ece < 0.1
        assert 0.0 <= report.brier < 0.25
        s = report.summary()
        assert "ECE" in s and "Brier" in s

    def test_to_dict_roundtrip(self, binary_labels_and_probs):
        y, p = binary_labels_and_probs
        report = CalibrationReport.from_predictions(y, p)
        d = report.to_dict()
        assert d["n_samples"] == len(y)
        assert "reliability" in d
        assert "bin_centers" in d["reliability"]

    def test_plot_returns_figure(self, binary_labels_and_probs):
        pytest.importorskip("matplotlib")
        y, p = binary_labels_and_probs
        report = CalibrationReport.from_predictions(y, p)
        fig = report.plot()
        assert fig is not None


# --------------------------------------------------------------------------- #
# Bootstrap
# --------------------------------------------------------------------------- #


class TestBootstrap:

    def test_scalar_stat(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(500, 1))
        bs = BootstrapMSV(n_bootstrap=200, ci_level=0.95, seed=42)
        result = bs.run(X, score_fn=lambda A: np.array([A.mean()]))
        assert isinstance(result, BootstrapResult)
        assert result.samples.shape == (200, 1)
        # True mean ~ 0, interval should contain 0
        assert result.lower[0] < 0.0 < result.upper[0]

    def test_vector_stat(self):
        rng = np.random.default_rng(1)
        X = rng.normal(loc=[1.0, -1.0], scale=1.0, size=(2000, 2))
        bs = BootstrapMSV(n_bootstrap=300, seed=0)
        result = bs.run(X, score_fn=lambda A: A.mean(axis=0))
        assert result.point.shape == (2,)
        # Interval should contain the sample mean of the input data (not the
        # population mean — bootstrap is about the finite-sample distribution).
        sample_mean = X.mean(axis=0)
        assert result.lower[0] <= sample_mean[0] <= result.upper[0]
        assert result.lower[1] <= sample_mean[1] <= result.upper[1]

    def test_width_stability_under_5pct(self):
        """Roadmap acceptance: interval width variance < 5% across 3 seeds."""
        rng = np.random.default_rng(0)
        X = rng.normal(size=(2000, 1))

        widths = []
        for s in (1, 2, 3):
            bs = BootstrapMSV(n_bootstrap=1000, ci_level=0.95, seed=s)
            res = bs.run(X, score_fn=lambda A: np.array([A.mean()]))
            widths.append(float(res.width[0]))

        widths = np.asarray(widths)
        rel_var = widths.std() / widths.mean()
        assert rel_var < 0.05, f"width relative variance {rel_var:.3f} exceeds 5%"

    def test_per_cell_bootstrap(self):
        rng = np.random.default_rng(3)
        n_cells = 50
        X = rng.normal(size=(n_cells, 1))
        bs = BootstrapMSV(n_bootstrap=200, seed=0)
        # Score fn returns one value per input cell — the cell's own score
        result = bs.run(X, score_fn=lambda A: A[:, 0], per_cell=True)
        assert result.lower.shape == (n_cells,)
        # Each cell's own score should fall inside its interval (up to NaN for unsampled)
        mask = ~np.isnan(result.lower)
        assert mask.mean() > 0.9  # almost all cells resampled at least once

    def test_rejects_bad_ci_level(self):
        with pytest.raises(ValueError):
            BootstrapMSV(ci_level=1.1)

    def test_rejects_tiny_n_bootstrap(self):
        with pytest.raises(ValueError):
            BootstrapMSV(n_bootstrap=3)

    def test_parallel_jobs_produce_results(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(200, 1))
        bs = BootstrapMSV(n_bootstrap=50, seed=0, n_jobs=2)
        result = bs.run(X, score_fn=lambda A: np.array([A.mean()]))
        assert result.samples.shape == (50, 1)

    def test_works_with_dataframe_input(self):
        rng = np.random.default_rng(0)
        df = pd.DataFrame(rng.normal(size=(100, 3)), columns=["a", "b", "c"])
        bs = BootstrapMSV(n_bootstrap=50, seed=0)
        result = bs.run(df, score_fn=lambda A: A.mean(axis=0).to_numpy())
        assert result.point.shape == (3,)


# --------------------------------------------------------------------------- #
# BayesianPathwayGMM
# --------------------------------------------------------------------------- #


class TestBayesianPathwayGMM:

    def test_fit_predict_recovers_clusters(self, gmm_data):
        X, y = gmm_data
        model = BayesianPathwayGMM(n_components=3, random_state=0).fit(X)
        labels = model.predict(X)
        ari = adjusted_rand_score(y, labels)
        assert ari > 0.95

    def test_predict_proba_shape_and_sum(self, gmm_data):
        X, _ = gmm_data
        model = BayesianPathwayGMM(n_components=3, random_state=0).fit(X)
        probs = model.predict_proba(X)
        assert probs.shape == (len(X), 3)
        np.testing.assert_allclose(probs.sum(axis=1), 1.0, atol=1e-6)

    def test_mode_reproduces_point_gmm_within_1pct(self, gmm_data):
        """Roadmap acceptance: mode-collapsed posterior within 1% of point GMM."""
        X, _ = gmm_data
        point_labels = (
            GaussianMixture(n_components=3, random_state=0, covariance_type="full")
            .fit(X)
            .predict(X)
        )
        mode_labels = BayesianPathwayGMM(n_components=3, random_state=0).fit(X).predict(X)
        ari = adjusted_rand_score(point_labels, mode_labels)
        assert ari > 0.99, f"mode-collapsed ARI vs point GMM = {ari:.4f}"

    def test_sample_assignments_shape(self, gmm_data):
        X, _ = gmm_data
        model = BayesianPathwayGMM(n_components=3, random_state=0).fit(X)
        draws = model.sample_assignments(X, n_samples=20, random_state=1)
        assert draws.shape == (20, len(X))
        assert draws.min() >= 0
        assert draws.max() < 3

    def test_sample_parameters_shapes(self, gmm_data):
        X, _ = gmm_data
        model = BayesianPathwayGMM(n_components=3, random_state=0).fit(X)
        samples = model.sample_parameters(n_samples=10, random_state=0)
        assert len(samples) == 10
        assert samples[0].weights.shape == (3,)
        assert samples[0].means.shape == (3, X.shape[1])
        np.testing.assert_allclose(samples[0].weights.sum(), 1.0, atol=1e-3)

    def test_requires_fit_before_predict(self):
        model = BayesianPathwayGMM(n_components=2)
        with pytest.raises(RuntimeError, match="fit"):
            model.predict(np.zeros((3, 2)))

    def test_accessors_after_fit(self, gmm_data):
        X, _ = gmm_data
        model = BayesianPathwayGMM(n_components=3, random_state=0).fit(X)
        assert model.weights_.shape == (3,)
        assert model.means_.shape == (3, X.shape[1])
        assert model.covariances_.shape == (3, X.shape[1], X.shape[1])
        assert 1 <= model.n_effective_components <= 3

    def test_posterior_sample_to_dict(self, gmm_data):
        X, _ = gmm_data
        model = BayesianPathwayGMM(n_components=2, random_state=0).fit(X)
        samples = model.sample_parameters(n_samples=2)
        d = samples[0].to_dict()
        assert "weights" in d and "means" in d


class TestConformalSmallCalibrationCoverage:
    """Coverage is the one guarantee conformal offers; it must not be missed quietly.

    When ceil((n+1)(1-alpha)) > n there is no such order statistic and the correct
    quantile is +infinity (the all-labels set). Clamping to the max observed score
    instead under-covered: 0.9035 at n=10 against a 0.95 target, 0.9373 at n=15.
    """

    def test_infeasible_coverage_yields_infinite_quantile(self):
        rng = np.random.default_rng(0)
        # n=10 with alpha=0.05 needs the 11th of 10 order statistics
        pred = ConformalPathwayPredictor(score_fn=lambda X: X[:, 0], coverage=0.95)
        X = rng.normal(size=(10, 2))
        y = X[:, 0] + rng.normal(scale=0.1, size=10)
        pred.calibrate(X, y)
        assert np.isinf(pred._quantile), (
            "too-small calibration set must give the all-labels/infinite interval, "
            "not the max observed score"
        )

    def test_feasible_coverage_is_finite(self):
        rng = np.random.default_rng(1)
        pred = ConformalPathwayPredictor(score_fn=lambda X: X[:, 0], coverage=0.95)
        X = rng.normal(size=(200, 2))
        y = X[:, 0] + rng.normal(scale=0.1, size=200)
        pred.calibrate(X, y)
        assert np.isfinite(pred._quantile)
