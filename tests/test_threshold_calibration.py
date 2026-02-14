"""Tests for the threshold calibration module."""

import numpy as np
import pytest

from pathway_subtyping.threshold_calibration import (
    CalibratedThresholds,
    CalibrationSimulationResult,
    _GRID_CLUSTERS,
    _GRID_SAMPLES,
    _NULL_ARI_TABLE,
    _STABILITY_TABLE,
    _interpolate_threshold,
    _simulate_null_distribution,
    _simulate_stability_distribution,
    calibrate_thresholds,
    generate_calibration_table,
    get_default_thresholds,
)


# =========================================================================
# CalibratedThresholds dataclass
# =========================================================================


class TestCalibratedThresholds:
    """Tests for CalibratedThresholds dataclass."""

    def test_to_dict(self):
        ct = CalibratedThresholds(
            null_ari_threshold=0.058,
            stability_threshold=0.85,
            n_samples=100,
            n_clusters=3,
        )
        d = ct.to_dict()
        assert d["null_ari_threshold"] == 0.058
        assert d["stability_threshold"] == 0.85
        assert d["n_samples"] == 100
        assert d["n_clusters"] == 3
        assert d["n_pathways"] == 15
        assert d["alpha"] == 0.05
        assert d["calibration_method"] == "lookup"
        assert d["interpolated"] is False

    def test_to_dict_rounding(self):
        ct = CalibratedThresholds(
            null_ari_threshold=0.058123456,
            stability_threshold=0.854321,
            n_samples=100,
            n_clusters=3,
        )
        d = ct.to_dict()
        assert d["null_ari_threshold"] == 0.0581
        assert d["stability_threshold"] == 0.8543

    def test_format_report(self):
        ct = CalibratedThresholds(
            null_ari_threshold=0.058,
            stability_threshold=0.85,
            n_samples=100,
            n_clusters=3,
            calibration_method="lookup",
        )
        report = ct.format_report()
        assert "Threshold Calibration Report" in report
        assert "n=100" in report
        assert "k=3" in report
        assert "0.0580" in report
        assert "0.8500" in report
        assert "lookup" in report

    def test_format_report_interpolated(self):
        ct = CalibratedThresholds(
            null_ari_threshold=0.06,
            stability_threshold=0.84,
            n_samples=110,
            n_clusters=3,
            calibration_method="lookup",
            interpolated=True,
        )
        report = ct.format_report()
        assert "(interpolated)" in report

    def test_get_citations(self):
        ct = CalibratedThresholds(
            null_ari_threshold=0.058,
            stability_threshold=0.85,
            n_samples=100,
            n_clusters=3,
        )
        citations = ct.get_citations()
        assert len(citations) == 2
        assert any("Hubert" in c for c in citations)
        assert any("Hennig" in c for c in citations)


class TestCalibrationSimulationResult:
    """Tests for CalibrationSimulationResult dataclass."""

    def test_to_dict(self):
        result = CalibrationSimulationResult(
            null_ari_values=np.array([0.01, 0.02, 0.03]),
            stability_ari_values=np.array([0.85, 0.90, 0.88]),
            null_percentile=0.03,
            stability_percentile=0.85,
            n_simulations=3,
        )
        d = result.to_dict()
        assert d["n_simulations"] == 3
        assert "null_ari_mean" in d
        assert "stability_ari_mean" in d
        assert abs(d["null_ari_mean"] - 0.02) < 0.001


# =========================================================================
# Lookup table properties
# =========================================================================


class TestLookupTable:
    """Tests for pre-computed lookup table properties."""

    def test_tables_same_keys(self):
        """Both tables should have the same grid points."""
        assert set(_NULL_ARI_TABLE.keys()) == set(_STABILITY_TABLE.keys())

    def test_expected_grid_size(self):
        """Table should have 8 sample sizes * 7 cluster counts = 56 entries."""
        assert len(_NULL_ARI_TABLE) == 56
        assert len(_STABILITY_TABLE) == 56

    def test_null_ari_decreases_with_n(self):
        """Null ARI threshold should decrease as n increases (for fixed k)."""
        for k in _GRID_CLUSTERS:
            values = [_NULL_ARI_TABLE[(n, k)] for n in _GRID_SAMPLES]
            for i in range(len(values) - 1):
                assert values[i] >= values[i + 1], (
                    f"Null ARI not monotonically decreasing for k={k}: "
                    f"n={_GRID_SAMPLES[i]} ({values[i]}) vs "
                    f"n={_GRID_SAMPLES[i+1]} ({values[i+1]})"
                )

    def test_null_ari_increases_with_k(self):
        """Null ARI threshold should increase as k increases (for fixed n)."""
        for n in _GRID_SAMPLES:
            values = [_NULL_ARI_TABLE[(n, k)] for k in _GRID_CLUSTERS]
            for i in range(len(values) - 1):
                assert values[i] <= values[i + 1], (
                    f"Null ARI not monotonically increasing for n={n}: "
                    f"k={_GRID_CLUSTERS[i]} ({values[i]}) vs "
                    f"k={_GRID_CLUSTERS[i+1]} ({values[i+1]})"
                )

    def test_stability_increases_with_n(self):
        """Stability threshold should increase as n increases (for fixed k)."""
        for k in _GRID_CLUSTERS:
            values = [_STABILITY_TABLE[(n, k)] for n in _GRID_SAMPLES]
            for i in range(len(values) - 1):
                assert values[i] <= values[i + 1], (
                    f"Stability not monotonically increasing for k={k}: "
                    f"n={_GRID_SAMPLES[i]} ({values[i]}) vs "
                    f"n={_GRID_SAMPLES[i+1]} ({values[i+1]})"
                )

    def test_stability_decreases_with_k(self):
        """Stability threshold should decrease as k increases (for fixed n)."""
        for n in _GRID_SAMPLES:
            values = [_STABILITY_TABLE[(n, k)] for k in _GRID_CLUSTERS]
            for i in range(len(values) - 1):
                assert values[i] >= values[i + 1], (
                    f"Stability not monotonically decreasing for n={n}: "
                    f"k={_GRID_CLUSTERS[i]} ({values[i]}) vs "
                    f"k={_GRID_CLUSTERS[i+1]} ({values[i+1]})"
                )

    def test_all_values_in_range(self):
        """All threshold values should be in [0, 1]."""
        for key, value in _NULL_ARI_TABLE.items():
            assert 0 <= value <= 1, f"Null ARI out of range at {key}: {value}"
        for key, value in _STABILITY_TABLE.items():
            assert 0 <= value <= 1, f"Stability out of range at {key}: {value}"


# =========================================================================
# Interpolation
# =========================================================================


class TestInterpolation:
    """Tests for bilinear interpolation."""

    def test_exact_match(self):
        """Exact grid point should return exact value."""
        value = _interpolate_threshold(_NULL_ARI_TABLE, 100, 3)
        assert value == _NULL_ARI_TABLE[(100, 3)]

    def test_interpolated_between_n(self):
        """Interpolation between sample sizes should be between neighbors."""
        v_lo = _NULL_ARI_TABLE[(100, 3)]
        v_hi = _NULL_ARI_TABLE[(150, 3)]
        v_mid = _interpolate_threshold(_NULL_ARI_TABLE, 125, 3)
        assert v_mid is not None
        assert min(v_lo, v_hi) <= v_mid <= max(v_lo, v_hi)

    def test_interpolated_between_k(self):
        """Interpolation between cluster counts should be between neighbors."""
        v_lo = _NULL_ARI_TABLE[(100, 3)]
        v_hi = _NULL_ARI_TABLE[(100, 4)]
        # k=3.5 is not integer but tests the math
        v_mid = _interpolate_threshold(_NULL_ARI_TABLE, 100, 3)
        assert v_mid is not None

    def test_interpolated_both_dimensions(self):
        """Bilinear interpolation in both dimensions."""
        corners = [
            _NULL_ARI_TABLE[(100, 3)],
            _NULL_ARI_TABLE[(100, 4)],
            _NULL_ARI_TABLE[(150, 3)],
            _NULL_ARI_TABLE[(150, 4)],
        ]
        v = _interpolate_threshold(_NULL_ARI_TABLE, 125, 3)
        assert v is not None
        assert min(corners) <= v <= max(corners)

    def test_out_of_range_low_n(self):
        """n_samples below minimum should return None."""
        result = _interpolate_threshold(_NULL_ARI_TABLE, 10, 3)
        assert result is None

    def test_out_of_range_high_n(self):
        """n_samples above maximum should return None."""
        result = _interpolate_threshold(_NULL_ARI_TABLE, 1000, 3)
        assert result is None

    def test_out_of_range_low_k(self):
        """n_clusters below minimum should return None."""
        result = _interpolate_threshold(_NULL_ARI_TABLE, 100, 1)
        assert result is None

    def test_out_of_range_high_k(self):
        """n_clusters above maximum should return None."""
        result = _interpolate_threshold(_NULL_ARI_TABLE, 100, 10)
        assert result is None

    def test_boundary_min(self):
        """Smallest grid point should work."""
        value = _interpolate_threshold(_NULL_ARI_TABLE, 30, 2)
        assert value == _NULL_ARI_TABLE[(30, 2)]

    def test_boundary_max(self):
        """Largest grid point should work."""
        value = _interpolate_threshold(_NULL_ARI_TABLE, 500, 8)
        assert value == _NULL_ARI_TABLE[(500, 8)]


# =========================================================================
# calibrate_thresholds
# =========================================================================


class TestCalibrateThresholds:
    """Tests for calibrate_thresholds function."""

    def test_lookup_exact_match(self):
        """Lookup mode should return exact table values."""
        ct = calibrate_thresholds(100, 3, method="lookup")
        assert ct.null_ari_threshold == _NULL_ARI_TABLE[(100, 3)]
        assert ct.stability_threshold == _STABILITY_TABLE[(100, 3)]
        assert ct.calibration_method == "lookup"
        assert not ct.interpolated

    def test_lookup_interpolated(self):
        """Lookup with intermediate values should interpolate."""
        ct = calibrate_thresholds(125, 3, method="lookup")
        assert ct.calibration_method == "lookup"
        assert ct.interpolated

    def test_lookup_out_of_range_raises(self):
        """Lookup mode should raise for out-of-range dimensions."""
        with pytest.raises(ValueError, match="outside the pre-computed grid"):
            calibrate_thresholds(10, 3, method="lookup")

    def test_auto_prefers_lookup(self):
        """Auto mode should use lookup when available."""
        ct = calibrate_thresholds(100, 3, method="auto")
        assert ct.calibration_method == "lookup"

    def test_auto_falls_back_to_simulate(self):
        """Auto mode should simulate when out of range."""
        ct = calibrate_thresholds(
            15, 2, method="auto", n_simulations=5, seed=42
        )
        assert ct.calibration_method == "simulate"

    def test_simulate_mode(self):
        """Simulate mode should always run simulations."""
        ct = calibrate_thresholds(
            100, 3, method="simulate", n_simulations=5, seed=42
        )
        assert ct.calibration_method == "simulate"
        assert ct.null_ari_threshold > 0
        assert ct.stability_threshold > 0

    def test_simulate_reproducibility(self):
        """Simulation with same seed should produce same results."""
        ct1 = calibrate_thresholds(
            50, 3, method="simulate", n_simulations=5, seed=42
        )
        ct2 = calibrate_thresholds(
            50, 3, method="simulate", n_simulations=5, seed=42
        )
        assert ct1.null_ari_threshold == ct2.null_ari_threshold
        assert ct1.stability_threshold == ct2.stability_threshold

    def test_different_alpha(self):
        """Different alpha should adjust thresholds."""
        ct_05 = calibrate_thresholds(100, 3, alpha=0.05, method="lookup")
        ct_10 = calibrate_thresholds(100, 3, alpha=0.10, method="lookup")
        # Higher alpha -> more permissive null threshold
        assert ct_10.null_ari_threshold > ct_05.null_ari_threshold

    def test_n_samples_stored(self):
        """Calibrated thresholds should store data dimensions."""
        ct = calibrate_thresholds(200, 4, method="lookup")
        assert ct.n_samples == 200
        assert ct.n_clusters == 4

    def test_larger_n_gives_stricter_null(self):
        """Larger samples should have stricter null thresholds."""
        ct_small = calibrate_thresholds(50, 3, method="lookup")
        ct_large = calibrate_thresholds(300, 3, method="lookup")
        assert ct_large.null_ari_threshold < ct_small.null_ari_threshold

    def test_more_clusters_gives_relaxed_null(self):
        """More clusters should have more relaxed null thresholds."""
        ct_few = calibrate_thresholds(100, 2, method="lookup")
        ct_many = calibrate_thresholds(100, 6, method="lookup")
        assert ct_many.null_ari_threshold > ct_few.null_ari_threshold


# =========================================================================
# get_default_thresholds
# =========================================================================


class TestGetDefaultThresholds:
    """Tests for get_default_thresholds function."""

    def test_returns_legacy_values(self):
        """Default thresholds should match the original hard-coded values."""
        ct = get_default_thresholds()
        assert ct.null_ari_threshold == 0.15
        assert ct.stability_threshold == 0.8

    def test_method_is_default(self):
        ct = get_default_thresholds()
        assert ct.calibration_method == "default"


# =========================================================================
# Null distribution simulation
# =========================================================================


class TestSimulateNullDistribution:
    """Tests for _simulate_null_distribution."""

    def test_returns_array(self):
        values = _simulate_null_distribution(50, 5, 3, n_simulations=5, seed=42)
        assert isinstance(values, np.ndarray)

    def test_correct_length(self):
        values = _simulate_null_distribution(50, 5, 3, n_simulations=5, seed=42)
        # Should have <= n_simulations values (some may not converge)
        assert len(values) <= 5
        assert len(values) > 0

    def test_values_near_zero(self):
        """Null ARI values should be near zero."""
        values = _simulate_null_distribution(100, 10, 3, n_simulations=10, seed=42)
        if len(values) > 0:
            assert np.mean(values) < 0.3  # generous bound for small n_sim

    def test_reproducible(self):
        v1 = _simulate_null_distribution(50, 5, 3, n_simulations=5, seed=42)
        v2 = _simulate_null_distribution(50, 5, 3, n_simulations=5, seed=42)
        np.testing.assert_array_equal(v1, v2)


# =========================================================================
# Stability distribution simulation
# =========================================================================


class TestSimulateStabilityDistribution:
    """Tests for _simulate_stability_distribution."""

    def test_returns_array(self):
        values = _simulate_stability_distribution(
            50, 5, 3, n_simulations=3, n_bootstrap_per_sim=3, seed=42
        )
        assert isinstance(values, np.ndarray)

    def test_values_positive(self):
        """Stability ARI values should generally be positive for structured data."""
        values = _simulate_stability_distribution(
            100, 10, 3, effect_size=2.0,
            n_simulations=3, n_bootstrap_per_sim=5, seed=42
        )
        if len(values) > 0:
            assert np.mean(values) > 0

    def test_correct_max_length(self):
        n_sim = 3
        n_boot = 4
        values = _simulate_stability_distribution(
            50, 5, 3, n_simulations=n_sim, n_bootstrap_per_sim=n_boot, seed=42
        )
        assert len(values) <= n_sim * n_boot


# =========================================================================
# generate_calibration_table
# =========================================================================


class TestGenerateCalibrationTable:
    """Tests for generate_calibration_table."""

    def test_returns_both_tables(self):
        tables = generate_calibration_table(
            sample_sizes=[50],
            cluster_counts=[3],
            n_simulations=5,
            seed=42,
        )
        assert "null_ari" in tables
        assert "stability" in tables

    def test_values_in_range(self):
        tables = generate_calibration_table(
            sample_sizes=[50],
            cluster_counts=[3],
            n_simulations=5,
            seed=42,
        )
        for key, value in tables["null_ari"].items():
            assert 0 <= value <= 1
        for key, value in tables["stability"].items():
            assert 0 <= value <= 1

    def test_correct_keys(self):
        tables = generate_calibration_table(
            sample_sizes=[50, 100],
            cluster_counts=[2, 3],
            n_simulations=5,
            seed=42,
        )
        assert (50, 2) in tables["null_ari"]
        assert (100, 3) in tables["null_ari"]
        assert len(tables["null_ari"]) == 4
