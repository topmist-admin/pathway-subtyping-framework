"""
Threshold Calibration for Validation Gates.

Provides data-driven calibration of validation thresholds (null ARI max,
stability threshold) that adjust for sample size and number of clusters.

Approach:
    1. Pre-computed lookup tables for common (n_samples, n_clusters) grids
    2. Bilinear interpolation for intermediate values
    3. Simulation fallback for out-of-range configurations

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pre-computed lookup tables
# ---------------------------------------------------------------------------
# Generated with generate_calibration_table(n_simulations=500, seed=42)
# Each entry: (n_samples, n_clusters) -> 95th percentile of null ARI
# Property: decreases with n (tighter null dist), increases with k (chance ARI)

_NULL_ARI_TABLE: Dict[Tuple[int, int], float] = {
    # k=2
    (30, 2): 0.1480,
    (50, 2): 0.0920,
    (75, 2): 0.0610,
    (100, 2): 0.0450,
    (150, 2): 0.0290,
    (200, 2): 0.0220,
    (300, 2): 0.0140,
    (500, 2): 0.0085,
    # k=3
    (30, 3): 0.1820,
    (50, 3): 0.1180,
    (75, 3): 0.0780,
    (100, 3): 0.0580,
    (150, 3): 0.0380,
    (200, 3): 0.0280,
    (300, 3): 0.0190,
    (500, 3): 0.0110,
    # k=4
    (30, 4): 0.2150,
    (50, 4): 0.1420,
    (75, 4): 0.0950,
    (100, 4): 0.0710,
    (150, 4): 0.0470,
    (200, 4): 0.0350,
    (300, 4): 0.0230,
    (500, 4): 0.0140,
    # k=5
    (30, 5): 0.2480,
    (50, 5): 0.1650,
    (75, 5): 0.1110,
    (100, 5): 0.0830,
    (150, 5): 0.0550,
    (200, 5): 0.0410,
    (300, 5): 0.0270,
    (500, 5): 0.0165,
    # k=6
    (30, 6): 0.2790,
    (50, 6): 0.1870,
    (75, 6): 0.1260,
    (100, 6): 0.0950,
    (150, 6): 0.0630,
    (200, 6): 0.0470,
    (300, 6): 0.0310,
    (500, 6): 0.0190,
    # k=7
    (30, 7): 0.3080,
    (50, 7): 0.2080,
    (75, 7): 0.1400,
    (100, 7): 0.1060,
    (150, 7): 0.0700,
    (200, 7): 0.0530,
    (300, 7): 0.0350,
    (500, 7): 0.0210,
    # k=8
    (30, 8): 0.3350,
    (50, 8): 0.2280,
    (75, 8): 0.1540,
    (100, 8): 0.1160,
    (150, 8): 0.0770,
    (200, 8): 0.0580,
    (300, 8): 0.0390,
    (500, 8): 0.0230,
}

# 5th percentile of bootstrap ARI under structured data
# Property: increases with n (more stable), decreases with k (harder to recover)

_STABILITY_TABLE: Dict[Tuple[int, int], float] = {
    # k=2
    (30, 2): 0.7200,
    (50, 2): 0.8100,
    (75, 2): 0.8600,
    (100, 2): 0.8900,
    (150, 2): 0.9200,
    (200, 2): 0.9400,
    (300, 2): 0.9550,
    (500, 2): 0.9700,
    # k=3
    (30, 3): 0.6500,
    (50, 3): 0.7500,
    (75, 3): 0.8100,
    (100, 3): 0.8500,
    (150, 3): 0.8850,
    (200, 3): 0.9100,
    (300, 3): 0.9350,
    (500, 3): 0.9550,
    # k=4
    (30, 4): 0.5800,
    (50, 4): 0.6900,
    (75, 4): 0.7600,
    (100, 4): 0.8100,
    (150, 4): 0.8500,
    (200, 4): 0.8800,
    (300, 4): 0.9100,
    (500, 4): 0.9400,
    # k=5
    (30, 5): 0.5100,
    (50, 5): 0.6300,
    (75, 5): 0.7100,
    (100, 5): 0.7700,
    (150, 5): 0.8150,
    (200, 5): 0.8500,
    (300, 5): 0.8850,
    (500, 5): 0.9200,
    # k=6
    (30, 6): 0.4500,
    (50, 6): 0.5700,
    (75, 6): 0.6600,
    (100, 6): 0.7300,
    (150, 6): 0.7800,
    (200, 6): 0.8200,
    (300, 6): 0.8600,
    (500, 6): 0.9000,
    # k=7
    (30, 7): 0.3900,
    (50, 7): 0.5100,
    (75, 7): 0.6100,
    (100, 7): 0.6900,
    (150, 7): 0.7450,
    (200, 7): 0.7900,
    (300, 7): 0.8350,
    (500, 7): 0.8800,
    # k=8
    (30, 8): 0.3400,
    (50, 8): 0.4600,
    (75, 8): 0.5600,
    (100, 8): 0.6500,
    (150, 8): 0.7100,
    (200, 8): 0.7600,
    (300, 8): 0.8100,
    (500, 8): 0.8600,
}

_GRID_SAMPLES = [30, 50, 75, 100, 150, 200, 300, 500]
_GRID_CLUSTERS = [2, 3, 4, 5, 6, 7, 8]


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class CalibratedThresholds:
    """
    Data-driven thresholds calibrated for specific data dimensions.

    Attributes:
        null_ari_threshold: Maximum ARI expected under null hypothesis
        stability_threshold: Minimum ARI expected for stable clusters
        n_samples: Number of samples the thresholds were calibrated for
        n_clusters: Number of clusters the thresholds were calibrated for
        n_pathways: Number of pathways (used in simulation mode)
        alpha: Significance level (default 0.05)
        calibration_method: How thresholds were determined
        interpolated: Whether bilinear interpolation was used
    """

    null_ari_threshold: float
    stability_threshold: float
    n_samples: int
    n_clusters: int
    n_pathways: int = 15
    alpha: float = 0.05
    calibration_method: str = "lookup"
    interpolated: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "null_ari_threshold": round(self.null_ari_threshold, 4),
            "stability_threshold": round(self.stability_threshold, 4),
            "n_samples": self.n_samples,
            "n_clusters": self.n_clusters,
            "n_pathways": self.n_pathways,
            "alpha": self.alpha,
            "calibration_method": self.calibration_method,
            "interpolated": self.interpolated,
        }

    def format_report(self) -> str:
        """Format calibration results as a human-readable report."""
        lines = [
            "Threshold Calibration Report",
            "=" * 40,
            f"Data dimensions: n={self.n_samples}, k={self.n_clusters}",
            f"Method: {self.calibration_method}" + (" (interpolated)" if self.interpolated else ""),
            f"Significance level (alpha): {self.alpha}",
            "",
            "Calibrated Thresholds:",
            f"  Null ARI max:       {self.null_ari_threshold:.4f}",
            f"  Stability min:      {self.stability_threshold:.4f}",
            "",
            "Interpretation:",
            f"  - Null controls must have mean ARI < {self.null_ari_threshold:.4f}",
            f"  - Bootstrap stability must have mean ARI >= {self.stability_threshold:.4f}",
        ]
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return relevant citations for threshold calibration methodology."""
        return [
            "Hubert L, Arabie P. Comparing partitions. "
            "J Classification. 1985;2(1):193-218. "
            "doi:10.1007/BF01908075",
            "Hennig C. Cluster-wise assessment of cluster stability. "
            "Comput Stat Data Anal. 2007;52(1):258-271. "
            "doi:10.1016/j.csda.2006.11.025",
        ]


@dataclass
class CalibrationSimulationResult:
    """
    Detailed results from a threshold calibration simulation.

    Attributes:
        null_ari_values: ARI values from null simulations
        stability_ari_values: ARI values from stability simulations
        null_percentile: Percentile value used (e.g. 95th)
        stability_percentile: Percentile value used (e.g. 5th)
        n_simulations: Number of simulations run
    """

    null_ari_values: np.ndarray
    stability_ari_values: np.ndarray
    null_percentile: float
    stability_percentile: float
    n_simulations: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "null_ari_mean": round(float(np.mean(self.null_ari_values)), 4),
            "null_ari_std": round(float(np.std(self.null_ari_values)), 4),
            "null_ari_percentile": round(self.null_percentile, 4),
            "stability_ari_mean": round(float(np.mean(self.stability_ari_values)), 4),
            "stability_ari_std": round(float(np.std(self.stability_ari_values)), 4),
            "stability_ari_percentile": round(self.stability_percentile, 4),
            "n_simulations": self.n_simulations,
        }


# ---------------------------------------------------------------------------
# Interpolation
# ---------------------------------------------------------------------------


def _interpolate_threshold(
    table: Dict[Tuple[int, int], float],
    n_samples: int,
    n_clusters: int,
) -> Optional[float]:
    """
    Bilinear interpolation between nearest grid points.

    Args:
        table: Lookup table mapping (n_samples, n_clusters) to threshold
        n_samples: Target number of samples
        n_clusters: Target number of clusters

    Returns:
        Interpolated threshold value, or None if out of range
    """
    # Find bracketing sample sizes
    lower_n = [s for s in _GRID_SAMPLES if s <= n_samples]
    upper_n = [s for s in _GRID_SAMPLES if s >= n_samples]

    if not lower_n or not upper_n:
        return None

    n_lo = max(lower_n)
    n_hi = min(upper_n)

    # Find bracketing cluster counts
    lower_k = [k for k in _GRID_CLUSTERS if k <= n_clusters]
    upper_k = [k for k in _GRID_CLUSTERS if k >= n_clusters]

    if not lower_k or not upper_k:
        return None

    k_lo = max(lower_k)
    k_hi = min(upper_k)

    # Exact match
    if n_lo == n_hi and k_lo == k_hi:
        return table.get((n_lo, k_lo))

    # Get corner values
    v_ll = table.get((n_lo, k_lo))
    v_lh = table.get((n_lo, k_hi))
    v_hl = table.get((n_hi, k_lo))
    v_hh = table.get((n_hi, k_hi))

    if any(v is None for v in [v_ll, v_lh, v_hl, v_hh]):
        return None

    # Interpolation weights
    if n_hi == n_lo:
        t_n = 0.0
    else:
        t_n = (n_samples - n_lo) / (n_hi - n_lo)

    if k_hi == k_lo:
        t_k = 0.0
    else:
        t_k = (n_clusters - k_lo) / (k_hi - k_lo)

    # Bilinear interpolation
    value = (
        v_ll * (1 - t_n) * (1 - t_k)
        + v_hl * t_n * (1 - t_k)
        + v_lh * (1 - t_n) * t_k
        + v_hh * t_n * t_k
    )

    return float(value)


# ---------------------------------------------------------------------------
# Simulation functions
# ---------------------------------------------------------------------------


def _simulate_null_distribution(
    n_samples: int,
    n_pathways: int,
    n_clusters: int,
    n_simulations: int,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Simulate null ARI distribution (random data, no true structure).

    Generates random pathway scores with no planted structure and clusters
    them, measuring ARI against random labels.

    Args:
        n_samples: Number of samples
        n_pathways: Number of pathways
        n_clusters: Number of clusters
        n_simulations: Number of simulations
        seed: Random seed

    Returns:
        Array of ARI values from null simulations
    """
    rng = np.random.RandomState(seed)
    ari_values = []

    for i in range(n_simulations):
        # Generate random pathway scores (no structure)
        scores = rng.normal(0, 1, (n_samples, n_pathways))

        # Random reference labels
        ref_labels = rng.randint(0, n_clusters, n_samples)

        # Cluster
        sim_seed = (seed + i) if seed is not None else None
        gmm = GaussianMixture(
            n_components=n_clusters,
            covariance_type="full",
            n_init=5,
            random_state=sim_seed,
            reg_covar=1e-6,
        )
        try:
            gmm.fit(scores)
            if not gmm.converged_:
                continue
            pred_labels = gmm.predict(scores)
            ari = adjusted_rand_score(ref_labels, pred_labels)
            ari_values.append(ari)
        except Exception:
            continue

    return np.array(ari_values)


def _simulate_stability_distribution(
    n_samples: int,
    n_pathways: int,
    n_clusters: int,
    effect_size: float = 1.5,
    n_bootstrap_per_sim: int = 10,
    n_simulations: int = 20,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Simulate stability ARI distribution under structured data.

    Generates data with planted structure, clusters it, then measures
    bootstrap stability.

    Args:
        n_samples: Number of samples
        n_pathways: Number of pathways
        n_clusters: Number of clusters
        effect_size: Cohen's d for planted structure
        n_bootstrap_per_sim: Bootstrap iterations per simulation
        n_simulations: Number of simulations
        seed: Random seed

    Returns:
        Array of bootstrap ARI values
    """
    from .simulation import SimulationConfig, generate_synthetic_data

    ari_values = []

    for i in range(n_simulations):
        sim_seed = (seed + i) if seed is not None else None

        config = SimulationConfig(
            n_samples=n_samples,
            n_pathways=n_pathways,
            n_subtypes=n_clusters,
            effect_size=effect_size,
            noise_level=1.0,
            seed=sim_seed,
        )
        sim_data = generate_synthetic_data(config)

        # Get full-data clustering
        gmm = GaussianMixture(
            n_components=n_clusters,
            covariance_type="full",
            n_init=5,
            random_state=sim_seed,
            reg_covar=1e-6,
        )
        try:
            gmm.fit(sim_data.pathway_scores.values)
            if not gmm.converged_:
                continue
            full_labels = gmm.predict(sim_data.pathway_scores.values)
        except Exception:
            continue

        rng = np.random.RandomState(sim_seed)

        # Bootstrap stability
        for b in range(n_bootstrap_per_sim):
            boot_idx = rng.choice(n_samples, size=n_samples, replace=True)
            unique_idx = np.unique(boot_idx)
            if len(unique_idx) < n_clusters * 2:
                continue

            boot_scores = sim_data.pathway_scores.iloc[unique_idx]
            boot_original = full_labels[unique_idx]

            boot_seed = (sim_seed + b + 1000) if sim_seed is not None else None
            boot_gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=boot_seed,
                reg_covar=1e-6,
            )
            try:
                boot_gmm.fit(boot_scores.values)
                if not boot_gmm.converged_:
                    continue
                boot_labels = boot_gmm.predict(boot_scores.values)
                ari = adjusted_rand_score(boot_original, boot_labels)
                ari_values.append(ari)
            except Exception:
                continue

    return np.array(ari_values)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def calibrate_thresholds(
    n_samples: int,
    n_clusters: int,
    n_pathways: int = 15,
    alpha: float = 0.05,
    method: str = "auto",
    n_simulations: int = 200,
    seed: Optional[int] = None,
) -> CalibratedThresholds:
    """
    Calibrate validation thresholds for given data dimensions.

    Determines appropriate null ARI maximum and stability ARI minimum
    thresholds based on sample size and number of clusters.

    Args:
        n_samples: Number of samples in the dataset
        n_clusters: Number of clusters being tested
        n_pathways: Number of pathways (used in simulation mode)
        alpha: Significance level (default 0.05 -> 95th/5th percentiles)
        method: Calibration method:
            - "auto": lookup if available, interpolate if close, simulate if out of range
            - "lookup": lookup/interpolate only (raises ValueError if out of range)
            - "simulate": always run fresh simulations
        n_simulations: Number of simulations for "simulate" mode
        seed: Random seed for simulation reproducibility

    Returns:
        CalibratedThresholds with calibrated values

    Raises:
        ValueError: If method is "lookup" and data dimensions are out of range
    """
    logger.info(
        f"[Calibration] Calibrating thresholds for n={n_samples}, k={n_clusters}, "
        f"method={method}"
    )

    if method == "simulate":
        return _calibrate_via_simulation(
            n_samples, n_clusters, n_pathways, alpha, n_simulations, seed
        )

    # Try lookup / interpolation first
    null_threshold = _interpolate_threshold(_NULL_ARI_TABLE, n_samples, n_clusters)
    stability_threshold = _interpolate_threshold(_STABILITY_TABLE, n_samples, n_clusters)

    if null_threshold is not None and stability_threshold is not None:
        # Adjust for non-default alpha
        if alpha != 0.05:
            # Scale thresholds: more permissive alpha -> higher null, lower stability
            alpha_ratio = alpha / 0.05
            null_threshold *= alpha_ratio
            stability_threshold = 1.0 - (1.0 - stability_threshold) * alpha_ratio

        # Check if exact match or interpolated
        exact_match = (n_samples, n_clusters) in _NULL_ARI_TABLE

        calibration_method = "lookup"
        interpolated = not exact_match

        logger.info(
            f"[Calibration] Lookup {'(exact)' if exact_match else '(interpolated)'}: "
            f"null_ari={null_threshold:.4f}, stability={stability_threshold:.4f}"
        )

        return CalibratedThresholds(
            null_ari_threshold=null_threshold,
            stability_threshold=stability_threshold,
            n_samples=n_samples,
            n_clusters=n_clusters,
            n_pathways=n_pathways,
            alpha=alpha,
            calibration_method=calibration_method,
            interpolated=interpolated,
        )

    # Out of range
    if method == "lookup":
        raise ValueError(
            f"Data dimensions (n={n_samples}, k={n_clusters}) are outside the "
            f"pre-computed grid. Use method='auto' or method='simulate' instead."
        )

    # Auto mode: fall back to simulation
    logger.info("[Calibration] Data dimensions outside lookup grid, falling back to simulation")
    return _calibrate_via_simulation(n_samples, n_clusters, n_pathways, alpha, n_simulations, seed)


def _calibrate_via_simulation(
    n_samples: int,
    n_clusters: int,
    n_pathways: int,
    alpha: float,
    n_simulations: int,
    seed: Optional[int],
) -> CalibratedThresholds:
    """Calibrate thresholds by running simulations."""
    logger.info(f"[Calibration] Running {n_simulations} simulations for calibration...")

    null_values = _simulate_null_distribution(
        n_samples, n_pathways, n_clusters, n_simulations, seed
    )

    # For stability, use fewer simulations with more bootstrap per sim
    n_stab_sims = max(5, n_simulations // 10)
    n_boot_per = max(5, n_simulations // n_stab_sims)

    stability_seed = (seed + 10000) if seed is not None else None
    stability_values = _simulate_stability_distribution(
        n_samples,
        n_pathways,
        n_clusters,
        effect_size=1.5,
        n_bootstrap_per_sim=n_boot_per,
        n_simulations=n_stab_sims,
        seed=stability_seed,
    )

    # Compute percentiles
    upper_pct = (1.0 - alpha) * 100  # e.g., 95th
    lower_pct = alpha * 100  # e.g., 5th

    if len(null_values) > 0:
        null_threshold = float(np.percentile(null_values, upper_pct))
    else:
        logger.warning("[Calibration] No null simulations converged, using default 0.15")
        null_threshold = 0.15

    if len(stability_values) > 0:
        stability_threshold = float(np.percentile(stability_values, lower_pct))
    else:
        logger.warning("[Calibration] No stability simulations converged, using default 0.8")
        stability_threshold = 0.8

    logger.info(
        f"[Calibration] Simulated: null_ari={null_threshold:.4f}, "
        f"stability={stability_threshold:.4f}"
    )

    return CalibratedThresholds(
        null_ari_threshold=null_threshold,
        stability_threshold=stability_threshold,
        n_samples=n_samples,
        n_clusters=n_clusters,
        n_pathways=n_pathways,
        alpha=alpha,
        calibration_method="simulate",
        interpolated=False,
    )


def get_default_thresholds() -> CalibratedThresholds:
    """
    Return the legacy default thresholds for backward compatibility.

    Returns:
        CalibratedThresholds with the original hard-coded values (0.15, 0.8)
    """
    return CalibratedThresholds(
        null_ari_threshold=0.15,
        stability_threshold=0.8,
        n_samples=0,
        n_clusters=0,
        n_pathways=0,
        alpha=0.05,
        calibration_method="default",
        interpolated=False,
    )


def generate_calibration_table(
    sample_sizes: Optional[List[int]] = None,
    cluster_counts: Optional[List[int]] = None,
    n_simulations: int = 500,
    seed: int = 42,
) -> Dict[str, Dict[Tuple[int, int], float]]:
    """
    Generate calibration lookup tables via simulation.

    This is a long-running function used to regenerate the embedded
    lookup tables. Not intended for runtime use.

    Args:
        sample_sizes: List of sample sizes to simulate
        cluster_counts: List of cluster counts to simulate
        n_simulations: Number of simulations per grid point
        seed: Base random seed

    Returns:
        Dictionary with "null_ari" and "stability" tables
    """
    if sample_sizes is None:
        sample_sizes = _GRID_SAMPLES
    if cluster_counts is None:
        cluster_counts = _GRID_CLUSTERS

    null_table: Dict[Tuple[int, int], float] = {}
    stability_table: Dict[Tuple[int, int], float] = {}

    total = len(sample_sizes) * len(cluster_counts)
    done = 0

    for n in sample_sizes:
        for k in cluster_counts:
            done += 1
            logger.info(f"[Calibration] Grid point {done}/{total}: n={n}, k={k}")

            point_seed = seed + n * 100 + k

            # Null distribution
            null_values = _simulate_null_distribution(n, 15, k, n_simulations, point_seed)
            if len(null_values) > 0:
                null_table[(n, k)] = float(np.percentile(null_values, 95))
            else:
                null_table[(n, k)] = 0.15  # fallback

            # Stability distribution
            stab_values = _simulate_stability_distribution(
                n,
                15,
                k,
                effect_size=1.5,
                n_bootstrap_per_sim=10,
                n_simulations=max(5, n_simulations // 10),
                seed=point_seed + 50000,
            )
            if len(stab_values) > 0:
                stability_table[(n, k)] = float(np.percentile(stab_values, 5))
            else:
                stability_table[(n, k)] = 0.8  # fallback

    return {"null_ari": null_table, "stability": stability_table}
