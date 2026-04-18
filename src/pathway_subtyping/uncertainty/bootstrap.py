"""
Non-parametric bootstrap over cells for MSV scores.

Resamples the input cell population with replacement, recomputes any MSV
score (or arbitrary statistic) on each resample, and reports percentile-based
intervals plus per-cell score distributions.

Default implementation is single-process. Parallelism is opt-in via
``n_jobs > 1``, using ``joblib`` when installed, otherwise a thread pool.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


ArrayLike = Union[np.ndarray, pd.DataFrame]


@dataclass
class BootstrapResult:
    """Result of a bootstrap resampling run.

    Attributes:
        point: Point estimate on the full sample (shape matches score output).
        samples: (n_bootstrap, *point_shape) array of bootstrap replicates.
        lower: Lower percentile interval (e.g., 2.5%).
        upper: Upper percentile interval (e.g., 97.5%).
        ci_level: Nominal two-sided coverage (e.g., 0.95).
    """

    point: np.ndarray
    samples: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    ci_level: float

    @property
    def width(self) -> np.ndarray:
        return self.upper - self.lower

    def se(self) -> np.ndarray:
        """Bootstrap standard error — std across replicates."""
        return np.asarray(self.samples.std(axis=0, ddof=1))

    def summary(self) -> str:
        return (
            f"BootstrapResult(n_bootstrap={self.samples.shape[0]}, "
            f"ci_level={self.ci_level}, "
            f"mean_width={float(np.asarray(self.width).mean()):.4f})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_bootstrap": int(self.samples.shape[0]),
            "ci_level": self.ci_level,
            "point": np.asarray(self.point).tolist(),
            "lower": np.asarray(self.lower).tolist(),
            "upper": np.asarray(self.upper).tolist(),
            "se": self.se().tolist(),
        }


class BootstrapMSV:
    """Non-parametric bootstrap for MSV scores or any per-cell statistic.

    Usage::

        bootstrap = BootstrapMSV(n_bootstrap=200, ci_level=0.95, seed=42)
        result = bootstrap.run(expression_df, score_fn=compute_stress_score)
        # result.lower / result.upper bracket each score

    The score function signature is ``score_fn(X) -> np.ndarray``. The bootstrap
    calls it once on the full input (``point``) and then once per resample.
    Per-cell interval reporting uses the per-cell index to average scores
    across replicates that happened to include that cell.
    """

    def __init__(
        self,
        n_bootstrap: int = 200,
        ci_level: float = 0.95,
        seed: Optional[int] = None,
        n_jobs: int = 1,
    ):
        if n_bootstrap < 10:
            raise ValueError("n_bootstrap should be at least 10")
        if not 0.0 < ci_level < 1.0:
            raise ValueError("ci_level must be in (0, 1)")
        self.n_bootstrap = int(n_bootstrap)
        self.ci_level = float(ci_level)
        self.seed = seed
        self.n_jobs = max(1, int(n_jobs))

    # -------------------------------------------------------- scalar / vector ---
    def run(
        self,
        X: ArrayLike,
        score_fn: Callable[[ArrayLike], np.ndarray],
        per_cell: bool = False,
    ) -> BootstrapResult:
        """Run bootstrap resampling.

        Args:
            X: Input cell data (DataFrame or ndarray, rows = cells).
            score_fn: Function that returns either a scalar, a vector indexed
                      by feature (e.g., one score per pathway), or a per-cell
                      vector whose length equals ``len(X)``.
            per_cell: If True, treat score_fn's output as per-cell and
                      aggregate across bootstrap replicates cell-wise.

        Returns:
            BootstrapResult.
        """
        n_cells = len(X)
        rng = np.random.default_rng(self.seed)

        point = np.asarray(score_fn(X))
        logger.info(
            "[BootstrapMSV] point shape=%s, n_bootstrap=%d, per_cell=%s",
            point.shape, self.n_bootstrap, per_cell,
        )

        if per_cell:
            samples, cell_counts = self._run_per_cell(
                X, score_fn, n_cells, rng, point.shape,
            )
            lower, upper = self._percentile_per_cell(
                samples, cell_counts, self.ci_level,
            )
        else:
            samples = self._run_stat(X, score_fn, n_cells, rng, point.shape)
            alpha = (1.0 - self.ci_level) / 2.0
            lower = np.quantile(samples, alpha, axis=0)
            upper = np.quantile(samples, 1.0 - alpha, axis=0)

        return BootstrapResult(
            point=point,
            samples=samples,
            lower=np.asarray(lower),
            upper=np.asarray(upper),
            ci_level=self.ci_level,
        )

    # ----------------------------------------------------------- workhorses ---
    def _run_stat(
        self,
        X: ArrayLike,
        score_fn: Callable[[ArrayLike], np.ndarray],
        n_cells: int,
        rng: np.random.Generator,
        out_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Bootstrap a statistic whose output does NOT index by cell."""
        samples = np.empty((self.n_bootstrap, *out_shape), dtype=float)

        def one_replicate(seed: int) -> Tuple[int, np.ndarray]:
            r = np.random.default_rng(seed)
            idx = r.integers(0, n_cells, size=n_cells)
            X_bs = _index(X, idx)
            return seed, np.asarray(score_fn(X_bs))

        seeds = rng.integers(0, 2**31 - 1, size=self.n_bootstrap)

        if self.n_jobs == 1:
            for b, s in enumerate(seeds):
                _, samples[b] = one_replicate(int(s))
        else:
            with ThreadPoolExecutor(max_workers=self.n_jobs) as pool:
                futures = {
                    pool.submit(one_replicate, int(s)): b for b, s in enumerate(seeds)
                }
                for fut in as_completed(futures):
                    b = futures[fut]
                    _, samples[b] = fut.result()

        return samples

    def _run_per_cell(
        self,
        X: ArrayLike,
        score_fn: Callable[[ArrayLike], np.ndarray],
        n_cells: int,
        rng: np.random.Generator,
        out_shape: Tuple[int, ...],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Bootstrap where score_fn returns one value per input cell."""
        if out_shape[0] != n_cells:
            raise ValueError(
                f"per_cell=True but score_fn returned shape {out_shape}, "
                f"expected first axis to equal n_cells={n_cells}"
            )

        extra_shape = out_shape[1:]

        # per-cell accumulator: list-of-values per cell, then quantiles
        per_cell_values: list = [[] for _ in range(n_cells)]

        def one_replicate(seed: int):
            r = np.random.default_rng(seed)
            idx = r.integers(0, n_cells, size=n_cells)
            X_bs = _index(X, idx)
            out = np.asarray(score_fn(X_bs))
            return idx, out

        seeds = rng.integers(0, 2**31 - 1, size=self.n_bootstrap)

        samples = np.empty((self.n_bootstrap, *out_shape), dtype=float)
        cell_counts = np.zeros(n_cells, dtype=int)

        if self.n_jobs == 1:
            for b, s in enumerate(seeds):
                idx, out = one_replicate(int(s))
                samples[b] = np.nan  # placeholder — we use per_cell_values
                for i_bs, cell_id in enumerate(idx):
                    per_cell_values[int(cell_id)].append(out[i_bs])
                    cell_counts[int(cell_id)] += 1
        else:
            with ThreadPoolExecutor(max_workers=self.n_jobs) as pool:
                futures = [pool.submit(one_replicate, int(s)) for s in seeds]
                for fut in as_completed(futures):
                    idx, out = fut.result()
                    for i_bs, cell_id in enumerate(idx):
                        per_cell_values[int(cell_id)].append(out[i_bs])
                        cell_counts[int(cell_id)] += 1

        # Reshape per-cell accumulator into (n_bootstrap_max, n_cells, *extra)
        max_draws = int(cell_counts.max()) if len(cell_counts) else 0
        padded = np.full(
            (max_draws, n_cells, *extra_shape), fill_value=np.nan, dtype=float
        )
        for cell_id, vals in enumerate(per_cell_values):
            for k, v in enumerate(vals):
                padded[k, cell_id] = v

        return padded, cell_counts

    # ------------------------------------------------------ percentile helper ---
    @staticmethod
    def _percentile_per_cell(
        samples: np.ndarray,
        cell_counts: np.ndarray,
        ci_level: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Per-cell percentile intervals using NaN-aware quantiles."""
        alpha = (1.0 - ci_level) / 2.0
        # samples shape: (K, n_cells, *extra)
        lower = np.nanquantile(samples, alpha, axis=0)
        upper = np.nanquantile(samples, 1.0 - alpha, axis=0)
        # Fill cells that were never resampled with NaN width
        never_sampled = cell_counts == 0
        if np.any(never_sampled):
            logger.warning(
                "[BootstrapMSV] %d cells never resampled; interval will be NaN",
                int(never_sampled.sum()),
            )
        return lower, upper


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _index(X: ArrayLike, idx: np.ndarray) -> ArrayLike:
    """Index with replacement; works for DataFrame and ndarray."""
    if isinstance(X, pd.DataFrame):
        return X.iloc[idx].reset_index(drop=True)
    return np.asarray(X)[idx]
