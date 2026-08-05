"""
Cross-platform harmonization benchmark.

Synthetic protocol that takes a single expression matrix, simulates several
"platforms" by applying structured distortions, and measures pathway-score
correlation between platforms before and after alignment. Used to establish
the roadmap acceptance criterion:

    Harmonized pathway-level rho across 10x vs Smart-seq2 on matched cortex
    samples exceeds 0.75 (pre-harmonization baseline: 0.55-0.65 typical).

Real-data benchmarking requires matched-sample cross-platform datasets
that PSF ships without. This module provides the synthetic ground-truth
framework and a CLI entry point so users with matched cohorts can drop in
their own data and reuse the protocol.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import zlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .align import CrossPlatformAligner
from .confidence import HarmonizationReport
from .uce import FallbackEmbedder

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Synthetic platform distortions
# --------------------------------------------------------------------------- #


def simulate_platform_distortion(
    scores: pd.DataFrame,
    platform: str,
    seed: Optional[int] = None,
    shift_scale: float = 0.3,
    noise_scale: float = 0.35,
    cell_covariate: Optional[np.ndarray] = None,
    coupling_scale: float = 0.8,
) -> pd.DataFrame:
    """Apply a deterministic "platform-style" distortion to pathway scores.

    The distortion has three components:

        1. A fixed additive per-pathway shift (platform-level detection bias).
        2. An embedding-dependent per-cell shift via a platform-specific
           coupling matrix ``A_p`` applied to ``cell_covariate``. This
           models the realistic case where platform bias interacts with cell
           biology (e.g., 10x dropout affects low-abundance marker genes in
           cell types more than others), producing rank-altering distortion
           that an embedding-anchored aligner can later remove.
        3. i.i.d. Gaussian noise (per-cell measurement variability).

    Args:
        scores: Reference pathway-score DataFrame (n_cells x n_pathways).
        platform: Platform label; drives the deterministic shift signature.
        seed: Random seed for the noise component.
        shift_scale: Magnitude of the per-platform additive shift.
        noise_scale: Stddev of the per-cell Gaussian noise.
        cell_covariate: Optional (n_cells, k) cell-level signal (e.g., PCs
            of the reference scores) that the platform-specific coupling
            is applied to. If None, only shift + noise are used.
        coupling_scale: Stddev of the platform-specific coupling matrix.

    Returns:
        Distorted DataFrame with the same shape and column names.
    """
    rng = np.random.default_rng(seed if seed is not None else 0)
    platform_rng = np.random.default_rng(zlib.crc32(str(platform).encode()) % (2**32))

    n_cells, n_pathways = scores.shape
    shift = platform_rng.standard_normal(n_pathways) * shift_scale
    distorted = scores.to_numpy() + shift

    if cell_covariate is not None:
        cov = np.asarray(cell_covariate)
        if cov.shape[0] != n_cells:
            raise ValueError(f"cell_covariate has {cov.shape[0]} rows; expected {n_cells}")
        coupling = platform_rng.standard_normal((cov.shape[1], n_pathways)) * coupling_scale
        distorted = distorted + cov @ coupling

    noise = rng.standard_normal(distorted.shape) * noise_scale
    distorted = distorted + noise
    return pd.DataFrame(distorted, index=scores.index, columns=scores.columns)


def _pairwise_rho(a: pd.DataFrame, b: pd.DataFrame) -> float:
    """Mean Spearman rho across matched pathway columns."""
    if a.shape != b.shape:
        raise ValueError("dataframes must have identical shape for rho")
    rhos: List[float] = []
    for col in a.columns:
        va, vb = a[col].to_numpy(), b[col].to_numpy()
        if np.std(va) == 0 or np.std(vb) == 0:
            continue
        rho, _ = spearmanr(va, vb)
        if not np.isnan(rho):
            rhos.append(float(rho))
    return float(np.mean(rhos)) if rhos else float("nan")


# --------------------------------------------------------------------------- #
# Benchmark harness
# --------------------------------------------------------------------------- #


@dataclass
class CrossPlatformBenchmark:
    """End-to-end benchmark harness.

    Usage::

        bench = CrossPlatformBenchmark(reference_scores=scores_df,
                                       platforms=["10x", "smartseq2"])
        result = bench.run(seed=0)
        # result["post_rho"] should exceed 0.75

    The benchmark holds the reference pathway-score matrix and a list of
    platform labels. Each ``run()`` call draws a fresh seed, distorts the
    reference matrix per platform, fits the embedder + aligner, and returns
    pre- and post-harmonization Spearman rho between platform pairs.
    """

    reference_scores: pd.DataFrame
    platforms: List[str]
    shift_scale: float = 0.3
    noise_scale: float = 0.2
    coupling_scale: float = 0.4
    covariate_dim: int = 4
    embedding_dim: int = 32
    _embedder: FallbackEmbedder = field(init=False)
    _covariate: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        if len(self.platforms) < 2:
            raise ValueError("benchmark needs at least 2 platforms")
        # Embedder is pretrained on the reference (undistorted) data. This
        # simulates UCE's behaviour: the upstream model has seen cells across
        # many platforms and provides a platform-invariant embedding. At run
        # time we project distorted data through this frozen embedder.
        self._embedder = FallbackEmbedder(embedding_dim=self.embedding_dim, seed=0).fit(
            self.reference_scores
        )
        # Biological covariate used to apply the platform-specific coupling:
        # top covariate_dim PCs of the reference scores, normalised to
        # unit stddev per dimension so coupling_scale has an interpretable
        # magnitude regardless of reference-data variance.
        covariate = (
            FallbackEmbedder(embedding_dim=self.covariate_dim, seed=1)
            .fit(self.reference_scores)
            .embed(self.reference_scores)
        )
        std = covariate.std(axis=0, keepdims=True)
        std[std == 0] = 1.0
        self._covariate = covariate / std

    # ---------------------------------------------------------------- run ---
    def run(self, seed: int = 0) -> Dict[str, Any]:
        """One benchmark trial.

        Returns a dict with ``pre_rho``, ``post_rho``, ``improvement``,
        per-platform drift, and the fitted alignment result.
        """
        distorted: Dict[str, pd.DataFrame] = {}
        for i, plat in enumerate(self.platforms):
            distorted[plat] = simulate_platform_distortion(
                self.reference_scores,
                platform=plat,
                seed=seed * len(self.platforms) + i,
                shift_scale=self.shift_scale,
                noise_scale=self.noise_scale,
                cell_covariate=self._covariate,
                coupling_scale=self.coupling_scale,
            )

        # Stack platforms for fitting; per-cell platform label
        all_scores = pd.concat(distorted.values(), axis=0).reset_index(drop=True)
        platform_labels = np.concatenate([np.full(len(df), plat) for plat, df in distorted.items()])

        # Platform-invariant reference embedding for each cell: UCE returns
        # a biology-only embedding that does NOT drift with platform, so
        # we use the frozen reference embedding here (tiled across each
        # platform's copy of the cells). When integrating real UCE, this
        # array is produced by UCE on the same cells regardless of platform.
        ref_embeddings = self._embedder.embed(self.reference_scores)
        embeddings = np.tile(ref_embeddings, (len(self.platforms), 1))

        aligner = CrossPlatformAligner()
        aligned = aligner.fit_transform(all_scores, platform_labels, embeddings)

        # Split aligned back per platform for rho computation
        counts = [len(df) for df in distorted.values()]
        boundaries = np.cumsum([0] + counts)
        per_plat_aligned = {
            plat: aligned.aligned_scores.iloc[boundaries[i] : boundaries[i + 1]].reset_index(
                drop=True
            )
            for i, plat in enumerate(distorted)
        }

        # Compute pre/post rho for every pair of platforms
        pre_rhos: List[float] = []
        post_rhos: List[float] = []
        pair_details: Dict[str, Dict[str, float]] = {}
        pnames = list(distorted.keys())
        for i in range(len(pnames)):
            for j in range(i + 1, len(pnames)):
                a, b = pnames[i], pnames[j]
                pre = _pairwise_rho(
                    distorted[a].reset_index(drop=True),
                    distorted[b].reset_index(drop=True),
                )
                post = _pairwise_rho(per_plat_aligned[a], per_plat_aligned[b])
                pre_rhos.append(pre)
                post_rhos.append(post)
                pair_details[f"{a}_vs_{b}"] = {"pre": pre, "post": post}

        pre_rho = float(np.mean(pre_rhos))
        post_rho = float(np.mean(post_rhos))

        report = HarmonizationReport.from_alignment(aligned)

        return {
            "seed": int(seed),
            "platforms": list(distorted.keys()),
            "pre_rho": pre_rho,
            "post_rho": post_rho,
            "improvement": post_rho - pre_rho,
            "pair_details": pair_details,
            "per_platform_drift": report.per_platform_drift,
            "mean_confidence": float(report.confidence.mean()),
            "alignment": aligned,
            "report": report,
        }

    # --------------------------------------------------------- multi-seed ---
    def run_many(self, n_seeds: int = 5) -> Dict[str, Any]:
        pre, post, imp = [], [], []
        for s in range(n_seeds):
            r = self.run(seed=s)
            pre.append(r["pre_rho"])
            post.append(r["post_rho"])
            imp.append(r["improvement"])
        return {
            "n_seeds": int(n_seeds),
            "pre_rho_mean": float(np.mean(pre)),
            "post_rho_mean": float(np.mean(post)),
            "improvement_mean": float(np.mean(imp)),
            "post_rho_min": float(np.min(post)),
            "post_rho_max": float(np.max(post)),
        }
