"""
Benchmark Comparison Module.

Compares the framework's pathway-based GMM subtyping approach against
established dimensionality reduction and clustering methods:
- NMF (Non-negative Matrix Factorization) on gene burdens
- PCA + K-means on pathway scores
- Direct K-means on gene burdens (no pathway aggregation)
- Random baseline (null model)

References:
- Lee DD, Seung HS. Learning the parts of objects by non-negative
  matrix factorization. Nature. 1999;401(6755):788-791.
- Jolliffe IT. Principal Component Analysis. 2nd ed.
  Springer; 2002.

Research use only. Not for clinical decision-making.
"""

import logging
import time
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF, PCA
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)

from .clustering import ClusteringAlgorithm, run_clustering

logger = logging.getLogger(__name__)


# =============================================================================
# ENUMS
# =============================================================================


class BenchmarkMethod(Enum):
    """
    Benchmark methods for comparison.

    PATHWAY_GMM: Framework's pathway-based GMM (reference method).
    NMF_CLUSTERING: NMF on gene burdens followed by K-means.
    PCA_KMEANS: PCA on pathway scores followed by K-means.
    GENE_KMEANS: Direct K-means on gene burdens (no pathway aggregation).
    RANDOM_BASELINE: Random label assignment (null model).
    """

    PATHWAY_GMM = "pathway_gmm"
    NMF_CLUSTERING = "nmf_clustering"
    PCA_KMEANS = "pca_kmeans"
    GENE_KMEANS = "gene_kmeans"
    RANDOM_BASELINE = "random_baseline"


# =============================================================================
# DATACLASSES
# =============================================================================


@dataclass
class BenchmarkResult:
    """
    Result from a single benchmark method.

    Attributes:
        method: Which benchmark method was used.
        predicted_labels: Cluster assignments.
        n_clusters_found: Number of unique clusters in output.
        ari: Adjusted Rand Index vs ground truth (None if no truth).
        nmi: Normalized Mutual Information vs ground truth.
        silhouette: Silhouette score of the clustering.
        runtime_seconds: Wall-clock execution time.
        converged: Whether the method converged successfully.
    """

    method: BenchmarkMethod
    predicted_labels: np.ndarray
    n_clusters_found: int
    ari: Optional[float] = None
    nmi: Optional[float] = None
    silhouette: float = 0.0
    runtime_seconds: float = 0.0
    converged: bool = True

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "method": self.method.value,
            "n_clusters_found": self.n_clusters_found,
            "silhouette": round(self.silhouette, 4),
            "runtime_seconds": round(self.runtime_seconds, 4),
            "converged": self.converged,
        }
        if self.ari is not None:
            result["ari"] = round(self.ari, 4)
        if self.nmi is not None:
            result["nmi"] = round(self.nmi, 4)
        return result


@dataclass
class BenchmarkComparisonResult:
    """
    Results from comparing multiple methods on the same data.

    Attributes:
        method_results: BenchmarkResult for each method.
        best_method: Method with highest ARI (or silhouette if no truth).
        ranking: Methods sorted by primary metric descending.
        n_samples: Number of samples compared.
        n_clusters: Number of clusters requested.
    """

    method_results: Dict[str, BenchmarkResult]
    best_method: str
    ranking: List[str]
    n_samples: int
    n_clusters: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_samples": self.n_samples,
            "n_clusters": self.n_clusters,
            "best_method": self.best_method,
            "ranking": self.ranking,
            "methods": {k: v.to_dict() for k, v in self.method_results.items()},
        }

    def format_report(self) -> str:
        lines = [
            "Benchmark Comparison Report",
            "=" * 40,
            f"Samples: {self.n_samples}, Clusters: {self.n_clusters}",
            f"Best method: {self.best_method}",
            "",
            "Ranking:",
        ]
        for i, method_name in enumerate(self.ranking, 1):
            r = self.method_results[method_name]
            ari_str = f"ARI={r.ari:.3f}" if r.ari is not None else "ARI=N/A"
            lines.append(
                f"  {i}. {method_name}: {ari_str}, "
                f"sil={r.silhouette:.3f}, "
                f"time={r.runtime_seconds:.3f}s"
            )
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        return [
            "Lee DD, Seung HS. Learning the parts of objects by "
            "non-negative matrix factorization. Nature. 1999;"
            "401(6755):788-791.",
            "Jolliffe IT. Principal Component Analysis. 2nd ed. " "Springer; 2002.",
            "Hubert L, Arabie P. Comparing partitions. " "J Classif. 1985;2(1):193-218.",
        ]


@dataclass
class BenchmarkSweepResult:
    """
    Results from sweeping benchmark across multiple conditions.

    Attributes:
        conditions: List of (effect_size, n_samples) tuples tested.
        results_per_condition: BenchmarkComparisonResult per condition.
        method_mean_ari: Mean ARI for each method across conditions.
        method_wins: Number of conditions each method ranked first.
    """

    conditions: List[Dict[str, float]]
    results_per_condition: List[BenchmarkComparisonResult]
    method_mean_ari: Dict[str, float]
    method_wins: Dict[str, int]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_conditions": len(self.conditions),
            "conditions": self.conditions,
            "method_mean_ari": {k: round(v, 4) for k, v in self.method_mean_ari.items()},
            "method_wins": self.method_wins,
        }

    def format_report(self) -> str:
        lines = [
            "Benchmark Sweep Report",
            "=" * 40,
            f"Conditions tested: {len(self.conditions)}",
            "",
            "Mean ARI across conditions:",
        ]
        sorted_methods = sorted(self.method_mean_ari, key=self.method_mean_ari.get, reverse=True)
        for method in sorted_methods:
            lines.append(
                f"  {method}: ARI={self.method_mean_ari[method]:.3f}, "
                f"wins={self.method_wins.get(method, 0)}"
            )
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        return [
            "Lee DD, Seung HS. Learning the parts of objects by "
            "non-negative matrix factorization. Nature. 1999;"
            "401(6755):788-791.",
            "Jolliffe IT. Principal Component Analysis. 2nd ed. " "Springer; 2002.",
            "Hubert L, Arabie P. Comparing partitions. " "J Classif. 1985;2(1):193-218.",
        ]


# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================


def run_single_benchmark(
    method: BenchmarkMethod,
    gene_burdens: pd.DataFrame,
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    true_labels: Optional[np.ndarray] = None,
    seed: Optional[int] = None,
) -> BenchmarkResult:
    """
    Run a single benchmark method and evaluate performance.

    Args:
        method: Which benchmark method to run.
        gene_burdens: Gene burden matrix (samples x genes).
        pathway_scores: Pathway score matrix (samples x pathways).
        n_clusters: Number of clusters to find.
        true_labels: Ground truth labels for ARI/NMI (optional).
        seed: Random seed for reproducibility.

    Returns:
        BenchmarkResult with labels, metrics, and timing.
    """
    logger.info("[Benchmark] Running %s with k=%d", method.value, n_clusters)

    start = time.time()
    converged = True

    if method == BenchmarkMethod.PATHWAY_GMM:
        labels, converged = _run_pathway_gmm(pathway_scores, n_clusters, seed)
    elif method == BenchmarkMethod.NMF_CLUSTERING:
        labels, converged = _run_nmf(gene_burdens, n_clusters, seed)
    elif method == BenchmarkMethod.PCA_KMEANS:
        labels, converged = _run_pca_kmeans(pathway_scores, n_clusters, seed)
    elif method == BenchmarkMethod.GENE_KMEANS:
        labels, converged = _run_gene_kmeans(gene_burdens, n_clusters, seed)
    elif method == BenchmarkMethod.RANDOM_BASELINE:
        labels, converged = _run_random(len(pathway_scores), n_clusters, seed)
    else:
        raise ValueError(f"Unknown benchmark method: {method}")

    elapsed = time.time() - start

    # Compute metrics
    ari = None
    nmi = None
    if true_labels is not None:
        ari = float(adjusted_rand_score(true_labels, labels))
        nmi = float(normalized_mutual_info_score(true_labels, labels))

    sil = 0.0
    n_unique = len(np.unique(labels))
    if n_unique >= 2:
        # Use pathway_scores as the feature space for silhouette
        try:
            sil = float(silhouette_score(pathway_scores.values, labels))
        except Exception:
            sil = 0.0

    logger.info(
        "[Benchmark] %s: ARI=%s, silhouette=%.3f, time=%.3fs",
        method.value,
        f"{ari:.3f}" if ari is not None else "N/A",
        sil,
        elapsed,
    )

    return BenchmarkResult(
        method=method,
        predicted_labels=labels,
        n_clusters_found=n_unique,
        ari=ari,
        nmi=nmi,
        silhouette=sil,
        runtime_seconds=elapsed,
        converged=converged,
    )


def run_benchmark_comparison(
    gene_burdens: pd.DataFrame,
    pathway_scores: pd.DataFrame,
    pathways: Dict[str, List[str]],
    true_labels: Optional[np.ndarray] = None,
    n_clusters: Optional[int] = None,
    methods: Optional[List[BenchmarkMethod]] = None,
    seed: Optional[int] = None,
) -> BenchmarkComparisonResult:
    """
    Compare multiple benchmark methods on the same data.

    Args:
        gene_burdens: Gene burden matrix (samples x genes).
        pathway_scores: Pathway score matrix (samples x pathways).
        pathways: Dictionary mapping pathway names to gene lists.
        true_labels: Ground truth labels for ARI/NMI (optional).
        n_clusters: Number of clusters. If None and true_labels given,
            uses the number of unique true labels.
        methods: Methods to benchmark. Default: all methods.
        seed: Random seed for reproducibility.

    Returns:
        BenchmarkComparisonResult with ranked method results.
    """
    if methods is None:
        methods = list(BenchmarkMethod)

    if n_clusters is None:
        if true_labels is not None:
            n_clusters = len(np.unique(true_labels))
        else:
            n_clusters = 3

    n_samples = len(pathway_scores)
    logger.info(
        "[Benchmark] Comparing %d methods on %d samples (k=%d)",
        len(methods),
        n_samples,
        n_clusters,
    )

    method_results = {}
    for method in methods:
        try:
            result = run_single_benchmark(
                method=method,
                gene_burdens=gene_burdens,
                pathway_scores=pathway_scores,
                n_clusters=n_clusters,
                true_labels=true_labels,
                seed=seed,
            )
            method_results[method.value] = result
        except Exception as e:
            logger.warning("[Benchmark] Method %s failed: %s", method.value, str(e))

    # Rank methods
    if true_labels is not None:
        # Rank by ARI
        ranking = sorted(
            method_results,
            key=lambda m: method_results[m].ari if method_results[m].ari is not None else -1.0,
            reverse=True,
        )
    else:
        # Rank by silhouette
        ranking = sorted(
            method_results,
            key=lambda m: method_results[m].silhouette,
            reverse=True,
        )

    best_method = ranking[0] if ranking else "none"

    logger.info("[Benchmark] Best method: %s", best_method)

    return BenchmarkComparisonResult(
        method_results=method_results,
        best_method=best_method,
        ranking=ranking,
        n_samples=n_samples,
        n_clusters=n_clusters,
    )


def run_benchmark_sweep(
    effect_sizes: Optional[List[float]] = None,
    sample_sizes: Optional[List[int]] = None,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    methods: Optional[List[BenchmarkMethod]] = None,
    seed: Optional[int] = None,
) -> BenchmarkSweepResult:
    """
    Sweep benchmark comparison across multiple conditions.

    Generates synthetic data at each (effect_size, n_samples) combination
    and runs all benchmark methods, producing publication-ready results.

    Args:
        effect_sizes: Effect sizes to test. Default: [0.25, 0.5, 1.0, 1.5, 2.0].
        sample_sizes: Sample sizes to test. Default: [100].
        n_pathways: Number of pathways per simulation.
        n_subtypes: Number of planted subtypes.
        methods: Methods to benchmark. Default: all methods.
        seed: Random seed for reproducibility.

    Returns:
        BenchmarkSweepResult with per-condition results and summary.
    """
    from .simulation import SimulationConfig, generate_synthetic_data

    if effect_sizes is None:
        effect_sizes = [0.25, 0.5, 1.0, 1.5, 2.0]
    if sample_sizes is None:
        sample_sizes = [100]
    if methods is None:
        methods = list(BenchmarkMethod)

    conditions = []
    results_per_condition = []

    # Collect ARI per method across conditions
    method_aris: Dict[str, List[float]] = {m.value: [] for m in methods}
    method_wins: Dict[str, int] = {m.value: 0 for m in methods}

    total = len(effect_sizes) * len(sample_sizes)
    logger.info("[Benchmark] Starting sweep: %d conditions", total)

    condition_idx = 0
    for n_samples in sample_sizes:
        for effect in effect_sizes:
            condition_idx += 1
            logger.info(
                "[Benchmark] Condition %d/%d: effect=%.2f, n=%d",
                condition_idx,
                total,
                effect,
                n_samples,
            )

            condition = {"effect_size": effect, "n_samples": n_samples}
            conditions.append(condition)

            config = SimulationConfig(
                n_samples=n_samples,
                n_pathways=n_pathways,
                n_subtypes=n_subtypes,
                effect_size=effect,
                noise_level=1.0,
                seed=seed + condition_idx if seed is not None else None,
            )
            sim_data = generate_synthetic_data(config)

            comp_result = run_benchmark_comparison(
                gene_burdens=sim_data.gene_burdens,
                pathway_scores=sim_data.pathway_scores,
                pathways=sim_data.pathways,
                true_labels=sim_data.true_labels,
                n_clusters=n_subtypes,
                methods=methods,
                seed=seed,
            )
            results_per_condition.append(comp_result)

            # Track ARI per method
            for m in methods:
                r = comp_result.method_results.get(m.value)
                if r is not None and r.ari is not None:
                    method_aris[m.value].append(r.ari)

            # Track wins
            if comp_result.best_method in method_wins:
                method_wins[comp_result.best_method] += 1

    # Compute mean ARI per method
    method_mean_ari = {}
    for m_name, aris in method_aris.items():
        method_mean_ari[m_name] = float(np.mean(aris)) if aris else 0.0

    logger.info("[Benchmark] Sweep complete")

    return BenchmarkSweepResult(
        conditions=conditions,
        results_per_condition=results_per_condition,
        method_mean_ari=method_mean_ari,
        method_wins=method_wins,
    )


# =============================================================================
# PRIVATE HELPERS
# =============================================================================


def _run_pathway_gmm(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> tuple:
    """Run framework's pathway-based GMM clustering."""
    result = run_clustering(
        data=pathway_scores.values,
        n_clusters=n_clusters,
        algorithm=ClusteringAlgorithm.GMM,
        seed=seed,
    )
    return result.labels, result.converged


def _run_nmf(
    gene_burdens: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> tuple:
    """Run NMF on gene burdens followed by K-means on latent factors."""
    data = gene_burdens.values.copy()

    # NMF requires non-negative input; clip any negatives
    data = np.clip(data, 0, None)

    # Replace zeros to avoid issues with NMF
    data[data == 0] = 1e-10

    try:
        nmf = NMF(
            n_components=n_clusters,
            init="nndsvda",
            random_state=seed,
            max_iter=500,
        )
        W = nmf.fit_transform(data)

        kmeans = KMeans(
            n_clusters=n_clusters,
            n_init=10,
            random_state=seed,
        )
        labels = kmeans.fit_predict(W)
        return labels, True
    except Exception as e:
        logger.warning("[Benchmark] NMF failed: %s, falling back to K-means", str(e))
        kmeans = KMeans(
            n_clusters=n_clusters,
            n_init=10,
            random_state=seed,
        )
        labels = kmeans.fit_predict(data)
        return labels, False


def _run_pca_kmeans(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> tuple:
    """Run PCA on pathway scores followed by K-means."""
    data = pathway_scores.values

    n_components = min(n_clusters, data.shape[1])
    pca = PCA(n_components=n_components, random_state=seed)
    transformed = pca.fit_transform(data)

    kmeans = KMeans(
        n_clusters=n_clusters,
        n_init=10,
        random_state=seed,
    )
    labels = kmeans.fit_predict(transformed)
    return labels, True


def _run_gene_kmeans(
    gene_burdens: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> tuple:
    """Run K-means directly on gene burdens (no pathway aggregation)."""
    kmeans = KMeans(
        n_clusters=n_clusters,
        n_init=10,
        random_state=seed,
    )
    labels = kmeans.fit_predict(gene_burdens.values)
    return labels, True


def _run_random(
    n_samples: int,
    n_clusters: int,
    seed: Optional[int] = None,
) -> tuple:
    """Generate random cluster assignments (null model)."""
    rng = np.random.RandomState(seed)
    labels = rng.randint(0, n_clusters, size=n_samples)
    return labels, True
