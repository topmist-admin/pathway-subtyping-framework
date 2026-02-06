"""
Clustering Module for the Pathway Subtyping Framework.

Implements multiple clustering algorithms for robustness analysis:
- Gaussian Mixture Models (GMM) - default
- K-means
- Hierarchical (agglomerative)
- Spectral clustering

Also includes model selection and cross-validation.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import KFold

logger = logging.getLogger(__name__)


class ClusteringAlgorithm(Enum):
    """Available clustering algorithms."""

    GMM = "gmm"
    KMEANS = "kmeans"
    HIERARCHICAL = "hierarchical"
    SPECTRAL = "spectral"


@dataclass
class ClusteringResult:
    """
    Result from clustering analysis.

    Attributes:
        labels: Cluster assignments for each sample
        n_clusters: Number of clusters found/used
        algorithm: Algorithm used
        silhouette: Silhouette score
        calinski_harabasz: Calinski-Harabasz index
        davies_bouldin: Davies-Bouldin index
        bic: BIC (for GMM only)
        converged: Whether algorithm converged (for GMM)
    """

    labels: np.ndarray
    n_clusters: int
    algorithm: ClusteringAlgorithm
    silhouette: float = 0.0
    calinski_harabasz: float = 0.0
    davies_bouldin: float = 0.0
    bic: Optional[float] = None
    converged: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "algorithm": self.algorithm.value,
            "n_clusters": self.n_clusters,
            "metrics": {
                "silhouette": round(self.silhouette, 4),
                "calinski_harabasz": round(self.calinski_harabasz, 4),
                "davies_bouldin": round(self.davies_bouldin, 4),
                "bic": round(self.bic, 2) if self.bic else None,
            },
            "converged": self.converged,
        }


def run_clustering(
    data: np.ndarray,
    n_clusters: int,
    algorithm: ClusteringAlgorithm = ClusteringAlgorithm.GMM,
    seed: Optional[int] = None,
    **kwargs,
) -> ClusteringResult:
    """
    Run clustering with specified algorithm.

    Args:
        data: Feature matrix (samples x features)
        n_clusters: Number of clusters
        algorithm: Clustering algorithm to use
        seed: Random seed
        **kwargs: Additional algorithm-specific parameters

    Returns:
        ClusteringResult with labels and metrics
    """
    if algorithm == ClusteringAlgorithm.GMM:
        return _run_gmm(data, n_clusters, seed, **kwargs)
    elif algorithm == ClusteringAlgorithm.KMEANS:
        return _run_kmeans(data, n_clusters, seed, **kwargs)
    elif algorithm == ClusteringAlgorithm.HIERARCHICAL:
        return _run_hierarchical(data, n_clusters, **kwargs)
    elif algorithm == ClusteringAlgorithm.SPECTRAL:
        return _run_spectral(data, n_clusters, seed, **kwargs)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")


def _run_gmm(
    data: np.ndarray,
    n_clusters: int,
    seed: Optional[int] = None,
    covariance_type: str = "full",
    n_init: int = 10,
    reg_covar: float = 1e-6,
) -> ClusteringResult:
    """Run Gaussian Mixture Model clustering."""
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type=covariance_type,
        n_init=n_init,
        random_state=seed,
        reg_covar=reg_covar,
    )

    gmm.fit(data)
    labels = gmm.predict(data)

    # Compute metrics
    metrics = _compute_cluster_metrics(data, labels)

    return ClusteringResult(
        labels=labels,
        n_clusters=n_clusters,
        algorithm=ClusteringAlgorithm.GMM,
        silhouette=metrics["silhouette"],
        calinski_harabasz=metrics["calinski_harabasz"],
        davies_bouldin=metrics["davies_bouldin"],
        bic=gmm.bic(data),
        converged=gmm.converged_,
    )


def _run_kmeans(
    data: np.ndarray,
    n_clusters: int,
    seed: Optional[int] = None,
    n_init: int = 10,
    max_iter: int = 300,
) -> ClusteringResult:
    """Run K-means clustering."""
    kmeans = KMeans(
        n_clusters=n_clusters,
        n_init=n_init,
        max_iter=max_iter,
        random_state=seed,
    )

    labels = kmeans.fit_predict(data)
    metrics = _compute_cluster_metrics(data, labels)

    return ClusteringResult(
        labels=labels,
        n_clusters=n_clusters,
        algorithm=ClusteringAlgorithm.KMEANS,
        silhouette=metrics["silhouette"],
        calinski_harabasz=metrics["calinski_harabasz"],
        davies_bouldin=metrics["davies_bouldin"],
    )


def _run_hierarchical(
    data: np.ndarray,
    n_clusters: int,
    linkage_method: str = "ward",
) -> ClusteringResult:
    """Run hierarchical agglomerative clustering."""
    clustering = AgglomerativeClustering(
        n_clusters=n_clusters,
        linkage=linkage_method,
    )

    labels = clustering.fit_predict(data)
    metrics = _compute_cluster_metrics(data, labels)

    return ClusteringResult(
        labels=labels,
        n_clusters=n_clusters,
        algorithm=ClusteringAlgorithm.HIERARCHICAL,
        silhouette=metrics["silhouette"],
        calinski_harabasz=metrics["calinski_harabasz"],
        davies_bouldin=metrics["davies_bouldin"],
    )


def _run_spectral(
    data: np.ndarray,
    n_clusters: int,
    seed: Optional[int] = None,
    affinity: str = "rbf",
) -> ClusteringResult:
    """Run spectral clustering."""
    spectral = SpectralClustering(
        n_clusters=n_clusters,
        affinity=affinity,
        random_state=seed,
        n_init=10,
    )

    labels = spectral.fit_predict(data)
    metrics = _compute_cluster_metrics(data, labels)

    return ClusteringResult(
        labels=labels,
        n_clusters=n_clusters,
        algorithm=ClusteringAlgorithm.SPECTRAL,
        silhouette=metrics["silhouette"],
        calinski_harabasz=metrics["calinski_harabasz"],
        davies_bouldin=metrics["davies_bouldin"],
    )


def _compute_cluster_metrics(
    data: np.ndarray,
    labels: np.ndarray,
) -> Dict[str, float]:
    """Compute clustering quality metrics."""
    n_unique = len(np.unique(labels))

    if n_unique < 2 or n_unique >= len(labels):
        return {
            "silhouette": 0.0,
            "calinski_harabasz": 0.0,
            "davies_bouldin": float("inf"),
        }

    try:
        sil = silhouette_score(data, labels)
    except Exception:
        sil = 0.0

    try:
        ch = calinski_harabasz_score(data, labels)
    except Exception:
        ch = 0.0

    try:
        db = davies_bouldin_score(data, labels)
    except Exception:
        db = float("inf")

    return {
        "silhouette": sil,
        "calinski_harabasz": ch,
        "davies_bouldin": db,
    }


@dataclass
class ModelSelectionResult:
    """
    Result of cluster number selection.

    Attributes:
        optimal_k: Recommended number of clusters
        k_range: Range of k values tested
        bic_values: BIC for each k (GMM only)
        silhouette_values: Silhouette score for each k
        method: Method used for selection
    """

    optimal_k: int
    k_range: List[int]
    bic_values: Dict[int, float] = field(default_factory=dict)
    silhouette_values: Dict[int, float] = field(default_factory=dict)
    method: str = "bic"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "optimal_k": self.optimal_k,
            "k_range": self.k_range,
            "method": self.method,
            "bic_values": {str(k): round(v, 2) for k, v in self.bic_values.items()},
            "silhouette_values": {str(k): round(v, 4) for k, v in self.silhouette_values.items()},
        }


def select_n_clusters(
    data: np.ndarray,
    k_range: List[int] = None,
    method: str = "bic",
    seed: Optional[int] = None,
) -> ModelSelectionResult:
    """
    Select optimal number of clusters.

    Args:
        data: Feature matrix (samples x features)
        k_range: Range of k values to test (default: [2, 8])
        method: Selection method ("bic" or "silhouette")
        seed: Random seed

    Returns:
        ModelSelectionResult with optimal k and scores
    """
    if k_range is None:
        k_range = list(range(2, 9))

    bic_values = {}
    silhouette_values = {}

    for k in k_range:
        result = run_clustering(data, k, ClusteringAlgorithm.GMM, seed)

        if result.bic is not None:
            bic_values[k] = result.bic
        silhouette_values[k] = result.silhouette

    # Select optimal k
    if method == "bic" and bic_values:
        optimal_k = min(bic_values, key=bic_values.get)
    elif silhouette_values:
        optimal_k = max(silhouette_values, key=silhouette_values.get)
    else:
        optimal_k = k_range[0]

    return ModelSelectionResult(
        optimal_k=optimal_k,
        k_range=k_range,
        bic_values=bic_values,
        silhouette_values=silhouette_values,
        method=method,
    )


@dataclass
class CrossValidationResult:
    """
    Result of cross-validation analysis.

    Attributes:
        mean_ari: Mean ARI across folds
        std_ari: Standard deviation of ARI
        fold_aris: ARI for each fold
        n_folds: Number of folds
    """

    mean_ari: float
    std_ari: float
    fold_aris: List[float]
    n_folds: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "mean_ari": round(self.mean_ari, 4),
            "std_ari": round(self.std_ari, 4),
            "n_folds": self.n_folds,
            "fold_aris": [round(a, 4) for a in self.fold_aris],
        }


def cross_validate_clustering(
    data: np.ndarray,
    n_clusters: int,
    algorithm: ClusteringAlgorithm = ClusteringAlgorithm.GMM,
    n_folds: int = 5,
    seed: Optional[int] = None,
) -> CrossValidationResult:
    """
    Cross-validate clustering stability.

    Splits data into folds, clusters each fold, and measures consistency.

    Args:
        data: Feature matrix (samples x features)
        n_clusters: Number of clusters
        algorithm: Clustering algorithm
        n_folds: Number of CV folds
        seed: Random seed

    Returns:
        CrossValidationResult with stability metrics
    """
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)

    # First, get reference clustering on full data
    full_result = run_clustering(data, n_clusters, algorithm, seed)
    reference_labels = full_result.labels

    fold_aris = []

    for train_idx, test_idx in kf.split(data):
        # Cluster on training fold
        train_data = data[train_idx]
        run_clustering(train_data, n_clusters, algorithm, seed)

        # Create cluster model and predict on test
        if algorithm == ClusteringAlgorithm.GMM:
            gmm = GaussianMixture(
                n_components=n_clusters,
                random_state=seed,
                reg_covar=1e-6,
            )
            gmm.fit(train_data)
            test_labels = gmm.predict(data[test_idx])
        else:
            # For other algorithms, re-cluster full data and compare
            full_fold_result = run_clustering(data, n_clusters, algorithm, seed)
            test_labels = full_fold_result.labels[test_idx]

        # Compare to reference labels on test set
        ref_test_labels = reference_labels[test_idx]
        ari = adjusted_rand_score(ref_test_labels, test_labels)
        fold_aris.append(ari)

    return CrossValidationResult(
        mean_ari=np.mean(fold_aris),
        std_ari=np.std(fold_aris),
        fold_aris=fold_aris,
        n_folds=n_folds,
    )


@dataclass
class AlgorithmComparisonResult:
    """
    Result of comparing multiple clustering algorithms.

    Attributes:
        results: ClusteringResult for each algorithm
        pairwise_ari: ARI between each pair of algorithms
        consensus_labels: Labels from most stable algorithm
    """

    results: Dict[str, ClusteringResult]
    pairwise_ari: Dict[str, float]
    consensus_labels: np.ndarray
    most_stable_algorithm: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "algorithms": {name: result.to_dict() for name, result in self.results.items()},
            "pairwise_ari": {k: round(v, 4) for k, v in self.pairwise_ari.items()},
            "most_stable_algorithm": self.most_stable_algorithm,
        }


def compare_algorithms(
    data: np.ndarray,
    n_clusters: int,
    algorithms: Optional[List[ClusteringAlgorithm]] = None,
    seed: Optional[int] = None,
) -> AlgorithmComparisonResult:
    """
    Compare multiple clustering algorithms.

    Args:
        data: Feature matrix (samples x features)
        n_clusters: Number of clusters
        algorithms: Algorithms to compare (default: all)
        seed: Random seed

    Returns:
        AlgorithmComparisonResult with comparison metrics
    """
    if algorithms is None:
        algorithms = [
            ClusteringAlgorithm.GMM,
            ClusteringAlgorithm.KMEANS,
            ClusteringAlgorithm.HIERARCHICAL,
            ClusteringAlgorithm.SPECTRAL,
        ]

    results = {}
    labels = {}

    for algo in algorithms:
        try:
            result = run_clustering(data, n_clusters, algo, seed)
            results[algo.value] = result
            labels[algo.value] = result.labels
        except Exception as e:
            logger.warning(f"Algorithm {algo.value} failed: {e}")

    # Compute pairwise ARI
    pairwise_ari = {}
    algo_names = list(labels.keys())

    for i, name1 in enumerate(algo_names):
        for name2 in algo_names[i + 1 :]:
            ari = adjusted_rand_score(labels[name1], labels[name2])
            pairwise_ari[f"{name1}_vs_{name2}"] = ari

    # Find most stable (highest average ARI with others)
    avg_ari = {}
    for name in algo_names:
        aris = []
        for key, ari in pairwise_ari.items():
            if name in key:
                aris.append(ari)
        avg_ari[name] = np.mean(aris) if aris else 0.0

    most_stable = max(avg_ari, key=avg_ari.get)

    return AlgorithmComparisonResult(
        results=results,
        pairwise_ari=pairwise_ari,
        consensus_labels=labels[most_stable],
        most_stable_algorithm=most_stable,
    )
