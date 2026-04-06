"""
Validation metrics for benchmark datasets.

Computes:
- Silhouette coefficient (cluster quality)
- Bootstrap ARI (stability across resampling)
- Adjusted Rand Index (ARI) vs ground truth
"""

import logging
import numpy as np
from typing import Callable, Optional
from sklearn.metrics import silhouette_score, adjusted_rand_score

logger = logging.getLogger(__name__)


def compute_silhouette(
    data: np.ndarray,
    labels: np.ndarray,
    sample_size: Optional[int] = None
) -> float:
    """
    Compute silhouette coefficient for clustering quality.

    Args:
        data: (n_samples, n_features) expression matrix
        labels: (n_samples,) cluster assignments
        sample_size: If set, subsample data for speed (silhouette on subset)

    Returns:
        Silhouette coefficient (-1 to 1; higher is better)
    """
    try:
        if sample_size and len(data) > sample_size:
            indices = np.random.choice(len(data), sample_size, replace=False)
            data_subset = data[indices]
            labels_subset = labels[indices]
        else:
            data_subset = data
            labels_subset = labels

        if len(np.unique(labels_subset)) < 2:
            return 0.0

        return float(silhouette_score(data_subset, labels_subset))

    except Exception as e:
        logger.warning(f"Could not compute silhouette: {e}")
        return 0.0


def bootstrap_ari_stability(
    data: np.ndarray,
    true_labels: np.ndarray,
    clustering_func: Callable[[np.ndarray], np.ndarray],
    n_iterations: int = 100,
    sample_fraction: float = 0.8,
    seed: Optional[int] = None
) -> np.ndarray:
    """
    Compute bootstrap stability: ARI between predicted and true labels
    across bootstrap resampled datasets.

    Tests robustness: how stable are cluster assignments when we resample
    the data?

    Args:
        data: (n_samples, n_features) expression matrix
        true_labels: (n_samples,) ground truth cluster labels
        clustering_func: Function(data) -> cluster_labels
        n_iterations: Number of bootstrap resamples
        sample_fraction: Fraction of data to resample per iteration (0.8 = 80%)
        seed: Random seed for reproducibility

    Returns:
        Array of ARI values (n_iterations,) between predicted and true labels
    """
    if seed is not None:
        np.random.seed(seed)

    ari_values = []
    n_samples = len(data)
    n_resample = int(n_samples * sample_fraction)

    logger.info(f"Computing bootstrap ARI stability ({n_iterations} iterations, {sample_fraction:.0%} resampling)...")

    for iteration in range(n_iterations):
        try:
            # Bootstrap resample with replacement
            indices = np.random.choice(n_samples, n_resample, replace=True)
            data_resampled = data[indices]
            true_labels_resampled = true_labels[indices]

            # Run clustering on resampled data
            predicted_labels = clustering_func(data_resampled)

            # Align labels if needed (relabel to match ground truth clusters)
            if len(np.unique(predicted_labels)) != len(np.unique(true_labels_resampled)):
                # Relabel to match cluster count
                predicted_labels = _relabel_clusters(predicted_labels, len(np.unique(true_labels_resampled)))

            # Compute ARI
            ari = adjusted_rand_score(true_labels_resampled, predicted_labels)
            ari_values.append(ari)

            if (iteration + 1) % 20 == 0:
                logger.debug(f"  Iteration {iteration + 1}/{n_iterations}: ARI = {ari:.4f}")

        except Exception as e:
            logger.warning(f"  Bootstrap iteration {iteration + 1} failed: {e}")
            ari_values.append(0.0)  # Assume 0 ARI if clustering fails

    ari_array = np.array(ari_values)
    logger.info(f"Bootstrap ARI: mean={ari_array.mean():.4f}, median={np.median(ari_array):.4f}, std={ari_array.std():.4f}")

    return ari_array


def _relabel_clusters(labels: np.ndarray, target_k: int) -> np.ndarray:
    """
    Relabel cluster assignments to target number of clusters.

    Simple approach: if more clusters than target, merge smallest clusters.
    If fewer clusters than target, split largest clusters.

    Args:
        labels: Current cluster assignments
        target_k: Target number of clusters

    Returns:
        Relabeled cluster assignments
    """
    unique_labels = np.unique(labels)
    current_k = len(unique_labels)

    if current_k == target_k:
        return labels

    if current_k > target_k:
        # Merge clusters: assign smallest clusters to nearest larger cluster
        label_counts = np.bincount(labels)
        sorted_indices = np.argsort(label_counts)

        # Keep the largest target_k clusters
        keep_labels = sorted(sorted_indices[-target_k:])

        # Remap
        relabeled = np.zeros_like(labels)
        for new_label, old_label in enumerate(keep_labels):
            relabeled[labels == old_label] = new_label

        # Merge remaining labels into closest cluster
        for old_label in sorted_indices[:-target_k]:
            closest = keep_labels[0]
            relabeled[labels == old_label] = np.where(keep_labels == closest)[0][0]

        return relabeled

    else:
        # Split clusters: duplicate labels (simple approach)
        relabeled = labels.copy()
        n_extra = target_k - current_k

        for i in range(n_extra):
            largest_label = np.argmax(np.bincount(relabeled))
            mask = relabeled == largest_label
            indices = np.where(mask)[0]
            split_point = len(indices) // 2
            relabeled[indices[split_point:]] = max(np.unique(relabeled)) + 1

        return relabeled


class BootstrapARI:
    """Callable wrapper for bootstrap ARI computation."""

    def __init__(
        self,
        true_labels: np.ndarray,
        clustering_func: Callable,
        n_iterations: int = 100,
        sample_fraction: float = 0.8,
        seed: Optional[int] = None
    ):
        """
        Args:
            true_labels: Ground truth cluster labels
            clustering_func: Function(data) -> cluster_labels
            n_iterations: Number of bootstrap iterations
            sample_fraction: Fraction of data per iteration
            seed: Random seed
        """
        self.true_labels = true_labels
        self.clustering_func = clustering_func
        self.n_iterations = n_iterations
        self.sample_fraction = sample_fraction
        self.seed = seed

    def __call__(self, data: np.ndarray) -> np.ndarray:
        """Compute bootstrap ARI for data."""
        return bootstrap_ari_stability(
            data,
            self.true_labels,
            self.clustering_func,
            n_iterations=self.n_iterations,
            sample_fraction=self.sample_fraction,
            seed=self.seed
        )
