"""
Sensitivity Analysis Module.

Implements methods to assess the robustness of subtype discovery
to analytical parameter choices:
- Systematic parameter variation across key settings
- Concordance measurement between parameter configurations
- Stability scoring for overall robustness
- Identification of parameters that most affect results

References:
- von Luxburg U. A tutorial on spectral clustering.
  Stat Comput. 2007;17(4):395-416.
- Hennig C. Cluster-wise assessment of cluster stability.
  Comput Stat Data Anal. 2007;52(1):258-271.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering

logger = logging.getLogger(__name__)


# =============================================================================
# ENUMS
# =============================================================================


class SensitivityParameter(Enum):
    """
    Parameters that can be varied in sensitivity analysis.

    CLUSTERING_ALGORITHM: Test GMM, K-means, Hierarchical, Spectral.
    N_CLUSTERS: Test a range of cluster counts.
    NORMALIZATION: Test different pathway score normalization methods.
    FEATURE_SUBSET: Test stability when dropping pathways one at a time.
    """

    CLUSTERING_ALGORITHM = "clustering_algorithm"
    N_CLUSTERS = "n_clusters"
    NORMALIZATION = "normalization"
    FEATURE_SUBSET = "feature_subset"


# =============================================================================
# DATACLASSES
# =============================================================================


@dataclass
class ParameterVariationResult:
    """
    Result from varying a single parameter.

    Attributes:
        parameter: Which parameter was varied.
        configurations: List of configuration values tested.
        labels_per_config: Cluster labels for each configuration.
        pairwise_ari: ARI between each pair of configurations.
        mean_ari: Mean pairwise ARI across all pairs.
        min_ari: Minimum pairwise ARI (worst-case agreement).
        reference_ari: ARI of each configuration vs. the reference (first config).
    """

    parameter: SensitivityParameter
    configurations: List[str]
    labels_per_config: Dict[str, np.ndarray]
    pairwise_ari: Dict[str, float]
    mean_ari: float
    min_ari: float
    reference_ari: Dict[str, float]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "parameter": self.parameter.value,
            "n_configurations": len(self.configurations),
            "configurations": self.configurations,
            "mean_pairwise_ari": round(self.mean_ari, 4),
            "min_pairwise_ari": round(self.min_ari, 4),
            "reference_ari": {
                k: round(v, 4) for k, v in self.reference_ari.items()
            },
        }


@dataclass
class SensitivityAnalysisResult:
    """
    Comprehensive sensitivity analysis result.

    Attributes:
        parameter_results: Results for each parameter varied.
        overall_stability: Mean ARI across all parameter variations.
        most_sensitive_parameter: Parameter with lowest mean ARI.
        least_sensitive_parameter: Parameter with highest mean ARI.
        is_robust: Whether overall stability exceeds threshold.
        robustness_threshold: Threshold for robustness (default 0.7).
    """

    parameter_results: Dict[str, ParameterVariationResult]
    overall_stability: float
    most_sensitive_parameter: str
    least_sensitive_parameter: str
    is_robust: bool
    robustness_threshold: float = 0.7

    def to_dict(self) -> Dict[str, Any]:
        return {
            "overall_stability": round(self.overall_stability, 4),
            "is_robust": self.is_robust,
            "robustness_threshold": self.robustness_threshold,
            "most_sensitive_parameter": self.most_sensitive_parameter,
            "least_sensitive_parameter": self.least_sensitive_parameter,
            "parameters": {
                k: v.to_dict() for k, v in self.parameter_results.items()
            },
        }

    def format_report(self) -> str:
        lines = [
            "Sensitivity Analysis Report",
            "=" * 40,
            f"Overall stability (mean ARI): {self.overall_stability:.3f}",
            f"Robust: {self.is_robust} (threshold: {self.robustness_threshold})",
            "",
            f"Most sensitive to:  {self.most_sensitive_parameter}",
            f"Least sensitive to: {self.least_sensitive_parameter}",
            "",
            "Per-parameter results:",
        ]
        for name, result in self.parameter_results.items():
            lines.append(
                f"  {name}: mean ARI = {result.mean_ari:.3f}, "
                f"min ARI = {result.min_ari:.3f} "
                f"({len(result.configurations)} configs)"
            )
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        return [
            "Hennig C. Cluster-wise assessment of cluster stability. "
            "Comput Stat Data Anal. 2007;52(1):258-271.",
            "von Luxburg U. A tutorial on spectral clustering. "
            "Stat Comput. 2007;17(4):395-416.",
        ]


# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================


def vary_clustering_algorithm(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> ParameterVariationResult:
    """
    Test sensitivity to clustering algorithm choice.

    Runs GMM, K-means, Hierarchical, and Spectral clustering
    on the same data and compares pairwise ARI.

    Args:
        pathway_scores: DataFrame of shape (n_samples, n_pathways).
        n_clusters: Number of clusters to use.
        seed: Random seed for reproducibility.

    Returns:
        ParameterVariationResult with algorithm comparison.
    """
    logger.info(
        "[Sensitivity] Varying clustering algorithm with k=%d", n_clusters
    )

    algorithms = {
        "GMM": lambda: GaussianMixture(
            n_components=n_clusters,
            covariance_type="full",
            n_init=5,
            reg_covar=1e-6,
            random_state=seed,
        ),
        "KMeans": lambda: KMeans(
            n_clusters=n_clusters, n_init=10, random_state=seed
        ),
        "Hierarchical": lambda: AgglomerativeClustering(
            n_clusters=n_clusters
        ),
    }

    # Only include spectral if n_samples > n_clusters
    n_samples = len(pathway_scores)
    if n_samples >= n_clusters + 2:
        algorithms["Spectral"] = lambda: SpectralClustering(
            n_clusters=n_clusters,
            random_state=seed,
            affinity="rbf",
            n_init=5,
        )

    data = pathway_scores.values
    labels_per_config = {}

    for name, make_model in algorithms.items():
        try:
            model = make_model()
            if name == "GMM":
                model.fit(data)
                labels = model.predict(data)
            else:
                labels = model.fit_predict(data)
            labels_per_config[name] = labels
        except Exception as e:
            logger.warning(
                "[Sensitivity] Algorithm %s failed: %s", name, str(e)
            )

    return _build_variation_result(
        SensitivityParameter.CLUSTERING_ALGORITHM,
        labels_per_config,
    )


def vary_n_clusters(
    pathway_scores: pd.DataFrame,
    cluster_range: Tuple[int, int] = (2, 6),
    seed: Optional[int] = None,
) -> ParameterVariationResult:
    """
    Test sensitivity to number of clusters.

    Runs GMM clustering with different k values and compares
    pairwise ARI between the resulting labelings.

    Args:
        pathway_scores: DataFrame of shape (n_samples, n_pathways).
        cluster_range: (min_k, max_k) inclusive.
        seed: Random seed for reproducibility.

    Returns:
        ParameterVariationResult with cluster count comparison.
    """
    min_k, max_k = cluster_range
    logger.info(
        "[Sensitivity] Varying n_clusters from %d to %d", min_k, max_k
    )

    data = pathway_scores.values
    labels_per_config = {}

    for k in range(min_k, max_k + 1):
        name = f"k={k}"
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                n_init=5,
                reg_covar=1e-6,
                random_state=seed,
            )
            gmm.fit(data)
            labels_per_config[name] = gmm.predict(data)
        except Exception as e:
            logger.warning(
                "[Sensitivity] GMM with k=%d failed: %s", k, str(e)
            )

    return _build_variation_result(
        SensitivityParameter.N_CLUSTERS,
        labels_per_config,
    )


def vary_feature_subset(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> ParameterVariationResult:
    """
    Test sensitivity to pathway inclusion via leave-one-out.

    For each pathway, removes it and re-clusters on the remaining
    pathways. Measures how much the results change.

    Args:
        pathway_scores: DataFrame of shape (n_samples, n_pathways).
        n_clusters: Number of clusters.
        seed: Random seed for reproducibility.

    Returns:
        ParameterVariationResult with feature subset comparison.
    """
    logger.info(
        "[Sensitivity] Leave-one-out feature sensitivity for %d features",
        len(pathway_scores.columns),
    )

    labels_per_config = {}

    # Full dataset baseline
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type="full",
        n_init=5,
        reg_covar=1e-6,
        random_state=seed,
    )
    gmm.fit(pathway_scores.values)
    labels_per_config["all_features"] = gmm.predict(pathway_scores.values)

    # Leave-one-out
    for col in pathway_scores.columns:
        subset = pathway_scores.drop(columns=[col])
        if len(subset.columns) < 2:
            continue
        name = f"without_{col}"
        try:
            gmm_sub = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                reg_covar=1e-6,
                random_state=seed,
            )
            gmm_sub.fit(subset.values)
            labels_per_config[name] = gmm_sub.predict(subset.values)
        except Exception as e:
            logger.warning(
                "[Sensitivity] Leave-one-out without %s failed: %s",
                col,
                str(e),
            )

    return _build_variation_result(
        SensitivityParameter.FEATURE_SUBSET,
        labels_per_config,
    )


def vary_normalization(
    pathway_scores_raw: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> ParameterVariationResult:
    """
    Test sensitivity to normalization method.

    Applies different normalization approaches to raw pathway
    scores and compares clustering results.

    Args:
        pathway_scores_raw: Un-normalized pathway scores.
        n_clusters: Number of clusters.
        seed: Random seed for reproducibility.

    Returns:
        ParameterVariationResult with normalization comparison.
    """
    logger.info("[Sensitivity] Varying normalization method")

    labels_per_config = {}
    normalizations = {
        "zscore": _normalize_zscore,
        "minmax": _normalize_minmax,
        "robust": _normalize_robust,
        "rank": _normalize_rank,
    }

    for name, norm_fn in normalizations.items():
        try:
            normalized = norm_fn(pathway_scores_raw)
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                reg_covar=1e-6,
                random_state=seed,
            )
            gmm.fit(normalized.values)
            labels_per_config[name] = gmm.predict(normalized.values)
        except Exception as e:
            logger.warning(
                "[Sensitivity] Normalization %s failed: %s", name, str(e)
            )

    return _build_variation_result(
        SensitivityParameter.NORMALIZATION,
        labels_per_config,
    )


def run_sensitivity_analysis(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
    robustness_threshold: float = 0.7,
    parameters: Optional[List[SensitivityParameter]] = None,
) -> SensitivityAnalysisResult:
    """
    Run comprehensive sensitivity analysis across multiple parameters.

    Tests how robust the clustering results are to changes in
    algorithm, cluster count, feature set, and normalization.

    Args:
        pathway_scores: DataFrame of shape (n_samples, n_pathways).
        n_clusters: Reference number of clusters.
        seed: Random seed for reproducibility.
        robustness_threshold: Minimum mean ARI for robustness.
        parameters: Which parameters to test. If None, tests all.

    Returns:
        SensitivityAnalysisResult with comprehensive stability metrics.
    """
    if parameters is None:
        parameters = list(SensitivityParameter)

    logger.info(
        "[Sensitivity] Running sensitivity analysis for %d parameters",
        len(parameters),
    )

    results = {}

    for param in parameters:
        if param == SensitivityParameter.CLUSTERING_ALGORITHM:
            results[param.value] = vary_clustering_algorithm(
                pathway_scores, n_clusters, seed
            )
        elif param == SensitivityParameter.N_CLUSTERS:
            min_k = max(2, n_clusters - 1)
            max_k = n_clusters + 2
            results[param.value] = vary_n_clusters(
                pathway_scores, (min_k, max_k), seed
            )
        elif param == SensitivityParameter.FEATURE_SUBSET:
            results[param.value] = vary_feature_subset(
                pathway_scores, n_clusters, seed
            )
        elif param == SensitivityParameter.NORMALIZATION:
            results[param.value] = vary_normalization(
                pathway_scores, n_clusters, seed
            )

    # Compute overall metrics
    mean_aris = {
        name: r.mean_ari for name, r in results.items()
    }

    if mean_aris:
        overall = np.mean(list(mean_aris.values()))
        most_sensitive = min(mean_aris, key=mean_aris.get)
        least_sensitive = max(mean_aris, key=mean_aris.get)
    else:
        overall = 0.0
        most_sensitive = "none"
        least_sensitive = "none"

    is_robust = overall >= robustness_threshold

    logger.info(
        "[Sensitivity] Overall stability: %.3f (robust: %s)",
        overall,
        is_robust,
    )

    return SensitivityAnalysisResult(
        parameter_results=results,
        overall_stability=float(overall),
        most_sensitive_parameter=most_sensitive,
        least_sensitive_parameter=least_sensitive,
        is_robust=is_robust,
        robustness_threshold=robustness_threshold,
    )


# =============================================================================
# PRIVATE HELPERS
# =============================================================================


def _build_variation_result(
    parameter: SensitivityParameter,
    labels_per_config: Dict[str, np.ndarray],
) -> ParameterVariationResult:
    """Build a ParameterVariationResult from config labels."""
    configs = list(labels_per_config.keys())

    # Compute pairwise ARI
    pairwise = {}
    ari_values = []
    for i, c1 in enumerate(configs):
        for c2 in configs[i + 1:]:
            key = f"{c1}_vs_{c2}"
            ari = adjusted_rand_score(
                labels_per_config[c1], labels_per_config[c2]
            )
            pairwise[key] = float(ari)
            ari_values.append(ari)

    mean_ari = float(np.mean(ari_values)) if ari_values else 0.0
    min_ari = float(np.min(ari_values)) if ari_values else 0.0

    # ARI vs reference (first config)
    reference_ari = {}
    if configs:
        ref_labels = labels_per_config[configs[0]]
        for c in configs:
            reference_ari[c] = float(
                adjusted_rand_score(ref_labels, labels_per_config[c])
            )

    return ParameterVariationResult(
        parameter=parameter,
        configurations=configs,
        labels_per_config=labels_per_config,
        pairwise_ari=pairwise,
        mean_ari=mean_ari,
        min_ari=min_ari,
        reference_ari=reference_ari,
    )


def _normalize_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score normalization (zero mean, unit variance)."""
    result = df.copy()
    for col in result.columns:
        std = result[col].std()
        if std > 0:
            result[col] = (result[col] - result[col].mean()) / std
    return result


def _normalize_minmax(df: pd.DataFrame) -> pd.DataFrame:
    """Min-max normalization to [0, 1]."""
    result = df.copy()
    for col in result.columns:
        col_min = result[col].min()
        col_max = result[col].max()
        rng = col_max - col_min
        if rng > 0:
            result[col] = (result[col] - col_min) / rng
    return result


def _normalize_robust(df: pd.DataFrame) -> pd.DataFrame:
    """Robust normalization using median and IQR."""
    result = df.copy()
    for col in result.columns:
        median = result[col].median()
        q1 = result[col].quantile(0.25)
        q3 = result[col].quantile(0.75)
        iqr = q3 - q1
        if iqr > 0:
            result[col] = (result[col] - median) / iqr
    return result


def _normalize_rank(df: pd.DataFrame) -> pd.DataFrame:
    """Rank-based normalization (percentile ranks)."""
    result = df.copy()
    for col in result.columns:
        result[col] = result[col].rank(pct=True)
    return result
