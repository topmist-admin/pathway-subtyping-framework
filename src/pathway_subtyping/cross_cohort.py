"""
Cross-cohort validation utilities.

Tools for comparing subtype definitions across different cohorts
to validate reproducibility of discovered molecular subtypes.
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


@dataclass
class CohortResult:
    """Results from a single cohort analysis."""

    name: str
    pathway_scores: pd.DataFrame
    cluster_labels: np.ndarray
    cluster_names: Dict[int, str]
    n_samples: int
    n_clusters: int


@dataclass
class CrossCohortResult:
    """Results from cross-cohort validation."""

    cohort_a: str
    cohort_b: str
    transfer_ari: float
    projection_ari: float
    pathway_correlation: float
    shared_subtypes: List[str]
    details: Dict[str, Any] = field(default_factory=dict)


def load_cohort_result(output_dir: str) -> CohortResult:
    """
    Load results from a completed pipeline run.

    Args:
        output_dir: Path to pipeline output directory

    Returns:
        CohortResult with loaded data
    """
    output_path = Path(output_dir)

    # Load pathway scores
    scores_path = output_path / "pathway_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"Pathway scores not found: {scores_path}")
    pathway_scores = pd.read_csv(scores_path, index_col=0)

    # Load cluster assignments
    assignments_path = output_path / "subtype_assignments.csv"
    if not assignments_path.exists():
        raise FileNotFoundError(f"Assignments not found: {assignments_path}")
    assignments = pd.read_csv(assignments_path)

    # Load report for metadata
    report_path = output_path / "report.json"
    if report_path.exists():
        with open(report_path) as f:
            report = json.load(f)
        name = report.get("pipeline_name", output_path.name)
    else:
        name = output_path.name

    # Build cluster name mapping
    cluster_names = dict(zip(assignments["cluster_id"], assignments["cluster_label"]))

    return CohortResult(
        name=name,
        pathway_scores=pathway_scores,
        cluster_labels=assignments["cluster_id"].values,
        cluster_names=cluster_names,
        n_samples=len(assignments),
        n_clusters=assignments["cluster_id"].nunique(),
    )


def compare_cohorts(
    cohort_a: CohortResult,
    cohort_b: CohortResult,
    seed: int = 42,
) -> CrossCohortResult:
    """
    Compare subtype definitions between two cohorts.

    Performs two validation approaches:
    1. Transfer Learning: Train GMM on cohort A, apply to cohort B
    2. Projection: Project cohort B samples into cohort A's PCA space

    Args:
        cohort_a: First cohort (reference)
        cohort_b: Second cohort (validation)
        seed: Random seed

    Returns:
        CrossCohortResult with comparison metrics
    """
    logger.info(f"Comparing {cohort_a.name} vs {cohort_b.name}")

    # Find common pathways
    common_pathways = list(
        set(cohort_a.pathway_scores.columns) & set(cohort_b.pathway_scores.columns)
    )

    if len(common_pathways) < 2:
        raise ValueError("Insufficient common pathways between cohorts")

    logger.info(f"Found {len(common_pathways)} common pathways")

    # Align pathway scores
    scores_a = cohort_a.pathway_scores[common_pathways]
    scores_b = cohort_b.pathway_scores[common_pathways]

    # 1. Transfer Learning: Train on A, predict on B
    transfer_ari = _transfer_learning_validation(
        scores_a,
        cohort_a.cluster_labels,
        scores_b,
        cohort_b.cluster_labels,
        cohort_a.n_clusters,
        seed,
    )

    # 2. Projection: Compare cluster means
    projection_ari = _projection_validation(
        scores_a, cohort_a.cluster_labels, scores_b, cohort_b.cluster_labels, seed
    )

    # 3. Pathway correlation
    pathway_corr = _pathway_correlation(scores_a, scores_b)

    # 4. Find shared subtypes (by name similarity)
    shared = _find_shared_subtypes(cohort_a.cluster_names, cohort_b.cluster_names)

    return CrossCohortResult(
        cohort_a=cohort_a.name,
        cohort_b=cohort_b.name,
        transfer_ari=transfer_ari,
        projection_ari=projection_ari,
        pathway_correlation=pathway_corr,
        shared_subtypes=shared,
        details={
            "common_pathways": len(common_pathways),
            "cohort_a_samples": cohort_a.n_samples,
            "cohort_b_samples": cohort_b.n_samples,
            "cohort_a_clusters": cohort_a.n_clusters,
            "cohort_b_clusters": cohort_b.n_clusters,
        },
    )


def _transfer_learning_validation(
    scores_a: pd.DataFrame,
    labels_a: np.ndarray,
    scores_b: pd.DataFrame,
    labels_b: np.ndarray,
    n_clusters: int,
    seed: int,
) -> float:
    """Train GMM on cohort A, predict on cohort B, compare to B's labels."""
    # Train on cohort A
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type="full",
        n_init=10,
        random_state=seed,
    )
    gmm.fit(scores_a.values)

    # Predict on cohort B
    predicted_b = gmm.predict(scores_b.values)

    # Compare to cohort B's actual labels
    ari = adjusted_rand_score(labels_b, predicted_b)
    return round(ari, 4)


def _projection_validation(
    scores_a: pd.DataFrame,
    labels_a: np.ndarray,
    scores_b: pd.DataFrame,
    labels_b: np.ndarray,
    seed: int,
) -> float:
    """Project both cohorts to shared space, compare cluster structure."""
    from sklearn.decomposition import PCA

    # Fit PCA on combined data
    combined = pd.concat([scores_a, scores_b], axis=0)
    pca = PCA(n_components=min(10, len(scores_a.columns)), random_state=seed)
    pca.fit(combined.values)

    # Transform both cohorts
    proj_a = pca.transform(scores_a.values)
    proj_b = pca.transform(scores_b.values)

    # Compute cluster centroids in cohort A
    unique_labels = np.unique(labels_a)
    centroids_a = {}
    for label in unique_labels:
        mask = labels_a == label
        centroids_a[label] = proj_a[mask].mean(axis=0)

    # Assign cohort B samples to nearest cohort A centroid
    assigned_b = []
    for sample in proj_b:
        distances = {
            label: np.linalg.norm(sample - centroid) for label, centroid in centroids_a.items()
        }
        assigned_b.append(min(distances, key=distances.get))

    # Compare assignments to cohort B's actual labels
    ari = adjusted_rand_score(labels_b, assigned_b)
    return round(ari, 4)


def _pathway_correlation(
    scores_a: pd.DataFrame,
    scores_b: pd.DataFrame,
) -> float:
    """Compute mean correlation of pathway variances between cohorts."""
    var_a = scores_a.var()
    var_b = scores_b.var()

    # Pearson correlation of pathway variances
    corr = var_a.corr(var_b)
    return round(corr, 4)


def _find_shared_subtypes(
    names_a: Dict[int, str],
    names_b: Dict[int, str],
) -> List[str]:
    """Find subtypes with matching names across cohorts."""
    labels_a = set(names_a.values())
    labels_b = set(names_b.values())
    return list(labels_a & labels_b)


def generate_cross_cohort_report(
    results: List[CrossCohortResult],
    output_path: str,
) -> None:
    """
    Generate a cross-cohort validation report.

    Args:
        results: List of CrossCohortResult objects
        output_path: Path to save report
    """
    lines = [
        "# Cross-Cohort Validation Report",
        "",
        "## Summary",
        "",
        "This report compares molecular subtype definitions across cohorts",
        "to assess reproducibility of discovered subtypes.",
        "",
        "## Validation Metrics",
        "",
        "| Cohort A | Cohort B | Transfer ARI | Projection ARI | Pathway Corr | Shared Subtypes |",
        "|----------|----------|--------------|----------------|--------------|-----------------|",
    ]

    for r in results:
        shared_str = ", ".join(r.shared_subtypes) if r.shared_subtypes else "None"
        lines.append(
            f"| {r.cohort_a} | {r.cohort_b} | {r.transfer_ari:.3f} | "
            f"{r.projection_ari:.3f} | {r.pathway_correlation:.3f} | {shared_str} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- **Transfer ARI > 0.5**: Good reproducibility of subtype structure",
            "- **Projection ARI > 0.5**: Consistent pathway-subtype relationships",
            "- **Pathway Correlation > 0.7**: Similar pathway importance across cohorts",
            "",
            "### Metric Definitions",
            "",
            "- **Transfer ARI**: Train GMM on cohort A, predict on cohort B, "
            "compare to cohort B's discovered labels",
            "- **Projection ARI**: Project both cohorts to shared PCA space, "
            "assign B samples to nearest A centroids",
            "- **Pathway Correlation**: Correlation of pathway variance between cohorts",
            "",
        ]
    )

    # Add detailed results
    lines.extend(
        [
            "## Detailed Results",
            "",
        ]
    )

    for i, r in enumerate(results, 1):
        lines.extend(
            [
                f"### {i}. {r.cohort_a} vs {r.cohort_b}",
                "",
                f"- Common pathways: {r.details.get('common_pathways', 'N/A')}",
                f"- Cohort A: {r.details.get('cohort_a_samples', 'N/A')} samples, "
                f"{r.details.get('cohort_a_clusters', 'N/A')} clusters",
                f"- Cohort B: {r.details.get('cohort_b_samples', 'N/A')} samples, "
                f"{r.details.get('cohort_b_clusters', 'N/A')} clusters",
                f"- Shared subtype labels: {', '.join(r.shared_subtypes) if r.shared_subtypes else 'None'}",
                "",
            ]
        )

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    logger.info(f"Saved cross-cohort report: {output_path}")


def batch_compare_cohorts(
    output_dirs: List[str],
    report_path: Optional[str] = None,
    seed: int = 42,
) -> List[CrossCohortResult]:
    """
    Compare all pairs of cohorts.

    Args:
        output_dirs: List of pipeline output directories
        report_path: Optional path to save report
        seed: Random seed

    Returns:
        List of CrossCohortResult for all pairs
    """
    # Load all cohorts
    cohorts = []
    for dir_path in output_dirs:
        try:
            cohort = load_cohort_result(dir_path)
            cohorts.append(cohort)
            logger.info(
                f"Loaded {cohort.name}: {cohort.n_samples} samples, {cohort.n_clusters} clusters"
            )
        except Exception as e:
            logger.warning(f"Could not load {dir_path}: {e}")

    if len(cohorts) < 2:
        raise ValueError("Need at least 2 cohorts for comparison")

    # Compare all pairs
    results = []
    for i, cohort_a in enumerate(cohorts):
        for cohort_b in cohorts[i + 1 :]:
            try:
                result = compare_cohorts(cohort_a, cohort_b, seed)
                results.append(result)
            except Exception as e:
                logger.warning(f"Could not compare {cohort_a.name} vs {cohort_b.name}: {e}")

    # Generate report if requested
    if report_path and results:
        generate_cross_cohort_report(results, report_path)

    return results
