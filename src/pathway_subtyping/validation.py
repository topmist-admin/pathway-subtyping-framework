"""
Validation Gates for the Pathway Subtyping Framework.

Implements mandatory negative controls and stability tests to validate
clustering results before accepting them as meaningful.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

from .utils.seed import get_rng

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Result from a single validation test."""

    name: str
    passed: bool
    metric_name: str
    metric_value: float
    threshold: float
    comparison: str  # ">" or "<" or "=="
    details: Dict[str, Any] = field(default_factory=dict)

    @property
    def status(self) -> str:
        return "PASS" if self.passed else "FAIL"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "metric": self.metric_name,
            "value": round(self.metric_value, 4),
            "threshold": self.threshold,
            "comparison": self.comparison,
            "details": self.details,
        }


@dataclass
class ValidationGatesResult:
    """Aggregated result from all validation gates."""

    results: List[ValidationResult]
    all_passed: bool
    summary: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "all_passed": self.all_passed,
            "summary": self.summary,
            "tests": [r.to_dict() for r in self.results],
        }


class ValidationGates:
    """
    Validation gates for clustering results.

    Implements:
    - Negative Control 1: Label shuffle (expect ARI ~0)
    - Negative Control 2: Random gene sets (expect no enrichment)
    - Stability Test: Bootstrap resampling (expect ARI >= threshold)
    """

    def __init__(
        self,
        seed: Optional[int] = 42,
        n_permutations: int = 100,
        n_bootstrap: int = 50,
        stability_threshold: float = 0.8,
        null_ari_max: float = 0.15,
    ):
        """
        Initialize validation gates.

        Args:
            seed: Random seed for reproducibility
            n_permutations: Number of permutations for null tests
            n_bootstrap: Number of bootstrap iterations for stability
            stability_threshold: Minimum ARI for stability test (default: 0.8)
            null_ari_max: Maximum expected ARI under null hypothesis (default: 0.15)
        """
        self.seed = seed
        self.n_permutations = n_permutations
        self.n_bootstrap = n_bootstrap
        self.stability_threshold = stability_threshold
        self.null_ari_max = null_ari_max
        self.rng = get_rng(seed, "validation")

    def run_all(
        self,
        pathway_scores: pd.DataFrame,
        cluster_labels: np.ndarray,
        pathways: Dict[str, List[str]],
        gene_burdens: pd.DataFrame,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationGatesResult:
        """
        Run all validation gates.

        Args:
            pathway_scores: DataFrame of pathway scores (samples x pathways)
            cluster_labels: Array of cluster assignments
            pathways: Dict mapping pathway names to gene lists
            gene_burdens: DataFrame of gene burdens (samples x genes)
            n_clusters: Number of clusters used
            gmm_seed: Seed for GMM clustering

        Returns:
            ValidationGatesResult with all test outcomes
        """
        logger.info("Running validation gates...")

        results = []

        # Negative Control 1: Label shuffle
        nc1_result = self.negative_control_label_shuffle(
            pathway_scores, cluster_labels, n_clusters, gmm_seed
        )
        results.append(nc1_result)
        logger.info(f"  - {nc1_result.name}: {nc1_result.status} (ARI = {nc1_result.metric_value:.3f})")

        # Negative Control 2: Random gene sets
        nc2_result = self.negative_control_random_gene_sets(
            gene_burdens, pathways, cluster_labels, n_clusters, gmm_seed
        )
        results.append(nc2_result)
        logger.info(f"  - {nc2_result.name}: {nc2_result.status} (ARI = {nc2_result.metric_value:.3f})")

        # Stability Test: Bootstrap
        stability_result = self.stability_test_bootstrap(
            pathway_scores, cluster_labels, n_clusters, gmm_seed
        )
        results.append(stability_result)
        logger.info(f"  - {stability_result.name}: {stability_result.status} (ARI = {stability_result.metric_value:.3f})")

        # Aggregate results
        all_passed = all(r.passed for r in results)
        n_passed = sum(1 for r in results if r.passed)
        n_total = len(results)

        if all_passed:
            summary = f"All {n_total} validation gates PASSED"
        else:
            failed = [r.name for r in results if not r.passed]
            summary = f"{n_passed}/{n_total} validation gates passed. FAILED: {', '.join(failed)}"

        return ValidationGatesResult(
            results=results,
            all_passed=all_passed,
            summary=summary,
        )

    def negative_control_label_shuffle(
        self,
        pathway_scores: pd.DataFrame,
        original_labels: np.ndarray,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationResult:
        """
        Negative Control 1: Label shuffle test.

        Shuffles cluster labels randomly and re-runs clustering. The clustering
        should NOT recover the shuffled labels (ARI should be ~0).

        If clustering CAN recover shuffled labels, it suggests the clusters
        are an artifact of the method rather than real structure in the data.
        """
        ari_values = []

        for i in range(self.n_permutations):
            # Shuffle labels
            shuffled_labels = self.rng.permutation(original_labels)

            # Re-cluster
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=(gmm_seed + i) if gmm_seed else None,
            )
            gmm.fit(pathway_scores.values)
            new_labels = gmm.predict(pathway_scores.values)

            # Compare to shuffled (should NOT match)
            ari = adjusted_rand_score(shuffled_labels, new_labels)
            ari_values.append(ari)

        mean_ari = np.mean(ari_values)
        std_ari = np.std(ari_values)

        # PASS if mean ARI is close to 0 (clustering doesn't recover random labels)
        passed = mean_ari < self.null_ari_max

        return ValidationResult(
            name="Negative Control 1: Label Shuffle",
            passed=passed,
            metric_name="mean_null_ARI",
            metric_value=mean_ari,
            threshold=self.null_ari_max,
            comparison="<",
            details={
                "std_ari": round(std_ari, 4),
                "n_permutations": self.n_permutations,
                "interpretation": "Clustering should NOT recover shuffled labels",
            },
        )

    def negative_control_random_gene_sets(
        self,
        gene_burdens: pd.DataFrame,
        real_pathways: Dict[str, List[str]],
        original_labels: np.ndarray,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationResult:
        """
        Negative Control 2: Random gene sets test.

        Replaces biological pathway definitions with random gene sets of the
        same sizes. Clustering on random gene sets should NOT recover the
        original cluster structure (ARI should be ~0).

        If clustering CAN recover structure with random gene sets, it suggests
        the results are not driven by meaningful biological pathways.
        """
        all_genes = list(gene_burdens.columns)
        ari_values = []

        for i in range(self.n_permutations):
            # Create random pathways with same sizes
            random_pathway_scores = {}

            for pathway_name, pathway_genes in real_pathways.items():
                # How many genes from this pathway are in our burden data?
                n_genes_needed = min(len(pathway_genes), len(all_genes))
                if n_genes_needed < 2:
                    continue

                # Sample random genes
                random_genes = list(self.rng.choice(all_genes, size=n_genes_needed, replace=False))

                # Compute score for random pathway
                common_genes = [g for g in random_genes if g in gene_burdens.columns]
                if len(common_genes) >= 2:
                    random_pathway_scores[pathway_name] = gene_burdens[common_genes].mean(axis=1)

            if len(random_pathway_scores) < 2:
                continue

            # Build random pathway score matrix
            random_scores_df = pd.DataFrame(random_pathway_scores)
            random_scores_df = (random_scores_df - random_scores_df.mean()) / random_scores_df.std()
            random_scores_df = random_scores_df.fillna(0)

            # Cluster on random pathway scores
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=(gmm_seed + i + 1000) if gmm_seed else None,
            )
            gmm.fit(random_scores_df.values)
            random_labels = gmm.predict(random_scores_df.values)

            # Compare to original clustering
            ari = adjusted_rand_score(original_labels, random_labels)
            ari_values.append(ari)

        if not ari_values:
            # Edge case: couldn't run any permutations
            return ValidationResult(
                name="Negative Control 2: Random Gene Sets",
                passed=False,
                metric_name="mean_random_ARI",
                metric_value=1.0,
                threshold=self.null_ari_max,
                comparison="<",
                details={"error": "Could not compute random gene set permutations"},
            )

        mean_ari = np.mean(ari_values)
        std_ari = np.std(ari_values)

        # PASS if mean ARI is close to 0 (random gene sets don't replicate structure)
        passed = mean_ari < self.null_ari_max

        return ValidationResult(
            name="Negative Control 2: Random Gene Sets",
            passed=passed,
            metric_name="mean_random_ARI",
            metric_value=mean_ari,
            threshold=self.null_ari_max,
            comparison="<",
            details={
                "std_ari": round(std_ari, 4),
                "n_permutations": len(ari_values),
                "interpretation": "Random gene sets should NOT replicate cluster structure",
            },
        )

    def stability_test_bootstrap(
        self,
        pathway_scores: pd.DataFrame,
        original_labels: np.ndarray,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationResult:
        """
        Stability Test: Bootstrap resampling.

        Resamples the data with replacement and re-runs clustering. Stable
        clusters should be recovered across bootstrap samples (high ARI).

        If clustering is unstable across bootstrap samples, it suggests the
        clusters are not robust features of the data.
        """
        n_samples = len(pathway_scores)
        ari_values = []

        for i in range(self.n_bootstrap):
            # Bootstrap sample (with replacement)
            bootstrap_idx = self.rng.choice(n_samples, size=n_samples, replace=True)

            # Get unique indices for clustering (avoid duplicate samples confusing GMM)
            unique_idx = np.unique(bootstrap_idx)
            if len(unique_idx) < n_clusters * 2:
                continue  # Not enough unique samples

            bootstrap_scores = pathway_scores.iloc[unique_idx]
            bootstrap_original_labels = original_labels[unique_idx]

            # Re-cluster on bootstrap sample
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=(gmm_seed + i + 2000) if gmm_seed else None,
            )
            try:
                gmm.fit(bootstrap_scores.values)
                bootstrap_labels = gmm.predict(bootstrap_scores.values)

                # Compare to original labels for these samples
                ari = adjusted_rand_score(bootstrap_original_labels, bootstrap_labels)
                ari_values.append(ari)
            except Exception:
                continue  # Skip failed fits

        if not ari_values:
            return ValidationResult(
                name="Stability Test: Bootstrap",
                passed=False,
                metric_name="mean_bootstrap_ARI",
                metric_value=0.0,
                threshold=self.stability_threshold,
                comparison=">=",
                details={"error": "Could not compute bootstrap stability"},
            )

        mean_ari = np.mean(ari_values)
        std_ari = np.std(ari_values)

        # PASS if mean ARI >= threshold (clustering is stable)
        passed = mean_ari >= self.stability_threshold

        return ValidationResult(
            name="Stability Test: Bootstrap",
            passed=passed,
            metric_name="mean_bootstrap_ARI",
            metric_value=mean_ari,
            threshold=self.stability_threshold,
            comparison=">=",
            details={
                "std_ari": round(std_ari, 4),
                "n_bootstrap": len(ari_values),
                "interpretation": "Clustering should be stable across bootstrap samples",
            },
        )


def format_validation_report(result: ValidationGatesResult) -> str:
    """Format validation results for display in report."""
    lines = [
        "## Validation Gates",
        "",
        f"**Overall Status:** {'PASS' if result.all_passed else 'FAIL'}",
        "",
        f"{result.summary}",
        "",
        "| Test | Status | Metric | Value | Threshold |",
        "|------|--------|--------|-------|-----------|",
    ]

    for test in result.results:
        status_icon = "✓" if test.passed else "✗"
        lines.append(
            f"| {test.name} | {status_icon} {test.status} | {test.metric_name} | "
            f"{test.metric_value:.3f} | {test.comparison} {test.threshold} |"
        )

    lines.extend([
        "",
        "### Interpretation",
        "",
        "- **Negative Control 1 (Label Shuffle)**: Tests if clustering can recover "
        "randomly shuffled labels. PASS means clustering does NOT find spurious patterns.",
        "",
        "- **Negative Control 2 (Random Gene Sets)**: Tests if clusters are driven by "
        "biological pathways vs. random gene groupings. PASS means biological pathways matter.",
        "",
        "- **Stability Test (Bootstrap)**: Tests if clusters are robust to resampling. "
        "PASS means clusters are stable features of the data.",
        "",
    ])

    return "\n".join(lines)
