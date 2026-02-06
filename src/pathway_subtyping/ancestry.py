"""
Ancestry and Population Stratification Correction Module.

Implements methods to detect and correct for population structure
in pathway-based subtype analysis:
- PCA computation from genotype data for ancestry inference
- Regression-based correction of pathway scores
- Ancestry independence testing for cluster validation
- Stratified analysis within ancestry groups

References:
- Price AL, et al. Principal components analysis corrects for
  stratification in genome-wide association studies.
  Nat Genet. 2006;38(8):904-909.
- Patterson N, et al. Population structure and eigenanalysis.
  PLoS Genet. 2006;2(12):e190.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


# =============================================================================
# ENUMS
# =============================================================================


class AncestryMethod(Enum):
    """
    Methods for ancestry correction in pathway scores.

    REGRESS_OUT: Remove ancestry-correlated variance by regressing
        pathway scores on ancestry PCs and taking residuals.
        Simple and widely used (Price et al., 2006).

    COVARIATE_AWARE: Include ancestry PCs as covariates during
        adjustment. Currently uses the same residualization as
        REGRESS_OUT; reserved for future multivariate mixed-model
        approaches.

    STRATIFIED: Perform clustering independently within each
        ancestry group, then compare subtype definitions across
        groups for concordance. Use with stratified_analysis().
    """

    REGRESS_OUT = "regress_out"
    COVARIATE_AWARE = "covariate_aware"
    STRATIFIED = "stratified"


# =============================================================================
# DATACLASSES
# =============================================================================


@dataclass
class AncestryPCs:
    """
    Principal components derived from genotype data for ancestry inference.

    Attributes:
        components: DataFrame of ancestry PCs (samples x n_components)
        explained_variance_ratio: Variance explained by each PC
        n_components: Number of PCs computed
        sample_ids: List of sample identifiers
    """

    components: pd.DataFrame
    explained_variance_ratio: np.ndarray
    n_components: int
    sample_ids: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_components": self.n_components,
            "n_samples": len(self.sample_ids),
            "explained_variance_ratio": [
                round(float(v), 4) for v in self.explained_variance_ratio
            ],
            "total_variance_explained": round(
                float(np.sum(self.explained_variance_ratio)), 4
            ),
        }


@dataclass
class AncestryAdjustmentResult:
    """
    Result of ancestry adjustment on pathway scores.

    Attributes:
        adjusted_scores: DataFrame of corrected pathway scores (samples x pathways)
        method: AncestryMethod used
        r_squared_per_pathway: R^2 of ancestry PCs for each pathway
        n_pcs_used: Number of ancestry PCs used in correction
        highly_confounded_pathways: Pathways with R^2 > threshold
    """

    adjusted_scores: pd.DataFrame
    method: AncestryMethod
    r_squared_per_pathway: Dict[str, float]
    n_pcs_used: int
    highly_confounded_pathways: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        r2_values = list(self.r_squared_per_pathway.values())
        return {
            "method": self.method.value,
            "n_pcs_used": self.n_pcs_used,
            "n_pathways_adjusted": len(self.r_squared_per_pathway),
            "n_highly_confounded": len(self.highly_confounded_pathways),
            "highly_confounded_pathways": self.highly_confounded_pathways,
            "mean_r_squared": round(float(np.mean(r2_values)), 4) if r2_values else 0.0,
            "max_r_squared": round(float(np.max(r2_values)), 4) if r2_values else 0.0,
            "r_squared_per_pathway": {
                k: round(v, 4) for k, v in self.r_squared_per_pathway.items()
            },
        }

    def format_report(self) -> str:
        """Format ancestry adjustment results for markdown report."""
        lines = [
            "## Ancestry Correction",
            "",
            f"**Method:** {self.method.value}",
            f"**Ancestry PCs used:** {self.n_pcs_used}",
            f"**Highly confounded pathways (R^2 > 0.1):** {len(self.highly_confounded_pathways)}",
            "",
        ]

        if self.highly_confounded_pathways:
            lines.extend([
                "### Pathways with Significant Ancestry Confounding",
                "",
                "| Pathway | R^2 (ancestry) |",
                "|---------|---------------|",
            ])
            for pathway in self.highly_confounded_pathways:
                r2 = self.r_squared_per_pathway[pathway]
                lines.append(f"| {pathway} | {r2:.4f} |")
            lines.append("")

        lines.extend([
            "### Interpretation",
            "",
            "R^2 indicates the proportion of pathway score variance explained by ancestry PCs.",
            "High R^2 values (> 0.1) suggest ancestry confounding that could create spurious subtypes.",
            "After correction, these pathways should have reduced ancestry-driven variance.",
        ])

        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return literature citations for the correction method."""
        return [
            "Price AL, et al. Principal components analysis corrects for "
            "stratification in genome-wide association studies. "
            "Nat Genet. 2006;38(8):904-909.",
            "Patterson N, et al. Population structure and eigenanalysis. "
            "PLoS Genet. 2006;2(12):e190.",
        ]


@dataclass
class AncestryStratificationReport:
    """
    Report on ancestry distribution relative to discovered clusters.

    Attributes:
        ancestry_cluster_correlation: Per-PC per-cluster mean values
        independence_test_pvalues: Per-PC Kruskal-Wallis p-values
        overall_independence_passed: Whether clusters are independent of ancestry
        significance_threshold: Alpha level used for testing
    """

    ancestry_cluster_correlation: Dict[str, Dict[str, float]]
    independence_test_pvalues: Dict[str, float]
    overall_independence_passed: bool
    significance_threshold: float = 0.05

    def to_dict(self) -> Dict[str, Any]:
        return {
            "overall_independence_passed": self.overall_independence_passed,
            "significance_threshold": self.significance_threshold,
            "independence_test_pvalues": {
                k: round(v, 6) for k, v in self.independence_test_pvalues.items()
            },
            "n_pcs_tested": len(self.independence_test_pvalues),
            "n_significant_pcs": sum(
                1 for v in self.independence_test_pvalues.values()
                if v < self.significance_threshold
            ),
        }

    def format_report(self) -> str:
        """Format ancestry stratification report for markdown."""
        status = "PASS" if self.overall_independence_passed else "FAIL"
        lines = [
            "## Ancestry Independence Test",
            "",
            f"**Status:** {status}",
            f"**Threshold:** p > {self.significance_threshold} for all PCs (Bonferroni-corrected)",
            "",
            "| Ancestry PC | Kruskal-Wallis p-value | Significant? |",
            "|-------------|----------------------|--------------|",
        ]

        for pc, pval in self.independence_test_pvalues.items():
            sig = "Yes" if pval < self.significance_threshold else "No"
            lines.append(f"| {pc} | {pval:.6f} | {sig} |")

        return "\n".join(lines)


# =============================================================================
# CORE FUNCTIONS
# =============================================================================


def compute_ancestry_pcs(
    genotype_matrix: pd.DataFrame,
    n_components: int = 10,
    seed: Optional[int] = None,
) -> AncestryPCs:
    """
    Compute principal components from genotype data for ancestry inference.

    Performs PCA on the genotype matrix (samples x variants) to extract
    top ancestry-informative principal components. These PCs capture
    population structure and can be used to correct pathway scores.

    Args:
        genotype_matrix: DataFrame of genotype dosages (samples x variants).
            Values should be 0, 1, or 2 (allele counts).
        n_components: Number of principal components to compute (default: 10).
            Price et al. (2006) recommend 10 PCs for most populations.
        seed: Random seed for reproducibility.

    Returns:
        AncestryPCs containing components and variance explained.

    Reference:
        Price AL, et al. Principal components analysis corrects for
        stratification in genome-wide association studies.
        Nat Genet. 2006;38(8):904-909.
    """
    logger.info(
        f"[Ancestry] Computing {n_components} ancestry PCs from "
        f"{genotype_matrix.shape[0]} samples x {genotype_matrix.shape[1]} variants"
    )

    # Cap n_components at min(n_samples, n_variants)
    max_components = min(genotype_matrix.shape[0], genotype_matrix.shape[1])
    n_components = min(n_components, max_components)

    if n_components < 1:
        raise ValueError(
            "Cannot compute ancestry PCs: need at least 1 sample and 1 variant"
        )

    # Standardize genotype matrix (mean-center, unit variance per variant)
    # Standard practice for genotype PCA (Patterson et al., 2006)
    geno_values = genotype_matrix.values.astype(float)
    means = np.nanmean(geno_values, axis=0)
    stds = np.nanstd(geno_values, axis=0)
    stds[stds == 0] = 1.0  # Avoid division by zero for monomorphic variants
    standardized = (geno_values - means) / stds

    # Replace any remaining NaN with 0 (missing genotypes)
    standardized = np.nan_to_num(standardized, nan=0.0)

    # Run PCA
    pca = PCA(n_components=n_components, random_state=seed)
    pc_values = pca.fit_transform(standardized)

    # Build DataFrame
    pc_columns = [f"PC{i+1}" for i in range(n_components)]
    components_df = pd.DataFrame(
        pc_values,
        index=genotype_matrix.index,
        columns=pc_columns,
    )

    logger.info(
        f"[Ancestry] Top {n_components} PCs explain "
        f"{np.sum(pca.explained_variance_ratio_):.1%} of genotype variance"
    )

    return AncestryPCs(
        components=components_df,
        explained_variance_ratio=pca.explained_variance_ratio_,
        n_components=n_components,
        sample_ids=list(genotype_matrix.index),
    )


def adjust_pathway_scores(
    pathway_scores: pd.DataFrame,
    ancestry_pcs: AncestryPCs,
    method: AncestryMethod = AncestryMethod.REGRESS_OUT,
    n_pcs: Optional[int] = None,
    confounding_threshold: float = 0.1,
) -> AncestryAdjustmentResult:
    """
    Adjust pathway scores for ancestry/population stratification.

    Removes ancestry-correlated variance from pathway scores to prevent
    spurious subtypes driven by population structure rather than biology.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways).
        ancestry_pcs: AncestryPCs from compute_ancestry_pcs().
        method: Correction method (default: REGRESS_OUT).
        n_pcs: Number of PCs to use (default: all available).
        confounding_threshold: R^2 threshold for flagging confounded pathways.

    Returns:
        AncestryAdjustmentResult with corrected scores and diagnostics.

    Raises:
        ValueError: If method is STRATIFIED (use stratified_analysis() instead).

    Reference:
        Price AL, et al. Nat Genet. 2006;38(8):904-909.
    """
    if method == AncestryMethod.STRATIFIED:
        raise ValueError(
            "AncestryMethod.STRATIFIED is not a score adjustment method. "
            "Use stratified_analysis() instead."
        )

    logger.info(f"[Ancestry] Adjusting pathway scores using method: {method.value}")

    # Align samples between pathway scores and ancestry PCs
    common_samples = sorted(
        set(pathway_scores.index) & set(ancestry_pcs.components.index)
    )

    if not common_samples:
        raise ValueError(
            "No common samples between pathway scores and ancestry PCs. "
            "Check that sample IDs match between data sources."
        )

    if len(common_samples) < len(pathway_scores):
        logger.warning(
            f"[Ancestry] {len(pathway_scores) - len(common_samples)} samples "
            f"missing ancestry PCs; using {len(common_samples)} common samples"
        )

    scores_aligned = pathway_scores.loc[common_samples]
    pcs_aligned = ancestry_pcs.components.loc[common_samples]

    # Select number of PCs
    if n_pcs is None:
        n_pcs = ancestry_pcs.n_components
    n_pcs = min(n_pcs, ancestry_pcs.n_components)
    pc_matrix = pcs_aligned.iloc[:, :n_pcs].values

    if method == AncestryMethod.REGRESS_OUT:
        adjusted, r2_per_pathway = _regress_out_ancestry(scores_aligned, pc_matrix)
    elif method == AncestryMethod.COVARIATE_AWARE:
        adjusted, r2_per_pathway = _covariate_aware_adjustment(scores_aligned, pc_matrix)
    else:
        raise ValueError(f"Unknown ancestry method: {method}")

    # Identify highly confounded pathways
    highly_confounded = [
        pathway for pathway, r2 in r2_per_pathway.items()
        if r2 > confounding_threshold
    ]

    if highly_confounded:
        logger.warning(
            f"[Ancestry] {len(highly_confounded)} pathway(s) have R^2 > {confounding_threshold} "
            f"with ancestry PCs: {highly_confounded[:5]}"
        )

    return AncestryAdjustmentResult(
        adjusted_scores=adjusted,
        method=method,
        r_squared_per_pathway=r2_per_pathway,
        n_pcs_used=n_pcs,
        highly_confounded_pathways=highly_confounded,
    )


def check_ancestry_independence(
    cluster_labels: np.ndarray,
    ancestry_pcs: AncestryPCs,
    significance_threshold: float = 0.05,
    n_pcs_to_test: Optional[int] = None,
) -> AncestryStratificationReport:
    """
    Check whether discovered clusters are independent of ancestry PCs.

    For each ancestry PC, runs a Kruskal-Wallis H-test to determine if
    PC values differ significantly across clusters. If clusters correlate
    with ancestry, the subtypes may reflect population structure rather
    than biological differences.

    Uses Bonferroni correction across PCs tested.

    Args:
        cluster_labels: Array of cluster assignments.
        ancestry_pcs: AncestryPCs from compute_ancestry_pcs().
        significance_threshold: P-value threshold before Bonferroni (default 0.05).
        n_pcs_to_test: Number of top PCs to test (default: min(5, n_components)).

    Returns:
        AncestryStratificationReport with independence test results.
    """
    logger.info("[Ancestry] Testing cluster independence from ancestry PCs...")

    if n_pcs_to_test is None:
        n_pcs_to_test = min(5, ancestry_pcs.n_components)

    pvalues = {}
    correlations = {}
    unique_clusters = np.unique(cluster_labels)

    for i in range(n_pcs_to_test):
        pc_name = f"PC{i+1}"
        pc_values = ancestry_pcs.components.iloc[:, i].values

        # Kruskal-Wallis H-test: are PC values different across clusters?
        groups = [pc_values[cluster_labels == c] for c in unique_clusters]
        groups = [g for g in groups if len(g) > 0]

        if len(groups) < 2:
            pvalues[pc_name] = 1.0
            correlations[pc_name] = {}
            continue

        stat, pval = stats.kruskal(*groups)
        pvalues[pc_name] = float(pval)

        # Per-cluster mean for reporting
        cluster_means = {}
        for c in unique_clusters:
            mask = cluster_labels == c
            if np.any(mask):
                cluster_means[str(c)] = round(float(np.mean(pc_values[mask])), 4)
        correlations[pc_name] = cluster_means

    # Bonferroni correction: pass if NO PC is significant after correction
    n_tests = len(pvalues)
    corrected_threshold = significance_threshold / n_tests if n_tests > 0 else significance_threshold
    overall_passed = all(p > corrected_threshold for p in pvalues.values())

    status = "PASS" if overall_passed else "FAIL"
    logger.info(
        f"[Ancestry] Independence test: {status} "
        f"(tested {n_pcs_to_test} PCs, Bonferroni threshold={corrected_threshold:.4f})"
    )

    return AncestryStratificationReport(
        ancestry_cluster_correlation=correlations,
        independence_test_pvalues=pvalues,
        overall_independence_passed=overall_passed,
        significance_threshold=significance_threshold,
    )


def stratified_analysis(
    pathway_scores: pd.DataFrame,
    ancestry_groups: np.ndarray,
    n_clusters: int,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Perform clustering independently within each ancestry group.

    Clusters samples within each group separately, then evaluates
    cross-group concordance to determine if subtypes replicate.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways).
        ancestry_groups: Array of ancestry group labels per sample.
        n_clusters: Number of clusters to fit per group.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with per-group results and cross-group concordance.
    """
    logger.info("[Ancestry] Running stratified analysis across ancestry groups...")

    unique_groups = np.unique(ancestry_groups)
    per_group_results = {}
    all_labels = {}

    for group in unique_groups:
        mask = ancestry_groups == group
        group_scores = pathway_scores.loc[mask]

        min_samples = n_clusters * 3
        if len(group_scores) < min_samples:
            logger.warning(
                f"[Ancestry] Group '{group}' has only {len(group_scores)} samples, "
                f"skipping (need at least {min_samples})"
            )
            continue

        gmm = GaussianMixture(
            n_components=n_clusters,
            covariance_type="full",
            n_init=10,
            random_state=seed,
            reg_covar=1e-6,
        )
        gmm.fit(group_scores.values)
        labels = gmm.predict(group_scores.values)

        per_group_results[str(group)] = {
            "n_samples": len(group_scores),
            "converged": bool(gmm.converged_),
            "bic": float(gmm.bic(group_scores.values)),
        }
        all_labels[str(group)] = (group_scores.index.tolist(), labels)

    # Cross-group concordance
    concordance = _compute_cross_group_concordance(
        pathway_scores, all_labels
    )

    return {
        "n_groups": len(unique_groups),
        "groups_analyzed": len(per_group_results),
        "per_group_results": per_group_results,
        "cross_group_concordance": concordance,
    }


def compute_ancestry_correlation(
    pathway_scores: pd.DataFrame,
    ancestry_pcs: AncestryPCs,
    n_pcs: Optional[int] = None,
) -> pd.DataFrame:
    """
    Compute Pearson correlation between each pathway score and each ancestry PC.

    Useful for identifying which pathways are most confounded by population
    structure before and after correction.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways).
        ancestry_pcs: AncestryPCs from compute_ancestry_pcs().
        n_pcs: Number of PCs to correlate (default: min(5, n_components)).

    Returns:
        DataFrame of correlations (pathways x PCs).
    """
    if n_pcs is None:
        n_pcs = min(5, ancestry_pcs.n_components)

    common_samples = sorted(
        set(pathway_scores.index) & set(ancestry_pcs.components.index)
    )

    if len(common_samples) < 3:
        raise ValueError(
            "Need at least 3 common samples to compute correlations"
        )

    scores_aligned = pathway_scores.loc[common_samples]
    pcs_aligned = ancestry_pcs.components.loc[common_samples].iloc[:, :n_pcs]

    correlations = {}
    for pathway in scores_aligned.columns:
        pathway_corr = {}
        for pc_col in pcs_aligned.columns:
            r, _ = stats.pearsonr(
                scores_aligned[pathway].values,
                pcs_aligned[pc_col].values,
            )
            pathway_corr[pc_col] = round(r, 4)
        correlations[pathway] = pathway_corr

    return pd.DataFrame(correlations).T


# =============================================================================
# PRIVATE HELPERS
# =============================================================================


def _regress_out_ancestry(
    pathway_scores: pd.DataFrame,
    pc_matrix: np.ndarray,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Regress out ancestry PCs from pathway scores.

    For each pathway, fits: score ~ PC1 + PC2 + ... + PCn
    and returns the residuals as corrected scores.
    """
    adjusted_data = {}
    r2_per_pathway = {}

    for pathway in pathway_scores.columns:
        y = pathway_scores[pathway].values
        model = LinearRegression()
        model.fit(pc_matrix, y)
        predicted = model.predict(pc_matrix)
        residuals = y - predicted

        # R^2 = proportion of variance explained by ancestry
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        adjusted_data[pathway] = residuals
        r2_per_pathway[pathway] = max(0.0, r_squared)

    adjusted_df = pd.DataFrame(adjusted_data, index=pathway_scores.index)

    # Re-normalize (z-score) after adjustment
    means = adjusted_df.mean()
    stds = adjusted_df.std()
    stds = stds.replace(0, 1e-10)
    adjusted_df = (adjusted_df - means) / stds

    return adjusted_df, r2_per_pathway


def _covariate_aware_adjustment(
    pathway_scores: pd.DataFrame,
    pc_matrix: np.ndarray,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Covariate-aware adjustment for pathway scores.

    Currently uses the same residualization as REGRESS_OUT.
    Reserved for future multivariate mixed-model approaches.
    """
    return _regress_out_ancestry(pathway_scores, pc_matrix)


def _compute_cross_group_concordance(
    pathway_scores: pd.DataFrame,
    all_labels: Dict[str, Tuple[List[str], np.ndarray]],
) -> Dict[str, Any]:
    """Compute concordance of cluster centroids across ancestry groups."""
    if len(all_labels) < 2:
        return {"concordance_score": 0.0, "sufficient_groups": False}

    group_names = list(all_labels.keys())
    centroids_per_group = {}

    for group_name in group_names:
        sample_ids, labels = all_labels[group_name]
        group_scores = pathway_scores.loc[sample_ids]
        centroids = []
        for c in np.unique(labels):
            mask = labels == c
            centroids.append(group_scores.iloc[mask].mean().values)
        centroids_per_group[group_name] = np.array(centroids)

    # Pairwise centroid correlations between groups
    correlations = []
    for i, g1 in enumerate(group_names):
        for g2 in group_names[i + 1:]:
            c1 = centroids_per_group[g1]
            c2 = centroids_per_group[g2]
            if len(c1) > 0 and len(c2) > 0:
                for centroid_a in c1:
                    best_corr = max(
                        float(np.corrcoef(centroid_a, centroid_b)[0, 1])
                        for centroid_b in c2
                    )
                    correlations.append(best_corr)

    mean_concordance = float(np.mean(correlations)) if correlations else 0.0

    return {
        "concordance_score": round(mean_concordance, 4),
        "sufficient_groups": True,
        "n_pairwise_comparisons": len(correlations),
    }
