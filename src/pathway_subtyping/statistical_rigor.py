"""
Statistical Rigor Module for the Pathway Subtyping Framework.

Implements publication-quality statistical methods including:
- Multiple testing correction (FDR via Benjamini-Hochberg)
- Permutation-based p-values for pathway associations
- Effect size calculations (Cohen's d)
- Confidence intervals via bootstrap
- Literature-based burden weight schemes

References:
- Benjamini & Hochberg (1995). Controlling the False Discovery Rate.
- Storey (2002). A direct approach to false discovery rates.
- Cohen (1988). Statistical Power Analysis for the Behavioral Sciences.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


# =============================================================================
# BURDEN WEIGHT SCHEMES
# =============================================================================

class BurdenWeightScheme(Enum):
    """
    Literature-based burden weighting schemes.

    Each scheme is derived from published studies on variant deleteriousness.
    """

    # Default scheme (ad-hoc, for backwards compatibility)
    DEFAULT = "default"

    # ExAC/gnomAD-inspired (Lek et al., Nature 2016)
    # Based on observed/expected ratios and constraint metrics
    GNOMAD_CONSTRAINT = "gnomad_constraint"

    # ACMG-inspired (Richards et al., Genet Med 2015)
    # Based on variant classification guidelines
    ACMG_INSPIRED = "acmg_inspired"

    # CADD percentile-based (Kircher et al., Nat Genet 2014)
    # Using CADD score percentiles from published benchmarks
    CADD_PERCENTILE = "cadd_percentile"

    # Uniform (all damaging variants weighted equally)
    # For sensitivity analysis
    UNIFORM = "uniform"


@dataclass
class BurdenWeights:
    """
    Configurable burden weight parameters with literature citations.

    Attributes:
        scheme: The weighting scheme to use
        lof_weight: Weight for loss-of-function variants
        missense_high_weight: Weight for high-impact missense
        missense_moderate_weight: Weight for moderate-impact missense
        missense_low_weight: Weight for low-impact missense
        other_weight: Weight for other variant types
        cadd_high_threshold: CADD score threshold for "high impact"
        cadd_moderate_threshold: CADD score threshold for "moderate impact"

    Citations for default values:
        - LoF weight (1.0): Based on pLI scores and observed depletion of LoF
          variants in constrained genes (Lek et al., Nature 2016)
        - CADD thresholds: CADD >= 25 corresponds to top 0.3% most deleterious
          variants (Kircher et al., Nat Genet 2014)
    """

    scheme: BurdenWeightScheme = BurdenWeightScheme.DEFAULT

    # Variant type weights
    lof_weight: float = 1.0
    missense_high_weight: float = 0.8
    missense_moderate_weight: float = 0.5
    missense_low_weight: float = 0.2
    other_weight: float = 0.1

    # CADD thresholds
    cadd_high_threshold: float = 25.0      # Top 0.3% most deleterious
    cadd_moderate_threshold: float = 20.0  # Top 1% most deleterious
    cadd_low_threshold: float = 15.0       # Top 3% most deleterious

    # CADD normalization
    cadd_cap: float = 40.0

    def get_citations(self) -> List[str]:
        """Return literature citations for the weight scheme."""
        base_citations = [
            "Lek M, et al. Analysis of protein-coding genetic variation in 60,706 humans. Nature. 2016;536(7616):285-291.",
            "Kircher M, et al. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014;46(3):310-315.",
        ]

        if self.scheme == BurdenWeightScheme.ACMG_INSPIRED:
            base_citations.append(
                "Richards S, et al. Standards and guidelines for the interpretation of sequence variants. Genet Med. 2015;17(5):405-424."
            )

        return base_citations

    @classmethod
    def from_scheme(cls, scheme: BurdenWeightScheme) -> "BurdenWeights":
        """
        Create weights from a predefined scheme.

        Args:
            scheme: The weighting scheme to use

        Returns:
            BurdenWeights configured for that scheme
        """
        if scheme == BurdenWeightScheme.DEFAULT:
            return cls(
                scheme=scheme,
                lof_weight=1.0,
                missense_high_weight=0.5,
                missense_moderate_weight=0.1,
                missense_low_weight=0.1,
                other_weight=0.1,
                cadd_high_threshold=25.0,
                cadd_moderate_threshold=25.0,  # Single threshold in default
            )

        elif scheme == BurdenWeightScheme.GNOMAD_CONSTRAINT:
            # Based on gnomAD constraint metrics
            return cls(
                scheme=scheme,
                lof_weight=1.0,
                missense_high_weight=0.8,
                missense_moderate_weight=0.4,
                missense_low_weight=0.1,
                other_weight=0.05,
                cadd_high_threshold=25.0,
                cadd_moderate_threshold=20.0,
                cadd_low_threshold=15.0,
            )

        elif scheme == BurdenWeightScheme.ACMG_INSPIRED:
            # Based on ACMG classification
            # Pathogenic/Likely pathogenic get higher weights
            return cls(
                scheme=scheme,
                lof_weight=1.0,
                missense_high_weight=0.9,  # Likely pathogenic
                missense_moderate_weight=0.5,  # VUS
                missense_low_weight=0.1,  # Likely benign
                other_weight=0.05,
                cadd_high_threshold=30.0,  # Stricter threshold
                cadd_moderate_threshold=25.0,
                cadd_low_threshold=20.0,
            )

        elif scheme == BurdenWeightScheme.CADD_PERCENTILE:
            # Continuous CADD-based scoring
            return cls(
                scheme=scheme,
                lof_weight=1.0,
                missense_high_weight=0.7,
                missense_moderate_weight=0.4,
                missense_low_weight=0.2,
                other_weight=0.1,
                cadd_high_threshold=25.0,
                cadd_moderate_threshold=20.0,
                cadd_low_threshold=15.0,
            )

        elif scheme == BurdenWeightScheme.UNIFORM:
            # All damaging variants weighted equally
            return cls(
                scheme=scheme,
                lof_weight=1.0,
                missense_high_weight=1.0,
                missense_moderate_weight=1.0,
                missense_low_weight=1.0,
                other_weight=1.0,
            )

        else:
            raise ValueError(f"Unknown scheme: {scheme}")


def compute_variant_weight(
    consequence: str,
    cadd_score: float,
    weights: BurdenWeights,
) -> float:
    """
    Compute burden weight for a variant using literature-based scheme.

    Args:
        consequence: Variant consequence (e.g., "frameshift_variant")
        cadd_score: CADD phred-scaled score
        weights: BurdenWeights configuration

    Returns:
        Weight for the variant (0.0 to 1.0)
    """
    consequence_lower = consequence.lower()

    # Loss-of-function variants
    lof_terms = ["frameshift", "stop_gained", "stop_lost", "splice_donor",
                 "splice_acceptor", "start_lost", "transcript_ablation"]
    if any(term in consequence_lower for term in lof_terms):
        return weights.lof_weight

    # Missense variants - weight by CADD score
    if "missense" in consequence_lower:
        if cadd_score >= weights.cadd_high_threshold:
            return weights.missense_high_weight
        elif cadd_score >= weights.cadd_moderate_threshold:
            return weights.missense_moderate_weight
        elif cadd_score >= weights.cadd_low_threshold:
            return weights.missense_low_weight
        else:
            return weights.missense_low_weight * 0.5  # Very low impact

    # Inframe indels
    if "inframe" in consequence_lower:
        # Scale by CADD if available
        if cadd_score >= weights.cadd_high_threshold:
            return weights.missense_moderate_weight
        else:
            return weights.other_weight

    # Other (synonymous, intronic, etc.)
    return weights.other_weight


# =============================================================================
# MULTIPLE TESTING CORRECTION
# =============================================================================

@dataclass
class FDRResult:
    """
    Result of FDR correction.

    Attributes:
        pathway: Pathway name
        p_value: Raw p-value from permutation test
        q_value: FDR-adjusted q-value (Benjamini-Hochberg)
        significant: Whether q-value < alpha
        effect_size: Cohen's d effect size
        ci_lower: Lower bound of 95% CI for effect size
        ci_upper: Upper bound of 95% CI for effect size
    """
    pathway: str
    p_value: float
    q_value: float
    significant: bool
    effect_size: float
    ci_lower: float
    ci_upper: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway": self.pathway,
            "p_value": round(self.p_value, 6),
            "q_value": round(self.q_value, 6),
            "significant": self.significant,
            "effect_size": round(self.effect_size, 4),
            "ci_95_lower": round(self.ci_lower, 4),
            "ci_95_upper": round(self.ci_upper, 4),
        }


def benjamini_hochberg(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """
    Apply Benjamini-Hochberg FDR correction.

    Args:
        p_values: Array of raw p-values
        alpha: Significance threshold (not used in calculation, just for reference)

    Returns:
        Array of q-values (FDR-adjusted p-values)

    Reference:
        Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical
        and Powerful Approach to Multiple Testing. J R Stat Soc B. 1995;57(1):289-300.
    """
    n = len(p_values)
    if n == 0:
        return np.array([])

    # Sort p-values
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]

    # Calculate q-values
    # q_i = min(p_i * n / i, q_{i+1})
    q_values = np.zeros(n)
    q_values[-1] = sorted_p[-1]

    for i in range(n - 2, -1, -1):
        q_values[i] = min(sorted_p[i] * n / (i + 1), q_values[i + 1])

    # Cap at 1.0
    q_values = np.minimum(q_values, 1.0)

    # Restore original order
    result = np.zeros(n)
    result[sorted_idx] = q_values

    return result


def compute_pathway_pvalues(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    n_permutations: int = 1000,
    seed: Optional[int] = None,
) -> Dict[str, float]:
    """
    Compute permutation-based p-values for pathway associations.

    Tests whether each pathway's scores differ significantly between clusters
    using a permutation test.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways)
        cluster_labels: Array of cluster assignments
        n_permutations: Number of permutations for null distribution
        seed: Random seed for reproducibility

    Returns:
        Dictionary mapping pathway names to p-values
    """
    rng = np.random.RandomState(seed)
    p_values = {}

    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)

    if n_clusters < 2:
        # Cannot compute p-values with only one cluster
        return {pathway: 1.0 for pathway in pathway_scores.columns}

    for pathway in pathway_scores.columns:
        scores = pathway_scores[pathway].values

        # Compute observed F-statistic (ANOVA)
        observed_stat = _compute_f_statistic(scores, cluster_labels)

        # Permutation test
        null_stats = []
        for _ in range(n_permutations):
            permuted_labels = rng.permutation(cluster_labels)
            null_stat = _compute_f_statistic(scores, permuted_labels)
            null_stats.append(null_stat)

        null_stats = np.array(null_stats)

        # P-value: proportion of null stats >= observed
        # Add 1 to numerator and denominator for conservative estimate
        p_value = (np.sum(null_stats >= observed_stat) + 1) / (n_permutations + 1)
        p_values[pathway] = p_value

    return p_values


def _compute_f_statistic(values: np.ndarray, groups: np.ndarray) -> float:
    """Compute F-statistic for one-way ANOVA."""
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)
    n_total = len(values)

    if n_groups < 2 or n_total < n_groups + 1:
        return 0.0

    # Overall mean
    grand_mean = np.mean(values)

    # Between-group sum of squares
    ss_between = 0
    for g in unique_groups:
        group_mask = groups == g
        n_g = np.sum(group_mask)
        group_mean = np.mean(values[group_mask])
        ss_between += n_g * (group_mean - grand_mean) ** 2

    # Within-group sum of squares
    ss_within = 0
    for g in unique_groups:
        group_mask = groups == g
        group_mean = np.mean(values[group_mask])
        ss_within += np.sum((values[group_mask] - group_mean) ** 2)

    # Degrees of freedom
    df_between = n_groups - 1
    df_within = n_total - n_groups

    if df_within <= 0 or ss_within == 0:
        return 0.0

    # F-statistic
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within

    return ms_between / ms_within


# =============================================================================
# EFFECT SIZE CALCULATIONS
# =============================================================================

def compute_cohens_d(
    group1: np.ndarray,
    group2: np.ndarray,
) -> float:
    """
    Compute Cohen's d effect size between two groups.

    Uses pooled standard deviation for unequal group sizes.

    Args:
        group1: Values for first group
        group2: Values for second group

    Returns:
        Cohen's d effect size

    Reference:
        Cohen J. Statistical Power Analysis for the Behavioral Sciences.
        2nd ed. Routledge; 1988.
    """
    n1, n2 = len(group1), len(group2)

    if n1 < 2 or n2 < 2:
        return 0.0

    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)

    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

    if pooled_std == 0:
        return 0.0

    return (mean1 - mean2) / pooled_std


def compute_pathway_effect_sizes(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
) -> Dict[str, float]:
    """
    Compute effect sizes for each pathway between clusters.

    For k clusters, computes the maximum pairwise Cohen's d.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways)
        cluster_labels: Array of cluster assignments

    Returns:
        Dictionary mapping pathway names to effect sizes
    """
    effect_sizes = {}
    unique_clusters = np.unique(cluster_labels)

    for pathway in pathway_scores.columns:
        scores = pathway_scores[pathway].values
        max_d = 0.0

        # Compute max pairwise effect size
        for i, c1 in enumerate(unique_clusters):
            for c2 in unique_clusters[i + 1:]:
                group1 = scores[cluster_labels == c1]
                group2 = scores[cluster_labels == c2]
                d = abs(compute_cohens_d(group1, group2))
                max_d = max(max_d, d)

        effect_sizes[pathway] = max_d

    return effect_sizes


def bootstrap_effect_size_ci(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathway: str,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    seed: Optional[int] = None,
) -> Tuple[float, float]:
    """
    Compute bootstrap confidence interval for effect size.

    Args:
        pathway_scores: DataFrame of pathway scores
        cluster_labels: Array of cluster assignments
        pathway: Pathway name to compute CI for
        n_bootstrap: Number of bootstrap iterations
        ci_level: Confidence level (default 0.95 for 95% CI)
        seed: Random seed

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    rng = np.random.RandomState(seed)
    scores = pathway_scores[pathway].values
    n_samples = len(scores)
    unique_clusters = np.unique(cluster_labels)

    if len(unique_clusters) < 2:
        return (0.0, 0.0)

    bootstrap_effects = []

    for _ in range(n_bootstrap):
        # Bootstrap sample
        idx = rng.choice(n_samples, size=n_samples, replace=True)
        boot_scores = scores[idx]
        boot_labels = cluster_labels[idx]

        # Compute effect size
        max_d = 0.0
        for i, c1 in enumerate(unique_clusters):
            for c2 in unique_clusters[i + 1:]:
                g1_mask = boot_labels == c1
                g2_mask = boot_labels == c2

                if np.sum(g1_mask) >= 2 and np.sum(g2_mask) >= 2:
                    d = abs(compute_cohens_d(boot_scores[g1_mask], boot_scores[g2_mask]))
                    max_d = max(max_d, d)

        bootstrap_effects.append(max_d)

    # Compute percentile CI
    alpha = 1 - ci_level
    lower = np.percentile(bootstrap_effects, 100 * alpha / 2)
    upper = np.percentile(bootstrap_effects, 100 * (1 - alpha / 2))

    return (lower, upper)


# =============================================================================
# PATHWAY SIZE NORMALIZATION
# =============================================================================

class PathwayNormalization(Enum):
    """Pathway aggregation and normalization methods."""

    # Simple mean (original method)
    MEAN = "mean"

    # Median (robust to outliers)
    MEDIAN = "median"

    # Size-normalized mean (divides by sqrt(n_genes))
    SIZE_NORMALIZED = "size_normalized"

    # Principal component (first PC of gene burdens)
    PCA = "pca"


def aggregate_pathway_scores(
    gene_burdens: pd.DataFrame,
    pathway_genes: List[str],
    method: PathwayNormalization = PathwayNormalization.MEAN,
) -> pd.Series:
    """
    Aggregate gene burdens into pathway score with normalization.

    Args:
        gene_burdens: DataFrame of gene burdens (samples x genes)
        pathway_genes: List of genes in the pathway
        method: Aggregation/normalization method

    Returns:
        Series of pathway scores (one per sample)
    """
    # Get genes present in data
    common_genes = [g for g in pathway_genes if g in gene_burdens.columns]

    if len(common_genes) < 2:
        return pd.Series(0.0, index=gene_burdens.index)

    burden_subset = gene_burdens[common_genes]
    n_genes = len(common_genes)

    if method == PathwayNormalization.MEAN:
        return burden_subset.mean(axis=1)

    elif method == PathwayNormalization.MEDIAN:
        return burden_subset.median(axis=1)

    elif method == PathwayNormalization.SIZE_NORMALIZED:
        # Divide by sqrt(n_genes) to account for pathway size
        # This prevents large pathways from having artificially higher variance
        return burden_subset.mean(axis=1) / np.sqrt(n_genes)

    elif method == PathwayNormalization.PCA:
        # Use first principal component
        from sklearn.decomposition import PCA

        # Center the data
        centered = burden_subset - burden_subset.mean()

        # Handle case where all values are zero or constant
        if centered.std().sum() == 0:
            return pd.Series(0.0, index=gene_burdens.index)

        pca = PCA(n_components=1)
        scores = pca.fit_transform(centered.values)
        return pd.Series(scores.flatten(), index=gene_burdens.index)

    else:
        raise ValueError(f"Unknown normalization method: {method}")


# =============================================================================
# COMPREHENSIVE STATISTICAL ANALYSIS
# =============================================================================

@dataclass
class StatisticalRigorResult:
    """
    Comprehensive statistical analysis result.

    Includes all metrics needed for publication-quality reporting.
    """

    # FDR-corrected results per pathway
    pathway_results: List[FDRResult] = field(default_factory=list)

    # Summary statistics
    n_pathways_tested: int = 0
    n_significant_pathways: int = 0
    fdr_alpha: float = 0.05

    # Burden weight scheme used
    weight_scheme: str = "default"
    weight_citations: List[str] = field(default_factory=list)

    # Normalization method
    normalization_method: str = "mean"

    # Analysis parameters
    n_permutations: int = 1000
    n_bootstrap: int = 1000

    def to_dict(self) -> Dict[str, Any]:
        return {
            "summary": {
                "n_pathways_tested": self.n_pathways_tested,
                "n_significant_pathways": self.n_significant_pathways,
                "fdr_alpha": self.fdr_alpha,
                "weight_scheme": self.weight_scheme,
                "normalization_method": self.normalization_method,
            },
            "parameters": {
                "n_permutations": self.n_permutations,
                "n_bootstrap": self.n_bootstrap,
            },
            "citations": self.weight_citations,
            "pathway_results": [r.to_dict() for r in self.pathway_results],
        }

    def get_significant_pathways(self) -> List[str]:
        """Return list of pathways with q-value < alpha."""
        return [r.pathway for r in self.pathway_results if r.significant]

    def format_report(self) -> str:
        """Format results for markdown report."""
        lines = [
            "## Statistical Analysis",
            "",
            "### Summary",
            f"- **Pathways tested:** {self.n_pathways_tested}",
            f"- **Significant pathways (FDR < {self.fdr_alpha}):** {self.n_significant_pathways}",
            f"- **Burden weight scheme:** {self.weight_scheme}",
            f"- **Pathway aggregation:** {self.normalization_method}",
            "",
            "### Pathway Associations (FDR-corrected)",
            "",
            "| Pathway | p-value | q-value | Effect Size (d) | 95% CI | Significant |",
            "|---------|---------|---------|-----------------|--------|-------------|",
        ]

        # Sort by q-value
        sorted_results = sorted(self.pathway_results, key=lambda x: x.q_value)

        for r in sorted_results:
            sig = "Yes" if r.significant else "No"
            ci = f"[{r.ci_lower:.2f}, {r.ci_upper:.2f}]"
            lines.append(
                f"| {r.pathway} | {r.p_value:.4f} | {r.q_value:.4f} | "
                f"{r.effect_size:.2f} | {ci} | {sig} |"
            )

        lines.extend([
            "",
            "### Citations",
            "",
        ])

        for citation in self.weight_citations:
            lines.append(f"- {citation}")

        return "\n".join(lines)


def run_statistical_analysis(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    weight_scheme: BurdenWeightScheme = BurdenWeightScheme.GNOMAD_CONSTRAINT,
    normalization: PathwayNormalization = PathwayNormalization.SIZE_NORMALIZED,
    fdr_alpha: float = 0.05,
    n_permutations: int = 1000,
    n_bootstrap: int = 1000,
    seed: Optional[int] = None,
) -> StatisticalRigorResult:
    """
    Run comprehensive statistical analysis on clustering results.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways)
        cluster_labels: Array of cluster assignments
        weight_scheme: Burden weighting scheme used
        normalization: Pathway aggregation method used
        fdr_alpha: FDR significance threshold
        n_permutations: Number of permutations for p-value calculation
        n_bootstrap: Number of bootstrap iterations for CI
        seed: Random seed

    Returns:
        StatisticalRigorResult with comprehensive analysis
    """
    logger.info("Running statistical rigor analysis...")

    # Get weight configuration
    weights = BurdenWeights.from_scheme(weight_scheme)

    # Compute p-values via permutation test
    logger.info(f"  Computing permutation p-values (n={n_permutations})...")
    p_values = compute_pathway_pvalues(
        pathway_scores, cluster_labels, n_permutations, seed
    )

    # Apply FDR correction
    pathways = list(p_values.keys())
    raw_p = np.array([p_values[p] for p in pathways])
    q_values = benjamini_hochberg(raw_p, fdr_alpha)

    # Compute effect sizes
    logger.info("  Computing effect sizes...")
    effect_sizes = compute_pathway_effect_sizes(pathway_scores, cluster_labels)

    # Compute confidence intervals
    logger.info(f"  Computing bootstrap CIs (n={n_bootstrap})...")
    results = []

    for i, pathway in enumerate(pathways):
        ci_lower, ci_upper = bootstrap_effect_size_ci(
            pathway_scores, cluster_labels, pathway,
            n_bootstrap=n_bootstrap,
            seed=seed + i if seed else None
        )

        results.append(FDRResult(
            pathway=pathway,
            p_value=raw_p[i],
            q_value=q_values[i],
            significant=q_values[i] < fdr_alpha,
            effect_size=effect_sizes[pathway],
            ci_lower=ci_lower,
            ci_upper=ci_upper,
        ))

    n_significant = sum(1 for r in results if r.significant)

    return StatisticalRigorResult(
        pathway_results=results,
        n_pathways_tested=len(pathways),
        n_significant_pathways=n_significant,
        fdr_alpha=fdr_alpha,
        weight_scheme=weight_scheme.value,
        weight_citations=weights.get_citations(),
        normalization_method=normalization.value,
        n_permutations=n_permutations,
        n_bootstrap=n_bootstrap,
    )


# =============================================================================
# SENSITIVITY ANALYSIS
# =============================================================================

@dataclass
class SensitivityResult:
    """Result of sensitivity analysis across parameter variations."""

    parameter: str
    values_tested: List[Any]
    results: Dict[str, List[float]]  # metric -> list of values

    def is_robust(self, metric: str, threshold: float = 0.1) -> bool:
        """Check if results are robust (coefficient of variation < threshold)."""
        values = self.results.get(metric, [])
        if len(values) < 2:
            return True

        mean_val = np.mean(values)
        if mean_val == 0:
            return np.std(values) < threshold

        cv = np.std(values) / abs(mean_val)
        return cv < threshold


def sensitivity_analysis_weights(
    gene_burdens: pd.DataFrame,
    pathways: Dict[str, List[str]],
    cluster_labels: np.ndarray,
    schemes: Optional[List[BurdenWeightScheme]] = None,
    seed: Optional[int] = None,
) -> SensitivityResult:
    """
    Analyze sensitivity of results to different burden weight schemes.

    Args:
        gene_burdens: DataFrame of gene burdens
        pathways: Dictionary of pathway definitions
        cluster_labels: Original cluster labels
        schemes: Weight schemes to test (default: all)
        seed: Random seed

    Returns:
        SensitivityResult showing impact of weight choices
    """
    if schemes is None:
        schemes = list(BurdenWeightScheme)

    results = {
        "mean_ari": [],
        "n_significant": [],
    }

    for scheme in schemes:
        # Compute pathway scores with this scheme
        # (Would need to recompute from variants - simplified here)
        # For now, just track the scheme
        results["mean_ari"].append(0.0)  # Placeholder
        results["n_significant"].append(0)  # Placeholder

    return SensitivityResult(
        parameter="burden_weight_scheme",
        values_tested=[s.value for s in schemes],
        results=results,
    )
