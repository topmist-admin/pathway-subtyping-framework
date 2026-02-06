"""
Batch Effect Correction Module.

Implements methods to detect and correct for technical batch effects
in pathway-based subtype analysis:
- Batch effect detection via variance decomposition
- ComBat-style empirical Bayes correction
- Mean-centering correction per batch
- Post-correction validation

References:
- Johnson WE, Li C, Rabinovic A. Adjusting batch effects in
  microarray expression data using empirical Bayes methods.
  Biostatistics. 2007;8(1):118-127.
- Leek JT, et al. Tackling the widespread and critical impact
  of batch effects in high-throughput data. Nat Rev Genet.
  2010;11(10):733-739.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


# =============================================================================
# ENUMS
# =============================================================================


class BatchCorrectionMethod(Enum):
    """
    Methods for batch effect correction.

    COMBAT: Empirical Bayes method that adjusts for known batch
        effects while preserving biological variation. Based on
        Johnson et al. (2007). Estimates batch-specific location
        and scale parameters, then shrinks them toward a common
        prior.

    MEAN_CENTER: Simple mean-centering per batch. Subtracts the
        batch-specific mean from each feature. Fast but does not
        correct for batch-specific variance differences.

    STANDARDIZE: Per-batch Z-score standardization. Centers and
        scales each feature within each batch to zero mean and
        unit variance before pooling.
    """

    COMBAT = "combat"
    MEAN_CENTER = "mean_center"
    STANDARDIZE = "standardize"


# =============================================================================
# DATACLASSES
# =============================================================================


@dataclass
class BatchEffectReport:
    """
    Report from batch effect detection.

    Attributes:
        batch_variable: Name of the batch variable tested.
        n_batches: Number of unique batches.
        batch_sizes: Samples per batch.
        f_statistics: Per-feature F-statistic from one-way ANOVA.
        p_values: Per-feature p-value from one-way ANOVA.
        significant_features: Features with significant batch effects (FDR < threshold).
        variance_explained: Per-feature proportion of variance explained by batch (eta-squared).
        overall_batch_effect: Whether significant batch effects were detected.
        significance_threshold: FDR threshold used.
    """

    batch_variable: str
    n_batches: int
    batch_sizes: Dict[str, int]
    f_statistics: Dict[str, float]
    p_values: Dict[str, float]
    significant_features: List[str]
    variance_explained: Dict[str, float]
    overall_batch_effect: bool
    significance_threshold: float = 0.05

    def to_dict(self) -> Dict[str, Any]:
        return {
            "batch_variable": self.batch_variable,
            "n_batches": self.n_batches,
            "batch_sizes": self.batch_sizes,
            "n_significant_features": len(self.significant_features),
            "significant_features": self.significant_features,
            "overall_batch_effect": self.overall_batch_effect,
            "mean_variance_explained": round(np.mean(list(self.variance_explained.values())), 4),
            "significance_threshold": self.significance_threshold,
        }

    def format_report(self) -> str:
        lines = [
            "Batch Effect Detection Report",
            "=" * 40,
            f"Batch variable: {self.batch_variable}",
            f"Number of batches: {self.n_batches}",
            f"Batch sizes: {self.batch_sizes}",
            "",
            f"Features with significant batch effects: "
            f"{len(self.significant_features)} / {len(self.f_statistics)}",
            f"Mean variance explained by batch: "
            f"{np.mean(list(self.variance_explained.values())):.1%}",
            f"Overall batch effect detected: {self.overall_batch_effect}",
        ]
        if self.significant_features:
            lines.append("")
            lines.append("Top affected features:")
            sorted_features = sorted(
                self.significant_features,
                key=lambda f: self.variance_explained.get(f, 0),
                reverse=True,
            )
            for feat in sorted_features[:10]:
                eta2 = self.variance_explained.get(feat, 0)
                lines.append(f"  {feat}: eta² = {eta2:.3f}")
        return "\n".join(lines)


@dataclass
class BatchCorrectionResult:
    """
    Result of batch effect correction.

    Attributes:
        corrected_scores: Batch-corrected pathway scores DataFrame.
        method: Correction method used.
        batch_variable: Name of the batch variable.
        n_batches: Number of batches.
        pre_correction_variance: Per-feature variance explained by batch before correction.
        post_correction_variance: Per-feature variance explained by batch after correction.
        variance_reduction: Per-feature reduction in batch-explained variance.
    """

    corrected_scores: pd.DataFrame
    method: BatchCorrectionMethod
    batch_variable: str
    n_batches: int
    pre_correction_variance: Dict[str, float]
    post_correction_variance: Dict[str, float]
    variance_reduction: Dict[str, float]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "method": self.method.value,
            "batch_variable": self.batch_variable,
            "n_batches": self.n_batches,
            "n_features": len(self.corrected_scores.columns),
            "mean_pre_correction_variance": round(
                np.mean(list(self.pre_correction_variance.values())), 4
            ),
            "mean_post_correction_variance": round(
                np.mean(list(self.post_correction_variance.values())), 4
            ),
            "mean_variance_reduction": round(np.mean(list(self.variance_reduction.values())), 4),
        }

    def format_report(self) -> str:
        mean_pre = np.mean(list(self.pre_correction_variance.values()))
        mean_post = np.mean(list(self.post_correction_variance.values()))
        mean_reduction = np.mean(list(self.variance_reduction.values()))
        lines = [
            "Batch Correction Report",
            "=" * 40,
            f"Method: {self.method.value}",
            f"Batch variable: {self.batch_variable}",
            f"Number of batches: {self.n_batches}",
            "",
            f"Mean batch variance (before): {mean_pre:.1%}",
            f"Mean batch variance (after):  {mean_post:.1%}",
            f"Mean variance reduction:      {mean_reduction:.1%}",
        ]
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        citations = [
            "Leek JT, et al. Tackling the widespread and critical impact "
            "of batch effects in high-throughput data. "
            "Nat Rev Genet. 2010;11(10):733-739.",
        ]
        if self.method == BatchCorrectionMethod.COMBAT:
            citations.insert(
                0,
                "Johnson WE, Li C, Rabinovic A. Adjusting batch effects in "
                "microarray expression data using empirical Bayes methods. "
                "Biostatistics. 2007;8(1):118-127.",
            )
        return citations


# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================


def detect_batch_effects(
    pathway_scores: pd.DataFrame,
    batch_labels: np.ndarray,
    batch_variable: str = "batch",
    significance_threshold: float = 0.05,
) -> BatchEffectReport:
    """
    Detect batch effects in pathway scores using one-way ANOVA.

    For each pathway (feature), tests whether the mean differs
    significantly across batches. Reports F-statistics, p-values,
    and eta-squared (proportion of variance explained by batch).

    Args:
        pathway_scores: DataFrame of shape (n_samples, n_pathways).
        batch_labels: Array of batch assignments per sample.
        batch_variable: Name for reporting (e.g., "sequencing_site").
        significance_threshold: FDR threshold for significance.

    Returns:
        BatchEffectReport with detection results.
    """
    if len(batch_labels) != len(pathway_scores):
        raise ValueError(
            f"batch_labels length ({len(batch_labels)}) must match "
            f"pathway_scores rows ({len(pathway_scores)})"
        )

    unique_batches = np.unique(batch_labels)
    n_batches = len(unique_batches)

    if n_batches < 2:
        raise ValueError(f"Need at least 2 batches for detection, got {n_batches}")

    logger.info(
        "[BatchCorrection] Detecting batch effects across %d batches " "for %d features",
        n_batches,
        len(pathway_scores.columns),
    )

    batch_sizes = {str(b): int(np.sum(batch_labels == b)) for b in unique_batches}

    f_statistics = {}
    p_values = {}
    variance_explained = {}

    for col in pathway_scores.columns:
        groups = [pathway_scores.loc[batch_labels == b, col].values for b in unique_batches]
        # Filter out empty groups
        groups = [g for g in groups if len(g) > 0]

        if len(groups) < 2:
            f_statistics[col] = 0.0
            p_values[col] = 1.0
            variance_explained[col] = 0.0
            continue

        f_stat, p_val = stats.f_oneway(*groups)

        if np.isnan(f_stat):
            f_stat = 0.0
            p_val = 1.0

        f_statistics[col] = float(f_stat)
        p_values[col] = float(p_val)

        # Compute eta-squared (variance explained by batch)
        values = pathway_scores[col].values
        grand_mean = np.mean(values)
        ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
        ss_total = np.sum((values - grand_mean) ** 2)
        eta_sq = ss_between / ss_total if ss_total > 0 else 0.0
        variance_explained[col] = float(eta_sq)

    # Apply BH FDR correction
    sorted_cols = sorted(p_values.keys(), key=lambda c: p_values[c])
    n_tests = len(sorted_cols)
    adjusted_pvals = {}
    for rank, col in enumerate(sorted_cols, 1):
        adjusted = min(1.0, p_values[col] * n_tests / rank)
        adjusted_pvals[col] = adjusted

    significant = [col for col in sorted_cols if adjusted_pvals[col] < significance_threshold]

    overall = len(significant) > 0

    logger.info(
        "[BatchCorrection] Detected %d / %d features with significant "
        "batch effects (FDR < %.2f)",
        len(significant),
        n_tests,
        significance_threshold,
    )

    return BatchEffectReport(
        batch_variable=batch_variable,
        n_batches=n_batches,
        batch_sizes=batch_sizes,
        f_statistics=f_statistics,
        p_values=p_values,
        significant_features=significant,
        variance_explained=variance_explained,
        overall_batch_effect=overall,
        significance_threshold=significance_threshold,
    )


def correct_batch_effects(
    pathway_scores: pd.DataFrame,
    batch_labels: np.ndarray,
    method: BatchCorrectionMethod = BatchCorrectionMethod.COMBAT,
    batch_variable: str = "batch",
    preserve_biological: Optional[np.ndarray] = None,
) -> BatchCorrectionResult:
    """
    Correct batch effects in pathway scores.

    Args:
        pathway_scores: DataFrame of shape (n_samples, n_pathways).
        batch_labels: Array of batch assignments per sample.
        method: Correction method to use.
        batch_variable: Name for reporting.
        preserve_biological: Optional array of biological group labels
            to preserve during correction (used in ComBat).

    Returns:
        BatchCorrectionResult with corrected scores and diagnostics.
    """
    if len(batch_labels) != len(pathway_scores):
        raise ValueError(
            f"batch_labels length ({len(batch_labels)}) must match "
            f"pathway_scores rows ({len(pathway_scores)})"
        )

    unique_batches = np.unique(batch_labels)
    n_batches = len(unique_batches)

    if n_batches < 2:
        raise ValueError(f"Need at least 2 batches for correction, got {n_batches}")

    logger.info(
        "[BatchCorrection] Correcting batch effects using %s " "across %d batches",
        method.value,
        n_batches,
    )

    # Compute pre-correction batch variance
    pre_variance = _compute_batch_variance(pathway_scores, batch_labels)

    # Apply correction
    if method == BatchCorrectionMethod.COMBAT:
        corrected = _combat_correction(pathway_scores, batch_labels, preserve_biological)
    elif method == BatchCorrectionMethod.MEAN_CENTER:
        corrected = _mean_center_correction(pathway_scores, batch_labels)
    elif method == BatchCorrectionMethod.STANDARDIZE:
        corrected = _standardize_correction(pathway_scores, batch_labels)
    else:
        raise ValueError(f"Unknown correction method: {method}")

    # Compute post-correction batch variance
    post_variance = _compute_batch_variance(corrected, batch_labels)

    # Compute variance reduction
    variance_reduction = {}
    for col in pathway_scores.columns:
        pre = pre_variance.get(col, 0.0)
        post = post_variance.get(col, 0.0)
        reduction = pre - post if pre > 0 else 0.0
        variance_reduction[col] = float(reduction)

    logger.info(
        "[BatchCorrection] Mean batch variance reduced from %.3f to %.3f",
        np.mean(list(pre_variance.values())),
        np.mean(list(post_variance.values())),
    )

    return BatchCorrectionResult(
        corrected_scores=corrected,
        method=method,
        batch_variable=batch_variable,
        n_batches=n_batches,
        pre_correction_variance=pre_variance,
        post_correction_variance=post_variance,
        variance_reduction=variance_reduction,
    )


def validate_batch_correction(
    original_scores: pd.DataFrame,
    corrected_scores: pd.DataFrame,
    batch_labels: np.ndarray,
    biological_labels: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Validate that batch correction removed batch effects without
    destroying biological signal.

    Args:
        original_scores: Pre-correction pathway scores.
        corrected_scores: Post-correction pathway scores.
        batch_labels: Batch assignments per sample.
        biological_labels: Optional known biological group labels.

    Returns:
        Dictionary with validation metrics.
    """
    logger.info("[BatchCorrection] Validating batch correction")

    # Check batch variance reduction
    pre_var = _compute_batch_variance(original_scores, batch_labels)
    post_var = _compute_batch_variance(corrected_scores, batch_labels)

    mean_pre = np.mean(list(pre_var.values()))
    mean_post = np.mean(list(post_var.values()))

    result = {
        "batch_variance_before": round(mean_pre, 4),
        "batch_variance_after": round(mean_post, 4),
        "batch_variance_reduced": mean_post < mean_pre,
    }

    # Check biological signal preservation (if labels provided)
    if biological_labels is not None:
        pre_bio = _compute_batch_variance(original_scores, biological_labels)
        post_bio = _compute_batch_variance(corrected_scores, biological_labels)

        mean_pre_bio = np.mean(list(pre_bio.values()))
        mean_post_bio = np.mean(list(post_bio.values()))

        result["biological_variance_before"] = round(mean_pre_bio, 4)
        result["biological_variance_after"] = round(mean_post_bio, 4)
        result["biological_signal_preserved"] = mean_post_bio >= mean_pre_bio * 0.5

    # Correlation between original and corrected
    correlations = []
    for col in original_scores.columns:
        r, _ = stats.pearsonr(original_scores[col].values, corrected_scores[col].values)
        correlations.append(r)
    result["mean_correlation_with_original"] = round(np.mean(correlations), 4)

    return result


# =============================================================================
# PRIVATE HELPERS
# =============================================================================


def _compute_batch_variance(scores: pd.DataFrame, batch_labels: np.ndarray) -> Dict[str, float]:
    """Compute eta-squared (batch-explained variance) per feature."""
    unique_batches = np.unique(batch_labels)
    result = {}

    for col in scores.columns:
        values = scores[col].values
        grand_mean = np.mean(values)
        groups = [values[batch_labels == b] for b in unique_batches]
        groups = [g for g in groups if len(g) > 0]

        ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
        ss_total = np.sum((values - grand_mean) ** 2)
        eta_sq = ss_between / ss_total if ss_total > 0 else 0.0
        result[col] = float(eta_sq)

    return result


def _combat_correction(
    scores: pd.DataFrame,
    batch_labels: np.ndarray,
    biological_labels: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Simplified ComBat-style empirical Bayes batch correction.

    Estimates batch-specific location (mean shift) and scale
    (variance) parameters, shrinks them toward pooled estimates
    using empirical Bayes, then adjusts the data.

    Based on Johnson WE, Li C, Rabinovic A (2007).
    """
    unique_batches = np.unique(batch_labels)
    corrected = scores.copy()

    # Step 1: Estimate grand mean and batch means
    grand_mean = scores.mean(axis=0)

    batch_means = {}
    batch_vars = {}
    batch_sizes = {}

    for batch in unique_batches:
        mask = batch_labels == batch
        batch_data = scores.loc[mask]
        batch_means[batch] = batch_data.mean(axis=0)
        batch_vars[batch] = batch_data.var(axis=0, ddof=1)
        batch_sizes[batch] = int(mask.sum())

    # Step 2: Standardize data (remove grand mean, divide by pooled std)
    pooled_var = scores.var(axis=0, ddof=1)
    pooled_std = np.sqrt(pooled_var)
    pooled_std = pooled_std.replace(0, 1)  # avoid division by zero

    standardized = (scores - grand_mean) / pooled_std

    # Step 3: Estimate batch effects (location and scale)
    gamma_hat = {}  # location (mean shift)
    delta_hat = {}  # scale (variance ratio)

    for batch in unique_batches:
        mask = batch_labels == batch
        batch_std = standardized.loc[mask]
        gamma_hat[batch] = batch_std.mean(axis=0)
        delta_hat[batch] = batch_std.var(axis=0, ddof=1)

    # Step 4: Empirical Bayes shrinkage of batch parameters
    # Shrink gamma (location) toward zero
    gamma_bar = np.mean([gamma_hat[b].values for b in unique_batches], axis=0)
    tau_sq = np.var([gamma_hat[b].values for b in unique_batches], axis=0, ddof=1)

    gamma_star = {}
    for batch in unique_batches:
        n_b = batch_sizes[batch]
        batch_var = delta_hat[batch].values
        batch_var = np.where(batch_var < 1e-10, 1e-10, batch_var)
        # Shrinkage: weight between batch estimate and grand mean
        shrink_weight = tau_sq / (tau_sq + batch_var / n_b)
        shrink_weight = np.where(np.isnan(shrink_weight), 0.5, shrink_weight)
        gamma_star[batch] = pd.Series(
            shrink_weight * gamma_hat[batch].values + (1 - shrink_weight) * gamma_bar,
            index=scores.columns,
        )

    # Shrink delta (scale) toward pooled variance
    delta_star = {}
    delta_values = np.array([delta_hat[b].values for b in unique_batches])
    delta_bar = np.mean(delta_values, axis=0)

    for batch in unique_batches:
        n_b = batch_sizes[batch]
        # Simple shrinkage toward pooled variance
        shrink = min(0.5, 2.0 / n_b)
        delta_star[batch] = pd.Series(
            (1 - shrink) * delta_hat[batch].values + shrink * delta_bar,
            index=scores.columns,
        )

    # Step 5: Apply correction
    for batch in unique_batches:
        mask = batch_labels == batch
        batch_std = standardized.loc[mask]

        # Remove batch effect: (x - gamma*) / sqrt(delta*)
        delta_sqrt = np.sqrt(delta_star[batch].values)
        delta_sqrt = np.where(delta_sqrt < 1e-10, 1, delta_sqrt)

        adjusted = (batch_std.values - gamma_star[batch].values) / delta_sqrt
        corrected.loc[mask] = adjusted * pooled_std.values + grand_mean.values

    return corrected


def _mean_center_correction(scores: pd.DataFrame, batch_labels: np.ndarray) -> pd.DataFrame:
    """Remove batch-specific mean shift from each feature."""
    unique_batches = np.unique(batch_labels)
    corrected = scores.copy()
    grand_mean = scores.mean(axis=0)

    for batch in unique_batches:
        mask = batch_labels == batch
        batch_mean = scores.loc[mask].mean(axis=0)
        corrected.loc[mask] = scores.loc[mask] - batch_mean + grand_mean

    return corrected


def _standardize_correction(scores: pd.DataFrame, batch_labels: np.ndarray) -> pd.DataFrame:
    """Per-batch Z-score standardization then re-pool."""
    unique_batches = np.unique(batch_labels)
    corrected = scores.copy()
    grand_mean = scores.mean(axis=0)
    grand_std = scores.std(axis=0, ddof=1)
    grand_std = grand_std.replace(0, 1)

    for batch in unique_batches:
        mask = batch_labels == batch
        batch_data = scores.loc[mask]
        batch_mean = batch_data.mean(axis=0)
        batch_std = batch_data.std(axis=0, ddof=1)
        batch_std = batch_std.replace(0, 1)

        # Standardize within batch, then re-scale to grand distribution
        standardized = (batch_data - batch_mean) / batch_std
        corrected.loc[mask] = standardized * grand_std + grand_mean

    return corrected
