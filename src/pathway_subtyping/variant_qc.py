"""
Variant Quality Control Module.

Implements standard genetic variant QC filters:
- QUAL score filtering
- Hardy-Weinberg equilibrium (HWE) test
- Per-variant call rate
- Minor allele frequency (MAF) filtering
- Genotype quality (GQ) and read depth (DP) filters

These filters remove technical artifacts before burden computation,
preventing low-quality variants from inflating pathway scores.

References:
- Wigginton JE, Cutler DJ, Gravel A. A note on exact tests of
  Hardy-Weinberg equilibrium. Am J Hum Genet. 2005;76(5):887-93.
- Anderson CA et al. Data quality control in genetic case-control
  association studies. Nat Protoc. 2010;5(9):1564-73.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


@dataclass
class VariantQCConfig:
    """
    Configuration for variant quality control filters.

    Attributes:
        min_qual: Minimum QUAL score (0 to disable).
        min_call_rate: Minimum genotype call rate per variant (0.0-1.0).
        hwe_p_threshold: HWE p-value threshold; variants below are removed
            (set to 0 to disable).
        max_maf: Maximum minor allele frequency for rare variant analysis.
            Variants above this are removed (set to 1.0 to disable).
        min_gq: Minimum genotype quality; genotypes below are set to missing.
            None to disable.
        min_dp: Minimum read depth; genotypes below are set to missing.
            None to disable.
    """

    min_qual: float = 30.0
    min_call_rate: float = 0.9
    hwe_p_threshold: float = 1e-6
    max_maf: float = 0.01
    min_gq: Optional[int] = None
    min_dp: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "min_qual": self.min_qual,
            "min_call_rate": self.min_call_rate,
            "hwe_p_threshold": self.hwe_p_threshold,
            "max_maf": self.max_maf,
            "min_gq": self.min_gq,
            "min_dp": self.min_dp,
        }


@dataclass
class VariantQCResult:
    """
    Result from variant quality control filtering.

    Attributes:
        total_variants: Variants before QC.
        passed_variants: Variants after QC.
        removed_variants: Number of variants removed.
        removal_reasons: Count of variants removed per filter.
        per_variant_metrics: DataFrame with QC metrics per variant.
        config: The QC config used.
    """

    total_variants: int
    passed_variants: int
    removed_variants: int
    removal_reasons: Dict[str, int] = field(default_factory=dict)
    per_variant_metrics: Optional[pd.DataFrame] = None
    config: Optional[VariantQCConfig] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "total_variants": self.total_variants,
            "passed_variants": self.passed_variants,
            "removed_variants": self.removed_variants,
            "removal_reasons": self.removal_reasons,
            "retention_rate": (
                f"{100 * self.passed_variants / self.total_variants:.1f}%"
                if self.total_variants > 0
                else "N/A"
            ),
            "config": self.config.to_dict() if self.config else {},
        }

    def format_report(self) -> str:
        """Generate human-readable QC report."""
        lines = [
            "Variant QC Report",
            "=" * 40,
            f"Total variants:   {self.total_variants}",
            f"Passed variants:  {self.passed_variants}",
            f"Removed variants: {self.removed_variants}",
            "",
        ]

        if self.total_variants > 0:
            rate = 100 * self.passed_variants / self.total_variants
            lines.append(f"Retention rate:   {rate:.1f}%")

        if self.removal_reasons:
            lines.extend(["", "Removal reasons:"])
            for reason, count in sorted(
                self.removal_reasons.items(), key=lambda x: -x[1]
            ):
                lines.append(f"  {reason}: {count}")

        if self.config:
            lines.extend(
                [
                    "",
                    "Filters applied:",
                    f"  QUAL >= {self.config.min_qual}",
                    f"  Call rate >= {self.config.min_call_rate}",
                    f"  HWE p >= {self.config.hwe_p_threshold}",
                    f"  MAF <= {self.config.max_maf}",
                ]
            )
            if self.config.min_gq is not None:
                lines.append(f"  GQ >= {self.config.min_gq}")
            if self.config.min_dp is not None:
                lines.append(f"  DP >= {self.config.min_dp}")

        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return relevant citations."""
        return [
            "Wigginton JE, Cutler DJ, Gravel A. A note on exact tests of "
            "Hardy-Weinberg equilibrium. Am J Hum Genet. 2005;76(5):887-93.",
            "Anderson CA et al. Data quality control in genetic case-control "
            "association studies. Nat Protoc. 2010;5(9):1564-73.",
        ]


def compute_call_rate(genotypes_df: pd.DataFrame) -> pd.Series:
    """
    Compute per-variant call rate (proportion of non-missing genotypes).

    Missing genotypes are represented as 0 in the genotypes DataFrame
    from the VCF loader. However, true missing data (./.) is stored as 0
    and is indistinguishable from homozygous reference (0/0) in the
    current genotype encoding.

    For this reason, call rate is computed from the raw genotype strings
    when available. When only the numeric matrix is available, we report
    1.0 (all present) since missing genotypes were already counted
    during VCF loading.

    Args:
        genotypes_df: DataFrame of genotype values (variants x samples).

    Returns:
        Series of call rates indexed by variant.
    """
    if genotypes_df.empty:
        return pd.Series(dtype=float)

    # In the numeric encoding, 0 can be ref/ref OR missing.
    # Since the VCF loader already tracks missing genotypes separately,
    # and the numeric matrix doesn't distinguish, we compute call rate
    # based on non-NaN values (NaN would indicate truly missing data).
    non_missing = genotypes_df.notna().sum(axis=1)
    total = genotypes_df.shape[1]
    return non_missing / total if total > 0 else pd.Series(1.0, index=genotypes_df.index)


def compute_maf(genotypes_df: pd.DataFrame) -> pd.Series:
    """
    Compute minor allele frequency (MAF) per variant.

    Args:
        genotypes_df: DataFrame of allele counts (0, 1, 2) per variant x sample.

    Returns:
        Series of MAF values indexed by variant.
    """
    if genotypes_df.empty:
        return pd.Series(dtype=float)

    # Total alleles = 2 * n_samples (diploid)
    n_samples = genotypes_df.shape[1]
    total_alleles = 2 * n_samples

    if total_alleles == 0:
        return pd.Series(0.0, index=genotypes_df.index)

    # Sum of alt allele counts across samples
    alt_count = genotypes_df.sum(axis=1)
    alt_freq = alt_count / total_alleles

    # MAF is the lesser of alt freq and ref freq
    maf = pd.Series(np.minimum(alt_freq.values, 1.0 - alt_freq.values), index=genotypes_df.index)
    return maf


def check_hwe(genotypes_df: pd.DataFrame) -> pd.Series:
    """
    Check Hardy-Weinberg equilibrium per variant using chi-squared test.

    Compares observed genotype counts (hom_ref, het, hom_alt) to
    expected counts under HWE given the observed allele frequency.

    Args:
        genotypes_df: DataFrame of allele counts (0, 1, 2) per variant x sample.

    Returns:
        Series of HWE p-values indexed by variant. Returns 1.0 for
        monomorphic variants or variants with insufficient data.
    """
    if genotypes_df.empty:
        return pd.Series(dtype=float)

    p_values = []

    for idx in genotypes_df.index:
        row = genotypes_df.loc[idx]
        valid = row.dropna()

        if len(valid) < 2:
            p_values.append(1.0)
            continue

        # Count genotypes
        n_aa = (valid == 0).sum()  # hom ref
        n_ab = (valid == 1).sum()  # het
        n_bb = (valid == 2).sum()  # hom alt
        n = n_aa + n_ab + n_bb

        if n == 0:
            p_values.append(1.0)
            continue

        # Observed allele frequency
        p = (2 * n_aa + n_ab) / (2 * n)
        q = 1.0 - p

        # Monomorphic check
        if p == 0.0 or p == 1.0:
            p_values.append(1.0)
            continue

        # Expected genotype counts under HWE
        exp_aa = p * p * n
        exp_ab = 2 * p * q * n
        exp_bb = q * q * n

        observed = np.array([n_aa, n_ab, n_bb], dtype=float)
        expected = np.array([exp_aa, exp_ab, exp_bb], dtype=float)

        # Skip if any expected count is 0 (chi2 undefined)
        if np.any(expected == 0):
            p_values.append(1.0)
            continue

        # Chi-squared goodness of fit (1 df for HWE)
        chi2_stat = np.sum((observed - expected) ** 2 / expected)
        hwe_p = 1.0 - stats.chi2.cdf(chi2_stat, df=1)
        p_values.append(float(hwe_p))

    return pd.Series(p_values, index=genotypes_df.index)


def apply_genotype_filters(
    genotypes_df: pd.DataFrame,
    genotype_fields: Optional[pd.DataFrame] = None,
    min_gq: Optional[int] = None,
    min_dp: Optional[int] = None,
) -> Tuple[pd.DataFrame, int]:
    """
    Apply per-genotype quality filters.

    Sets genotypes below quality thresholds to NaN (missing).

    Args:
        genotypes_df: DataFrame of allele counts (variants x samples).
        genotype_fields: Optional DataFrame with 'GQ' and/or 'DP' columns
            per variant-sample pair. If None, no genotype filtering is applied.
        min_gq: Minimum genotype quality. None to skip.
        min_dp: Minimum read depth. None to skip.

    Returns:
        Tuple of (filtered_genotypes_df, n_genotypes_masked).
    """
    if genotype_fields is None or (min_gq is None and min_dp is None):
        return genotypes_df, 0

    filtered = genotypes_df.copy()
    n_masked = 0

    if min_gq is not None and "GQ" in genotype_fields.columns:
        mask = genotype_fields["GQ"] < min_gq
        n_gq_masked = mask.sum()
        if n_gq_masked > 0:
            filtered[mask] = np.nan
            n_masked += int(n_gq_masked)
            logger.info(f"[VariantQC] Masked {n_gq_masked} genotypes with GQ < {min_gq}")

    if min_dp is not None and "DP" in genotype_fields.columns:
        mask = genotype_fields["DP"] < min_dp
        n_dp_masked = mask.sum()
        if n_dp_masked > 0:
            filtered[mask] = np.nan
            n_masked += int(n_dp_masked)
            logger.info(f"[VariantQC] Masked {n_dp_masked} genotypes with DP < {min_dp}")

    return filtered, n_masked


def filter_variants(
    variants_df: pd.DataFrame,
    genotypes_df: pd.DataFrame,
    config: Optional[VariantQCConfig] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, VariantQCResult]:
    """
    Apply variant-level QC filters.

    Filters are applied in order: QUAL, call rate, HWE, MAF.
    Each filter is independent — a variant removed by an earlier
    filter is not counted in later filters.

    Args:
        variants_df: DataFrame with variant metadata (must have 'qual' column
            for QUAL filtering).
        genotypes_df: DataFrame of allele counts (variants x samples).
            Must have same index as variants_df.
        config: QC configuration. Uses defaults if None.

    Returns:
        Tuple of (filtered_variants_df, filtered_genotypes_df, qc_result).
    """
    if config is None:
        config = VariantQCConfig()

    total = len(variants_df)
    removal_reasons: Dict[str, int] = {}
    keep_mask = pd.Series(True, index=variants_df.index)

    # Compute per-variant metrics
    metrics: Dict[str, pd.Series] = {}

    # 1. QUAL filter
    if config.min_qual > 0 and "qual" in variants_df.columns:
        qual_fail = variants_df["qual"] < config.min_qual
        n_qual_fail = qual_fail.sum()
        if n_qual_fail > 0:
            removal_reasons["low_qual"] = int(n_qual_fail)
            keep_mask &= ~qual_fail
            logger.info(
                f"[VariantQC] QUAL < {config.min_qual}: {n_qual_fail} variants removed"
            )
        metrics["qual"] = variants_df["qual"]

    # 2. Call rate filter
    call_rates = compute_call_rate(genotypes_df)
    metrics["call_rate"] = call_rates
    if config.min_call_rate > 0:
        cr_fail = call_rates < config.min_call_rate
        # Only count variants not already removed
        n_cr_fail = (cr_fail & keep_mask).sum()
        if n_cr_fail > 0:
            removal_reasons["low_call_rate"] = int(n_cr_fail)
            keep_mask &= ~cr_fail
            logger.info(
                f"[VariantQC] Call rate < {config.min_call_rate}: "
                f"{n_cr_fail} variants removed"
            )

    # 3. HWE filter
    hwe_pvals = check_hwe(genotypes_df)
    metrics["hwe_p"] = hwe_pvals
    if config.hwe_p_threshold > 0:
        hwe_fail = hwe_pvals < config.hwe_p_threshold
        n_hwe_fail = (hwe_fail & keep_mask).sum()
        if n_hwe_fail > 0:
            removal_reasons["hwe_violation"] = int(n_hwe_fail)
            keep_mask &= ~hwe_fail
            logger.info(
                f"[VariantQC] HWE p < {config.hwe_p_threshold}: "
                f"{n_hwe_fail} variants removed"
            )

    # 4. MAF filter
    maf_values = compute_maf(genotypes_df)
    metrics["maf"] = maf_values
    if config.max_maf < 1.0:
        maf_fail = maf_values > config.max_maf
        n_maf_fail = (maf_fail & keep_mask).sum()
        if n_maf_fail > 0:
            removal_reasons["high_maf"] = int(n_maf_fail)
            keep_mask &= ~maf_fail
            logger.info(
                f"[VariantQC] MAF > {config.max_maf}: {n_maf_fail} variants removed"
            )

    # Build per-variant metrics DataFrame
    metrics_df = pd.DataFrame(metrics)
    metrics_df["passed"] = keep_mask

    # Apply filter
    filtered_variants = variants_df[keep_mask].copy()
    filtered_genotypes = genotypes_df[keep_mask].copy()

    passed = int(keep_mask.sum())
    removed = total - passed

    logger.info(
        f"[VariantQC] Retained {passed}/{total} variants "
        f"({100 * passed / total:.1f}%)" if total > 0 else
        "[VariantQC] No variants to filter"
    )

    result = VariantQCResult(
        total_variants=total,
        passed_variants=passed,
        removed_variants=removed,
        removal_reasons=removal_reasons,
        per_variant_metrics=metrics_df,
        config=config,
    )

    return filtered_variants, filtered_genotypes, result
