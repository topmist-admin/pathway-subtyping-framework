"""
Bulk Deconvolution Integration

Estimates cell-type proportions from bulk RNA-seq data using a
single-cell reference expression profile. Uses non-negative least
squares (NNLS) to solve: bulk ≈ reference.T @ proportions.

The estimated proportions can be combined with pathway scores to
create subtypes that are both pathway-defined and cell-type-aware.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from tqdm import tqdm

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class DeconvolutionMethod(Enum):
    """Deconvolution algorithms for bulk expression."""

    NNLS = "nnls"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class DeconvolutionQualityReport:
    """Quality metrics for a deconvolution run."""

    n_samples: int = 0
    n_genes_bulk: int = 0
    n_genes_reference: int = 0
    n_genes_shared: int = 0
    n_cell_types: int = 0
    gene_coverage: float = 0.0
    proportion_stats: Dict[str, Dict[str, float]] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    is_usable: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_samples": self.n_samples,
            "n_genes_bulk": self.n_genes_bulk,
            "n_genes_reference": self.n_genes_reference,
            "n_genes_shared": self.n_genes_shared,
            "n_cell_types": self.n_cell_types,
            "gene_coverage": round(float(self.gene_coverage), 4),
            "proportion_stats": {
                ct: {k: round(float(v), 4) for k, v in stats.items()}
                for ct, stats in self.proportion_stats.items()
            },
            "warnings": list(self.warnings),
            "is_usable": self.is_usable,
        }


@dataclass
class DeconvolutionResult:
    """Result of bulk deconvolution."""

    cell_type_proportions: pd.DataFrame
    method: DeconvolutionMethod
    quality_report: DeconvolutionQualityReport
    reference_cell_types: List[str]
    n_samples: int = 0
    n_cell_types: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cell_type_proportions": {
                "index": list(self.cell_type_proportions.index),
                "columns": list(self.cell_type_proportions.columns),
                "shape": list(self.cell_type_proportions.shape),
            },
            "method": self.method.value,
            "quality_report": self.quality_report.to_dict(),
            "reference_cell_types": list(self.reference_cell_types),
            "n_samples": self.n_samples,
            "n_cell_types": self.n_cell_types,
        }

    def format_report(self) -> str:
        """Human-readable deconvolution report."""
        qr = self.quality_report
        lines = [
            "# Bulk Deconvolution Report",
            "",
            f"**Method:** {self.method.value.upper()}",
            f"**Samples:** {self.n_samples}",
            f"**Cell types:** {self.n_cell_types} ({', '.join(self.reference_cell_types)})",
            "",
            "## Gene Overlap",
            f"- Bulk genes: {qr.n_genes_bulk}",
            f"- Reference genes: {qr.n_genes_reference}",
            f"- Shared genes: {qr.n_genes_shared} ({qr.gene_coverage:.1%})",
            "",
            "## Estimated Proportions (mean across samples)",
        ]

        for ct in self.reference_cell_types:
            stats = qr.proportion_stats.get(ct, {})
            mean_val = stats.get("mean", 0.0)
            min_val = stats.get("min", 0.0)
            max_val = stats.get("max", 0.0)
            lines.append(f"- {ct}: {mean_val:.3f} (range: {min_val:.3f}–{max_val:.3f})")

        if qr.warnings:
            lines.extend(["", "## Warnings"])
            for w in qr.warnings:
                lines.append(f"- {w}")

        lines.extend(
            [
                "",
                "## Interpretation",
                "",
                "Cell-type proportions are estimated from bulk expression using a "
                "single-cell reference profile. Proportions sum to 1.0 per sample. "
                "Use `combine_features()` to merge with pathway scores for "
                "cell-type-aware subtype discovery.",
            ]
        )

        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        return [
            "Lawson CL, Hanson RJ. Solving Least Squares Problems. " "Prentice-Hall, 1974.",
            "Newman AM, et al. Robust enumeration of cell subsets from "
            "tissue expression profiles. Nat Methods. 2015;12(5):453-457.",
            "Wang X, et al. Bulk tissue cell type deconvolution with "
            "multi-subject single-cell expression reference. "
            "Nat Commun. 2019;10(1):380.",
        ]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _nnls_deconvolve(bulk_vector: np.ndarray, reference_matrix: np.ndarray) -> np.ndarray:
    """
    Estimate cell-type proportions for a single sample via NNLS.

    Solves: bulk ≈ reference_matrix.T @ proportions, s.t. proportions >= 0.

    Args:
        bulk_vector: Gene expression for one sample (n_genes,).
        reference_matrix: Reference profiles (n_cell_types x n_genes).

    Returns:
        Proportion vector (n_cell_types,), normalized to sum=1.
    """
    # NNLS expects: minimize ||A @ x - b||, s.t. x >= 0
    # A = reference_matrix.T (n_genes x n_cell_types)
    # b = bulk_vector (n_genes,)
    A = reference_matrix.T
    proportions, _ = nnls(A, bulk_vector)

    total = proportions.sum()
    if total > 0:
        proportions = proportions / total
    else:
        # Uniform if NNLS finds no signal
        n_types = len(proportions)
        proportions = np.ones(n_types) / n_types

    return proportions


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def build_reference_profile(
    reference_expr: pd.DataFrame,
    cell_type_labels: np.ndarray,
    min_cells_per_type: int = 5,
) -> pd.DataFrame:
    """
    Build cell-type mean expression profiles from single-cell reference.

    Aggregates expression per cell type (mean), excluding types with
    fewer than min_cells_per_type cells.

    Args:
        reference_expr: Single-cell expression matrix (cells x genes).
        cell_type_labels: Array of cell-type labels for each cell.
        min_cells_per_type: Minimum cells required per type.

    Returns:
        DataFrame of mean expression profiles (cell_types x genes).

    Raises:
        ValueError: If inputs are empty or mismatched.
    """
    if reference_expr is None or reference_expr.empty:
        raise ValueError("[Deconvolution] reference_expr is empty or None")

    labels = np.asarray(cell_type_labels)
    if len(labels) != len(reference_expr):
        raise ValueError(
            f"[Deconvolution] cell_type_labels length ({len(labels)}) "
            f"does not match reference_expr rows ({len(reference_expr)})"
        )

    unique_types = sorted(set(labels))
    profiles = {}
    excluded = []

    for ct in unique_types:
        mask = labels == ct
        n_cells = int(mask.sum())

        if n_cells < min_cells_per_type:
            excluded.append(str(ct))
            logger.info(
                f"[Deconvolution] Excluding cell type '{ct}': "
                f"{n_cells} cells < min {min_cells_per_type}"
            )
            continue

        profiles[str(ct)] = reference_expr.values[mask].mean(axis=0)

    if not profiles:
        raise ValueError(
            "[Deconvolution] No cell types passed min_cells_per_type filter "
            f"(threshold={min_cells_per_type})"
        )

    result = pd.DataFrame(profiles, index=reference_expr.columns).T

    if excluded:
        logger.info(f"[Deconvolution] Excluded {len(excluded)} cell types: {excluded}")

    logger.info(
        f"[Deconvolution] Built reference: {len(result)} cell types x " f"{result.shape[1]} genes"
    )

    return result


def deconvolve_bulk(
    bulk_expr: pd.DataFrame,
    reference_profile: pd.DataFrame,
    method: DeconvolutionMethod = DeconvolutionMethod.NNLS,
    min_genes: int = 50,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> DeconvolutionResult:
    """
    Estimate cell-type proportions from bulk RNA-seq expression.

    For each sample, solves bulk ≈ reference.T @ proportions using NNLS,
    then normalizes proportions to sum to 1.

    Args:
        bulk_expr: Bulk expression matrix (samples x genes).
        reference_profile: Cell-type mean expression (cell_types x genes),
            from build_reference_profile().
        method: Deconvolution algorithm (default: NNLS).
        min_genes: Minimum shared genes required.
        seed: Random seed (reserved for future stochastic methods).
        show_progress: Show tqdm progress bar.

    Returns:
        DeconvolutionResult with cell_type_proportions DataFrame.

    Raises:
        ValueError: If inputs are invalid or insufficient gene overlap.
    """
    if bulk_expr is None or bulk_expr.empty:
        raise ValueError("[Deconvolution] bulk_expr is empty or None")

    if reference_profile is None or reference_profile.empty:
        raise ValueError("[Deconvolution] reference_profile is empty or None")

    # Find shared genes
    bulk_genes = set(bulk_expr.columns)
    ref_genes = set(reference_profile.columns)
    shared_genes = sorted(bulk_genes & ref_genes)

    if len(shared_genes) < min_genes:
        raise ValueError(
            f"[Deconvolution] Only {len(shared_genes)} shared genes between bulk and "
            f"reference (minimum {min_genes} required). "
            f"Bulk has {len(bulk_genes)} genes, reference has {len(ref_genes)} genes."
        )

    logger.info(
        f"[Deconvolution] Using {len(shared_genes)} shared genes "
        f"(bulk: {len(bulk_genes)}, reference: {len(ref_genes)})"
    )

    # Align to shared genes
    bulk_aligned = bulk_expr[shared_genes].values
    ref_aligned = reference_profile[shared_genes].values
    cell_types = list(reference_profile.index)
    sample_ids = list(bulk_expr.index)

    # Per-sample deconvolution
    n_samples = len(sample_ids)
    n_cell_types = len(cell_types)
    proportions = np.zeros((n_samples, n_cell_types))

    iterator = range(n_samples)
    if show_progress:
        iterator = tqdm(iterator, desc="Deconvolving", unit="sample")

    for i in iterator:
        proportions[i] = _nnls_deconvolve(bulk_aligned[i], ref_aligned)

    # Build proportions DataFrame
    proportions_df = pd.DataFrame(
        proportions,
        index=sample_ids,
        columns=cell_types,
    )

    # Build quality report
    warnings = []
    gene_coverage = len(shared_genes) / max(len(ref_genes), 1)

    if gene_coverage < 0.3:
        warnings.append(
            f"Low gene overlap: {len(shared_genes)}/{len(ref_genes)} "
            f"reference genes found in bulk ({gene_coverage:.1%})"
        )

    # Per-cell-type proportion statistics
    proportion_stats: Dict[str, Dict[str, float]] = {}
    for ct in cell_types:
        col = proportions_df[ct].values
        proportion_stats[ct] = {
            "mean": float(col.mean()),
            "min": float(col.min()),
            "max": float(col.max()),
            "std": float(col.std()),
        }

    # Check for cell types with near-zero proportions everywhere
    for ct in cell_types:
        if proportion_stats[ct]["max"] < 0.01:
            warnings.append(
                f"Cell type '{ct}' has near-zero estimated proportions "
                f"(max={proportion_stats[ct]['max']:.4f})"
            )

    quality_report = DeconvolutionQualityReport(
        n_samples=n_samples,
        n_genes_bulk=len(bulk_genes),
        n_genes_reference=len(ref_genes),
        n_genes_shared=len(shared_genes),
        n_cell_types=n_cell_types,
        gene_coverage=gene_coverage,
        proportion_stats=proportion_stats,
        warnings=warnings,
        is_usable=len(shared_genes) >= min_genes and n_samples >= 2,
    )

    logger.info(
        f"[Deconvolution] Estimated proportions for {n_samples} samples "
        f"across {n_cell_types} cell types"
    )

    return DeconvolutionResult(
        cell_type_proportions=proportions_df,
        method=method,
        quality_report=quality_report,
        reference_cell_types=cell_types,
        n_samples=n_samples,
        n_cell_types=n_cell_types,
    )


def combine_features(
    pathway_scores: pd.DataFrame,
    cell_type_proportions: pd.DataFrame,
    proportion_weight: float = 0.5,
) -> pd.DataFrame:
    """
    Combine pathway scores with cell-type proportions into a unified feature matrix.

    Both inputs are Z-normalized independently, then scaled by their
    respective weights before concatenation.

    Args:
        pathway_scores: Pathway scores (samples x pathways).
        cell_type_proportions: Cell-type proportions (samples x cell_types).
        proportion_weight: Weight for proportions (0-1). Pathway scores
            get weight (1 - proportion_weight). Default: 0.5.

    Returns:
        Combined DataFrame (shared_samples x [pathways + cell_type_features]).

    Raises:
        ValueError: If no shared samples or invalid weight.
    """
    if not 0.0 <= proportion_weight <= 1.0:
        raise ValueError(
            f"[Deconvolution] proportion_weight must be in [0, 1], got {proportion_weight}"
        )

    # Find shared samples
    shared_samples = sorted(set(pathway_scores.index) & set(cell_type_proportions.index))

    if len(shared_samples) == 0:
        raise ValueError(
            "[Deconvolution] No shared samples between pathway_scores and " "cell_type_proportions"
        )

    if len(shared_samples) < len(pathway_scores):
        logger.warning(
            f"[Deconvolution] Using {len(shared_samples)}/{len(pathway_scores)} "
            f"samples (intersection with cell-type proportions)"
        )

    # Align to shared samples
    ps_aligned = pathway_scores.loc[shared_samples].copy()
    ct_aligned = cell_type_proportions.loc[shared_samples].copy()

    # Z-normalize pathway scores
    ps_mean = ps_aligned.mean(axis=0)
    ps_std = ps_aligned.std(axis=0)
    ps_std = ps_std.replace(0, 1)  # avoid division by zero
    ps_z = (ps_aligned - ps_mean) / ps_std

    # Z-normalize cell-type proportions
    ct_mean = ct_aligned.mean(axis=0)
    ct_std = ct_aligned.std(axis=0)
    ct_std = ct_std.replace(0, 1)
    ct_z = (ct_aligned - ct_mean) / ct_std

    # Apply weights
    pathway_weight = 1.0 - proportion_weight
    ps_weighted = ps_z * pathway_weight
    ct_weighted = ct_z * proportion_weight

    # Rename cell-type columns with prefix
    ct_weighted.columns = [f"CELLTYPE:{ct}_proportion" for ct in ct_weighted.columns]

    # Concatenate
    combined = pd.concat([ps_weighted, ct_weighted], axis=1)

    logger.info(
        f"[Deconvolution] Combined features: {combined.shape[0]} samples x "
        f"{combined.shape[1]} features "
        f"({ps_aligned.shape[1]} pathways + {ct_aligned.shape[1]} cell types, "
        f"weight={proportion_weight:.2f})"
    )

    return combined


def generate_synthetic_bulk(
    reference_profile: pd.DataFrame,
    n_samples: int = 100,
    n_subtypes: int = 3,
    noise_level: float = 0.1,
    seed: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """
    Generate synthetic bulk expression from known cell-type proportions.

    Each subtype has a distinct cell-type composition pattern. Bulk
    expression is computed as proportions @ reference_profile + noise.

    Args:
        reference_profile: Cell-type mean expression (cell_types x genes).
        n_samples: Number of bulk samples to generate.
        n_subtypes: Number of distinct subtypes (composition patterns).
        noise_level: Standard deviation of Gaussian noise.
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (bulk_expr, true_proportions, subtype_labels):
            - bulk_expr: DataFrame (samples x genes)
            - true_proportions: DataFrame (samples x cell_types)
            - subtype_labels: Array of subtype assignments
    """
    rng = np.random.RandomState(seed)
    cell_types = list(reference_profile.index)
    genes = list(reference_profile.columns)
    n_cell_types = len(cell_types)

    if n_subtypes > n_cell_types:
        raise ValueError(
            f"[Deconvolution] n_subtypes ({n_subtypes}) cannot exceed "
            f"n_cell_types ({n_cell_types})"
        )

    # Generate subtype-specific proportion templates
    # Each subtype has one dominant cell type
    samples_per_subtype = n_samples // n_subtypes
    subtype_labels = np.repeat(np.arange(n_subtypes), samples_per_subtype)
    # Handle remainder
    remainder = n_samples - len(subtype_labels)
    if remainder > 0:
        subtype_labels = np.concatenate([subtype_labels, np.arange(remainder)])

    proportions = np.zeros((n_samples, n_cell_types))

    for s in range(n_subtypes):
        mask = subtype_labels == s

        # Dominant cell type for this subtype
        template = np.ones(n_cell_types) * 0.05
        template[s % n_cell_types] = 0.6
        # Distribute remaining among others
        remaining = 1.0 - template.sum()
        non_dominant = [i for i in range(n_cell_types) if i != s % n_cell_types]
        for idx in non_dominant:
            template[idx] += remaining / len(non_dominant)

        # Add per-sample variation
        for i in np.where(mask)[0]:
            noise = rng.dirichlet(template * 20)
            proportions[i] = noise

    # Generate bulk expression: proportions @ reference + noise
    ref_matrix = reference_profile.values  # cell_types x genes
    bulk_matrix = proportions @ ref_matrix

    # Add Gaussian noise
    if noise_level > 0:
        noise = rng.normal(0, noise_level * bulk_matrix.std(), bulk_matrix.shape)
        bulk_matrix = np.maximum(bulk_matrix + noise, 0)  # ensure non-negative

    sample_ids = [f"SAMPLE_{i:04d}" for i in range(n_samples)]

    bulk_df = pd.DataFrame(bulk_matrix, index=sample_ids, columns=genes)
    proportions_df = pd.DataFrame(proportions, index=sample_ids, columns=cell_types)

    logger.info(
        f"[Deconvolution] Generated synthetic bulk: {n_samples} samples, "
        f"{n_subtypes} subtypes, {n_cell_types} cell types"
    )

    return bulk_df, proportions_df, subtype_labels
