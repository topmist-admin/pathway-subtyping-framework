"""
Multi-Omic Pathway Score Fusion

Fuses pathway scores from multiple modalities (VCF variant burden,
bulk RNA-seq expression, single-cell scRNA-seq) into a combined
feature matrix for unified subtype discovery.

Fusion strategies:
    - concatenate: Column-bind pathway scores with modality prefixes;
      preserves all information; handles partial sample overlap
    - weighted_average: Weighted mean of shared pathway scores;
      same dimensionality as single-modality output
    - intersection_only: Restrict to samples and pathways present
      in all modalities

This is the Phase 4 integration module. Advanced fusion methods
(MOFA+, kernel-based SNF) are planned for Phase 5.

Research use only. Not for clinical decision-making.
"""

import itertools
import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class ModalityType(Enum):
    """Types of omic modalities supported for fusion."""

    VCF = "vcf"
    EXPRESSION = "expression"
    SINGLE_CELL = "single_cell"


class FusionStrategy(Enum):
    """Strategies for fusing multi-omic pathway scores."""

    CONCATENATE = "concatenate"
    WEIGHTED_AVERAGE = "weighted_average"
    INTERSECTION_ONLY = "intersection_only"


class MissingStrategy(Enum):
    """How to handle samples missing from some modalities."""

    IMPUTE_ZERO = "impute_zero"
    IMPUTE_MEAN = "impute_mean"
    DROP = "drop"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class ModalityInput:
    """Wrapped input from a single omic modality.

    Attributes:
        modality_type: Which omic type this represents.
        pathway_scores: Z-normalized pathway scores (samples x pathways).
        gene_data: Gene-level data (gene_burdens or gene_expression).
                   Optional; used for downstream characterization.
        quality_report: Quality report from the modality's scoring step.
        label: Human-readable label for this modality (e.g. "WES", "RNA-seq").
    """

    modality_type: ModalityType
    pathway_scores: pd.DataFrame
    gene_data: Optional[pd.DataFrame] = None
    quality_report: Optional[Any] = None
    label: Optional[str] = None

    @property
    def sample_ids(self) -> List[str]:
        return list(self.pathway_scores.index)

    @property
    def pathway_names(self) -> List[str]:
        return list(self.pathway_scores.columns)

    @property
    def n_samples(self) -> int:
        return len(self.pathway_scores)

    @property
    def n_pathways(self) -> int:
        return len(self.pathway_scores.columns)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "modality_type": self.modality_type.value,
            "label": self.label or self.modality_type.value,
            "n_samples": self.n_samples,
            "n_pathways": self.n_pathways,
            "sample_ids": self.sample_ids[:10],
            "pathway_names": self.pathway_names[:10],
        }


@dataclass
class SampleOverlapStats:
    """Statistics about sample overlap across modalities.

    Attributes:
        n_modalities: Number of modalities.
        samples_per_modality: Dict of label -> list of sample IDs.
        n_samples_per_modality: Dict of label -> count.
        shared_samples: Samples present in ALL modalities.
        n_shared: Number of shared samples.
        union_samples: All unique samples across modalities.
        n_union: Total unique samples.
        pairwise_overlap: Dict of "label_a_vs_label_b" -> overlap count.
        overlap_fraction: Fraction of union that is shared across all.
    """

    n_modalities: int = 0
    samples_per_modality: Dict[str, List[str]] = field(default_factory=dict)
    n_samples_per_modality: Dict[str, int] = field(default_factory=dict)
    shared_samples: List[str] = field(default_factory=list)
    n_shared: int = 0
    union_samples: List[str] = field(default_factory=list)
    n_union: int = 0
    pairwise_overlap: Dict[str, int] = field(default_factory=dict)
    overlap_fraction: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_modalities": self.n_modalities,
            "n_samples_per_modality": dict(self.n_samples_per_modality),
            "n_shared_samples": self.n_shared,
            "n_union_samples": self.n_union,
            "overlap_fraction": round(self.overlap_fraction, 4),
            "pairwise_overlap": dict(self.pairwise_overlap),
        }


@dataclass
class MultiOmicQualityReport:
    """Quality report for multi-omic fusion.

    Attributes:
        per_modality_reports: Quality report dict per modality label.
        sample_overlap: Sample overlap statistics.
        pathway_overlap_count: Pathways shared across ALL modalities.
        pathway_union_count: Unique pathways across all modalities.
        pathway_overlap_names: Names of shared pathways.
        warnings: List of warning messages.
        is_usable: Whether the fused data is usable.
    """

    per_modality_reports: Dict[str, Any] = field(default_factory=dict)
    sample_overlap: Optional[SampleOverlapStats] = None
    pathway_overlap_count: int = 0
    pathway_union_count: int = 0
    pathway_overlap_names: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    is_usable: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "per_modality_reports": dict(self.per_modality_reports),
            "sample_overlap": self.sample_overlap.to_dict() if self.sample_overlap else {},
            "pathway_overlap_count": self.pathway_overlap_count,
            "pathway_union_count": self.pathway_union_count,
            "pathway_overlap_names": self.pathway_overlap_names[:20],
            "warnings": list(self.warnings),
            "is_usable": self.is_usable,
        }


@dataclass
class MultiOmicFusionResult:
    """Result of multi-omic pathway score fusion.

    Attributes:
        fused_pathway_scores: Combined pathway scores (samples x pathways).
            Z-normalized. Compatible with run_clustering(), ValidationGates,
            and characterize_subtypes().
        fused_gene_data: Combined gene-level data (union of gene columns).
        per_modality_scores: Original per-modality scores, keyed by label.
        strategy: Fusion strategy used.
        quality_report: Multi-omic quality report.
        sample_overlap: Sample overlap statistics.
        modality_weights: Weights used (relevant for WEIGHTED_AVERAGE).
        n_modalities: Number of modalities fused.
        modality_labels: Labels for each modality.
        cross_modality_correlations: Pearson correlations between modality
            pairs for shared pathways.
    """

    fused_pathway_scores: pd.DataFrame = field(default_factory=pd.DataFrame)
    fused_gene_data: Optional[pd.DataFrame] = None
    per_modality_scores: Dict[str, pd.DataFrame] = field(default_factory=dict)
    strategy: FusionStrategy = FusionStrategy.CONCATENATE
    quality_report: Optional[MultiOmicQualityReport] = None
    sample_overlap: Optional[SampleOverlapStats] = None
    modality_weights: Dict[str, float] = field(default_factory=dict)
    n_modalities: int = 0
    modality_labels: List[str] = field(default_factory=list)
    cross_modality_correlations: Optional[pd.DataFrame] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "strategy": self.strategy.value,
            "n_modalities": self.n_modalities,
            "modality_labels": list(self.modality_labels),
            "fused_shape": {
                "n_samples": len(self.fused_pathway_scores),
                "n_pathways": len(self.fused_pathway_scores.columns),
            },
            "modality_weights": {k: round(v, 4) for k, v in self.modality_weights.items()},
            "sample_overlap": self.sample_overlap.to_dict() if self.sample_overlap else {},
            "quality_report": self.quality_report.to_dict() if self.quality_report else {},
            "has_cross_modality_correlations": self.cross_modality_correlations is not None,
        }

    def format_report(self) -> str:
        """Generate a human-readable multi-omic fusion report."""
        lines = [
            "## Multi-Omic Fusion Report",
            "",
            f"**Strategy:** {self.strategy.value}",
            f"**Modalities:** {', '.join(self.modality_labels)}",
            f"**Fused samples:** {len(self.fused_pathway_scores)}",
            f"**Fused pathways:** {len(self.fused_pathway_scores.columns)}",
            "",
            "### Modality Summary",
            "",
            "| Modality | Samples | Pathways |",
            "|----------|---------|----------|",
        ]
        for label, scores in self.per_modality_scores.items():
            lines.append(f"| {label} | {len(scores)} | {len(scores.columns)} |")

        if self.sample_overlap:
            lines.extend(
                [
                    "",
                    "### Sample Overlap",
                    "",
                    f"- Shared across all: {self.sample_overlap.n_shared}",
                    f"- Union (any modality): {self.sample_overlap.n_union}",
                    f"- Overlap fraction: {self.sample_overlap.overlap_fraction:.1%}",
                ]
            )

        if self.modality_weights:
            lines.extend(["", "### Modality Weights", ""])
            for label, weight in self.modality_weights.items():
                lines.append(f"- {label}: {weight:.3f}")

        qr = self.quality_report
        if qr and qr.warnings:
            lines.extend(["", "### Warnings", ""])
            for w in qr.warnings:
                lines.append(f"- {w}")

        lines.append("")
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return method-specific citations."""
        citations = [
            "Ritchie MD, et al. Methods of integrating data to uncover "
            "genotype-phenotype interactions. Nat Rev Genet. "
            "2015;16(2):85-97.",
        ]
        if self.strategy == FusionStrategy.WEIGHTED_AVERAGE:
            citations.append(
                "Speicher NK, Pfeifer N. Integrating different data types "
                "by regularized unsupervised multiple kernel learning with "
                "application to cancer subtype discovery. "
                "Bioinformatics. 2015;31(12):i268-i275."
            )
        return citations


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _z_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score normalize each column of a DataFrame.

    Zero-variance columns are left as zeros.
    """
    means = df.mean()
    stds = df.std()
    stds = stds.replace(0, 1.0)  # avoid division by zero
    return (df - means) / stds


def _resolve_labels(modalities: List[ModalityInput]) -> List[str]:
    """Assign unique labels to modalities, using defaults if not provided."""
    labels: List[str] = []
    counts: Dict[str, int] = {}
    for m in modalities:
        base = m.label or m.modality_type.value
        if base in counts:
            counts[base] += 1
            labels.append(f"{base}_{counts[base]}")
        else:
            counts[base] = 1
            labels.append(base)
    return labels


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def prepare_modality(
    modality_type: ModalityType,
    pathway_scores: pd.DataFrame,
    gene_data: Optional[pd.DataFrame] = None,
    quality_report: Optional[Any] = None,
    label: Optional[str] = None,
) -> ModalityInput:
    """Wrap a single modality's outputs into a ModalityInput.

    Validates the pathway_scores DataFrame and creates a ModalityInput
    ready for fusion.

    Args:
        modality_type: Which omic modality this represents.
        pathway_scores: Z-normalized pathway scores (samples x pathways).
        gene_data: Gene-level data for characterization (optional).
        quality_report: Quality report from the modality's scoring step.
        label: Human-readable label (defaults to modality_type.value).

    Returns:
        ModalityInput ready for fusion.

    Raises:
        ValueError: If pathway_scores is empty or invalid.
    """
    if pathway_scores is None or pathway_scores.empty:
        raise ValueError(f"[MultiOmic] pathway_scores for {modality_type.value} is empty or None")

    if pathway_scores.shape[0] < 2:
        raise ValueError(
            f"[MultiOmic] pathway_scores for {modality_type.value} has "
            f"{pathway_scores.shape[0]} sample(s); need at least 2"
        )

    if pathway_scores.shape[1] < 1:
        raise ValueError(f"[MultiOmic] pathway_scores for {modality_type.value} has no pathways")

    # Ensure string index
    pathway_scores = pathway_scores.copy()
    pathway_scores.index = pathway_scores.index.astype(str)

    result_label = label or modality_type.value
    logger.info(
        f"[MultiOmic] Prepared {result_label} modality: "
        f"{pathway_scores.shape[0]} samples x {pathway_scores.shape[1]} pathways"
    )

    return ModalityInput(
        modality_type=modality_type,
        pathway_scores=pathway_scores,
        gene_data=gene_data,
        quality_report=quality_report,
        label=result_label,
    )


def compute_sample_overlap(
    modalities: List[ModalityInput],
) -> SampleOverlapStats:
    """Analyze sample overlap across modalities.

    Args:
        modalities: List of ModalityInput objects.

    Returns:
        SampleOverlapStats with overlap analysis.
    """
    labels = _resolve_labels(modalities)
    sample_sets: Dict[str, set] = {}
    samples_per_modality: Dict[str, List[str]] = {}
    n_samples_per_modality: Dict[str, int] = {}

    for label, m in zip(labels, modalities):
        ids = set(m.sample_ids)
        sample_sets[label] = ids
        samples_per_modality[label] = m.sample_ids
        n_samples_per_modality[label] = len(ids)

    # Intersection (shared across all)
    shared = set.intersection(*sample_sets.values()) if sample_sets else set()
    # Union (any modality)
    union = set.union(*sample_sets.values()) if sample_sets else set()

    # Pairwise overlap
    pairwise: Dict[str, int] = {}
    for (la, sa), (lb, sb) in itertools.combinations(sample_sets.items(), 2):
        key = f"{la}_vs_{lb}"
        pairwise[key] = len(sa & sb)

    overlap_fraction = len(shared) / len(union) if union else 0.0

    if overlap_fraction < 0.5 and len(modalities) >= 2:
        logger.warning(
            f"[MultiOmic] Low sample overlap: {len(shared)}/{len(union)} "
            f"({overlap_fraction:.1%}) shared across all modalities"
        )

    return SampleOverlapStats(
        n_modalities=len(modalities),
        samples_per_modality=samples_per_modality,
        n_samples_per_modality=n_samples_per_modality,
        shared_samples=sorted(shared),
        n_shared=len(shared),
        union_samples=sorted(union),
        n_union=len(union),
        pairwise_overlap=pairwise,
        overlap_fraction=overlap_fraction,
    )


def fuse_modalities(
    modalities: List[ModalityInput],
    strategy: FusionStrategy = FusionStrategy.CONCATENATE,
    weights: Optional[Dict[str, float]] = None,
    missing_strategy: MissingStrategy = MissingStrategy.IMPUTE_ZERO,
    renormalize: bool = True,
    seed: Optional[int] = None,
) -> MultiOmicFusionResult:
    """Fuse pathway scores from multiple omic modalities.

    Aligns samples across modalities, applies the chosen fusion strategy,
    and produces a single fused_pathway_scores DataFrame compatible with
    downstream clustering, validation, and characterization.

    Args:
        modalities: List of ModalityInput objects (at least 2).
        strategy: Fusion strategy to use.
        weights: Modality weights for WEIGHTED_AVERAGE strategy.
            Dict mapping modality label -> weight. Must sum to 1.0.
            If None, equal weights are used.
        missing_strategy: How to handle samples absent from some modalities.
            Ignored for INTERSECTION_ONLY strategy.
        renormalize: Whether to Z-normalize the fused output.
        seed: Random seed for reproducibility.

    Returns:
        MultiOmicFusionResult with fused pathway scores.

    Raises:
        ValueError: If fewer than 2 modalities, invalid weights, etc.
    """
    if len(modalities) < 2:
        raise ValueError(
            f"[MultiOmic] Need at least 2 modalities for fusion, got {len(modalities)}"
        )

    labels = _resolve_labels(modalities)
    # Assign resolved labels back to modalities for consistency
    for m, label in zip(modalities, labels):
        m.label = label

    logger.info(
        f"[MultiOmic] Fusing {len(modalities)} modalities "
        f"({', '.join(labels)}) using {strategy.value}"
    )

    # Compute sample overlap
    overlap = compute_sample_overlap(modalities)

    # Compute pathway overlap
    pathway_sets = [set(m.pathway_names) for m in modalities]
    shared_pathways = sorted(set.intersection(*pathway_sets)) if pathway_sets else []
    union_pathways = sorted(set.union(*pathway_sets)) if pathway_sets else []

    # Build quality report
    warnings: List[str] = []
    per_modality_reports: Dict[str, Any] = {}
    for m in modalities:
        per_modality_reports[m.label] = {
            "modality_type": m.modality_type.value,
            "n_samples": m.n_samples,
            "n_pathways": m.n_pathways,
        }

    if overlap.n_shared == 0 and strategy != FusionStrategy.CONCATENATE:
        raise ValueError(
            "[MultiOmic] No shared samples across modalities. "
            "Use CONCATENATE strategy with imputation, or check sample IDs."
        )

    if len(shared_pathways) == 0 and strategy in (
        FusionStrategy.WEIGHTED_AVERAGE,
        FusionStrategy.INTERSECTION_ONLY,
    ):
        raise ValueError(
            "[MultiOmic] No shared pathways across modalities. "
            "Use CONCATENATE strategy or ensure modalities share the same GMT."
        )

    if overlap.overlap_fraction < 0.5:
        warnings.append(
            f"Low sample overlap: {overlap.n_shared}/{overlap.n_union} "
            f"({overlap.overlap_fraction:.1%})"
        )

    if len(shared_pathways) < len(union_pathways) * 0.5:
        warnings.append(
            f"Low pathway overlap: {len(shared_pathways)}/{len(union_pathways)} "
            f"pathways shared across all modalities"
        )

    # Dispatch to strategy
    if strategy == FusionStrategy.CONCATENATE:
        fused, effective_weights = _fuse_concatenate(modalities, missing_strategy, renormalize)
    elif strategy == FusionStrategy.WEIGHTED_AVERAGE:
        fused, effective_weights = _fuse_weighted_average(
            modalities, weights, missing_strategy, renormalize
        )
    elif strategy == FusionStrategy.INTERSECTION_ONLY:
        fused, effective_weights = _fuse_intersection_only(modalities, renormalize)
    else:
        raise ValueError(f"[MultiOmic] Unknown fusion strategy: {strategy}")

    # Fuse gene-level data
    fused_gene_data = _fuse_gene_data(modalities, list(fused.index))

    quality_report = MultiOmicQualityReport(
        per_modality_reports=per_modality_reports,
        sample_overlap=overlap,
        pathway_overlap_count=len(shared_pathways),
        pathway_union_count=len(union_pathways),
        pathway_overlap_names=shared_pathways,
        warnings=warnings,
        is_usable=len(fused) >= 4 and len(fused.columns) >= 2,
    )

    per_modality_scores = {m.label: m.pathway_scores for m in modalities}

    logger.info(
        f"[MultiOmic] Fused result: {fused.shape[0]} samples x " f"{fused.shape[1]} pathways"
    )

    return MultiOmicFusionResult(
        fused_pathway_scores=fused,
        fused_gene_data=fused_gene_data,
        per_modality_scores=per_modality_scores,
        strategy=strategy,
        quality_report=quality_report,
        sample_overlap=overlap,
        modality_weights=effective_weights,
        n_modalities=len(modalities),
        modality_labels=labels,
    )


def correlation_analysis(
    modalities: List[ModalityInput],
    method: str = "pearson",
) -> pd.DataFrame:
    """Compute cross-modality correlations for shared pathways.

    For each pair of modalities, computes the correlation between
    pathway scores for pathways present in both, using only samples
    present in both.

    Args:
        modalities: List of ModalityInput objects.
        method: Correlation method ("pearson" or "spearman").

    Returns:
        DataFrame with columns [modality_a, modality_b, pathway,
        correlation, p_value, n_samples].
    """
    if method not in ("pearson", "spearman"):
        raise ValueError(f"[MultiOmic] method must be 'pearson' or 'spearman', got '{method}'")

    labels = _resolve_labels(modalities)
    results: List[Dict[str, Any]] = []

    for (i, m_a), (j, m_b) in itertools.combinations(enumerate(modalities), 2):
        label_a, label_b = labels[i], labels[j]

        shared_pathways = sorted(set(m_a.pathway_names) & set(m_b.pathway_names))
        shared_samples = sorted(set(m_a.sample_ids) & set(m_b.sample_ids))

        if len(shared_samples) < 3 or len(shared_pathways) == 0:
            logger.warning(
                f"[MultiOmic] Skipping correlation for {label_a} vs {label_b}: "
                f"{len(shared_samples)} shared samples, "
                f"{len(shared_pathways)} shared pathways"
            )
            continue

        scores_a = m_a.pathway_scores.loc[shared_samples, shared_pathways]
        scores_b = m_b.pathway_scores.loc[shared_samples, shared_pathways]

        corr_fn = stats.pearsonr if method == "pearson" else stats.spearmanr

        for pathway in shared_pathways:
            a_vals = scores_a[pathway].values
            b_vals = scores_b[pathway].values

            # Skip if either has zero variance
            if np.std(a_vals) == 0 or np.std(b_vals) == 0:
                continue

            corr, pval = corr_fn(a_vals, b_vals)
            results.append(
                {
                    "modality_a": label_a,
                    "modality_b": label_b,
                    "pathway": pathway,
                    "correlation": float(corr),
                    "p_value": float(pval),
                    "n_samples": len(shared_samples),
                }
            )

    if results:
        df = pd.DataFrame(results)
        mean_corr = df["correlation"].mean()
        logger.info(
            f"[MultiOmic] Cross-modality correlation: "
            f"mean r = {mean_corr:.3f} across {len(results)} pathway pairs"
        )
        return df

    logger.warning("[MultiOmic] No pathway pairs available for correlation analysis")
    return pd.DataFrame(
        columns=["modality_a", "modality_b", "pathway", "correlation", "p_value", "n_samples"]
    )


# ---------------------------------------------------------------------------
# Strategy implementations
# ---------------------------------------------------------------------------


def _fuse_concatenate(
    modalities: List[ModalityInput],
    missing_strategy: MissingStrategy,
    renormalize: bool,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Concatenate modality scores with prefixed column names."""
    # Determine sample universe
    if missing_strategy == MissingStrategy.DROP:
        sample_sets = [set(m.sample_ids) for m in modalities]
        sample_universe = sorted(set.intersection(*sample_sets))
    else:
        sample_sets = [set(m.sample_ids) for m in modalities]
        sample_universe = sorted(set.union(*sample_sets))

    if len(sample_universe) == 0:
        raise ValueError("[MultiOmic] No samples available after alignment")

    frames: List[pd.DataFrame] = []
    for m in modalities:
        prefixed = m.pathway_scores.copy()
        prefixed.columns = [f"{m.label}:{col}" for col in prefixed.columns]
        aligned = prefixed.reindex(sample_universe)

        if missing_strategy == MissingStrategy.IMPUTE_ZERO:
            aligned = aligned.fillna(0.0)
        elif missing_strategy == MissingStrategy.IMPUTE_MEAN:
            aligned = aligned.fillna(aligned.mean())

        frames.append(aligned)

    fused = pd.concat(frames, axis=1)

    if renormalize:
        fused = _z_normalize(fused)

    # Equal implicit weight for concatenation
    weights = {m.label: 1.0 / len(modalities) for m in modalities}
    return fused, weights


def _fuse_weighted_average(
    modalities: List[ModalityInput],
    weights: Optional[Dict[str, float]],
    missing_strategy: MissingStrategy,
    renormalize: bool,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Weighted average of shared pathway scores."""
    # Identify shared pathways
    pathway_sets = [set(m.pathway_names) for m in modalities]
    shared_pathways = sorted(set.intersection(*pathway_sets))

    if len(shared_pathways) < 2:
        raise ValueError(
            f"[MultiOmic] Only {len(shared_pathways)} shared pathway(s) for "
            f"WEIGHTED_AVERAGE; need at least 2"
        )

    # Warn about dropped pathways
    all_pathways = set.union(*pathway_sets)
    n_dropped = len(all_pathways) - len(shared_pathways)
    if n_dropped > 0:
        logger.warning(
            f"[MultiOmic] WEIGHTED_AVERAGE drops {n_dropped} modality-specific "
            f"pathways; keeping {len(shared_pathways)} shared pathways"
        )

    # Determine sample universe
    if missing_strategy == MissingStrategy.DROP:
        sample_sets = [set(m.sample_ids) for m in modalities]
        sample_universe = sorted(set.intersection(*sample_sets))
    else:
        sample_sets = [set(m.sample_ids) for m in modalities]
        sample_universe = sorted(set.union(*sample_sets))

    if len(sample_universe) == 0:
        raise ValueError("[MultiOmic] No samples available after alignment")

    # Resolve weights
    if weights is None:
        effective_weights = {m.label: 1.0 / len(modalities) for m in modalities}
    else:
        effective_weights = dict(weights)
        # Validate weights
        for m in modalities:
            if m.label not in effective_weights:
                raise ValueError(
                    f"[MultiOmic] Weight missing for modality '{m.label}'. "
                    f"Provided weights: {list(effective_weights.keys())}"
                )
        weight_sum = sum(effective_weights.values())
        if abs(weight_sum - 1.0) > 0.01:
            raise ValueError(f"[MultiOmic] Weights must sum to 1.0, got {weight_sum:.4f}")

    # Compute weighted average
    fused = pd.DataFrame(0.0, index=sample_universe, columns=shared_pathways)

    for m in modalities:
        w = effective_weights[m.label]
        aligned = m.pathway_scores[shared_pathways].reindex(sample_universe)

        if missing_strategy == MissingStrategy.IMPUTE_ZERO:
            aligned = aligned.fillna(0.0)
        elif missing_strategy == MissingStrategy.IMPUTE_MEAN:
            aligned = aligned.fillna(aligned.mean())

        fused = fused + w * aligned

    if renormalize:
        fused = _z_normalize(fused)

    return fused, effective_weights


def _fuse_intersection_only(
    modalities: List[ModalityInput],
    renormalize: bool,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Average scores restricted to shared samples and pathways."""
    # Shared samples
    sample_sets = [set(m.sample_ids) for m in modalities]
    shared_samples = sorted(set.intersection(*sample_sets))

    if len(shared_samples) < 2:
        raise ValueError(
            f"[MultiOmic] Only {len(shared_samples)} shared sample(s) for "
            f"INTERSECTION_ONLY; need at least 2"
        )

    # Shared pathways
    pathway_sets = [set(m.pathway_names) for m in modalities]
    shared_pathways = sorted(set.intersection(*pathway_sets))

    if len(shared_pathways) < 2:
        raise ValueError(
            f"[MultiOmic] Only {len(shared_pathways)} shared pathway(s) for "
            f"INTERSECTION_ONLY; need at least 2"
        )

    # Average across modalities
    n = len(modalities)
    fused = pd.DataFrame(0.0, index=shared_samples, columns=shared_pathways)
    for m in modalities:
        fused = fused + m.pathway_scores.loc[shared_samples, shared_pathways]
    fused = fused / n

    if renormalize:
        fused = _z_normalize(fused)

    weights = {m.label: 1.0 / n for m in modalities}
    return fused, weights


def _fuse_gene_data(
    modalities: List[ModalityInput],
    sample_universe: List[str],
) -> Optional[pd.DataFrame]:
    """Fuse gene-level data across modalities."""
    frames = []
    for m in modalities:
        if m.gene_data is not None and not m.gene_data.empty:
            # Prefix gene columns with modality label to avoid collisions
            prefixed = m.gene_data.copy()
            prefixed.index = prefixed.index.astype(str)
            prefixed.columns = [f"{m.label}:{col}" for col in prefixed.columns]
            frames.append(prefixed)

    if not frames:
        return None

    combined = pd.concat(frames, axis=1)
    aligned = combined.reindex(sample_universe).fillna(0.0)
    return aligned
