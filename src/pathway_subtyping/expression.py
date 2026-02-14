"""
Bulk RNA-seq Expression Pathway Scoring

Computes pathway-level scores from gene expression matrices using
ssGSEA, GSVA, or mean-Z methods. Produces the same pathway_scores
DataFrame format as variant-based scoring, enabling unified
downstream clustering and validation.

Methods:
    - mean_z: Z-score normalize genes, average per pathway (fast, simple)
    - ssgsea: Single-sample GSEA with rank-based enrichment (recommended)
    - gsva: Simplified Gene Set Variation Analysis (empirical CDF + KS)

For publication-grade GSVA, precompute scores using the R GSVA package
and import the resulting matrix with input_type=LOG2.
"""

import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


class ExpressionScoringMethod(Enum):
    """Pathway scoring methods for expression data."""

    MEAN_Z = "mean_z"
    SSGSEA = "ssgsea"
    GSVA = "gsva"


class ExpressionInputType(Enum):
    """Type of expression data."""

    COUNTS = "counts"
    TPM = "tpm"
    FPKM = "fpkm"
    LOG2 = "log2"


@dataclass
class ExpressionDataQualityReport:
    """Report on expression data quality after loading."""

    n_samples: int = 0
    n_genes: int = 0
    n_genes_before_filter: int = 0
    n_zero_genes: int = 0
    n_low_variance_genes: int = 0
    n_pathways_covered: int = 0
    n_pathways_total: int = 0
    mean_pathway_gene_coverage: float = 0.0
    input_type: str = "unknown"
    orientation_detected: str = ""
    was_transposed: bool = False
    warnings: List[str] = field(default_factory=list)
    is_usable: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_samples": self.n_samples,
            "n_genes": self.n_genes,
            "n_genes_before_filter": self.n_genes_before_filter,
            "n_zero_genes": self.n_zero_genes,
            "n_low_variance_genes": self.n_low_variance_genes,
            "n_pathways_covered": self.n_pathways_covered,
            "n_pathways_total": self.n_pathways_total,
            "mean_pathway_gene_coverage": round(
                float(self.mean_pathway_gene_coverage), 4
            ),
            "input_type": self.input_type,
            "orientation_detected": self.orientation_detected,
            "was_transposed": self.was_transposed,
            "warnings": list(self.warnings),
            "is_usable": self.is_usable,
        }


@dataclass
class ExpressionScoringResult:
    """Result of pathway scoring from expression data."""

    pathway_scores: pd.DataFrame
    gene_expression: pd.DataFrame
    method: ExpressionScoringMethod
    quality_report: ExpressionDataQualityReport
    n_pathways_scored: int = 0
    n_pathways_skipped: int = 0
    skipped_pathways: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "method": self.method.value,
            "n_samples": len(self.pathway_scores),
            "n_pathways_scored": self.n_pathways_scored,
            "n_pathways_skipped": self.n_pathways_skipped,
            "skipped_pathways": list(self.skipped_pathways),
            "pathway_names": list(self.pathway_scores.columns),
            "quality_report": self.quality_report.to_dict(),
        }

    def format_report(self) -> str:
        """Generate a human-readable report."""
        lines = [
            "## Expression Pathway Scoring Report",
            "",
            f"**Method:** {self.method.value}",
            f"**Samples:** {len(self.pathway_scores)}",
            f"**Genes:** {len(self.gene_expression.columns)}",
            f"**Pathways scored:** {self.n_pathways_scored}",
            f"**Pathways skipped:** {self.n_pathways_skipped}",
            "",
        ]
        if self.skipped_pathways:
            lines.append("### Skipped Pathways")
            lines.append("")
            for pw in self.skipped_pathways[:10]:
                lines.append(f"- {pw}")
            if len(self.skipped_pathways) > 10:
                lines.append(
                    f"- ... and {len(self.skipped_pathways) - 10} more"
                )
            lines.append("")
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return method-specific citations."""
        citations = []
        if self.method == ExpressionScoringMethod.SSGSEA:
            citations.append(
                "Barbie DA, et al. Systematic RNA interference reveals "
                "that oncogenic KRAS-driven cancers require TBK1. "
                "Nature. 2009;462(7269):108-112."
            )
        elif self.method == ExpressionScoringMethod.GSVA:
            citations.append(
                "Hanzelmann S, Castelo R, Guinney J. GSVA: gene set "
                "variation analysis for microarray and RNA-seq data. "
                "BMC Bioinformatics. 2013;14:7."
            )
        elif self.method == ExpressionScoringMethod.MEAN_Z:
            citations.append(
                "Lee E, et al. Inferring pathway activity toward precise "
                "disease classification. PLoS Comput Biol. "
                "2008;4(11):e1000217."
            )
        return citations


# ---------------------------------------------------------------------------
# Orientation detection
# ---------------------------------------------------------------------------

_GENE_SYMBOL_PATTERN = re.compile(r"^[A-Z][A-Z0-9]{1,14}$")


def _looks_like_gene_symbols(values: list, threshold: float = 0.5) -> bool:
    """Check whether a list of values looks like gene symbols."""
    if not values:
        return False
    str_values = [str(v) for v in values]
    matches = sum(1 for v in str_values if _GENE_SYMBOL_PATTERN.match(v))
    return (matches / len(str_values)) >= threshold


def _detect_orientation(
    df: pd.DataFrame,
    gene_column: Optional[str] = None,
) -> str:
    """
    Detect whether genes are rows or columns.

    Returns 'genes_as_columns' or 'genes_as_rows'.
    """
    if gene_column is not None:
        if gene_column in df.columns:
            return "genes_as_rows"
        logger.warning(
            f"[Expression] Specified gene_column '{gene_column}' not found "
            f"in columns. Falling back to auto-detection."
        )

    index_values = list(df.index[:100])
    col_values = list(df.columns[:100])

    index_looks_like_genes = _looks_like_gene_symbols(index_values)
    cols_look_like_genes = _looks_like_gene_symbols(col_values)

    if index_looks_like_genes and not cols_look_like_genes:
        return "genes_as_rows"
    if cols_look_like_genes and not index_looks_like_genes:
        return "genes_as_columns"

    # Heuristic: expression matrices typically have more genes than samples
    if df.shape[0] > df.shape[1] * 2:
        return "genes_as_rows"

    return "genes_as_columns"


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------


def _preprocess_expression(
    gene_expression: pd.DataFrame,
    input_type: ExpressionInputType,
) -> pd.DataFrame:
    """
    Preprocess expression data based on input type.

    - COUNTS: log2(x + 1) transformation
    - TPM/FPKM: log2(x + 1) if values suggest non-log scale (max > 20)
    - LOG2: no transformation
    """
    df = gene_expression.copy()

    if input_type == ExpressionInputType.LOG2:
        return df

    if input_type == ExpressionInputType.COUNTS:
        logger.info("[Expression] Applying log2(x + 1) transformation to counts")
        df = np.log2(df + 1)
        return df

    # TPM or FPKM: check if already log-transformed
    max_val = df.values.max()
    if max_val > 20:
        logger.info(
            f"[Expression] Values appear non-log-scaled (max={max_val:.1f}). "
            f"Applying log2(x + 1) transformation."
        )
        df = np.log2(df + 1)
    else:
        logger.info(
            f"[Expression] Values appear already log-scaled (max={max_val:.1f}). "
            f"No transformation applied."
        )

    return df


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_expression_matrix(
    path: str,
    input_type: ExpressionInputType = ExpressionInputType.TPM,
    gene_column: Optional[str] = None,
    sample_column: Optional[str] = None,
    min_genes_per_sample: int = 100,
    min_samples_per_gene: int = 3,
) -> Tuple[pd.DataFrame, ExpressionDataQualityReport]:
    """
    Load and validate an expression matrix from CSV/TSV.

    Auto-detects orientation (genes as rows vs columns).
    Applies basic QC filtering (zero-expression genes, low-variance).
    Handles counts by log2(x+1) transformation.

    Args:
        path: Path to CSV or TSV expression file.
        input_type: Type of expression values.
        gene_column: Column containing gene symbols (if genes are rows).
        sample_column: Column containing sample IDs (if genes are columns).
        min_genes_per_sample: Minimum genes required per sample.
        min_samples_per_gene: Minimum non-zero samples per gene.

    Returns:
        Tuple of (gene_expression DataFrame [samples x genes], quality report).

    Raises:
        FileNotFoundError: If file does not exist.
        ValueError: If data is empty or unusable.
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"Expression file not found: {path}")

    report = ExpressionDataQualityReport(input_type=input_type.value)

    # Read file (detect separator)
    sep = "\t" if file_path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(file_path, sep=sep, index_col=0)

    if df.empty:
        raise ValueError(f"Expression file is empty: {path}")

    # Detect orientation
    orientation = _detect_orientation(df, gene_column=gene_column)
    report.orientation_detected = orientation

    if orientation == "genes_as_rows":
        if gene_column and gene_column in df.columns:
            df = df.set_index(gene_column)
        df = df.T
        report.was_transposed = True
        logger.info(
            "[Expression] Detected genes as rows — transposed to "
            f"({df.shape[0]} samples x {df.shape[1]} genes)"
        )
    else:
        logger.info(
            f"[Expression] Detected genes as columns "
            f"({df.shape[0]} samples x {df.shape[1]} genes)"
        )

    # Ensure numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    report.n_genes_before_filter = len(df.columns)

    # Preprocess (log transform if needed)
    df = _preprocess_expression(df, input_type)

    # QC: remove all-zero genes
    gene_nonzero = (df != 0).sum(axis=0)
    zero_genes = gene_nonzero[gene_nonzero == 0].index
    report.n_zero_genes = len(zero_genes)
    if len(zero_genes) > 0:
        df = df.drop(columns=zero_genes)
        logger.info(f"[Expression] Removed {len(zero_genes)} all-zero genes")

    # QC: remove low-occurrence genes
    low_genes = gene_nonzero[
        (gene_nonzero > 0) & (gene_nonzero < min_samples_per_gene)
    ].index
    low_genes = [g for g in low_genes if g in df.columns]
    if len(low_genes) > 0:
        df = df.drop(columns=low_genes)
        logger.info(
            f"[Expression] Removed {len(low_genes)} genes expressed in "
            f"<{min_samples_per_gene} samples"
        )

    # QC: remove zero-variance genes
    gene_vars = df.var()
    low_var_genes = gene_vars[gene_vars == 0].index
    report.n_low_variance_genes = len(low_var_genes)
    if len(low_var_genes) > 0:
        df = df.drop(columns=low_var_genes)
        logger.info(
            f"[Expression] Removed {len(low_var_genes)} zero-variance genes"
        )

    report.n_samples = len(df)
    report.n_genes = len(df.columns)

    # Usability checks
    if report.n_genes < min_genes_per_sample:
        report.is_usable = False
        report.warnings.append(
            f"Only {report.n_genes} genes after filtering "
            f"(minimum {min_genes_per_sample} required)"
        )

    if report.n_samples < 2:
        report.is_usable = False
        report.warnings.append(
            f"Only {report.n_samples} sample(s). Need at least 2."
        )

    logger.info(
        f"[Expression] Loaded {report.n_samples} samples x "
        f"{report.n_genes} genes (from {report.n_genes_before_filter} raw)"
    )

    return df, report


# ---------------------------------------------------------------------------
# Scoring: Mean-Z
# ---------------------------------------------------------------------------


def _score_mean_z(
    gene_expression: pd.DataFrame,
    pathways: Dict[str, List[str]],
    min_genes: int = 2,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Mean-Z pathway scoring.

    1. Z-score normalize each gene across samples.
    2. For each pathway: mean of member gene Z-scores.

    Returns:
        (pathway_scores DataFrame, list of skipped pathway names)
    """
    # Z-score normalize genes (columns)
    means = gene_expression.mean(axis=0)
    stds = gene_expression.std(axis=0)
    stds = stds.replace(0, 1e-10)
    z_expression = (gene_expression - means) / stds

    scores = {}
    skipped = []

    for pathway_name, pathway_genes in pathways.items():
        common = [g for g in pathway_genes if g in z_expression.columns]
        if len(common) < min_genes:
            skipped.append(pathway_name)
            continue
        scores[pathway_name] = z_expression[common].mean(axis=1)

    if not scores:
        return pd.DataFrame(index=gene_expression.index), skipped

    return pd.DataFrame(scores, index=gene_expression.index), skipped


# ---------------------------------------------------------------------------
# Scoring: ssGSEA
# ---------------------------------------------------------------------------


def _score_ssgsea(
    gene_expression: pd.DataFrame,
    pathways: Dict[str, List[str]],
    alpha: float = 0.25,
    min_genes: int = 2,
    show_progress: bool = True,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Single-sample GSEA scoring.

    For each sample:
        1. Rank genes by expression.
        2. For each pathway, compute weighted enrichment score:
           - Walk along ranked list.
           - Step up by |rank|^alpha for pathway genes.
           - Step down by 1/(n_total - n_pathway) for non-pathway genes.
           - Score = sum of running enrichment statistic.

    Args:
        gene_expression: DataFrame (samples x genes), preprocessed.
        pathways: Pathway name -> gene list.
        alpha: Weighting exponent for ranks (default 0.25).
        min_genes: Skip pathways with fewer overlapping genes.

    Returns:
        (pathway_scores DataFrame, list of skipped pathway names)
    """
    n_samples, n_genes = gene_expression.shape
    all_genes = gene_expression.columns.tolist()

    # Pre-compute ranks for all samples (higher expression = higher rank)
    # rank() gives 1-based ranks; ascending=True means smallest gets rank 1
    ranks = gene_expression.rank(axis=1, method="average", ascending=True)
    ranks_array = ranks.values  # (n_samples, n_genes)

    scores = {}
    skipped = []

    for pathway_name, pathway_genes in tqdm(
        pathways.items(),
        desc="ssGSEA",
        total=len(pathways),
        disable=not show_progress,
    ):
        common = [g for g in pathway_genes if g in all_genes]
        if len(common) < min_genes:
            skipped.append(pathway_name)
            continue

        # Boolean mask for pathway genes
        gene_in_pathway = np.array([g in common for g in all_genes])
        n_hit = gene_in_pathway.sum()
        n_miss = n_genes - n_hit

        if n_miss == 0:
            # All genes in pathway — score is 0
            scores[pathway_name] = np.zeros(n_samples)
            continue

        # Compute per-sample enrichment scores (vectorized across samples)
        sample_scores = np.zeros(n_samples)

        for i in range(n_samples):
            sample_ranks = ranks_array[i]

            # Sort genes by rank (descending — highest expressed first)
            sort_idx = np.argsort(-sample_ranks)
            sorted_in_pathway = gene_in_pathway[sort_idx]
            sorted_ranks = sample_ranks[sort_idx]

            # Weighted step-up for pathway genes
            weights = np.abs(sorted_ranks) ** alpha
            hit_weights = np.where(sorted_in_pathway, weights, 0.0)
            hit_sum = hit_weights.sum()
            if hit_sum == 0:
                sample_scores[i] = 0.0
                continue

            # Running enrichment statistic
            p_hit = np.cumsum(hit_weights) / hit_sum
            p_miss = np.cumsum(~sorted_in_pathway) / n_miss

            # Enrichment score: sum of (p_hit - p_miss)
            es = (p_hit - p_miss).sum()
            sample_scores[i] = es

        scores[pathway_name] = sample_scores

    if not scores:
        return pd.DataFrame(index=gene_expression.index), skipped

    return pd.DataFrame(scores, index=gene_expression.index), skipped


# ---------------------------------------------------------------------------
# Scoring: GSVA (simplified)
# ---------------------------------------------------------------------------


def _score_gsva(
    gene_expression: pd.DataFrame,
    pathways: Dict[str, List[str]],
    min_genes: int = 2,
    show_progress: bool = True,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Simplified GSVA scoring.

    1. For each gene, compute empirical CDF across samples.
    2. Rank CDF values within each sample.
    3. For each pathway, compute KS-like statistic: pathway genes
       vs. non-pathway genes within each sample.

    This is a simplified implementation. For publication-grade GSVA,
    use the R GSVA package and import precomputed scores.

    Returns:
        (pathway_scores DataFrame, list of skipped pathway names)
    """
    n_samples, n_genes = gene_expression.shape
    all_genes = gene_expression.columns.tolist()

    # Step 1: Transform to empirical CDF per gene (column-wise ranking)
    # Each gene's values are replaced by their rank percentile across samples
    ecdf = gene_expression.rank(axis=0, method="average", pct=True)
    ecdf_array = ecdf.values  # (n_samples, n_genes)

    scores = {}
    skipped = []

    for pathway_name, pathway_genes in tqdm(
        pathways.items(),
        desc="GSVA",
        total=len(pathways),
        disable=not show_progress,
    ):
        common = [g for g in pathway_genes if g in all_genes]
        if len(common) < min_genes:
            skipped.append(pathway_name)
            continue

        gene_in_pathway = np.array([g in common for g in all_genes])
        n_hit = gene_in_pathway.sum()
        n_miss = n_genes - n_hit

        if n_miss == 0:
            scores[pathway_name] = np.zeros(n_samples)
            continue

        sample_scores = np.zeros(n_samples)

        for i in range(n_samples):
            sample_ecdf = ecdf_array[i]

            # Sort by ECDF value (descending)
            sort_idx = np.argsort(-sample_ecdf)
            sorted_in_pathway = gene_in_pathway[sort_idx]

            # KS-like statistic
            p_hit = np.cumsum(sorted_in_pathway) / n_hit
            p_miss = np.cumsum(~sorted_in_pathway) / n_miss
            diff = p_hit - p_miss

            # Max deviation (signed: positive = enriched, negative = depleted)
            pos_max = diff.max()
            neg_min = diff.min()
            if abs(pos_max) >= abs(neg_min):
                sample_scores[i] = pos_max
            else:
                sample_scores[i] = neg_min

        scores[pathway_name] = sample_scores

    if not scores:
        return pd.DataFrame(index=gene_expression.index), skipped

    return pd.DataFrame(scores, index=gene_expression.index), skipped


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def score_pathways_from_expression(
    gene_expression: pd.DataFrame,
    pathways: Dict[str, List[str]],
    method: ExpressionScoringMethod = ExpressionScoringMethod.SSGSEA,
    min_genes_per_pathway: int = 2,
    alpha: float = 0.25,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> ExpressionScoringResult:
    """
    Compute pathway scores from a gene expression matrix.

    This is the main entry point for expression-based pathway scoring.
    Dispatches to the appropriate scoring method and Z-normalizes output.

    Args:
        gene_expression: DataFrame (n_samples x n_genes), preprocessed.
        pathways: Dict mapping pathway names to gene lists (GMT format).
        method: Scoring method to use.
        min_genes_per_pathway: Skip pathways with fewer genes in data.
        alpha: Weight parameter for ssGSEA (ignored for other methods).
        seed: Random seed for reproducibility.
        show_progress: Show tqdm progress bar for ssGSEA/GSVA scoring.

    Returns:
        ExpressionScoringResult with Z-normalized pathway scores.
    """
    if seed is not None:
        np.random.seed(seed)

    logger.info(
        f"[Expression] Scoring {len(pathways)} pathways via {method.value} "
        f"({gene_expression.shape[0]} samples, "
        f"{gene_expression.shape[1]} genes)"
    )

    # Build quality report stub (pathway coverage)
    report = ExpressionDataQualityReport(
        n_samples=len(gene_expression),
        n_genes=len(gene_expression.columns),
        n_pathways_total=len(pathways),
    )

    # Compute pathway gene coverage
    available_genes = set(gene_expression.columns)
    coverages = []
    n_covered = 0
    for pw_genes in pathways.values():
        overlap = len(set(pw_genes) & available_genes)
        total = len(pw_genes)
        if total > 0:
            coverages.append(overlap / total)
            if overlap >= min_genes_per_pathway:
                n_covered += 1
    report.n_pathways_covered = n_covered
    report.mean_pathway_gene_coverage = (
        float(np.mean(coverages)) if coverages else 0.0
    )

    # Dispatch to scoring method
    if method == ExpressionScoringMethod.MEAN_Z:
        raw_scores, skipped = _score_mean_z(
            gene_expression, pathways, min_genes=min_genes_per_pathway
        )
    elif method == ExpressionScoringMethod.SSGSEA:
        raw_scores, skipped = _score_ssgsea(
            gene_expression,
            pathways,
            alpha=alpha,
            min_genes=min_genes_per_pathway,
            show_progress=show_progress,
        )
    elif method == ExpressionScoringMethod.GSVA:
        raw_scores, skipped = _score_gsva(
            gene_expression, pathways, min_genes=min_genes_per_pathway,
            show_progress=show_progress,
        )
    else:
        raise ValueError(f"Unknown scoring method: {method}")

    n_scored = len(raw_scores.columns)
    n_skipped = len(skipped)

    if n_scored == 0:
        logger.warning(
            "[Expression] No pathways scored. Check that pathway gene names "
            "match expression matrix column names."
        )
        return ExpressionScoringResult(
            pathway_scores=raw_scores,
            gene_expression=gene_expression,
            method=method,
            quality_report=report,
            n_pathways_scored=0,
            n_pathways_skipped=n_skipped,
            skipped_pathways=skipped,
        )

    # Remove zero-variance pathways
    pathway_stds = raw_scores.std()
    zero_var = pathway_stds[pathway_stds == 0].index.tolist()
    if zero_var:
        raw_scores = raw_scores.drop(columns=zero_var)
        skipped.extend(zero_var)
        n_skipped += len(zero_var)
        n_scored -= len(zero_var)
        logger.info(
            f"[Expression] Removed {len(zero_var)} zero-variance pathway(s)"
        )

    # Z-score normalize pathway scores
    if n_scored > 0:
        means = raw_scores.mean()
        stds = raw_scores.std().replace(0, 1e-10)
        pathway_scores = (raw_scores - means) / stds
    else:
        pathway_scores = raw_scores

    logger.info(
        f"[Expression] Scored {n_scored} pathways, skipped {n_skipped} "
        f"(mean gene coverage: {report.mean_pathway_gene_coverage:.1%})"
    )

    return ExpressionScoringResult(
        pathway_scores=pathway_scores,
        gene_expression=gene_expression,
        method=method,
        quality_report=report,
        n_pathways_scored=n_scored,
        n_pathways_skipped=n_skipped,
        skipped_pathways=skipped,
    )
