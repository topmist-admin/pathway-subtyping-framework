"""
Subtype Characterization Module for the Pathway Subtyping Framework.

Provides tools to interpret and describe discovered subtypes:
- Pathway enrichment analysis (Kruskal-Wallis + FDR correction)
- Gene-level contribution scores (Cohen's d effect sizes)
- Subtype comparison heatmaps (publication-ready)
- Export to standard formats (CSV, Excel)

References:
- Kruskal & Wallis (1952). Use of ranks in one-criterion variance analysis.
  JASA, 47(260), 583-621.
- Cohen (1988). Statistical Power Analysis for the Behavioral Sciences.
  2nd ed. Lawrence Erlbaum Associates.
- Benjamini & Hochberg (1995). Controlling the False Discovery Rate.
  JRSS-B, 57(1), 289-300.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from .statistical_rigor import benjamini_hochberg

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================


@dataclass
class PathwayEnrichment:
    """
    Per-pathway enrichment result for one subtype.

    Attributes:
        pathway: Pathway name
        mean_score: Mean pathway score in this subtype
        overall_mean: Mean pathway score across all samples
        fold_change: Ratio of subtype mean to overall mean
        effect_size: Cohen's d (subtype vs rest)
        p_value: Kruskal-Wallis p-value (across all subtypes)
        q_value: FDR-corrected p-value
        significant: Whether q_value < alpha
    """

    pathway: str
    mean_score: float
    overall_mean: float
    fold_change: float
    effect_size: float
    p_value: float
    q_value: float
    significant: bool

    def to_dict(self) -> Dict[str, Any]:
        fc = self.fold_change
        if not np.isfinite(fc):
            fc = None
        else:
            fc = round(float(fc), 4)
        return {
            "pathway": self.pathway,
            "mean_score": round(float(self.mean_score), 4),
            "overall_mean": round(float(self.overall_mean), 4),
            "fold_change": fc,
            "effect_size": round(float(self.effect_size), 4),
            "p_value": round(float(self.p_value), 6),
            "q_value": round(float(self.q_value), 6),
            "significant": bool(self.significant),
        }


@dataclass
class GeneContribution:
    """
    Per-gene contribution to a subtype's pathway profile.

    Attributes:
        gene: Gene name
        pathway: Pathway the gene belongs to
        mean_burden: Mean gene burden in this subtype
        overall_mean: Mean gene burden across all samples
        fold_change: Ratio of subtype mean to overall mean
        effect_size: Cohen's d (subtype vs rest)
    """

    gene: str
    pathway: str
    mean_burden: float
    overall_mean: float
    fold_change: float
    effect_size: float

    def to_dict(self) -> Dict[str, Any]:
        fc = self.fold_change
        if not np.isfinite(fc):
            fc = None
        else:
            fc = round(float(fc), 4)
        return {
            "gene": self.gene,
            "pathway": self.pathway,
            "mean_burden": round(float(self.mean_burden), 4),
            "overall_mean": round(float(self.overall_mean), 4),
            "fold_change": fc,
            "effect_size": round(float(self.effect_size), 4),
        }


@dataclass
class SubtypeProfile:
    """
    Complete characterization of one subtype.

    Attributes:
        subtype_id: Cluster ID
        subtype_label: Human-readable label
        n_samples: Number of samples in this subtype
        fraction: Fraction of total samples
        mean_confidence: Mean GMM assignment confidence
        enriched_pathways: Pathway enrichment results sorted by effect size
        top_genes: Top contributing genes sorted by effect size
        pathway_score_means: Mean pathway score for each pathway
    """

    subtype_id: int
    subtype_label: str
    n_samples: int
    fraction: float
    mean_confidence: float
    enriched_pathways: List[PathwayEnrichment] = field(default_factory=list)
    top_genes: List[GeneContribution] = field(default_factory=list)
    pathway_score_means: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "subtype_id": int(self.subtype_id),
            "subtype_label": str(self.subtype_label),
            "n_samples": int(self.n_samples),
            "fraction": round(float(self.fraction), 4),
            "mean_confidence": round(float(self.mean_confidence), 4),
            "enriched_pathways": [p.to_dict() for p in self.enriched_pathways],
            "top_genes": [g.to_dict() for g in self.top_genes],
            "pathway_score_means": {
                k: round(float(v), 4) for k, v in self.pathway_score_means.items()
            },
        }


@dataclass
class CharacterizationResult:
    """
    Complete subtype characterization result.

    Attributes:
        subtype_profiles: One SubtypeProfile per discovered subtype
        n_subtypes: Number of subtypes
        n_samples: Total number of samples
        n_pathways: Number of pathways analyzed
        n_genes: Number of genes analyzed (0 if gene burdens not provided)
        fdr_alpha: FDR significance threshold used
    """

    subtype_profiles: List[SubtypeProfile] = field(default_factory=list)
    n_subtypes: int = 0
    n_samples: int = 0
    n_pathways: int = 0
    n_genes: int = 0
    fdr_alpha: float = 0.05

    def to_dict(self) -> Dict[str, Any]:
        return {
            "summary": {
                "n_subtypes": self.n_subtypes,
                "n_samples": self.n_samples,
                "n_pathways": self.n_pathways,
                "n_genes": self.n_genes,
                "fdr_alpha": self.fdr_alpha,
            },
            "subtype_profiles": [p.to_dict() for p in self.subtype_profiles],
        }

    def format_report(self) -> str:
        """Generate publication-ready markdown report."""
        lines = [
            "## Subtype Characterization",
            "",
            "### Summary",
            f"- **Subtypes discovered:** {self.n_subtypes}",
            f"- **Total samples:** {self.n_samples}",
            f"- **Pathways analyzed:** {self.n_pathways}",
            f"- **Genes analyzed:** {self.n_genes}",
            f"- **FDR threshold:** {self.fdr_alpha}",
            "",
        ]

        for profile in self.subtype_profiles:
            lines.extend(
                [
                    f"### Subtype {profile.subtype_id}: {profile.subtype_label}",
                    "",
                    f"- **Samples:** {profile.n_samples} " f"({profile.fraction * 100:.1f}%)",
                    f"- **Mean confidence:** {profile.mean_confidence:.3f}",
                    "",
                ]
            )

            # Enriched pathways table
            sig_pathways = [p for p in profile.enriched_pathways if p.significant]
            if sig_pathways:
                lines.extend(
                    [
                        "**Significantly enriched pathways:**",
                        "",
                        "| Pathway | Effect Size | Fold Change | q-value |",
                        "|---------|------------|-------------|---------|",
                    ]
                )
                for p in sig_pathways:
                    lines.append(
                        f"| {p.pathway} | {p.effect_size:.2f} | "
                        f"{p.fold_change:.2f} | {p.q_value:.4f} |"
                    )
                lines.append("")

            # Top genes table
            if profile.top_genes:
                lines.extend(
                    [
                        "**Top contributing genes:**",
                        "",
                        "| Gene | Pathway | Effect Size | Fold Change |",
                        "|------|---------|------------|-------------|",
                    ]
                )
                for g in profile.top_genes[:10]:
                    lines.append(
                        f"| {g.gene} | {g.pathway} | "
                        f"{g.effect_size:.2f} | {g.fold_change:.2f} |"
                    )
                lines.append("")

        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return literature citations for methods used."""
        return [
            "Kruskal WH, Wallis WA. Use of ranks in one-criterion variance "
            "analysis. JASA. 1952;47(260):583-621.",
            "Cohen J. Statistical Power Analysis for the Behavioral Sciences. "
            "2nd ed. Lawrence Erlbaum Associates; 1988.",
            "Benjamini Y, Hochberg Y. Controlling the false discovery rate: "
            "a practical and powerful approach to multiple testing. "
            "JRSS-B. 1995;57(1):289-300.",
        ]


# =============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# =============================================================================


def _cohens_d(group: np.ndarray, rest: np.ndarray) -> float:
    """
    Compute Cohen's d effect size between two groups.

    Args:
        group: Values for the subtype of interest
        rest: Values for all other samples

    Returns:
        Cohen's d (positive = group higher than rest)
    """
    n1, n2 = len(group), len(rest)
    if n1 < 2 or n2 < 2:
        return 0.0

    mean1, mean2 = np.mean(group), np.mean(rest)
    var1, var2 = np.var(group, ddof=1), np.var(rest, ddof=1)

    # Pooled standard deviation
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    pooled_sd = np.sqrt(pooled_var)

    if pooled_sd == 0:
        return 0.0

    return float((mean1 - mean2) / pooled_sd)


def pathway_enrichment_analysis(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    fdr_alpha: float = 0.05,
) -> Dict[int, List[PathwayEnrichment]]:
    """
    Run pathway enrichment analysis across all subtypes.

    For each pathway, performs a Kruskal-Wallis test across all subtypes,
    then computes Cohen's d effect size for each subtype vs the rest.
    P-values are FDR-corrected per subtype using Benjamini-Hochberg.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways)
        cluster_labels: Array of cluster assignments
        fdr_alpha: FDR significance threshold

    Returns:
        Dict mapping subtype ID to list of PathwayEnrichment results
    """
    logger.info("[Characterization] Running pathway enrichment analysis...")

    cluster_labels = np.asarray(cluster_labels)
    unique_clusters = np.unique(cluster_labels)
    pathways = list(pathway_scores.columns)

    if len(unique_clusters) < 2:
        logger.warning(
            "[Characterization] Only one cluster found — "
            "enrichment analysis requires at least 2 subtypes"
        )
        # Return single-cluster result with no significance
        cid = unique_clusters[0]
        enrichments = []
        for pw in pathways:
            vals = pathway_scores[pw].values
            mean_val = float(np.mean(vals))
            enrichments.append(
                PathwayEnrichment(
                    pathway=pw,
                    mean_score=mean_val,
                    overall_mean=mean_val,
                    fold_change=1.0,
                    effect_size=0.0,
                    p_value=1.0,
                    q_value=1.0,
                    significant=False,
                )
            )
        return {int(cid): enrichments}

    # Compute Kruskal-Wallis p-values for each pathway (across all clusters)
    kw_pvalues = {}
    for pw in pathways:
        groups = [pathway_scores[pw].values[cluster_labels == cid] for cid in unique_clusters]
        # Filter out groups with fewer than 2 samples
        valid_groups = [g for g in groups if len(g) >= 2]
        if len(valid_groups) < 2:
            kw_pvalues[pw] = 1.0
        else:
            try:
                _, p = stats.kruskal(*valid_groups)
                kw_pvalues[pw] = float(p) if not np.isnan(p) else 1.0
            except ValueError:
                kw_pvalues[pw] = 1.0

    # Per-subtype enrichment with effect sizes
    result: Dict[int, List[PathwayEnrichment]] = {}

    for cid in unique_clusters:
        mask = cluster_labels == cid
        rest_mask = ~mask

        p_values_list = []
        enrichments_raw = []

        for pw in pathways:
            vals = pathway_scores[pw].values
            group_vals = vals[mask]
            rest_vals = vals[rest_mask]

            mean_score = float(np.mean(group_vals))
            overall_mean = float(np.mean(vals))

            # Fold change (handle zero overall mean)
            if overall_mean != 0:
                fold_change = mean_score / overall_mean
            elif mean_score != 0:
                fold_change = float("inf")
            else:
                fold_change = 1.0

            effect_size = _cohens_d(group_vals, rest_vals)

            p_values_list.append(kw_pvalues[pw])
            enrichments_raw.append((pw, mean_score, overall_mean, fold_change, effect_size))

        # FDR correction
        raw_p = np.array(p_values_list)
        q_values = benjamini_hochberg(raw_p, fdr_alpha)

        enrichments = []
        for i, (pw, ms, om, fc, es) in enumerate(enrichments_raw):
            enrichments.append(
                PathwayEnrichment(
                    pathway=pw,
                    mean_score=ms,
                    overall_mean=om,
                    fold_change=fc,
                    effect_size=es,
                    p_value=float(raw_p[i]),
                    q_value=float(q_values[i]),
                    significant=bool(q_values[i] < fdr_alpha),
                )
            )

        # Sort by absolute effect size descending
        enrichments.sort(key=lambda x: abs(x.effect_size), reverse=True)
        result[int(cid)] = enrichments

    n_sig_total = sum(sum(1 for e in enrs if e.significant) for enrs in result.values())
    logger.info(
        f"[Characterization] Found {n_sig_total} significant "
        f"pathway-subtype associations (FDR < {fdr_alpha})"
    )

    return result


# =============================================================================
# GENE-LEVEL CONTRIBUTION SCORES
# =============================================================================


def gene_contribution_scores(
    gene_burdens: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathways: Dict[str, List[str]],
    top_n: int = 20,
) -> Dict[int, List[GeneContribution]]:
    """
    Identify top genes driving each subtype's pathway profile.

    For each gene in each pathway, computes Cohen's d effect size
    comparing the subtype to the rest. Returns the top N genes
    per subtype ranked by absolute effect size.

    Args:
        gene_burdens: DataFrame of gene burden scores (samples x genes)
        cluster_labels: Array of cluster assignments
        pathways: Dict mapping pathway names to gene lists
        top_n: Number of top genes to return per subtype

    Returns:
        Dict mapping subtype ID to list of GeneContribution results
    """
    logger.info("[Characterization] Computing gene contribution scores...")

    cluster_labels = np.asarray(cluster_labels)
    unique_clusters = np.unique(cluster_labels)
    available_genes = set(gene_burdens.columns)

    # Build gene-to-pathway mapping
    gene_to_pathway: Dict[str, str] = {}
    for pw_name, genes in pathways.items():
        for gene in genes:
            if gene in available_genes:
                gene_to_pathway[gene] = pw_name

    genes_in_pathways = list(gene_to_pathway.keys())
    if not genes_in_pathways:
        logger.warning("[Characterization] No pathway genes found in gene burden data")
        return {int(cid): [] for cid in unique_clusters}

    result: Dict[int, List[GeneContribution]] = {}

    for cid in unique_clusters:
        mask = cluster_labels == cid
        rest_mask = ~mask

        contributions = []
        for gene in genes_in_pathways:
            vals = gene_burdens[gene].values
            group_vals = vals[mask]
            rest_vals = vals[rest_mask]

            mean_burden = float(np.mean(group_vals))
            overall_mean = float(np.mean(vals))

            if overall_mean != 0:
                fold_change = mean_burden / overall_mean
            elif mean_burden != 0:
                fold_change = float("inf")
            else:
                fold_change = 1.0

            effect_size = _cohens_d(group_vals, rest_vals)

            contributions.append(
                GeneContribution(
                    gene=gene,
                    pathway=gene_to_pathway[gene],
                    mean_burden=mean_burden,
                    overall_mean=overall_mean,
                    fold_change=fold_change,
                    effect_size=effect_size,
                )
            )

        # Sort by absolute effect size and take top N
        contributions.sort(key=lambda x: abs(x.effect_size), reverse=True)
        result[int(cid)] = contributions[:top_n]

    logger.info(
        f"[Characterization] Analyzed {len(genes_in_pathways)} genes "
        f"across {len(pathways)} pathways"
    )

    return result


# =============================================================================
# MAIN CHARACTERIZATION FUNCTION
# =============================================================================


def characterize_subtypes(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    gene_burdens: Optional[pd.DataFrame] = None,
    pathways: Optional[Dict[str, List[str]]] = None,
    cluster_names: Optional[Dict[int, str]] = None,
    confidence_scores: Optional[np.ndarray] = None,
    fdr_alpha: float = 0.05,
    top_n_genes: int = 20,
    seed: Optional[int] = None,
) -> CharacterizationResult:
    """
    Characterize discovered subtypes with pathway enrichment and gene contributions.

    This is the main entry point for subtype characterization. It runs pathway
    enrichment analysis (Kruskal-Wallis + FDR correction) and optionally
    gene-level contribution scoring.

    Args:
        pathway_scores: DataFrame of pathway scores (samples x pathways)
        cluster_labels: Array of cluster assignments (one per sample)
        gene_burdens: Optional DataFrame of gene burdens (samples x genes)
        pathways: Optional dict mapping pathway names to gene lists
        cluster_names: Optional dict mapping cluster ID to human-readable name
        confidence_scores: Optional array of assignment confidence (e.g. GMM proba)
        fdr_alpha: FDR significance threshold
        top_n_genes: Number of top genes per subtype
        seed: Random seed (reserved for future stochastic extensions)

    Returns:
        CharacterizationResult with complete subtype profiles
    """
    logger.info("[Characterization] Starting subtype characterization...")

    cluster_labels = np.asarray(cluster_labels)
    n_samples = len(cluster_labels)

    # Guard: no samples
    if n_samples == 0:
        logger.warning("[Characterization] No samples — returning empty result")
        return CharacterizationResult(
            n_subtypes=0,
            n_samples=0,
            n_pathways=0,
            n_genes=0,
            fdr_alpha=fdr_alpha,
        )

    # Guard: empty pathway scores
    if pathway_scores is None or pathway_scores.empty:
        logger.warning("[Characterization] No pathway scores — returning empty result")
        return CharacterizationResult(
            n_subtypes=0,
            n_samples=n_samples,
            n_pathways=0,
            n_genes=0,
            fdr_alpha=fdr_alpha,
        )

    unique_clusters = np.unique(cluster_labels)

    if len(pathway_scores) != n_samples:
        raise ValueError(
            f"pathway_scores has {len(pathway_scores)} samples but "
            f"cluster_labels has {n_samples} samples"
        )

    # Default cluster names
    if cluster_names is None:
        cluster_names = {int(cid): f"Subtype_{cid}" for cid in unique_clusters}

    # Pathway enrichment analysis
    enrichment_results = pathway_enrichment_analysis(pathway_scores, cluster_labels, fdr_alpha)

    # Gene contribution scores (if data available)
    gene_results: Dict[int, List[GeneContribution]] = {}
    n_genes = 0
    if gene_burdens is not None and pathways is not None:
        gene_results = gene_contribution_scores(gene_burdens, cluster_labels, pathways, top_n_genes)
        n_genes = len(gene_burdens.columns)
    elif gene_burdens is not None and pathways is None:
        logger.info(
            "[Characterization] Gene burdens provided but no pathway "
            "definitions — skipping gene contribution analysis"
        )

    # Build subtype profiles
    profiles = []
    for cid in unique_clusters:
        mask = cluster_labels == cid
        cid_int = int(cid)

        # Mean confidence
        if confidence_scores is not None:
            mean_conf = float(np.mean(confidence_scores[mask]))
        else:
            mean_conf = 0.0

        # Pathway score means
        pw_means = {}
        for pw in pathway_scores.columns:
            pw_means[pw] = float(np.mean(pathway_scores[pw].values[mask]))

        profile = SubtypeProfile(
            subtype_id=cid_int,
            subtype_label=cluster_names.get(cid_int, f"Subtype_{cid_int}"),
            n_samples=int(np.sum(mask)),
            fraction=float(np.sum(mask)) / n_samples,
            mean_confidence=mean_conf,
            enriched_pathways=enrichment_results.get(cid_int, []),
            top_genes=gene_results.get(cid_int, []),
            pathway_score_means=pw_means,
        )
        profiles.append(profile)

    result = CharacterizationResult(
        subtype_profiles=profiles,
        n_subtypes=len(unique_clusters),
        n_samples=n_samples,
        n_pathways=len(pathway_scores.columns),
        n_genes=n_genes,
        fdr_alpha=fdr_alpha,
    )

    logger.info(
        f"[Characterization] Characterized {result.n_subtypes} subtypes "
        f"across {result.n_pathways} pathways"
    )

    return result


# =============================================================================
# VISUALIZATION — HEATMAPS
# =============================================================================


def _build_heatmap_matrix(
    result: CharacterizationResult,
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Build subtypes x pathways matrix of mean Z-scores.

    Returns:
        (matrix, row_labels, col_labels)
    """
    if not result.subtype_profiles:
        return np.array([[]]), [], []

    pathways = list(result.subtype_profiles[0].pathway_score_means.keys())
    row_labels = [f"{p.subtype_label} (n={p.n_samples})" for p in result.subtype_profiles]

    matrix = np.zeros((len(result.subtype_profiles), len(pathways)))
    for i, profile in enumerate(result.subtype_profiles):
        for j, pw in enumerate(pathways):
            matrix[i, j] = profile.pathway_score_means.get(pw, 0.0)

    return matrix, row_labels, pathways


def generate_subtype_heatmap(
    result: CharacterizationResult,
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 8),
) -> Optional[Any]:
    """
    Generate a subtypes x pathways heatmap of mean pathway Z-scores.

    Args:
        result: CharacterizationResult from characterize_subtypes()
        output_path: Optional file path to save the figure
        figsize: Figure size (width, height) in inches

    Returns:
        matplotlib Figure object, or None if matplotlib is unavailable
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning(
            "[Characterization] matplotlib not installed — " "skipping heatmap generation"
        )
        return None

    matrix, row_labels, col_labels = _build_heatmap_matrix(result)

    if matrix.size == 0:
        logger.warning("[Characterization] No data for heatmap")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    # Use diverging colormap centered at 0
    vmax = max(abs(matrix.min()), abs(matrix.max()), 0.01)
    im = ax.imshow(
        matrix,
        cmap="RdBu_r",
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
    )

    # Labels
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)

    # Annotate cells
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            color = "white" if abs(val) > vmax * 0.6 else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                color=color,
                fontsize=8,
            )

    ax.set_title("Subtype Pathway Profiles (Mean Z-scores)", fontsize=13)
    fig.colorbar(im, ax=ax, label="Mean Z-score", shrink=0.8)
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"[Characterization] Saved heatmap: {output_path}")
        plt.close(fig)

    return fig


def generate_gene_heatmap(
    result: CharacterizationResult,
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (14, 10),
    top_n: int = 15,
) -> Optional[Any]:
    """
    Generate a subtypes x top-genes heatmap of effect sizes.

    Args:
        result: CharacterizationResult from characterize_subtypes()
        output_path: Optional file path to save the figure
        figsize: Figure size (width, height) in inches
        top_n: Max number of genes to show per subtype

    Returns:
        matplotlib Figure object, or None if matplotlib is unavailable
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning(
            "[Characterization] matplotlib not installed — " "skipping gene heatmap generation"
        )
        return None

    # Collect all unique top genes across subtypes
    all_genes = []
    for profile in result.subtype_profiles:
        for g in profile.top_genes[:top_n]:
            if g.gene not in all_genes:
                all_genes.append(g.gene)

    if not all_genes:
        logger.warning("[Characterization] No gene contributions available for heatmap")
        return None

    row_labels = [f"{p.subtype_label} (n={p.n_samples})" for p in result.subtype_profiles]

    # Build effect size matrix
    matrix = np.zeros((len(result.subtype_profiles), len(all_genes)))
    for i, profile in enumerate(result.subtype_profiles):
        gene_map = {g.gene: g.effect_size for g in profile.top_genes}
        for j, gene in enumerate(all_genes):
            matrix[i, j] = gene_map.get(gene, 0.0)

    fig, ax = plt.subplots(figsize=figsize)

    vmax = max(abs(matrix.min()), abs(matrix.max()), 0.01)
    im = ax.imshow(
        matrix,
        cmap="RdBu_r",
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
    )

    ax.set_xticks(range(len(all_genes)))
    ax.set_xticklabels(all_genes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)

    # Annotate cells
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            color = "white" if abs(val) > vmax * 0.6 else "black"
            ax.text(
                j,
                i,
                f"{val:.1f}",
                ha="center",
                va="center",
                color=color,
                fontsize=7,
            )

    ax.set_title("Gene Contributions by Subtype (Cohen's d)", fontsize=13)
    fig.colorbar(im, ax=ax, label="Effect Size (Cohen's d)", shrink=0.8)
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"[Characterization] Saved gene heatmap: {output_path}")
        plt.close(fig)

    return fig


# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================


def export_characterization(
    result: CharacterizationResult,
    output_dir: str,
    formats: Optional[List[str]] = None,
) -> List[str]:
    """
    Export characterization results to CSV and/or Excel.

    Generates the following files:
    - subtype_summary.{csv,xlsx} — one row per subtype
    - pathway_enrichment.{csv,xlsx} — all enrichment results
    - gene_contributions.{csv,xlsx} — top gene contributions
    - pathway_scores_matrix.{csv,xlsx} — subtypes x pathways mean scores

    Args:
        result: CharacterizationResult from characterize_subtypes()
        output_dir: Directory to write output files
        formats: List of formats to export ("csv", "excel"). Defaults to ["csv"].

    Returns:
        List of file paths written
    """
    if formats is None:
        formats = ["csv"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    written_files: List[str] = []

    # 1. Subtype summary
    summary_rows = []
    for p in result.subtype_profiles:
        n_sig = sum(1 for e in p.enriched_pathways if e.significant)
        top_pw = p.enriched_pathways[0].pathway if p.enriched_pathways else "N/A"
        summary_rows.append(
            {
                "subtype_id": p.subtype_id,
                "subtype_label": p.subtype_label,
                "n_samples": p.n_samples,
                "fraction": round(p.fraction, 4),
                "mean_confidence": round(p.mean_confidence, 4),
                "n_significant_pathways": n_sig,
                "top_pathway": top_pw,
            }
        )
    summary_df = pd.DataFrame(summary_rows)

    # 2. Pathway enrichment
    enrichment_rows = []
    for p in result.subtype_profiles:
        for e in p.enriched_pathways:
            enrichment_rows.append(
                {
                    "subtype_id": p.subtype_id,
                    "subtype_label": p.subtype_label,
                    **e.to_dict(),
                }
            )
    enrichment_df = pd.DataFrame(enrichment_rows)

    # 3. Gene contributions
    gene_rows = []
    for p in result.subtype_profiles:
        for g in p.top_genes:
            gene_rows.append(
                {
                    "subtype_id": p.subtype_id,
                    "subtype_label": p.subtype_label,
                    **g.to_dict(),
                }
            )
    gene_df = pd.DataFrame(gene_rows)

    # 4. Pathway scores matrix
    matrix_data = {}
    for p in result.subtype_profiles:
        matrix_data[p.subtype_label] = {pw: round(v, 4) for pw, v in p.pathway_score_means.items()}
    matrix_df = pd.DataFrame(matrix_data).T
    matrix_df.index.name = "subtype"

    # Export CSV
    if "csv" in formats:
        for name, df in [
            ("subtype_summary", summary_df),
            ("pathway_enrichment", enrichment_df),
            ("gene_contributions", gene_df),
            ("pathway_scores_matrix", matrix_df),
        ]:
            path = output_path / f"{name}.csv"
            df.to_csv(path, index=(name == "pathway_scores_matrix"))
            written_files.append(str(path))

        logger.info(f"[Characterization] Exported CSV files to {output_dir}")

    # Export Excel
    if "excel" in formats:
        xlsx_path = output_path / "characterization_results.xlsx"
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
                summary_df.to_excel(writer, sheet_name="Summary", index=False)
                enrichment_df.to_excel(writer, sheet_name="Pathway Enrichment", index=False)
                if not gene_df.empty:
                    gene_df.to_excel(writer, sheet_name="Gene Contributions", index=False)
                matrix_df.to_excel(writer, sheet_name="Scores Matrix", index=True)
            written_files.append(str(xlsx_path))
            logger.info(f"[Characterization] Exported Excel: {xlsx_path}")
        except ImportError:
            logger.warning(
                "[Characterization] openpyxl not installed — "
                "skipping Excel export. Install with: pip install openpyxl"
            )

    return written_files
