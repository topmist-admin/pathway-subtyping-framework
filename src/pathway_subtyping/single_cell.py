"""
Single-Cell Pathway Scoring

Computes pathway-level scores from single-cell RNA-seq data using
pseudobulk aggregation or per-cell scoring. Produces the same
pathway_scores DataFrame format as bulk expression scoring, enabling
unified downstream clustering and validation.

Methods:
    - mean_z: Per-cell Z-score + mean per pathway (fast, for visualization)
    - pseudobulk_mean_z: Aggregate per cell-type, then mean-Z
    - pseudobulk_ssgsea: Aggregate per cell-type, then ssGSEA (recommended)
    - pseudobulk_gsva: Aggregate per cell-type, then GSVA

Requires: pip install pathway-subtyping[sc]
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class SingleCellScoringMethod(Enum):
    """Pathway scoring methods for single-cell data."""

    MEAN_Z = "mean_z"
    PSEUDOBULK_MEAN_Z = "pseudobulk_mean_z"
    PSEUDOBULK_SSGSEA = "pseudobulk_ssgsea"
    PSEUDOBULK_GSVA = "pseudobulk_gsva"


class SingleCellInputType(Enum):
    """Type of single-cell input data."""

    RAW_COUNTS = "raw_counts"
    LOG_NORMALIZED = "log_normalized"
    H5AD = "h5ad"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class SingleCellQualityReport:
    """Report on single-cell data quality after loading."""

    n_cells: int = 0
    n_genes: int = 0
    n_genes_after_filter: int = 0
    n_cell_types: int = 0
    cell_type_counts: Dict[str, int] = field(default_factory=dict)
    sparsity: float = 0.0
    median_genes_per_cell: float = 0.0
    median_umi_per_cell: float = 0.0
    n_pathways_covered: int = 0
    n_pathways_total: int = 0
    mean_pathway_gene_coverage: float = 0.0
    excluded_cell_types: List[str] = field(default_factory=list)
    input_type: str = "unknown"
    warnings: List[str] = field(default_factory=list)
    is_usable: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "n_genes_after_filter": self.n_genes_after_filter,
            "n_cell_types": self.n_cell_types,
            "cell_type_counts": dict(self.cell_type_counts),
            "sparsity": round(float(self.sparsity), 4),
            "median_genes_per_cell": round(float(self.median_genes_per_cell), 1),
            "median_umi_per_cell": round(float(self.median_umi_per_cell), 1),
            "n_pathways_covered": self.n_pathways_covered,
            "n_pathways_total": self.n_pathways_total,
            "mean_pathway_gene_coverage": round(float(self.mean_pathway_gene_coverage), 4),
            "excluded_cell_types": list(self.excluded_cell_types),
            "input_type": self.input_type,
            "warnings": list(self.warnings),
            "is_usable": self.is_usable,
        }


@dataclass
class SingleCellScoringResult:
    """Result of single-cell pathway scoring."""

    pathway_scores: pd.DataFrame
    per_cell_scores: Optional[pd.DataFrame]
    pseudobulk_expression: pd.DataFrame
    method: SingleCellScoringMethod
    quality_report: SingleCellQualityReport
    n_pathways_scored: int = 0
    n_pathways_skipped: int = 0
    skipped_pathways: List[str] = field(default_factory=list)
    cell_type_column: str = "cell_type"
    n_cells_per_type: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "method": self.method.value,
            "n_cell_types": len(self.pathway_scores),
            "n_pathways_scored": self.n_pathways_scored,
            "n_pathways_skipped": self.n_pathways_skipped,
            "skipped_pathways": list(self.skipped_pathways),
            "pathway_names": list(self.pathway_scores.columns),
            "cell_types": list(self.pathway_scores.index),
            "n_cells_per_type": dict(self.n_cells_per_type),
            "has_per_cell_scores": self.per_cell_scores is not None,
            "quality_report": self.quality_report.to_dict(),
        }

    def format_report(self) -> str:
        """Generate a human-readable report."""
        lines = [
            "## Single-Cell Pathway Scoring Report",
            "",
            f"**Method:** {self.method.value}",
            f"**Cells:** {self.quality_report.n_cells}",
            f"**Genes:** {self.quality_report.n_genes}",
            f"**Cell types:** {self.quality_report.n_cell_types}",
            f"**Sparsity:** {self.quality_report.sparsity:.1%}",
            f"**Pathways scored:** {self.n_pathways_scored}",
            f"**Pathways skipped:** {self.n_pathways_skipped}",
            "",
            "### Cell Type Summary",
            "",
            "| Cell Type | Cells |",
            "|-----------|-------|",
        ]
        for ct, n in sorted(self.n_cells_per_type.items()):
            lines.append(f"| {ct} | {n} |")
        if self.quality_report.excluded_cell_types:
            lines.append("")
            lines.append(
                f"**Excluded** (too few cells): "
                f"{', '.join(self.quality_report.excluded_cell_types)}"
            )
        if self.skipped_pathways:
            lines.append("")
            lines.append("### Skipped Pathways")
            lines.append("")
            for pw in self.skipped_pathways[:10]:
                lines.append(f"- {pw}")
            if len(self.skipped_pathways) > 10:
                lines.append(f"- ... and {len(self.skipped_pathways) - 10} more")
        lines.append("")
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return method-specific citations."""
        citations = [
            "Satija R, et al. Spatial reconstruction of single-cell "
            "gene expression data. Nat Biotechnol. 2015;33(5):495-502."
        ]
        method_val = self.method.value
        if "ssgsea" in method_val:
            citations.append(
                "Barbie DA, et al. Systematic RNA interference reveals "
                "that oncogenic KRAS-driven cancers require TBK1. "
                "Nature. 2009;462(7269):108-112."
            )
        elif "gsva" in method_val:
            citations.append(
                "Hanzelmann S, Castelo R, Guinney J. GSVA: gene set "
                "variation analysis for microarray and RNA-seq data. "
                "BMC Bioinformatics. 2013;14:7."
            )
        elif "mean_z" in method_val:
            citations.append(
                "Lee E, et al. Inferring pathway activity toward precise "
                "disease classification. PLoS Comput Biol. "
                "2008;4(11):e1000217."
            )
        return citations


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _ensure_anndata():
    """Lazy import anndata with helpful error message."""
    try:
        import anndata

        return anndata
    except ImportError:
        raise ImportError(
            "anndata is required for single-cell analysis. "
            "Install with: pip install pathway-subtyping[sc]"
        )


def _normalize_counts(
    X,
    scale_factor: float = 1e4,
) -> np.ndarray:
    """
    Log-normalize raw counts: log1p(x / total_counts * scale_factor).

    Handles both dense ndarray and scipy sparse matrices.
    Processes in row-chunks to limit memory for large sparse matrices.

    Args:
        X: Expression matrix (n_cells x n_genes), raw counts.
        scale_factor: Target sum for normalization (default 10,000).

    Returns:
        Dense ndarray of log-normalized values.
    """
    import scipy.sparse

    if scipy.sparse.issparse(X):
        n_cells = X.shape[0]
        chunk_size = 5000
        result = np.empty(X.shape, dtype=np.float64)
        for start in range(0, n_cells, chunk_size):
            end = min(start + chunk_size, n_cells)
            block = X[start:end].toarray().astype(np.float64)
            totals = block.sum(axis=1, keepdims=True)
            totals[totals == 0] = 1.0
            result[start:end] = np.log1p(block / totals * scale_factor)
        return result
    else:
        X = np.asarray(X, dtype=np.float64)
        totals = X.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1.0
        return np.log1p(X / totals * scale_factor)


def _compute_pseudobulk(
    adata,
    cell_type_column: str,
    min_cells_per_type: int = 10,
    normalized_X: Optional[np.ndarray] = None,
) -> Tuple[pd.DataFrame, Dict[str, int], List[str]]:
    """
    Aggregate expression per cell type (mean).

    Args:
        adata: AnnData object.
        cell_type_column: Column in adata.obs with cell type labels.
        min_cells_per_type: Exclude cell types with fewer cells.
        normalized_X: Pre-normalized expression matrix (if None, uses adata.X).

    Returns:
        Tuple of (pseudobulk DataFrame [cell_types x genes],
                  cell counts dict, excluded cell type names).
    """
    import scipy.sparse

    X = normalized_X if normalized_X is not None else adata.X
    cell_types = adata.obs[cell_type_column]
    unique_types = sorted(cell_types.unique())

    pseudobulk = {}
    cell_counts = {}
    excluded = []

    for ct in unique_types:
        mask = (cell_types == ct).values
        n_cells = int(mask.sum())
        cell_counts[str(ct)] = n_cells

        if n_cells < min_cells_per_type:
            excluded.append(str(ct))
            logger.info(
                f"[SingleCell] Excluding cell type '{ct}': "
                f"{n_cells} cells < min {min_cells_per_type}"
            )
            continue

        if scipy.sparse.issparse(X):
            ct_mean = np.asarray(X[mask].mean(axis=0)).flatten()
        else:
            ct_mean = np.asarray(X[mask].mean(axis=0)).flatten()

        pseudobulk[str(ct)] = ct_mean

    if not pseudobulk:
        logger.warning("[SingleCell] No cell types passed min_cells filter")
        gene_names = list(adata.var_names)
        return pd.DataFrame(columns=gene_names), cell_counts, excluded

    gene_names = list(adata.var_names)
    pseudobulk_df = pd.DataFrame(pseudobulk, index=gene_names).T
    # Shape: (n_cell_types, n_genes)

    logger.info(
        f"[SingleCell] Pseudobulk: {len(pseudobulk_df)} cell types, "
        f"{len(gene_names)} genes (excluded {len(excluded)} types)"
    )

    return pseudobulk_df, cell_counts, excluded


def _score_cells_mean_z(
    X,
    gene_names: List[str],
    pathways: Dict[str, List[str]],
    min_genes: int = 2,
    chunk_size: int = 5000,
    cell_ids: Optional[List[str]] = None,
    show_progress: bool = True,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Per-cell mean-Z pathway scoring with memory-efficient chunking.

    1. Compute gene means/stds across all cells.
    2. For each chunk of cells: Z-normalize, compute mean per pathway.

    Args:
        X: Expression matrix (n_cells x n_genes), dense or sparse.
        gene_names: Gene names matching columns of X.
        pathways: Pathway name -> gene list.
        min_genes: Skip pathways with fewer overlapping genes.
        chunk_size: Number of cells per processing chunk.
        cell_ids: Cell identifiers for DataFrame index.
        show_progress: Show tqdm progress bar.

    Returns:
        Tuple of (per-cell scores DataFrame [cells x pathways],
                  list of skipped pathway names).
    """
    import scipy.sparse

    n_cells, n_genes = X.shape

    # Compute gene-level stats across all cells (sparse-aware)
    if scipy.sparse.issparse(X):
        gene_means = np.asarray(X.mean(axis=0)).flatten()
        X_sq = X.copy()
        X_sq.data **= 2
        gene_vars = np.asarray(X_sq.mean(axis=0)).flatten() - gene_means**2
        gene_stds = np.sqrt(np.maximum(gene_vars, 0))
    else:
        X_arr = np.asarray(X)
        gene_means = X_arr.mean(axis=0)
        gene_stds = X_arr.std(axis=0)

    gene_stds[gene_stds == 0] = 1e-10

    # Build pathway index maps
    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}
    pathway_indices = {}
    skipped = []

    for pw_name, pw_genes in pathways.items():
        idxs = [gene_name_to_idx[g] for g in pw_genes if g in gene_name_to_idx]
        if len(idxs) < min_genes:
            skipped.append(pw_name)
            continue
        pathway_indices[pw_name] = np.array(idxs)

    if not pathway_indices:
        if cell_ids is None:
            cell_ids = [f"cell_{i}" for i in range(n_cells)]
        return pd.DataFrame(index=cell_ids), skipped

    pw_names = list(pathway_indices.keys())
    all_scores = np.empty((n_cells, len(pw_names)), dtype=np.float64)

    n_chunks = (n_cells + chunk_size - 1) // chunk_size
    for chunk_idx in tqdm(
        range(n_chunks),
        desc="Per-cell scoring",
        disable=not show_progress,
    ):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)
        if scipy.sparse.issparse(X):
            chunk_data = X[start:end].toarray().astype(np.float64)
        else:
            chunk_data = np.asarray(X[start:end], dtype=np.float64)

        z_chunk = (chunk_data - gene_means) / gene_stds

        for j, pw_name in enumerate(pw_names):
            idxs = pathway_indices[pw_name]
            all_scores[start:end, j] = z_chunk[:, idxs].mean(axis=1)

    if cell_ids is None:
        cell_ids = [f"cell_{i}" for i in range(n_cells)]

    scores_df = pd.DataFrame(all_scores, index=cell_ids, columns=pw_names)
    return scores_df, skipped


# ---------------------------------------------------------------------------
# Public: Load
# ---------------------------------------------------------------------------


def load_single_cell_data(
    path: str,
    cell_type_column: str = "cell_type",
    input_type: SingleCellInputType = SingleCellInputType.H5AD,
    min_cells_pct: float = 0.01,
    min_cells_per_type: int = 10,
    gene_column: Optional[str] = None,
    scale_factor: float = 1e4,
) -> Tuple:
    """
    Load and validate single-cell data from h5ad or CSV/TSV.

    For raw counts, applies log1p normalization. For h5ad files,
    reads directly with anndata. For CSV/TSV, wraps in AnnData.

    Args:
        path: Path to h5ad, CSV, or TSV file.
        cell_type_column: Column in obs with cell type annotations.
        input_type: Type of input data.
        min_cells_pct: Min fraction of cells a gene must be detected in.
        min_cells_per_type: Min cells per cell type.
        gene_column: Column name if genes are rows (CSV only).
        scale_factor: Scale factor for count normalization.

    Returns:
        Tuple of (AnnData with normalized X, SingleCellQualityReport).

    Raises:
        ImportError: If anndata is not installed.
        FileNotFoundError: If file does not exist.
        ValueError: If data is empty or cell_type_column is missing.
    """
    ad = _ensure_anndata()
    import scipy.sparse

    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"Single-cell data file not found: {path}")

    report = SingleCellQualityReport(input_type=input_type.value)

    # Load data
    if input_type == SingleCellInputType.H5AD or file_path.suffix == ".h5ad":
        adata = ad.read_h5ad(file_path)
        logger.info(f"[SingleCell] Loaded h5ad: {adata.n_obs} cells x {adata.n_vars} genes")
    else:
        # CSV/TSV
        sep = "\t" if file_path.suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(file_path, sep=sep, index_col=0)

        if df.empty:
            raise ValueError(f"Single-cell data file is empty: {path}")

        # Auto-detect orientation
        from .expression import _looks_like_gene_symbols

        col_genes = _looks_like_gene_symbols(list(df.columns[:100]))
        idx_genes = _looks_like_gene_symbols(list(df.index[:100]))

        if idx_genes and not col_genes:
            df = df.T
            logger.info("[SingleCell] Detected genes as rows, transposed")

        adata = ad.AnnData(
            X=df.values.astype(np.float64),
            obs=pd.DataFrame(index=df.index),
            var=pd.DataFrame(index=df.columns),
        )
        logger.info(f"[SingleCell] Loaded CSV/TSV: {adata.n_obs} cells x {adata.n_vars} genes")

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"Single-cell data is empty: {adata.n_obs} cells x {adata.n_vars} genes")

    # Validate cell type column
    if cell_type_column not in adata.obs.columns:
        raise ValueError(
            f"Cell type column '{cell_type_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    report.n_cells = adata.n_obs
    report.n_genes = adata.n_vars

    # Compute sparsity
    X = adata.X
    if scipy.sparse.issparse(X):
        report.sparsity = 1.0 - (X.nnz / (X.shape[0] * X.shape[1]))
        genes_per_cell = np.asarray((X > 0).sum(axis=1)).flatten()
        umi_per_cell = np.asarray(X.sum(axis=1)).flatten()
    else:
        X_arr = np.asarray(X)
        report.sparsity = 1.0 - (np.count_nonzero(X_arr) / X_arr.size) if X_arr.size > 0 else 0.0
        genes_per_cell = (X_arr > 0).sum(axis=1)
        umi_per_cell = X_arr.sum(axis=1)

    report.median_genes_per_cell = float(np.median(genes_per_cell))
    report.median_umi_per_cell = float(np.median(umi_per_cell))

    # Gene filtering: keep genes detected in >= min_cells_pct of cells
    n_min_cells = max(1, int(adata.n_obs * min_cells_pct))
    if scipy.sparse.issparse(X):
        gene_detection = np.asarray((X > 0).sum(axis=0)).flatten()
    else:
        gene_detection = (np.asarray(X) > 0).sum(axis=0)

    keep_mask = gene_detection >= n_min_cells
    n_kept = int(keep_mask.sum())

    if n_kept < adata.n_vars:
        adata = adata[:, keep_mask].copy()
        logger.info(
            f"[SingleCell] Gene filter: kept {n_kept}/{report.n_genes} genes "
            f"(detected in >= {min_cells_pct:.0%} of cells)"
        )

    report.n_genes_after_filter = adata.n_vars

    # Normalize if raw counts
    if input_type == SingleCellInputType.RAW_COUNTS:
        normalized = _normalize_counts(adata.X, scale_factor=scale_factor)
        adata.X = normalized
        logger.info("[SingleCell] Applied log1p normalization to raw counts")
    elif input_type == SingleCellInputType.H5AD:
        # Check if values look like counts (integers, max > 20)
        if scipy.sparse.issparse(adata.X):
            sample = adata.X[:100].toarray()
        else:
            sample = np.asarray(adata.X[:100])
        max_val = sample.max()
        is_integer = np.allclose(sample, np.round(sample), atol=1e-6)
        if is_integer and max_val > 20:
            normalized = _normalize_counts(adata.X, scale_factor=scale_factor)
            adata.X = normalized
            logger.info(
                f"[SingleCell] Detected raw counts (max={max_val:.0f}), "
                f"applied log1p normalization"
            )
        else:
            logger.info(f"[SingleCell] Values appear pre-normalized (max={max_val:.2f})")

    # Cell type counts
    cell_types = adata.obs[cell_type_column]
    for ct in sorted(cell_types.unique()):
        report.cell_type_counts[str(ct)] = int((cell_types == ct).sum())
    report.n_cell_types = len(report.cell_type_counts)

    # Usability checks
    if report.n_genes_after_filter < 50:
        report.is_usable = False
        report.warnings.append(
            f"Only {report.n_genes_after_filter} genes after filtering " f"(recommend >= 50)"
        )
    if report.n_cells < 10:
        report.is_usable = False
        report.warnings.append(f"Only {report.n_cells} cells. Need at least 10.")

    n_usable_types = sum(1 for n in report.cell_type_counts.values() if n >= min_cells_per_type)
    if n_usable_types < 2:
        report.warnings.append(
            f"Only {n_usable_types} cell types with >= {min_cells_per_type} cells. "
            f"Need at least 2 for clustering."
        )

    logger.info(
        f"[SingleCell] Loaded: {adata.n_obs} cells, {adata.n_vars} genes, "
        f"{report.n_cell_types} cell types, "
        f"sparsity={report.sparsity:.1%}"
    )

    return adata, report


# ---------------------------------------------------------------------------
# Public: Score
# ---------------------------------------------------------------------------


def score_single_cell_pathways(
    adata,
    pathways: Dict[str, List[str]],
    cell_type_column: str = "cell_type",
    method: SingleCellScoringMethod = SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
    min_genes_per_pathway: int = 2,
    min_cells_per_type: int = 10,
    compute_per_cell: bool = False,
    alpha: float = 0.25,
    chunk_size: int = 5000,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> SingleCellScoringResult:
    """
    Compute pathway scores from single-cell RNA-seq data.

    This is the main entry point. For pseudobulk methods, it aggregates
    expression per cell type, then scores pathways using the same algorithms
    as bulk expression scoring. For MEAN_Z, it scores per-cell directly.

    The primary output (pathway_scores) is a DataFrame with:
        index = cell type names
        columns = pathway names
    Compatible with run_clustering(), ValidationGates.run_all(),
    and characterize_subtypes().

    Args:
        adata: AnnData object (from load_single_cell_data or user-created).
        pathways: Dict mapping pathway names to gene lists.
        cell_type_column: Column in adata.obs with cell type labels.
        method: Scoring method to use.
        min_genes_per_pathway: Skip pathways with fewer genes in data.
        min_cells_per_type: Exclude cell types with fewer cells.
        compute_per_cell: Also compute per-cell scores (mean_z only).
        alpha: Weight parameter for ssGSEA.
        chunk_size: Cells per chunk for per-cell scoring.
        seed: Random seed for reproducibility.
        show_progress: Show progress bars.

    Returns:
        SingleCellScoringResult with Z-normalized pathway scores.

    Raises:
        ImportError: If anndata is not installed.
        ValueError: If cell_type_column not found or data is empty.
    """
    _ensure_anndata()
    import scipy.sparse

    if seed is not None:
        np.random.seed(seed)

    if cell_type_column not in adata.obs.columns:
        raise ValueError(
            f"Cell type column '{cell_type_column}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    logger.info(
        f"[SingleCell] Scoring {len(pathways)} pathways via {method.value} "
        f"({adata.n_obs} cells, {adata.n_vars} genes)"
    )

    X = adata.X
    gene_names = list(adata.var_names)
    cell_types_series = adata.obs[cell_type_column]

    # Compute pseudobulk
    pseudobulk_df, cell_counts, excluded_types = _compute_pseudobulk(
        adata, cell_type_column, min_cells_per_type
    )

    per_cell_scores = None
    skipped: List[str] = []

    if method == SingleCellScoringMethod.MEAN_Z:
        # Per-cell mean_z, then aggregate per cell-type for primary output
        per_cell_scores, skipped = _score_cells_mean_z(
            X,
            gene_names,
            pathways,
            min_genes=min_genes_per_pathway,
            chunk_size=chunk_size,
            cell_ids=list(adata.obs_names),
            show_progress=show_progress,
        )
        # Aggregate: mean per cell type
        if per_cell_scores.empty or pseudobulk_df.empty:
            pathway_scores = pd.DataFrame()
        else:
            temp = per_cell_scores.copy()
            temp["_ct"] = cell_types_series.values
            temp = temp[~temp["_ct"].isin(excluded_types)]
            pathway_scores = temp.groupby("_ct").mean()
            pathway_scores.index.name = None
    else:
        # Pseudobulk methods: reuse expression.py internals
        from .expression import _score_gsva, _score_mean_z, _score_ssgsea

        if pseudobulk_df.empty:
            pathway_scores = pd.DataFrame()
            skipped = list(pathways.keys())
        elif method == SingleCellScoringMethod.PSEUDOBULK_MEAN_Z:
            pathway_scores, skipped = _score_mean_z(
                pseudobulk_df, pathways, min_genes=min_genes_per_pathway
            )
        elif method == SingleCellScoringMethod.PSEUDOBULK_SSGSEA:
            pathway_scores, skipped = _score_ssgsea(
                pseudobulk_df,
                pathways,
                alpha=alpha,
                min_genes=min_genes_per_pathway,
                show_progress=show_progress,
            )
        elif method == SingleCellScoringMethod.PSEUDOBULK_GSVA:
            pathway_scores, skipped = _score_gsva(
                pseudobulk_df,
                pathways,
                min_genes=min_genes_per_pathway,
                show_progress=show_progress,
            )
        else:
            raise ValueError(f"Unknown scoring method: {method}")

        # Optionally compute per-cell scores
        if compute_per_cell:
            per_cell_scores, _ = _score_cells_mean_z(
                X,
                gene_names,
                pathways,
                min_genes=min_genes_per_pathway,
                chunk_size=chunk_size,
                cell_ids=list(adata.obs_names),
                show_progress=show_progress,
            )

    # Post-processing: remove zero-variance, Z-normalize
    n_scored = len(pathway_scores.columns) if not pathway_scores.empty else 0
    n_skipped = len(skipped)

    if n_scored > 0:
        pw_stds = pathway_scores.std()
        zero_var = pw_stds[pw_stds == 0].index.tolist()
        if zero_var:
            pathway_scores = pathway_scores.drop(columns=zero_var)
            skipped.extend(zero_var)
            n_skipped += len(zero_var)
            n_scored -= len(zero_var)
            logger.info(f"[SingleCell] Removed {len(zero_var)} zero-variance pathway(s)")

        if n_scored > 0:
            means = pathway_scores.mean()
            stds = pathway_scores.std().replace(0, 1e-10)
            pathway_scores = (pathway_scores - means) / stds

    # Build quality report
    available_genes = set(gene_names)
    coverages = []
    n_covered = 0
    for pw_genes in pathways.values():
        overlap = len(set(pw_genes) & available_genes)
        if len(pw_genes) > 0:
            coverages.append(overlap / len(pw_genes))
            if overlap >= min_genes_per_pathway:
                n_covered += 1

    if scipy.sparse.issparse(X):
        sparsity = 1.0 - (X.nnz / (X.shape[0] * X.shape[1]))
        genes_per_cell = np.asarray((X > 0).sum(axis=1)).flatten()
        umi_per_cell = np.asarray(X.sum(axis=1)).flatten()
    else:
        X_arr = np.asarray(X)
        sparsity = 1.0 - (np.count_nonzero(X_arr) / X_arr.size) if X_arr.size > 0 else 0.0
        genes_per_cell = (X_arr > 0).sum(axis=1)
        umi_per_cell = X_arr.sum(axis=1)

    n_cells_per_type = {ct: cell_counts[ct] for ct in pseudobulk_df.index if ct in cell_counts}

    report = SingleCellQualityReport(
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        n_genes_after_filter=adata.n_vars,
        n_cell_types=len(pseudobulk_df),
        cell_type_counts=cell_counts,
        sparsity=sparsity,
        median_genes_per_cell=float(np.median(genes_per_cell)),
        median_umi_per_cell=float(np.median(umi_per_cell)),
        n_pathways_covered=n_covered,
        n_pathways_total=len(pathways),
        mean_pathway_gene_coverage=(float(np.mean(coverages)) if coverages else 0.0),
        excluded_cell_types=excluded_types,
        input_type="single_cell",
        is_usable=n_scored > 0 and len(pseudobulk_df) >= 2,
    )

    if n_scored == 0:
        report.warnings.append(
            "No pathways scored. Check that pathway gene names " "match AnnData var_names."
        )
    if len(pseudobulk_df) < 2:
        report.warnings.append(
            f"Only {len(pseudobulk_df)} cell types after filtering. "
            f"Need at least 2 for clustering."
        )

    logger.info(
        f"[SingleCell] Scored {n_scored} pathways, skipped {n_skipped} "
        f"({len(pseudobulk_df)} cell types, "
        f"gene coverage: {report.mean_pathway_gene_coverage:.1%})"
    )

    return SingleCellScoringResult(
        pathway_scores=pathway_scores,
        per_cell_scores=per_cell_scores,
        pseudobulk_expression=pseudobulk_df,
        method=method,
        quality_report=report,
        n_pathways_scored=n_scored,
        n_pathways_skipped=n_skipped,
        skipped_pathways=skipped,
        cell_type_column=cell_type_column,
        n_cells_per_type=n_cells_per_type,
    )
