"""
Performance optimization utilities for large cohort processing.

Provides memory-efficient and parallelized implementations for
processing large datasets (10,000+ samples).
"""

import gzip
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Generator, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..data_quality import parse_genotype

logger = logging.getLogger(__name__)


@dataclass
class ProcessingStats:
    """Statistics from chunked processing."""

    total_samples: int
    total_variants: int
    chunks_processed: int
    peak_memory_mb: float
    processing_time_seconds: float


def chunked_vcf_reader(
    vcf_path: str,
    chunk_size: int = 1000,
    sample_subset: Optional[List[str]] = None,
) -> Generator[Tuple[pd.DataFrame, pd.DataFrame, List[str]], None, None]:
    """
    Memory-efficient chunked VCF reader.

    Reads VCF file in chunks to avoid loading entire file into memory.
    Suitable for cohorts with 10,000+ samples.

    Features:
    - Gzip support (.vcf.gz files)
    - Multi-allelic variant expansion
    - Consistent genotype parsing using parse_genotype()

    Args:
        vcf_path: Path to VCF file (supports .vcf and .vcf.gz)
        chunk_size: Number of variants per chunk
        sample_subset: Optional list of sample IDs to include

    Yields:
        Tuple of (variants_df, genotypes_df, samples)
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    samples: List[str] = []
    sample_indices: Optional[List[int]] = None

    # Determine if gzipped
    is_gzipped = str(vcf_path).endswith(".gz")
    open_func = gzip.open if is_gzipped else open

    with open_func(vcf_path, "rt") as f:
        # Parse header
        for line in f:
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                all_samples = parts[9:]

                if sample_subset:
                    sample_indices = [i for i, s in enumerate(all_samples) if s in sample_subset]
                    samples = [all_samples[i] for i in sample_indices]
                else:
                    samples = all_samples
                break

        # Process variants in chunks
        variants_chunk: List[Dict[str, Any]] = []
        genotypes_chunk: List[Dict[str, int]] = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            # Parse variant
            chrom, pos, vid, ref, alt, qual, filt, info, fmt = parts[:9]

            info_dict: Dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    info_dict[key] = val

            # Handle multi-allelic variants
            alts = alt.split(",")
            is_multi_allelic = len(alts) > 1

            # Parse genotypes for subset
            sample_gts = parts[9:]
            if sample_indices:
                sample_gts = [sample_gts[i] for i in sample_indices]

            # Expand multi-allelic variants into separate records
            for alt_idx, alt_allele in enumerate(alts):
                alt_num = alt_idx + 1  # ALT alleles are 1-indexed

                variant_id = f"{vid}_{alt_num}" if is_multi_allelic else vid

                variants_chunk.append(
                    {
                        "chrom": chrom,
                        "pos": int(pos),
                        "id": variant_id,
                        "ref": ref,
                        "alt": alt_allele,
                        "qual": float(qual) if qual != "." else 0,
                        "filter": filt,
                        "gene": info_dict.get("GENE", ""),
                        "consequence": info_dict.get("CONSEQUENCE", ""),
                        "cadd": _safe_float(info_dict.get("CADD"), 0.0),
                    }
                )

                # Parse genotypes using unified parse_genotype function
                gt_row: Dict[str, int] = {}
                for sample, gt_data in zip(samples, sample_gts):
                    # Use parse_genotype for consistent allele-specific counting
                    allele_count, is_valid = parse_genotype(gt_data, target_allele=alt_num)
                    gt_row[sample] = allele_count if is_valid else 0
                genotypes_chunk.append(gt_row)

            # Yield chunk when full
            if len(variants_chunk) >= chunk_size:
                yield (
                    pd.DataFrame(variants_chunk),
                    pd.DataFrame(genotypes_chunk),
                    samples,
                )
                variants_chunk = []
                genotypes_chunk = []

        # Yield remaining
        if variants_chunk:
            yield (
                pd.DataFrame(variants_chunk),
                pd.DataFrame(genotypes_chunk),
                samples,
            )


def _safe_float(value: Any, default: float = 0.0) -> float:
    """Safely convert a value to float."""
    if value is None or value == "." or value == "":
        return default
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def compute_gene_burdens_chunked(
    vcf_path: str,
    chunk_size: int = 1000,
    sample_subset: Optional[List[str]] = None,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> pd.DataFrame:
    """
    Compute gene burdens with chunked processing for memory efficiency.

    Features:
    - Memory-efficient chunked processing
    - Consequence-based CADD defaults for missing values
    - Multi-allelic variant support

    Args:
        vcf_path: Path to VCF file (supports .vcf and .vcf.gz)
        chunk_size: Variants per chunk
        sample_subset: Optional sample subset
        progress_callback: Optional callback(current, total) for progress

    Returns:
        Gene burden DataFrame (samples × genes)
    """
    burden_accumulator: Dict[str, Dict[str, float]] = {}
    total_variants = 0
    chunks_processed = 0

    for variants_df, genotypes_df, samples in chunked_vcf_reader(
        vcf_path, chunk_size, sample_subset
    ):
        total_variants += len(variants_df)
        chunks_processed += 1

        # Compute burdens for this chunk
        for idx, var in variants_df.iterrows():
            gene = var["gene"]
            if not gene:
                continue

            if gene not in burden_accumulator:
                burden_accumulator[gene] = {s: 0.0 for s in samples}

            consequence = str(var.get("consequence", "")).lower()

            # Compute weight based on consequence
            weight = 0.1
            if "frameshift" in consequence or "stop" in consequence:
                weight = 1.0
            elif "missense" in consequence:
                weight = 0.5 if var["cadd"] > 25 else 0.1

            # Handle missing CADD scores with consequence-based defaults
            # This prevents silent data loss when CADD is unavailable
            cadd_score = var["cadd"]
            if cadd_score <= 0:
                # Use consequence-based defaults (same as pipeline.py)
                if "frameshift" in consequence or "stop" in consequence:
                    cadd_score = 35.0  # High impact default
                elif "missense" in consequence:
                    cadd_score = 20.0  # Moderate impact default
                else:
                    cadd_score = 10.0  # Low impact default
            # Cap and normalize CADD score
            cadd_normalized = min(cadd_score, 40.0) / 40.0

            # Add burden for each sample
            for sample in samples:
                gt = genotypes_df.loc[idx, sample]
                if gt > 0:
                    burden_accumulator[gene][sample] += gt * weight * cadd_normalized

        if progress_callback:
            progress_callback(total_variants, -1)  # -1 = unknown total

    logger.info(f"Processed {total_variants} variants in {chunks_processed} chunks")
    return pd.DataFrame(burden_accumulator).fillna(0)


def parallel_pathway_scores(
    gene_burdens: pd.DataFrame,
    pathways: Dict[str, List[str]],
    n_workers: Optional[int] = None,
) -> pd.DataFrame:
    """
    Compute pathway scores in parallel.

    Features:
    - Parallel computation for large pathway sets
    - Zero-variance pathway handling (removed before normalization)
    - Robust Z-score normalization

    Args:
        gene_burdens: Gene burden DataFrame (samples × genes)
        pathways: Dict of pathway name → gene list
        n_workers: Number of parallel workers (None = CPU count)

    Returns:
        Pathway scores DataFrame (samples × pathways), Z-score normalized

    Raises:
        ValueError: If fewer than 2 pathways remain after filtering
    """
    if n_workers is None:
        n_workers = min(os.cpu_count() or 4, len(pathways))

    def compute_single_pathway(
        pathway_data: Tuple[str, List[str]],
    ) -> Tuple[str, Optional[pd.Series]]:
        pathway_name, pathway_genes = pathway_data
        common_genes = [g for g in pathway_genes if g in gene_burdens.columns]

        if len(common_genes) < 2:
            return pathway_name, None

        return pathway_name, gene_burdens[common_genes].mean(axis=1)

    pathway_items = list(pathways.items())

    # For small numbers of pathways, just use sequential processing
    if len(pathway_items) <= 4:
        results = [compute_single_pathway(item) for item in pathway_items]
    else:
        # Use thread pool for I/O bound operations
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(compute_single_pathway, pathway_items))

    # Assemble results
    pathway_scores = {}
    for name, scores in results:
        if scores is not None:
            pathway_scores[name] = scores

    scores_df = pd.DataFrame(pathway_scores)

    if scores_df.empty:
        raise ValueError("No valid pathways found (all have <2 genes in burden data)")

    # Check for and handle zero-variance pathways before normalization
    pathway_stds = scores_df.std()
    zero_variance_pathways = pathway_stds[pathway_stds == 0].index.tolist()

    if zero_variance_pathways:
        logger.warning(
            f"Removing {len(zero_variance_pathways)} zero-variance pathway(s): "
            f"{zero_variance_pathways[:5]}{'...' if len(zero_variance_pathways) > 5 else ''}"
        )
        scores_df = scores_df.drop(columns=zero_variance_pathways)

    if scores_df.empty or len(scores_df.columns) < 2:
        raise ValueError(
            f"Insufficient pathways after filtering: {len(scores_df.columns)} remaining. "
            f"Need at least 2 pathways with non-zero variance for clustering."
        )

    # Z-score normalize with numerical stability
    means = scores_df.mean()
    stds = scores_df.std()
    stds = stds.replace(0, 1e-10)  # Safety net (should not occur after filtering)
    scores_df = (scores_df - means) / stds

    return scores_df.fillna(0)


class ProgressTracker:
    """
    Progress tracking for long-running operations.

    Example:
        tracker = ProgressTracker(total=1000, desc="Processing variants")
        for i, variant in enumerate(variants):
            process(variant)
            tracker.update(i + 1)
        tracker.finish()
    """

    def __init__(
        self,
        total: Optional[int] = None,
        desc: str = "Processing",
        log_interval: int = 100,
    ):
        self.total = total
        self.desc = desc
        self.log_interval = log_interval
        self.current = 0

    def update(self, current: int) -> None:
        """Update progress."""
        self.current = current
        if current % self.log_interval == 0:
            if self.total:
                pct = 100 * current / self.total
                logger.info(f"{self.desc}: {current}/{self.total} ({pct:.1f}%)")
            else:
                logger.info(f"{self.desc}: {current}")

    def finish(self) -> None:
        """Mark as complete."""
        logger.info(f"{self.desc}: Complete ({self.current} processed)")


def estimate_memory_usage(
    n_samples: int,
    n_variants: int,
    n_genes: int,
    n_pathways: int,
) -> Dict[str, float]:
    """
    Estimate memory requirements for a given dataset size.

    Args:
        n_samples: Number of samples
        n_variants: Number of variants
        n_genes: Number of unique genes
        n_pathways: Number of pathways

    Returns:
        Dict with memory estimates in MB
    """
    # Rough estimates based on data types
    genotypes_mb = (n_samples * n_variants * 1) / (1024 * 1024)  # int8
    variants_mb = (n_variants * 200) / (1024 * 1024)  # ~200 bytes per variant
    burdens_mb = (n_samples * n_genes * 8) / (1024 * 1024)  # float64
    pathway_scores_mb = (n_samples * n_pathways * 8) / (1024 * 1024)  # float64

    total_mb = genotypes_mb + variants_mb + burdens_mb + pathway_scores_mb

    return {
        "genotypes_mb": round(genotypes_mb, 2),
        "variants_mb": round(variants_mb, 2),
        "gene_burdens_mb": round(burdens_mb, 2),
        "pathway_scores_mb": round(pathway_scores_mb, 2),
        "total_estimated_mb": round(total_mb, 2),
        "recommended_chunk_size": min(1000, max(100, 2000 // max(1, n_samples // 1000))),
    }


def downsample_cohort(
    phenotypes: pd.DataFrame,
    target_size: int,
    stratify_by: Optional[str] = None,
    seed: int = 42,
) -> List[str]:
    """
    Downsample cohort while optionally maintaining stratification.

    Args:
        phenotypes: Phenotype DataFrame with sample_id index
        target_size: Target number of samples
        stratify_by: Optional column to stratify by
        seed: Random seed

    Returns:
        List of selected sample IDs
    """
    rng = np.random.RandomState(seed)

    if stratify_by and stratify_by in phenotypes.columns:
        # Stratified sampling
        selected = []
        groups = phenotypes.groupby(stratify_by)
        n_groups = len(groups)
        samples_per_group = target_size // n_groups

        for name, group in groups:
            n_select = min(len(group), samples_per_group)
            selected.extend(rng.choice(group.index.tolist(), size=n_select, replace=False))

        # Fill remaining quota
        remaining = target_size - len(selected)
        if remaining > 0:
            available = [s for s in phenotypes.index if s not in selected]
            if available:
                selected.extend(
                    rng.choice(available, size=min(remaining, len(available)), replace=False)
                )

        return selected
    else:
        # Simple random sampling
        return rng.choice(
            phenotypes.index.tolist(),
            size=min(target_size, len(phenotypes)),
            replace=False,
        ).tolist()
