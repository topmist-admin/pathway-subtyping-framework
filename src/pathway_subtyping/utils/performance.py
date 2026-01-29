"""
Performance optimization utilities for large cohort processing.

Provides memory-efficient and parallelized implementations for
processing large datasets (10,000+ samples).
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Generator, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd

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

    Args:
        vcf_path: Path to VCF file
        chunk_size: Number of variants per chunk
        sample_subset: Optional list of sample IDs to include

    Yields:
        Tuple of (variants_df, genotypes_df, samples)
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    samples = []
    sample_indices = None

    with open(vcf_path, "r") as f:
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
        variants_chunk = []
        genotypes_chunk = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            # Parse variant
            chrom, pos, vid, ref, alt, qual, filt, info, fmt = parts[:9]

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    info_dict[key] = val

            variants_chunk.append(
                {
                    "chrom": chrom,
                    "pos": int(pos),
                    "id": vid,
                    "ref": ref,
                    "alt": alt,
                    "qual": float(qual) if qual != "." else 0,
                    "filter": filt,
                    "gene": info_dict.get("GENE", ""),
                    "consequence": info_dict.get("CONSEQUENCE", ""),
                    "cadd": float(info_dict.get("CADD", 0)),
                }
            )

            # Parse genotypes
            sample_gts = parts[9:]
            if sample_indices:
                sample_gts = [sample_gts[i] for i in sample_indices]

            gt_row = {}
            for sample, gt_data in zip(samples, sample_gts):
                gt = gt_data.split(":")[0]
                if gt in ("0/0", "0|0"):
                    gt_row[sample] = 0
                elif gt in ("0/1", "1/0", "0|1", "1|0"):
                    gt_row[sample] = 1
                elif gt in ("1/1", "1|1"):
                    gt_row[sample] = 2
                else:
                    gt_row[sample] = 0
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


def compute_gene_burdens_chunked(
    vcf_path: str,
    chunk_size: int = 1000,
    sample_subset: Optional[List[str]] = None,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> pd.DataFrame:
    """
    Compute gene burdens with chunked processing for memory efficiency.

    Args:
        vcf_path: Path to VCF file
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

            # Compute weight
            weight = 0.1
            if "frameshift" in var["consequence"] or "stop" in var["consequence"]:
                weight = 1.0
            elif "missense" in var["consequence"]:
                weight = 0.5 if var["cadd"] > 25 else 0.1

            # Add burden for each sample
            for sample in samples:
                gt = genotypes_df.loc[idx, sample]
                if gt > 0:
                    burden_accumulator[gene][sample] += gt * weight * (var["cadd"] / 40.0)

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

    Args:
        gene_burdens: Gene burden DataFrame (samples × genes)
        pathways: Dict of pathway name → gene list
        n_workers: Number of parallel workers (None = CPU count)

    Returns:
        Pathway scores DataFrame (samples × pathways)
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

    # Z-score normalize
    scores_df = (scores_df - scores_df.mean()) / scores_df.std()
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
