#!/usr/bin/env python3
"""
Run pathway-subtyping framework on TCGA benchmark datasets only.

Generates bootstrap_threshold_calibration_47datasets_tcga.csv with:
- dataset_id, dataset_name, domain, source
- silhouette, bootstrap_ari_5th_percentile
- n_samples, n_detected_clusters, n_true_clusters, status

For GEO datasets, use: run_47benchmarks_geo.py

This dataset is meant for independent verification and Zenodo publication.
"""

import logging
import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys
from typing import Optional, Dict, Tuple
import traceback

# Import pathway-subtyping framework
try:
    from pathway_subtyping.clustering import run_clustering, select_n_clusters, _compute_cluster_metrics
except ImportError as e:
    print(f"ERROR: Cannot import pathway-subtyping. Is it installed? {e}")
    sys.exit(1)

# Import benchmark helpers
try:
    from benchmark_data_loaders import BenchmarkDataLoaderFactory, DataLoaderError
    from benchmark_validation import compute_silhouette, bootstrap_ari_stability
except ImportError as e:
    print(f"ERROR: Cannot import benchmark helpers. Are they in scripts/ directory? {e}")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(name)s] %(message)s',
    handlers=[
        logging.FileHandler('benchmark_run.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Paths
REPO_ROOT = Path(__file__).parent.parent
MANIFEST_PATH = REPO_ROOT / 'data' / 'benchmarks' / 'benchmark_47datasets_manifest.csv'
RESULTS_DIR = REPO_ROOT / 'research-results' / 'benchmarks'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_CSV = RESULTS_DIR / 'bootstrap_threshold_calibration_47datasets_tcga.csv'


def load_benchmark_dataset(row: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load a single benchmark dataset using the factory loader.

    Args:
        row: Pandas Series from manifest CSV

    Returns:
        (gene_expression, ground_truth_labels)
    """
    try:
        return BenchmarkDataLoaderFactory.load(row)
    except DataLoaderError as e:
        logger.error(f"Failed to load {row['dataset_name']}: {e}")
        raise


def select_top_genes(gene_expression: np.ndarray, n_genes: int = 1000) -> np.ndarray:
    """
    Select top genes by variance for faster computation.

    Args:
        gene_expression: (n_samples, n_genes) array
        n_genes: Number of top genes to select (default 1000)

    Returns:
        Filtered gene expression array
    """
    if gene_expression.shape[1] <= n_genes:
        return gene_expression

    # Compute variance per gene
    gene_vars = np.var(gene_expression, axis=0)
    top_indices = np.argsort(gene_vars)[-n_genes:]

    logger.debug(f"        Selected {n_genes} genes (variance-based) from {gene_expression.shape[1]}")
    return gene_expression[:, top_indices]


def run_single_benchmark(
    gene_expression: np.ndarray,
    true_labels: np.ndarray,
    dataset_name: str,
    seed: int = 42,
    n_components_range: Tuple[int, int] = (2, 10),
    n_genes_selection: int = 1000
) -> Dict:
    """
    Run pathway-subtyping framework on a single dataset and extract metrics.

    Args:
        gene_expression: (n_samples, n_genes) array
        true_labels: (n_samples,) array of ground truth cluster labels
        dataset_name: Name for logging
        seed: Random seed for reproducibility
        n_components_range: Range of k values to test
        n_genes_selection: Number of top genes to select (default 1000)

    Returns:
        Dictionary with keys:
        - silhouette: Silhouette coefficient (-1 to 1)
        - bootstrap_ari_5th_percentile: 5th percentile of ARI over bootstrap samples
        - n_samples, n_detected_clusters, n_true_clusters
        - framework_metrics: Dict of other metrics
    """
    logger.info(f"    Running framework on {dataset_name}...")

    try:
        # Step 0: Select top genes by variance (for computational efficiency)
        logger.info(f"      - Selecting top {n_genes_selection} genes by variance...")
        gene_expression = select_top_genes(gene_expression, n_genes=n_genes_selection)
        logger.info(f"      - Using {gene_expression.shape[1]} genes for analysis")

        # Step 1: Select optimal number of clusters using BIC
        logger.info(f"      - Selecting k using BIC...")
        k_range = list(range(n_components_range[0], n_components_range[1] + 1))
        k_selection = select_n_clusters(
            data=gene_expression,
            k_range=k_range,
            method='bic',
            seed=seed
        )
        optimal_k = k_selection.optimal_k
        logger.info(f"      - Optimal k: {optimal_k}")

        # Step 2: Run clustering with selected k
        logger.info(f"      - Running clustering with k={optimal_k}...")
        clustering_result = run_clustering(
            data=gene_expression,
            n_clusters=optimal_k,
            seed=seed
        )

        clusters = clustering_result.labels

        # Compute silhouette score
        logger.info(f"      - Computing silhouette coefficient...")
        silhouette = compute_silhouette(gene_expression, clusters)

        # Create a clustering function for bootstrap
        def cluster_func(data):
            # Re-select k on resampled data
            k_sel = select_n_clusters(
                data=data,
                k_range=k_range,
                method='bic',
                seed=seed
            )
            result = run_clustering(
                data=data,
                n_clusters=k_sel.optimal_k,
                seed=seed
            )
            return result.labels

        # Compute bootstrap stability (ARI across bootstrap samples)
        logger.info(f"      - Running bootstrap stability analysis (50 replicates)...")
        bootstrap_ari = bootstrap_ari_stability(
            gene_expression,
            true_labels,
            cluster_func,
            n_iterations=50,
            sample_fraction=0.8,
            seed=seed
        )
        bootstrap_ari_5th = np.percentile(bootstrap_ari, 5)

        # Detected vs true cluster counts
        n_detected = len(np.unique(clusters))
        n_true = len(np.unique(true_labels))
        n_samples = gene_expression.shape[0]

        logger.info(f"      ✓ Silhouette: {silhouette:.4f}, Bootstrap ARI (5th pct): {bootstrap_ari_5th:.4f}")
        logger.info(f"      ✓ Detected clusters: {n_detected} (true: {n_true}, samples: {n_samples})")

        return {
            'silhouette': float(silhouette),
            'bootstrap_ari_5th_percentile': float(bootstrap_ari_5th),
            'n_samples': int(n_samples),
            'n_detected_clusters': int(n_detected),
            'n_true_clusters': int(n_true),
            'framework_metrics': {
                'bootstrap_ari_mean': float(np.mean(bootstrap_ari)),
                'bootstrap_ari_std': float(np.std(bootstrap_ari)),
                'bootstrap_ari_95th': float(np.percentile(bootstrap_ari, 95))
            },
            'status': 'PASS'
        }

    except Exception as e:
        logger.error(f"      ✗ ERROR: {e}")
        logger.debug(traceback.format_exc())

        return {
            'status': 'ERROR',
            'error_message': str(e),
            'n_samples': gene_expression.shape[0] if hasattr(gene_expression, 'shape') else None
        }


def main():
    """Main benchmark runner."""
    import argparse

    parser = argparse.ArgumentParser(description='Run pathway-subtyping benchmark on TCGA datasets only')
    parser.add_argument('--test', action='store_true', help='Test mode: run on first N datasets only')
    parser.add_argument('--n-datasets', type=int, default=3, help='Number of datasets to test (default: 3)')
    parser.add_argument('--domain', type=str, help='Filter datasets by domain (e.g., "Oncology")')
    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("PATHWAY-SUBTYPING TCGA-BENCHMARK CALIBRATION RUN")
    if args.test:
        logger.info(f"TEST MODE: Running on first {args.n_datasets} dataset(s)")
    logger.info("=" * 80)

    # Load manifest
    logger.info(f"Loading manifest from {MANIFEST_PATH}...")
    manifest_df = pd.read_csv(MANIFEST_PATH)
    logger.info(f"Loaded {len(manifest_df)} datasets")

    # Filter to TCGA datasets only
    manifest_df = manifest_df[manifest_df['source'] == 'GDC']
    logger.info(f"Filtered to TCGA (GDC source): {len(manifest_df)} dataset(s)")

    # Apply additional filters in order: domain → test limit
    if args.domain:
        manifest_df = manifest_df[manifest_df['domain'] == args.domain]
        logger.info(f"Filtered to domain '{args.domain}': {len(manifest_df)} dataset(s)")

    if args.test:
        manifest_df = manifest_df.head(args.n_datasets)
        logger.info(f"TEST MODE: Limited to first {len(manifest_df)} dataset(s)")

    results = []
    success_count = 0
    error_count = 0
    total_datasets = len(manifest_df)

    # Iterate through datasets
    for dataset_num, (idx, row) in enumerate(manifest_df.iterrows(), 1):
        dataset_id = row['dataset_id']
        dataset_name = row['dataset_name']
        domain = row['domain']
        source = row['source']

        logger.info(f"\n[{dataset_num:2d}/{total_datasets}] {dataset_name} ({domain}, {source})")
        logger.info("-" * 70)

        try:
            # Load dataset
            expression, labels = load_benchmark_dataset(row)

            logger.info(f"  Loaded: {expression.shape[0]} samples × {expression.shape[1]} genes")

            # Run framework
            framework_result = run_single_benchmark(
                gene_expression=expression,
                true_labels=labels,
                dataset_name=dataset_name,
                seed=42
            )

            # Compile result
            result = {
                'dataset_id': dataset_id,
                'dataset_name': dataset_name,
                'domain': domain,
                'source': source,
                'silhouette': framework_result.get('silhouette'),
                'bootstrap_ari_5th_percentile': framework_result.get('bootstrap_ari_5th_percentile'),
                'n_samples': framework_result.get('n_samples'),
                'n_detected_clusters': framework_result.get('n_detected_clusters'),
                'n_true_clusters': framework_result.get('n_true_clusters'),
                'status': framework_result.get('status', 'UNKNOWN'),
                'error_message': framework_result.get('error_message', '')
            }

            results.append(result)

            if framework_result.get('status') == 'PASS':
                success_count += 1
            else:
                error_count += 1

        except Exception as e:
            logger.error(f"  ✗ FATAL ERROR: {e}")
            logger.debug(traceback.format_exc())

            result = {
                'dataset_id': dataset_id,
                'dataset_name': dataset_name,
                'domain': domain,
                'source': source,
                'status': 'ERROR',
                'error_message': str(e)
            }
            results.append(result)
            error_count += 1

    # Save results
    logger.info("\n" + "=" * 80)
    logger.info("BENCHMARK RUN SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total datasets: {len(results)}")
    logger.info(f"Passed: {success_count}")
    logger.info(f"Failed: {error_count}")
    if len(results) > 0:
        logger.info(f"Success rate: {100*success_count/len(results):.1f}%")
    else:
        logger.warning("No datasets were processed")

    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_CSV, index=False)
    logger.info(f"\n✓ Results saved to: {RESULTS_CSV}")

    # Print summary statistics
    if success_count > 0:
        passed_df = results_df[results_df['status'] == 'PASS']
        logger.info("\nSilhouette Statistics (passed datasets only):")
        logger.info(f"  Mean: {passed_df['silhouette'].mean():.4f}")
        logger.info(f"  Median: {passed_df['silhouette'].median():.4f}")
        logger.info(f"  Std: {passed_df['silhouette'].std():.4f}")
        logger.info(f"  Range: [{passed_df['silhouette'].min():.4f}, {passed_df['silhouette'].max():.4f}]")

        logger.info("\nBootstrap ARI (5th percentile) Statistics:")
        logger.info(f"  Mean: {passed_df['bootstrap_ari_5th_percentile'].mean():.4f}")
        logger.info(f"  Median: {passed_df['bootstrap_ari_5th_percentile'].median():.4f}")
        logger.info(f"  Std: {passed_df['bootstrap_ari_5th_percentile'].std():.4f}")
        logger.info(f"  Range: [{passed_df['bootstrap_ari_5th_percentile'].min():.4f}, {passed_df['bootstrap_ari_5th_percentile'].max():.4f}]")

    logger.info("\n" + "=" * 80)
    logger.info("BENCHMARK RUN COMPLETE")
    logger.info("=" * 80)


if __name__ == '__main__':
    main()
