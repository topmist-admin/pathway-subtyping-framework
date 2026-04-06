#!/usr/bin/env python3
"""
Merge all benchmark result CSV files into complete 47-dataset file.

Combines results from:
1. bootstrap_threshold_calibration_15datasets_zenodo.csv (12 TCGA + 3 GEO test)
2. bootstrap_threshold_calibration_geo_remaining_32.csv (32 GEO with mixed status)
3. bootstrap_threshold_calibration_geo_final_10.csv (10 GEO final run - prioritized)

Output: bootstrap_threshold_calibration_47datasets_zenodo.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Paths
REPO_ROOT = Path(__file__).parent.parent
RESULTS_DIR = REPO_ROOT / 'research-results' / 'benchmarks'

# Input files
FILE_15_ZENODO = RESULTS_DIR / 'bootstrap_threshold_calibration_15datasets_zenodo.csv'
FILE_32_REMAINING = RESULTS_DIR / 'bootstrap_threshold_calibration_geo_remaining_32.csv'
FILE_10_FINAL = RESULTS_DIR / 'bootstrap_threshold_calibration_geo_final_10.csv'

# Output file
OUTPUT_FILE = RESULTS_DIR / 'bootstrap_threshold_calibration_47datasets_zenodo.csv'

def merge_results():
    """Merge all result files into complete 47-dataset CSV."""

    print("=" * 80)
    print("MERGING COMPLETE 47-DATASET BENCHMARK RESULTS")
    print("=" * 80)

    # Load 15-dataset file (12 TCGA + 3 GEO test)
    print(f"\n1. Loading {FILE_15_ZENODO.name}...")
    if not FILE_15_ZENODO.exists():
        print(f"ERROR: File not found: {FILE_15_ZENODO}")
        return False

    df_15 = pd.read_csv(FILE_15_ZENODO)
    print(f"   Loaded {len(df_15)} datasets")
    print(f"   Dataset IDs: {sorted(df_15['dataset_id'].unique())}")

    # Load 32-dataset file
    print(f"\n2. Loading {FILE_32_REMAINING.name}...")
    if not FILE_32_REMAINING.exists():
        print(f"ERROR: File not found: {FILE_32_REMAINING}")
        return False

    df_32 = pd.read_csv(FILE_32_REMAINING)
    print(f"   Loaded {len(df_32)} rows total")
    print(f"   Dataset IDs (all): {sorted(df_32['dataset_id'].unique())}")

    # Load final 10 file
    print(f"\n3. Loading {FILE_10_FINAL.name}...")
    if not FILE_10_FINAL.exists():
        print(f"ERROR: File not found: {FILE_10_FINAL}")
        return False

    df_10_final = pd.read_csv(FILE_10_FINAL)
    print(f"   Loaded {len(df_10_final)} datasets")
    print(f"   Dataset IDs: {sorted(df_10_final['dataset_id'].unique())}")

    # Merge: Start with 15 datasets
    print(f"\n4. Merging datasets...")

    # Identify overlapping dataset IDs between 32 and 10_final
    ids_32 = set(df_32['dataset_id'].unique())
    ids_10 = set(df_10_final['dataset_id'].unique())
    overlap = ids_32 & ids_10
    print(f"   Overlapping datasets between 32 and 10_final: {sorted(overlap)}")

    # Remove overlapping IDs from df_32 (we'll use df_10_final results instead)
    df_32_filtered = df_32[~df_32['dataset_id'].isin(overlap)].copy()
    print(f"   After removing overlaps from 32: {len(df_32_filtered)} datasets")

    # Combine: 15 + filtered_32 + 10_final
    df_all = pd.concat([df_15, df_32_filtered, df_10_final], ignore_index=True)

    # Sort by dataset_id for consistency
    df_all = df_all.sort_values('dataset_id', ignore_index=True)

    print(f"\n5. Final merged dataset:")
    print(f"   Total rows: {len(df_all)}")
    print(f"   Dataset ID range: {int(df_all['dataset_id'].min())} to {int(df_all['dataset_id'].max())}")

    # Count status
    pass_count = len(df_all[df_all['status'] == 'PASS'])
    error_count = len(df_all[df_all['status'] == 'ERROR'])
    print(f"   Total PASS: {pass_count}")
    print(f"   Total ERROR: {error_count}")
    print(f"   Success rate: {100*pass_count/len(df_all):.1f}%")

    # Save merged results
    print(f"\n6. Saving merged results to {OUTPUT_FILE.name}...")
    df_all.to_csv(OUTPUT_FILE, index=False)
    print(f"   ✓ Saved {len(df_all)} datasets")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS (ALL 47 DATASETS)")
    print("=" * 80)

    passed_df = df_all[df_all['status'] == 'PASS']
    error_df = df_all[df_all['status'] == 'ERROR']

    print(f"\nOverall Results:")
    print(f"  PASS: {len(passed_df)}")
    print(f"  ERROR: {len(error_df)}")
    print(f"  Success rate: {100*len(passed_df)/len(df_all):.1f}%")

    if len(passed_df) > 0:
        print(f"\nSilhouette Coefficient (PASS only, n={len(passed_df)}):")
        print(f"  Mean:   {passed_df['silhouette'].mean():.4f}")
        print(f"  Median: {passed_df['silhouette'].median():.4f}")
        print(f"  Std:    {passed_df['silhouette'].std():.4f}")
        print(f"  Min:    {passed_df['silhouette'].min():.4f}")
        print(f"  Max:    {passed_df['silhouette'].max():.4f}")

        print(f"\nBootstrap ARI (5th percentile, PASS only):")
        print(f"  Mean:   {passed_df['bootstrap_ari_5th_percentile'].mean():.4f}")
        print(f"  Median: {passed_df['bootstrap_ari_5th_percentile'].median():.4f}")
        print(f"  Std:    {passed_df['bootstrap_ari_5th_percentile'].std():.4f}")
        print(f"  Min:    {passed_df['bootstrap_ari_5th_percentile'].min():.4f}")
        print(f"  Max:    {passed_df['bootstrap_ari_5th_percentile'].max():.4f}")

        # Statistics by source
        print(f"\nResults by Source (PASS only):")
        for source in sorted(passed_df['source'].unique()):
            source_df = passed_df[passed_df['source'] == source]
            if len(source_df) > 0:
                print(f"\n  {source} (n={len(source_df)}):")
                print(f"    Silhouette mean: {source_df['silhouette'].mean():.4f}")
                print(f"    Bootstrap ARI mean: {source_df['bootstrap_ari_5th_percentile'].mean():.4f}")

        # Statistics by domain
        print(f"\nResults by Domain (PASS only):")
        for domain in sorted(passed_df['domain'].unique()):
            domain_df = passed_df[passed_df['domain'] == domain]
            if len(domain_df) > 0:
                print(f"  {domain}: {len(domain_df)} PASS")

    # Error summary
    if len(error_df) > 0:
        print(f"\nError Summary ({len(error_df)} datasets):")
        for idx, row in error_df.iterrows():
            error_msg = row['error_message'][:60] if pd.notna(row['error_message']) else "Unknown"
            print(f"  ID {int(row['dataset_id'])}: {row['dataset_name']} - {error_msg}...")

    print("\n" + "=" * 80)
    print(f"✓ MERGE COMPLETE!")
    print(f"  Total datasets: {len(df_all)}")
    print(f"  Passed: {len(passed_df)}")
    print(f"  Failed: {len(error_df)}")
    print(f"  Output: {OUTPUT_FILE}")
    print("=" * 80)

    return True

if __name__ == '__main__':
    success = merge_results()
    sys.exit(0 if success else 1)
