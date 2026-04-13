#!/usr/bin/env python3
"""
Merge all benchmark result CSV files into a single 47-dataset file.

Combines results from:
1. bootstrap_threshold_calibration_15datasets_zenodo.csv (12 TCGA + 3 GEO test)
2. bootstrap_threshold_calibration_geo_remaining_32.csv (22 GEO with PASS status)
3. bootstrap_threshold_calibration_geo_final_10.csv (remaining 10 GEO) [optional - may not exist]

Output: bootstrap_threshold_calibration_47datasets_final_zenodo.csv
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
OUTPUT_FILE = RESULTS_DIR / 'bootstrap_threshold_calibration_47datasets_final_zenodo.csv'

def merge_results():
    """Merge all result files into single 47-dataset CSV."""

    print("=" * 80)
    print("MERGING 47-DATASET BENCHMARK RESULTS")
    print("=" * 80)

    # Load 15-dataset file (12 TCGA + 3 GEO test)
    print(f"\n1. Loading {FILE_15_ZENODO.name}...")
    if not FILE_15_ZENODO.exists():
        print(f"ERROR: File not found: {FILE_15_ZENODO}")
        return False

    df_15 = pd.read_csv(FILE_15_ZENODO)
    print(f"   Loaded {len(df_15)} datasets")
    print(f"   Dataset IDs: {sorted(df_15['dataset_id'].unique())}")

    # Load 32-dataset file and filter to PASS entries
    print(f"\n2. Loading {FILE_32_REMAINING.name}...")
    if not FILE_32_REMAINING.exists():
        print(f"ERROR: File not found: {FILE_32_REMAINING}")
        return False

    df_32 = pd.read_csv(FILE_32_REMAINING)
    print(f"   Loaded {len(df_32)} rows total")

    # Filter to PASS status only
    df_32_pass = df_32[df_32['status'] == 'PASS'].copy()
    print(f"   Filtered to {len(df_32_pass)} PASS entries")
    print(f"   Dataset IDs (PASS): {sorted(df_32_pass['dataset_id'].unique())}")
    print(f"   Dataset IDs (ERROR): {sorted(df_32[df_32['status'] == 'ERROR']['dataset_id'].unique())}")

    # Try to load final 10 file if it exists
    df_10_final = None
    if FILE_10_FINAL.exists():
        print(f"\n3. Loading {FILE_10_FINAL.name}...")
        df_10_final = pd.read_csv(FILE_10_FINAL)
        print(f"   Loaded {len(df_10_final)} datasets")

        # Filter to PASS only
        df_10_pass = df_10_final[df_10_final['status'] == 'PASS'].copy()
        print(f"   Filtered to {len(df_10_pass)} PASS entries")

        # Use final 10 results to replace/supplement the 32 results
        # This prioritizes the final 10 run results if available
        if len(df_10_pass) > 0:
            print(f"   Using final 10 results for datasets: {sorted(df_10_pass['dataset_id'].unique())}")
            # Replace any overlapping dataset IDs with the final 10 results
            df_32_pass = df_32_pass[~df_32_pass['dataset_id'].isin(df_10_pass['dataset_id'])]
            df_32_pass = pd.concat([df_32_pass, df_10_pass], ignore_index=True)
            print(f"   Merged dataset count: {len(df_32_pass)}")
    else:
        print(f"\n3. {FILE_10_FINAL.name} not found (validation may still be running)")

    # Combine 15 + 32_pass results
    print(f"\n4. Merging all results...")
    df_all = pd.concat([df_15, df_32_pass], ignore_index=True)

    # Sort by dataset_id for consistency
    df_all = df_all.sort_values('dataset_id', ignore_index=True)

    print(f"   Total datasets: {len(df_all)}")
    print(f"   Dataset ID range: {int(df_all['dataset_id'].min())} to {int(df_all['dataset_id'].max())}")
    print(f"   Total PASS: {len(df_all[df_all['status'] == 'PASS'])}")
    print(f"   Total ERROR: {len(df_all[df_all['status'] == 'ERROR'])}")

    # Save merged results
    print(f"\n5. Saving merged results to {OUTPUT_FILE.name}...")
    df_all.to_csv(OUTPUT_FILE, index=False)
    print(f"   ✓ Saved {len(df_all)} datasets")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    passed_df = df_all[df_all['status'] == 'PASS']
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
        print(f"\nResults by Source:")
        for source in sorted(df_all['source'].unique()):
            source_df = passed_df[passed_df['source'] == source]
            if len(source_df) > 0:
                print(f"\n  {source} (n={len(source_df)}):")
                print(f"    Silhouette mean: {source_df['silhouette'].mean():.4f}")
                print(f"    Bootstrap ARI mean: {source_df['bootstrap_ari_5th_percentile'].mean():.4f}")

    print("\n" + "=" * 80)
    print(f"✓ Merge complete! Results saved to:")
    print(f"  {OUTPUT_FILE}")
    print("=" * 80)

    return True

if __name__ == '__main__':
    success = merge_results()
    sys.exit(0 if success else 1)
