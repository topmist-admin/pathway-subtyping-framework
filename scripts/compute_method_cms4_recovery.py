#!/usr/bin/env python3
"""
Compute CMS4 recovery rates for alternative methods on TCGA-COAD.

Runs PCA+k-means, gene-level k-means, and NMF clustering on TCGA-COAD,
then validates each against CMS classification using Fisher exact test.

Usage:
    cd /Users/rohitchauhan/Downloads/AI-Genetic-Research/pathway-subtyping-framework
    source pathwayenv/bin/activate
    python scripts/compute_method_cms4_recovery.py

Output:
    - method_cms4_recovery_comparison.csv (results table)
    - Updated results_summary.json with per-method CMS4 rates
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

try:
    from pathway_subtyping.benchmark import (
        BenchmarkMethod,
        run_single_benchmark,
    )
except ImportError:
    print("ERROR: Cannot import pathway_subtyping. Please activate the framework venv:")
    print("  source pathwayenv/bin/activate")
    sys.exit(1)

# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_DIR = Path(
    "/Users/rohitchauhan/Downloads/AI-Genetic-Research/pathway-subtyping-framework/examples/notebooks/research-results/tcga"
)
RESULTS_FILE = OUTPUT_DIR / "results_summary.json"
EXPRESSION_FILE = OUTPUT_DIR / "gene_expression_processed.csv"
PATHWAY_SCORES_FILE = OUTPUT_DIR / "pathway_scores.csv"
CMS_FILE = OUTPUT_DIR / "cms_ntp_results.csv"
METADATA_FILE = OUTPUT_DIR / "sample_metadata_with_subtypes.csv"

print("=" * 80)
print("COMPUTING CMS4 RECOVERY RATES FOR ALTERNATIVE METHODS")
print("=" * 80)

# =============================================================================
# STEP 1: Verify all files exist
# =============================================================================

print("\n[1/5] Checking for required data files...")

required_files = [
    (EXPRESSION_FILE, "Gene expression data"),
    (PATHWAY_SCORES_FILE, "Pathway scores"),
    (CMS_FILE, "CMS predictions"),
    (METADATA_FILE, "Sample metadata"),
    (RESULTS_FILE, "Existing results"),
]

for fpath, desc in required_files:
    if not fpath.exists():
        print(f"  ✗ MISSING: {desc}")
        print(f"    File: {fpath}")
        print("\nERROR: Please ensure NB17 (tcga_cancer_validation.ipynb) has been run.")
        sys.exit(1)
    print(f"  ✓ {desc}")

# =============================================================================
# STEP 2: Load data
# =============================================================================

print("\n[2/5] Loading data...")

expression_data = pd.read_csv(EXPRESSION_FILE, index_col=0)
pathway_scores = pd.read_csv(PATHWAY_SCORES_FILE, index_col=0)
cms_data = pd.read_csv(CMS_FILE, index_col="sample_id")
metadata = pd.read_csv(METADATA_FILE, index_col="sample_id")

print(f"  ✓ Expression: {expression_data.shape}")
print(f"  ✓ Pathway scores: {pathway_scores.shape}")
print(f"  ✓ CMS predictions: {cms_data.shape}")
print(f"  ✓ Sample metadata: {metadata.shape}")

# Align samples
common_samples = expression_data.index.intersection(cms_data.index).intersection(
    pathway_scores.index
)
print(f"  ✓ Common samples across all datasets: {len(common_samples)}")

expression_data = expression_data.loc[common_samples]
pathway_scores = pathway_scores.loc[common_samples]
cms_data = cms_data.loc[common_samples]

# Load existing results
with open(RESULTS_FILE) as f:
    results = json.load(f)

print(f"  ✓ Loaded existing results from {RESULTS_FILE}")

# =============================================================================
# STEP 3: Extract pathway_gmm clusters and get baseline
# =============================================================================

print("\n[3/5] Extracting pathway_gmm baseline...")

pathway_gmm_clusters = metadata.loc[common_samples, "subtype"].values
pathway_gmm_silhouette = results["benchmark"]["silhouette_scores"]["pathway_gmm"]
pathway_gmm_cms4_recovery = results["external_comparison"]["cms_ntp_validation"]["key_alignments"][
    "subtype_0_to_CMS4_pct"
]
pathway_gmm_fisher_or = results["external_comparison"]["cms_ntp_validation"]["key_alignments"][
    "fisher_S0_CMS4_OR"
]
pathway_gmm_fisher_p = results["external_comparison"]["cms_ntp_validation"]["key_alignments"][
    "fisher_S0_CMS4_p"
]

print(
    f"  ✓ pathway_gmm: silhouette={pathway_gmm_silhouette:.4f}, CMS4 recovery={pathway_gmm_cms4_recovery:.1%}"
)

# =============================================================================
# STEP 4: Run alternative methods and compute CMS4 recovery
# =============================================================================

print("\n[4/5] Running alternative methods and computing CMS4 recovery rates...")

methods_to_test = [
    BenchmarkMethod.PCA_KMEANS,
    BenchmarkMethod.GENE_KMEANS,
    BenchmarkMethod.NMF_CLUSTERING,
]

comparison_results = {
    "pathway_gmm": {
        "silhouette": pathway_gmm_silhouette,
        "cms4_recovery": pathway_gmm_cms4_recovery,
        "fisher_or": pathway_gmm_fisher_or,
        "fisher_p": pathway_gmm_fisher_p,
    }
}

# Get CMS4 classification
is_cms4_vector = (cms_data.loc[common_samples, "cms_prediction"] == "CMS4").values

for method in methods_to_test:
    method_name = method.value
    print(f"\n  {method_name}:")

    try:
        # Run clustering
        result = run_single_benchmark(
            method=method,
            gene_burdens=expression_data,
            pathway_scores=pathway_scores,
            n_clusters=3,  # Match pathway_gmm's k=3
            true_labels=None,
            seed=42,
        )

        clusters = result.predicted_labels
        silhouette = result.silhouette

        print(f"    Silhouette: {silhouette:.4f}")

        # Find best cluster for CMS4 recovery (highest recovery rate + p-value)
        unique_clusters = np.unique(clusters)
        best_cms4_rate = 0
        best_fisher_or = 0
        best_fisher_p = 1.0
        best_cluster_id = None

        for cluster_id in unique_clusters:
            is_in_cluster = clusters == cluster_id

            # Build 2x2 contingency table: [in_cluster & cms4, in_cluster & ~cms4]
            #                               [~in_cluster & cms4, ~in_cluster & ~cms4]
            table_2x2 = [
                [(is_in_cluster & is_cms4_vector).sum(), (is_in_cluster & ~is_cms4_vector).sum()],
                [(~is_in_cluster & is_cms4_vector).sum(), (~is_in_cluster & ~is_cms4_vector).sum()],
            ]

            # Fisher exact test
            odds_ratio, p_value = fisher_exact(table_2x2)

            # CMS4 recovery rate
            cms4_count = (is_in_cluster & is_cms4_vector).sum()
            total_in_cluster = is_in_cluster.sum()
            cms4_rate = cms4_count / total_in_cluster if total_in_cluster > 0 else 0

            # Keep track of best
            if cms4_rate > best_cms4_rate:
                best_cms4_rate = cms4_rate
                best_fisher_or = odds_ratio
                best_fisher_p = p_value
                best_cluster_id = cluster_id

            print(
                f"    Cluster {cluster_id}: {cms4_rate:.1%} CMS4 (n={total_in_cluster}, "
                f"OR={odds_ratio:.2f}, p={p_value:.2e})"
            )

        comparison_results[method_name] = {
            "silhouette": silhouette,
            "cms4_recovery": best_cms4_rate,
            "fisher_or": best_fisher_or,
            "fisher_p": best_fisher_p,
            "best_cluster": int(best_cluster_id),
        }

        print(
            f"    → Best cluster: Cluster {best_cluster_id} with {best_cms4_rate:.1%} CMS4 recovery"
        )

    except Exception as e:
        print(f"    ✗ ERROR: {e}")
        import traceback

        traceback.print_exc()
        continue

# =============================================================================
# STEP 5: Compile and output results
# =============================================================================

print("\n[5/5] Compiling results...")

# Build comparison dataframe
comparison_df = pd.DataFrame(
    [
        {
            "Method": method_name,
            "Silhouette": data["silhouette"],
            "CMS4 Recovery (%)": data["cms4_recovery"] * 100,
            "Fisher OR": data["fisher_or"],
            "Fisher p-value": data["fisher_p"],
        }
        for method_name, data in comparison_results.items()
    ]
)

# Sort: pathway_gmm first, then others by CMS4 recovery descending
comparison_df = pd.concat(
    [
        comparison_df[comparison_df["Method"] == "pathway_gmm"],
        comparison_df[comparison_df["Method"] != "pathway_gmm"].sort_values(
            "CMS4 Recovery (%)", ascending=False
        ),
    ],
    ignore_index=True,
)

# Save to CSV
output_csv = OUTPUT_DIR / "method_cms4_recovery_comparison.csv"
comparison_df.to_csv(output_csv, index=False)
print(f"  ✓ Saved comparison table to: {output_csv}")

# Update results_summary.json
results["benchmark"]["cms4_recovery_rates"] = {
    row["Method"]: {
        "silhouette": float(row["Silhouette"]),
        "cms4_recovery": float(row["CMS4 Recovery (%)"]) / 100,
        "fisher_or": float(row["Fisher OR"]),
        "fisher_p": float(row["Fisher p-value"]),
    }
    for _, row in comparison_df.iterrows()
}

with open(RESULTS_FILE, "w") as f:
    json.dump(results, f, indent=2)

print(f"  ✓ Updated {RESULTS_FILE}")

# =============================================================================
# DISPLAY RESULTS
# =============================================================================

print("\n" + "=" * 80)
print("RESULTS FOR MANUSCRIPT")
print("=" * 80)

print("\n### Result 1b: Method Comparison on TCGA-COAD (n=452)\n")
print(comparison_df.to_string(index=False))

print("\n" + "=" * 80)
print("MARKDOWN TABLE FOR MANUSCRIPT")
print("=" * 80)

print("\n| Method | Silhouette | CMS4 Recovery | Fisher OR | p-value |")
print("|--------|-----------|---------------|-----------|---------|")
for _, row in comparison_df.iterrows():
    method = row["Method"]
    sil = row["Silhouette"]
    cms4 = row["CMS4 Recovery (%)"]
    or_val = row["Fisher OR"]
    p_val = row["Fisher p-value"]

    # Format method name
    method_display = {
        "pathway_gmm": "**Pathway+GMM (Framework)**",
        "pca_kmeans": "PCA+k-means",
        "gene_kmeans": "Gene-level k-means",
        "nmf_clustering": "NMF",
    }.get(method, method)

    p_str = f"{p_val:.1e}" if p_val > 0 else "n.s."
    print(f"| {method_display:.<40} | {sil:.3f} | {cms4:.0f}% | {or_val:.1f} | {p_str} |")

print("\n" + "=" * 80)
print("✓ COMPLETE - Results ready for manuscript integration")
print("=" * 80)
