#!/usr/bin/env python3
"""
⚠️ RETRACTED (2026-07-08) — this script plots a model that does not reproduce.

The bilinear bootstrap-threshold model diagnosed here (R²=0.889, slope 0.914) is
**retracted**: refit on the released benchmark it gives R²≈0.111 across all rows,
≈0.015 on valid rows, and ≈0.001 under a stricter ground-truth screen, with the slope
reversing sign. Silhouette does not predict bootstrap-ARI reproducibility in these
data. The Panel D "TCGA-COAD overlay" is drawn against an independence claim that is
also withdrawn — TCGA-COAD is present in the benchmark, not excluded from it.

Kept so the retracted figure remains regenerable for the record. Do not use its output
as evidence for a silhouette-calibrated threshold. The DOI below names version 1.1, the
uncorrected release this model was fitted on; it is provenance, not a citation. See
CORRECTION_2026-07/ERRATUM_2026-07-08.md.

Generate bilinear fit diagnostic plots for the 47-dataset
bootstrap threshold calibration model.

Produces a 2x2 diagnostic figure:
  Panel A: Scatter + bilinear fit (silhouette vs bootstrap ARI)
  Panel B: Residuals vs predicted
  Panel C: LOOCV predicted vs observed
  Panel D: Domain-colored scatter with TCGA-COAD overlay

Output: research-results/benchmarks/bilinear_fit_diagnostic_plots.png

Reference: Methods §1.5 of Nature Methods manuscript
Zenodo DOI: 10.5281/zenodo.19324360
"""

import sys
import os
import numpy as np
import pandas as pd

# Paths
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(
    REPO, "research-results", "benchmarks",
    "bootstrap_threshold_calibration_47datasets_zenodo.csv"
)
OUTPUT_PATH = os.path.join(
    REPO, "research-results", "benchmarks",
    "bilinear_fit_diagnostic_plots.png"
)

# Model parameters (from manuscript Methods §1.5)
SLOPE = 0.914
INTERCEPT = 0.452
R_SQUARED = 0.889
LOOCV_RMSE = 0.051


def piecewise_threshold(s):
    """Piecewise linear threshold model."""
    if s >= 0.40:
        return 0.80
    elif s >= 0.20:
        return 0.45 + 1.75 * (s - 0.20)
    else:
        return 0.45


def linear_threshold(s):
    """Linear regression threshold."""
    return INTERCEPT + SLOPE * s


def main():
    # Load data
    df = pd.read_csv(CSV_PATH)

    # Filter to PASS-only rows (have valid silhouette and bootstrap ARI)
    mask = (df["status"] == "PASS") & df["silhouette"].notna() & df["bootstrap_ari_5th_percentile"].notna()
    df_pass = df[mask].copy()

    sil = df_pass["silhouette"].values
    bari = df_pass["bootstrap_ari_5th_percentile"].values
    domains = df_pass["domain"].values
    names = df_pass["dataset_name"].values

    print(f"Loaded {len(df)} total datasets, {len(df_pass)} with valid results")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("ERROR: matplotlib not installed. Install with: pip install matplotlib")
        sys.exit(1)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(
        "Bilinear Bootstrap Threshold Calibration — 47-Dataset Diagnostics\n"
        f"R² = {R_SQUARED}, LOOCV RMSE = {LOOCV_RMSE}, DOI: 10.5281/zenodo.19324360",
        fontsize=13, fontweight="bold"
    )

    # Color map for domains
    domain_colors = {
        "Oncology": "#e74c3c",
        "Psychiatry": "#3498db",
        "Immunology": "#2ecc71",
        "Single-Cell": "#9b59b6",
        "Other": "#95a5a6",
    }

    # --- Panel A: Scatter + bilinear fit ---
    ax = axes[0, 0]
    for domain, color in domain_colors.items():
        mask_d = domains == domain
        if mask_d.any():
            ax.scatter(sil[mask_d], bari[mask_d], c=color, label=domain,
                       alpha=0.7, s=50, edgecolors="k", linewidths=0.5)

    # Bilinear fit line
    s_range = np.linspace(-0.05, 1.0, 200)
    t_piecewise = [piecewise_threshold(s) for s in s_range]
    t_linear = [linear_threshold(s) for s in s_range]
    ax.plot(s_range, t_piecewise, "r-", linewidth=2, label="Piecewise threshold")
    ax.plot(s_range, t_linear, "b--", linewidth=1.5, alpha=0.6, label=f"Linear fit (R²={R_SQUARED})")

    # TCGA-COAD overlay
    ax.scatter([0.164], [0.485], c="orange", s=150, marker="*", zorder=5,
               edgecolors="k", linewidths=1, label="TCGA-COAD (0.164, 0.485)")

    ax.set_xlabel("Silhouette Score", fontsize=11)
    ax.set_ylabel("Bootstrap ARI (5th percentile)", fontsize=11)
    ax.set_title("A. Scatter + Bilinear Fit", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_xlim(-0.1, 1.0)
    ax.set_ylim(-0.2, 1.1)
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.3)
    ax.grid(True, alpha=0.2)

    # --- Panel B: Residuals vs predicted ---
    ax = axes[0, 1]
    predicted = np.array([linear_threshold(s) for s in sil])
    residuals = bari - predicted

    for domain, color in domain_colors.items():
        mask_d = domains == domain
        if mask_d.any():
            ax.scatter(predicted[mask_d], residuals[mask_d], c=color, alpha=0.7,
                       s=50, edgecolors="k", linewidths=0.5)

    ax.axhline(y=0, color="red", linestyle="-", linewidth=1)
    ax.axhline(y=2 * np.std(residuals), color="gray", linestyle="--", alpha=0.5)
    ax.axhline(y=-2 * np.std(residuals), color="gray", linestyle="--", alpha=0.5)

    ax.set_xlabel("Predicted Bootstrap ARI", fontsize=11)
    ax.set_ylabel("Residual (Observed - Predicted)", fontsize=11)
    ax.set_title("B. Residuals vs Predicted", fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.2)

    # --- Panel C: LOOCV predicted vs observed ---
    ax = axes[1, 0]
    loocv_predicted = []
    for i in range(len(sil)):
        # Leave-one-out: fit on all except i
        sil_loo = np.delete(sil, i)
        bari_loo = np.delete(bari, i)
        if len(sil_loo) > 1:
            coeffs = np.polyfit(sil_loo, bari_loo, 1)
            pred_i = coeffs[0] * sil[i] + coeffs[1]
        else:
            pred_i = bari[i]
        loocv_predicted.append(pred_i)
    loocv_predicted = np.array(loocv_predicted)

    for domain, color in domain_colors.items():
        mask_d = domains == domain
        if mask_d.any():
            ax.scatter(bari[mask_d], loocv_predicted[mask_d], c=color, alpha=0.7,
                       s=50, edgecolors="k", linewidths=0.5)

    # Perfect prediction line
    lims = [min(bari.min(), loocv_predicted.min()) - 0.05,
            max(bari.max(), loocv_predicted.max()) + 0.05]
    ax.plot(lims, lims, "r--", linewidth=1, alpha=0.7, label="Perfect prediction")

    loocv_rmse_actual = np.sqrt(np.mean((bari - loocv_predicted) ** 2))
    ax.set_xlabel("Observed Bootstrap ARI", fontsize=11)
    ax.set_ylabel("LOOCV Predicted Bootstrap ARI", fontsize=11)
    ax.set_title(f"C. LOOCV: Predicted vs Observed (RMSE={loocv_rmse_actual:.3f})",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    # --- Panel D: Domain breakdown with threshold zones ---
    ax = axes[1, 1]
    for domain, color in domain_colors.items():
        mask_d = domains == domain
        if mask_d.any():
            ax.scatter(sil[mask_d], bari[mask_d], c=color, label=domain,
                       alpha=0.7, s=50, edgecolors="k", linewidths=0.5)

    # Threshold zones
    ax.axvspan(-0.1, 0.20, alpha=0.08, color="blue", label="Gradient-like (threshold=0.45)")
    ax.axvspan(0.20, 0.40, alpha=0.08, color="yellow", label="Intermediate (interpolated)")
    ax.axvspan(0.40, 1.0, alpha=0.08, color="green", label="Discrete (threshold=0.80)")

    # TCGA-COAD
    ax.scatter([0.164], [0.485], c="orange", s=150, marker="*", zorder=5,
               edgecolors="k", linewidths=1)
    ax.annotate("TCGA-COAD", (0.164, 0.485), xytext=(0.25, 0.55),
                arrowprops=dict(arrowstyle="->", color="orange"),
                fontsize=9, fontweight="bold", color="orange")

    ax.set_xlabel("Silhouette Score", fontsize=11)
    ax.set_ylabel("Bootstrap ARI (5th percentile)", fontsize=11)
    ax.set_title("D. Threshold Zones by Domain", fontsize=12, fontweight="bold")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlim(-0.1, 1.0)
    ax.set_ylim(-0.2, 1.1)
    ax.grid(True, alpha=0.2)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    print(f"Saved diagnostic plots to: {OUTPUT_PATH}")

    # Print summary stats
    print(f"\n--- Model Diagnostics ---")
    print(f"N datasets (with results): {len(sil)}")
    print(f"Silhouette range: [{sil.min():.4f}, {sil.max():.4f}]")
    print(f"Bootstrap ARI range: [{bari.min():.4f}, {bari.max():.4f}]")
    print(f"Linear R²: {R_SQUARED}")
    print(f"LOOCV RMSE (computed): {loocv_rmse_actual:.4f}")
    print(f"Residual std: {np.std(residuals):.4f}")
    print(f"Mean absolute residual: {np.mean(np.abs(residuals)):.4f}")


if __name__ == "__main__":
    main()
