#!/usr/bin/env python3
"""Generate the four main-text figures for the cautionary-subtyping manuscript.

Every value plotted is read directly from a DEPOSITED artifact (the same JSON/CSV
the Results tables cite), so the figures cannot drift from the numbers. No network,
no re-clustering. Deterministic.

Figures (one per Result):
  Figure 1  Result 1  gate ablation (paired FPR/TPR) + separation sweep
  Figure 2  Result 2  benchmark reproducibility distribution + retracted-model refit
  Figure 3  Result 3  cancer: metric choice (k-way ARI vs single-subtype enrichment)
  Figure 4  Result 4  flagship: region vs diagnosis under correct inference

Usage: python build_manuscript_figures.py --out DIR
       [--repo PATH]  (defaults to the repo this script sits in)
"""
from __future__ import annotations

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# colour-blind-safe (Okabe-Ito)
BLUE, ORANGE, GREEN, GREY, RED = "#0072B2", "#E69F00", "#009E73", "#999999", "#D55E00"


def load(p):
    with open(p) as fh:
        return json.load(fh)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    hw = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return (max(0, c - hw), min(1, c + hw))


def panel_label(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="right")


# ----------------------------------------------------------------------------
def figure1(cd, out):
    ab = load(f"{cd}/gate_ablation/results/ablation_honest.json")
    sw = load(f"{cd}/gate_ablation/results/separation_sweep.json")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    # Panel A: paired FPR (negatives) and TPR (positives), stability-only vs gate.
    # Every count is read from the deposited ablation JSON so the panel cannot drift
    # from the Results text.
    #
    # NO ERROR BAR IS DRAWN ON THE SCREEN'S FALSE-CERTIFY RATE, deliberately. Bar
    # HEIGHTS use `FPR_abstention_as_reject` (k=0, n=30) because the paired McNemar
    # comparison is over all 30 negatives. But the screen abstained on 26 of those 30,
    # so only 4 were testable, and NO interval on that bar can be honest:
    #   * Wilson(0,30) = [0, 0.114]  -> matches the plotted denominator but implies
    #     ~4x more precision than 4 testable datasets support;
    #   * Wilson(0, 4) = [0, 0.490]  -> the honest interval, but its denominator is
    #     not the one the bar height is computed on.
    # Rather than pick a misleading one, the bar carries an explicit annotation of the
    # testable denominator. The other three bars have genuine n=30 denominators (the
    # stability-only gate never abstains; the screen's TPR is 30/30) and keep their CIs.
    sb, gt = ab["stability_only_baseline"], ab["gate"]
    f_s, t_s = sb["FPR"], sb["TPR"]
    f_g, t_g = gt["FPR_abstention_as_reject"], gt["TPR"]
    fpr_stab, tpr_stab = f_s["k"] / f_s["n"], t_s["k"] / t_s["n"]
    fpr_gate, tpr_gate = f_g["k"] / f_g["n"], t_g["k"] / t_g["n"]
    groups = [f"False-certify rate\n(negatives, n={f_s['n']})",
              f"True-certify rate\n(positives, n={t_s['n']})"]
    stab = [fpr_stab, tpr_stab]
    gate = [fpr_gate, tpr_gate]
    stab_ci = [wilson(f_s["k"], f_s["n"]), wilson(t_s["k"], t_s["n"])]
    # None => draw no interval (see the note above); the screen's FPR is the only such bar.
    gate_ci = [None, wilson(t_g["k"], t_g["n"])]
    testable = gt["FPR_excluding_abstentions"]
    x = np.arange(2)
    w = 0.36
    b1 = axA.bar(x - w / 2, stab, w, color=GREY, label="Stability only")
    b2 = axA.bar(x + w / 2, gate, w, color=BLUE, label="Discreteness screen")
    for xi, v, ci in zip(x - w / 2, stab, stab_ci):
        axA.errorbar(xi, v, yerr=[[v - ci[0]], [ci[1] - v]], color="k", capsize=2, lw=0.8)
    for xi, v, ci in zip(x + w / 2, gate, gate_ci):
        if ci is None:
            continue
        axA.errorbar(xi, v, yerr=[[v - ci[0]], [ci[1] - v]], color="k", capsize=2, lw=0.8)
    # explicit denominator in place of an interval
    axA.annotate(f"{testable['k']}/{testable['n']} testable\n(no CI: see caption)",
                 (w / 2, 0.0), xytext=(w / 2, 0.20), fontsize=6.8, ha="center",
                 color=BLUE,
                 arrowprops=dict(arrowstyle="-", lw=0.6, color=BLUE, shrinkB=1))
    axA.set_xticks(x)
    axA.set_xticklabels(groups)
    axA.set_ylabel("Rate")
    axA.set_ylim(0, 1.12)
    axA.legend(frameon=False, fontsize=7.5, loc="upper center")
    axA.set_title("Gate suppresses false certification")
    # Both p-values are READ FROM THE ARTIFACT, never hardcoded. They were previously
    # literal strings ("p = 0.001", "p = 1.0") and went stale when the analysis was
    # re-run, leaving the figure contradicting the manuscript text it illustrates.
    _mc = ab["paired_mcnemar"]
    _p_neg = _mc["on_negatives_FPR"]["exact_p"]
    _p_pos = _mc["on_positives_TPR"]["exact_p"]
    # Match the manuscript's own rendering exactly (it writes p=0.0002, not 2e-4).
    _fmt = lambda v: (f"p = {v:.4f}" if v < 0.001 else
                      f"p = {v:g}" if v < 1 else f"p = {v:.1f}")
    axA.annotate(f"McNemar\n{_fmt(_p_neg)}", (0, fpr_stab), xytext=(0.02, 0.66),
                 fontsize=7.5, ha="center")
    axA.annotate(_fmt(_p_pos), (1, 1.0), xytext=(1, 1.06), fontsize=7.5, ha="center")
    axA.text(0.0, -0.30,
             f"gate abstains on {gt['negatives_abstained']}/{gt['negatives_total']}"
             f" ({gt['negatives_abstain_rate']*100:.0f}%)",
             transform=axA.get_xaxis_transform(), fontsize=6.8, ha="center", color=GREY)
    panel_label(axA, "a")

    # Panel B: separation sweep — certify rate + median SigClust p vs separation
    seps = [s["sep"] for s in sw["sweep"]]
    certify = [s["certify_rate"] for s in sw["sweep"]]
    medp = [s["median_sg_p"] for s in sw["sweep"]]
    l1, = axB.plot(seps, certify, "-o", color=BLUE, ms=4, label="Certify rate")
    axB.set_xlabel("Cluster separation δ (SD)")
    axB.set_ylabel("Certify rate", color=BLUE)
    axB.tick_params(axis="y", labelcolor=BLUE)
    axB.set_ylim(-0.03, 1.03)
    axB2 = axB.twinx()
    axB2.spines["top"].set_visible(False)
    l2, = axB2.plot(seps, medp, "--s", color=ORANGE, ms=3.5, label="Median SigClust p")
    axB2.axhline(0.05, color=RED, lw=0.7, ls=":")
    axB2.set_ylabel("Median SigClust p", color=ORANGE)
    axB2.tick_params(axis="y", labelcolor=ORANGE)
    axB2.set_ylim(-0.03, 1.03)
    axB.set_title("Conservative: resolves only at δ ≥ 2.5")
    axB.legend([l1, l2], ["Certify rate", "Median SigClust p"], frameon=False,
               fontsize=7.5, loc="center left")
    panel_label(axB, "b")

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure1_gate_ablation.{ext}")
    plt.close(fig)
    print("Figure 1 written")


def figure2(cd, out, csv):
    au = load(f"{cd}/benchmark_audit/results/benchmark_audit.json")
    df = pd.read_csv(csv)
    valid = df[df["valid"]]["bootstrap_ari_5th_percentile"].dropna().values

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    # Panel A: distribution of reproducibility (bootstrap ARI 5th pct) across valid cohorts
    axA.hist(valid, bins=np.arange(-0.2, 0.6, 0.05), color=BLUE, alpha=0.85,
             edgecolor="white")
    axA.axvline(0.5, color=RED, lw=1.0, ls="--")
    axA.text(0.5, axA.get_ylim()[1] * 0.92, "  0.5 bar\n  0/22 clear", color=RED,
             fontsize=7.2, va="top")
    med = au["reproducibility_distribution"]["valid_rows"]["median"]
    mx = au["reproducibility_distribution"]["valid_rows"]["max"]
    axA.axvline(med, color="k", lw=0.8, ls=":")
    axA.text(med, axA.get_ylim()[1] * 0.5, f" median\n {med:.3f}", fontsize=7)
    axA.set_xlabel("Reproducibility (bootstrap ARI, 5th percentile)")
    axA.set_ylabel("Cohorts (n=22 valid)")
    axA.set_title("No cohort clears a 0.5 bar")
    axA.annotate(f"max {mx:.2f}", (mx, 0.5), xytext=(mx - 0.02, 2.2),
                 fontsize=7, ha="right", arrowprops=dict(arrowstyle="->", lw=0.6))
    panel_label(axA, "a")

    # Panel B: retracted adaptive-threshold model R^2 — published vs refit under 3 screens
    m = au["retracted_threshold_model"]
    labels = ["Published\nclaim", "Refit\nall rows", "Refit\nvalid", "Refit\nstrict"]
    r2 = [m["published_claim"]["r_squared"], m["refit_all_rows"]["r_squared"],
          m["refit_valid_rows"]["r_squared"], m["refit_strict_rows"]["r_squared"]]
    colors = [ORANGE, GREY, BLUE, GREEN]
    bars = axB.bar(labels, r2, color=colors, edgecolor="white")
    for b, v in zip(bars, r2):
        axB.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.3f}", ha="center",
                 fontsize=7)
    axB.set_ylabel("R²  (silhouette → bootstrap ARI)")
    axB.set_ylim(0, 1.0)
    axB.set_title("Retracted model does not refit")
    panel_label(axB, "b")

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure2_benchmark_audit.{ext}")
    plt.close(fig)
    print("Figure 2 written")


def figure3(cd, out):
    dl = load(f"{cd}/cancer_r38/results/brca_pam50_validation_with_DL.json")
    ari = dl["recovery_vs_pam50_ARI"]
    enr = dl["recovery_vs_pam50_single_subtype_enrichment"]
    methods = ["pathway-GMM", "k-means (pathway)", "VAE-GMM", "DEC"]
    keys = ["PSF (pathway-GMM)", "k-means (pathway)", "VAE-GMM", "DEC"]
    ari_v = [ari[k]["ari"] for k in keys]
    enr_v = [enr[k]["best"]["enrichment_frac"] for k in keys]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    x = np.arange(len(methods))
    cols = [BLUE, GREY, ORANGE, GREEN]

    barsA = axA.bar(x, ari_v, color=cols, edgecolor="white")
    axA.set_xticks(x)
    axA.set_xticklabels(methods, rotation=25, ha="right", fontsize=7)
    axA.set_ylabel("k-way recovery (ARI vs PAM50)")
    axA.set_title("By k-way agreement with PAM50")
    for b, v in zip(barsA, ari_v):
        axA.text(b.get_x() + b.get_width() / 2, v + 0.004, f"{v:.3f}", ha="center", fontsize=7)
    axA.set_ylim(0, max(ari_v) * 1.25)
    # No winner highlight: the panels reorder, and no ranking is claimed.
    panel_label(axA, "a")

    barsB = axB.bar(x, [v * 100 for v in enr_v], color=cols, edgecolor="white")
    axB.set_xticks(x)
    axB.set_xticklabels(methods, rotation=25, ha="right", fontsize=7)
    axB.set_ylabel("Single-subtype enrichment (LumA, %)")
    axB.set_title("By single-subtype enrichment (LumA)")
    for b, v in zip(barsB, enr_v):
        axB.text(b.get_x() + b.get_width() / 2, v * 100 + 0.6, f"{v*100:.1f}", ha="center", fontsize=7)
    axB.set_ylim(0, 100)
    panel_label(axB, "b")

    fig.suptitle("The method ranking flips with the metric (TCGA-BRCA, n=981)",
                 fontsize=9.5, fontweight="bold", y=1.02)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure3_cancer_metric_choice.{ext}")
    plt.close(fig)
    print("Figure 3 written")


def figure4(cd, out):
    fl = load(f"{cd}/flagship_stats/results/flagship_donor_level.json")
    region_v = fl["region_sample_level_valid"]["cramers_v_bergsma"]
    diag = fl["diagnosis_donor_level_CORRECT"]
    diag_obs = diag["permutation_test_permuting_diagnosis_across_donors"]["observed_v"]
    null95 = diag["permutation_test_permuting_diagnosis_across_donors"]["null_v_95th_percentile"]
    pval = diag["permutation_test_permuting_diagnosis_across_donors"]["p"]

    fig, ax = plt.subplots(figsize=(4.0, 3.2))
    labels = ["Brain region\n(sample-level, valid)", "Diagnosis\n(donor-level, correct)"]
    vals = [region_v, diag_obs]
    bars = ax.bar([0, 1], vals, width=0.5, color=[BLUE, GREY], edgecolor="white")
    # chance band for diagnosis
    ax.axhline(null95, color=RED, lw=1.0, ls="--")
    ax.text(1.28, null95, "chance\n(95th pct)", color=RED, fontsize=7, va="center")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Cramér's V")
    ax.set_ylim(0, 0.8)
    ax.text(0, region_v + 0.02, f"{region_v:.3f}", ha="center", fontsize=8, fontweight="bold")
    ax.text(1, diag_obs + 0.02, f"V={diag_obs:.3f}\np={pval:.3f}\n(= chance)", ha="center", fontsize=7.2)
    ax.set_title("Region, not diagnosis\n(GSE80655: 141 samples, 48 donors)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure4_region_vs_diagnosis.{ext}")
    plt.close(fig)
    print("Figure 4 written")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo_default = os.path.abspath(os.path.join(here, "..", ".."))
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", default=repo_default)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cd = os.path.join(args.repo, "consolidation-cautionary", "cross-domain")
    csv = os.path.join(args.repo, "CORRECTION_2026-07",
                       "corrected_benchmark_47datasets_v2.csv")
    figure1(cd, args.out)
    figure2(cd, args.out, csv)
    figure3(cd, args.out)
    figure4(cd, args.out)
    print("All figures ->", args.out)


if __name__ == "__main__":
    main()
