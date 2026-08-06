#!/usr/bin/env python3
"""SUPERSEDED 2026-08-06 — do not use for publication. Kept for provenance only.

This script renders the gate's false-positive rate as a bare `0.00`, which is the
**abstention-as-rejection** convention the analysis has since retracted: on the 30
negatives the screen removes all 13 stability-only false certifications by *abstaining*
(26/30) and none by rejecting, so a `0.00` bar credits the screen with rejections it
never made. The honest rendering, which shows the abstention decomposition and the
0/4 testable-negative denominator, is **Figure 1 of the manuscript**
(`consolidation-cautionary/scripts/build_manuscript_figures.py`, from
`ablation_honest.json`).

Its output, `results/gate_ablation_figure.{png,svg}`, was **withdrawn from the deposit
on 2026-08-06**: rendered 2026-07-23, it still printed the pre-seeding-fix `FPR 0.37`
and `TPR 0.97` against current values of 13/30 and 30/30, and labelled its arms
`pre-v0.8.0`/`v0.8.0` inside a v0.9.0 deposit. Because the numbers lived only as
rendered glyphs, no text sweep could see them.

If this figure is ever needed again, redesign panel A to plot three outcomes
(certify / abstain / reject) rather than a two-way TPR-vs-FPR split, and drop the
version-numbered arm labels.

Original purpose (reviewer R3.10): reads
outputs/gate_ablation/gate_ablation_results.json + gate_ablation_raw.csv and renders
  A. TPR vs FPR by gate subset.
  B. Per-condition certification rate — false positives concentrated on the 1-D
     continuum (the erratum failure mode).

Usage: python scripts/plot_gate_ablation.py [--out DIR]   (requires --i-know-this-is-superseded)
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# colorblind-safe (Okabe-Ito subset)
C_TPR, C_FPR = "#0072B2", "#D55E00"
COND_COLORS = {"discrete_k2": "#009E73", "discrete_k3": "#56B4E9",
               "single_gaussian": "#999999", "continuum_1d": "#D55E00"}
LABELS = {"stability_only": "stability\nonly\n(pre-v0.8.0)",
          "discreteness_only": "discreteness\nonly",
          "stability+discreteness": "stability +\ndiscreteness\n(v0.8.0)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("outputs/gate_ablation"))
    ap.add_argument("--i-know-this-is-superseded", action="store_true",
                    help="required: this figure uses the retracted abstention-as-rejection "
                         "convention (see module docstring)")
    args = ap.parse_args()
    if not getattr(args, "i_know_this_is_superseded"):
        raise SystemExit(
            "plot_gate_ablation.py is SUPERSEDED (2026-08-06) and will not run by default.\n"
            "It plots the gate's FPR as a bare 0.00, which is the abstention-as-rejection\n"
            "convention this analysis retracts: on the 30 negatives the screen removes all 13\n"
            "stability-only false certifications by abstaining (26/30) and none by rejecting.\n"
            "Use Figure 1 of the manuscript instead — consolidation-cautionary/scripts/\n"
            "build_manuscript_figures.py, which reads ablation_honest.json and shows the\n"
            "abstention decomposition. Pass --i-know-this-is-superseded to override.")
    d = json.loads((args.out / "gate_ablation_results.json").read_text())
    summ = d["summary"]
    policies = list(summ.keys())
    x = np.arange(len(policies))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: TPR / FPR grouped bars
    w = 0.38
    tpr = [summ[p]["TPR"] for p in policies]
    fpr = [summ[p]["FPR"] for p in policies]
    axA.bar(x - w / 2, tpr, w, label="TPR (recover real subtypes)", color=C_TPR)
    axA.bar(x + w / 2, fpr, w, label="FPR (false-certify no-structure)", color=C_FPR)
    for xi, (t, f) in enumerate(zip(tpr, fpr)):
        axA.text(xi - w / 2, t + 0.02, f"{t:.2f}", ha="center", fontsize=9)
        axA.text(xi + w / 2, f + 0.02, f"{f:.2f}", ha="center", fontsize=9)
    axA.set_xticks(x)
    axA.set_xticklabels([LABELS[p] for p in policies], fontsize=8.5)
    axA.set_ylabel("rate")
    axA.set_ylim(0, 1.12)
    axA.set_title("A  True- vs false-positive rate by gate subset", fontsize=11, loc="left")
    axA.legend(fontsize=8.5, loc="upper center")
    axA.spines[["top", "right"]].set_visible(False)

    # Panel B: per-condition certification rate
    conds = ["discrete_k2", "discrete_k3", "single_gaussian", "continuum_1d"]
    wB = 0.8 / len(conds)
    for j, c in enumerate(conds):
        vals = [summ[p]["cert_rate_by_condition"][c] for p in policies]
        axB.bar(x + (j - 1.5) * wB, vals, wB, label=c, color=COND_COLORS[c])
    axB.set_xticks(x)
    axB.set_xticklabels([LABELS[p] for p in policies], fontsize=8.5)
    axB.set_ylabel("certified-as-subtype rate")
    axB.set_ylim(0, 1.12)
    axB.set_title("B  Certification rate by data type", fontsize=11, loc="left")
    axB.legend(fontsize=8, title="synthetic condition", title_fontsize=8)
    axB.spines[["top", "right"]].set_visible(False)
    axB.annotate("continuum false-\npositives eliminated",
                 xy=(1.5, 0.02), xytext=(1.2, 0.55), fontsize=8,
                 arrowprops=dict(arrowstyle="->", color="#D55E00"))

    params = d.get("params", {})
    fig.suptitle(
        f"Validation-gate ablation (n={params.get('n','?')} samples/dataset, "
        f"{params.get('reps','?')} reps/condition) — the discreteness gate removes "
        f"continuum false-positives",
        fontsize=10.5, y=1.02)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(args.out / f"gate_ablation_figure.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out}/gate_ablation_figure.png (+ .svg)")


if __name__ == "__main__":
    main()
