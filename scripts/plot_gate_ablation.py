#!/usr/bin/env python3
"""Publication figure for the gate ablation study (reviewer R3.10).

Reads outputs/gate_ablation/gate_ablation_results.json + gate_ablation_raw.csv
and renders a two-panel figure:
  A. TPR vs FPR by gate subset — the discreteness gate removes the false-positive
     rate the stability-only gate incurs, at negligible true-positive cost.
  B. Per-condition certification rate — the false positives are concentrated on
     the 1-D continuum (the erratum failure mode); the discreteness gate zeroes it.

Usage: python scripts/plot_gate_ablation.py [--out DIR]
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
    args = ap.parse_args()
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
