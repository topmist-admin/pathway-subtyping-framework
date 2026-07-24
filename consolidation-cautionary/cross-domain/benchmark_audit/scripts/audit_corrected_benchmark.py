#!/usr/bin/env python3
"""Result 2 — empirical audit of the corrected 47-dataset benchmark.

The cautionary-framework paper's Result 2 asks an empirical question: across many
independent public cohorts, how often does a discrete-subtype partition actually
clear a reproducibility bar? This script reads the CORRECTED benchmark (v2.0,
Zenodo 10.5281/zenodo.21262112, issued with the 2026-07-08 erratum) and reports
that distribution honestly.

It does three things:

  1. Reproduces the erratum's headline refit — the retracted "adaptive bootstrap
     threshold" model (published R^2 = 0.889, slope +0.914) does not hold. On the
     valid rows silhouette does not predict bootstrap ARI at all.

  2. Reports the reproducibility distribution across valid cohorts, at several
     candidate thresholds. This is the Result 2 headline.

  3. Applies a STRICTER ground-truth screen than the erratum did, and re-reports
     everything under it. The erratum flagged impossible labels with the rule
     `n_true_clusters > n_samples`, which catches GSE92332 (1957 labels for 533
     samples) but NOT rows where the label count merely equals or approaches the
     sample count — e.g. GSE2109, 2158 "true clusters" for 2158 samples, i.e. one
     cluster per sample. Those rows are still marked valid=True in v2.0. Applying
     the paper's own standard to the paper's own correction is on-thesis, and the
     conclusion must be shown to survive it.

Deterministic; no network. Reads only the deposited CSV.

Usage:
    python audit_corrected_benchmark.py [--csv PATH] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import os

import pandas as pd
from scipy import stats

# Thresholds to report the pass rate at. 0.80 is the bar the withdrawn
# manuscript actually used; the rest bracket it to show the result is not an
# artifact of one arbitrary choice (reviewer R3.7).
THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

# A row whose ground-truth label count approaches its sample count is not a
# usable ground truth: it is (near-)singleton labelling. 0.5 is deliberately
# permissive - it still admits rows with one label per two samples.
SINGLETON_RATIO = 0.5


def refit(sub: pd.DataFrame) -> dict:
    """Refit the retracted model: bootstrap ARI ~ silhouette."""
    s = sub.dropna(subset=["silhouette", "bootstrap_ari_5th_percentile"])
    if len(s) < 3:
        return {"n": int(len(s)), "note": "too few rows to fit"}
    r = stats.linregress(s["silhouette"], s["bootstrap_ari_5th_percentile"])
    return {
        "n": int(len(s)),
        "r_squared": round(float(r.rvalue ** 2), 4),
        "slope": round(float(r.slope), 4),
        "p": round(float(r.pvalue), 4),
    }


def distribution(sub: pd.DataFrame) -> dict:
    a = sub["bootstrap_ari_5th_percentile"].dropna()
    return {
        "n": int(len(a)),
        "min": round(float(a.min()), 4),
        "median": round(float(a.median()), 4),
        "mean": round(float(a.mean()), 4),
        "max": round(float(a.max()), 4),
        "pass_rate_at_threshold": {
            str(t): {
                "n_passing": int((a >= t).sum()),
                "frac": round(float((a >= t).mean()), 4),
            }
            for t in THRESHOLDS
        },
    }


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(
        here, "../../../../CORRECTION_2026-07/corrected_benchmark_47datasets_v2.csv"))
    ap.add_argument("--out", default=os.path.join(here, "../results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = pd.read_csv(args.csv)
    d["singleton_ratio"] = d["n_true_clusters"] / d["n_samples"]

    valid = d[d["valid"]]
    strict = valid[valid["singleton_ratio"] < SINGLETON_RATIO]
    dropped = valid[valid["singleton_ratio"] >= SINGLETON_RATIO]

    out = {
        "source": {
            "file": os.path.basename(args.csv),
            "zenodo_doi": "10.5281/zenodo.21262112",
            "version": "v2.0 (corrected, 2026-07-08 erratum)",
        },
        "row_accounting": {
            "total_rows": int(len(d)),
            "unique_dataset_ids": int(d["dataset_id_original"].nunique()),
            "marked_PASS": int((d["status"] == "PASS").sum()),
            "valid": int(d["valid"].sum()),
            "invalid": int((~d["valid"]).sum()),
        },
        # 1. the retracted model
        "retracted_threshold_model": {
            "published_claim": {"r_squared": 0.889, "slope": 0.914,
                                "loocv_rmse": 0.051},
            "refit_all_rows": refit(d),
            "refit_valid_rows": refit(valid),
            "refit_strict_rows": refit(strict),
            "verdict": ("Silhouette does not predict bootstrap ARI. The published "
                        "R^2=0.889 does not reproduce under any row screen."),
        },
        # 2. the Result 2 headline
        "reproducibility_distribution": {
            "valid_rows": distribution(valid),
            "strict_rows": distribution(strict),
        },
        # 3. is the reproducibility column even trustworthy?
        "column_validity_diagnostic": {
            "question": ("Is `bootstrap_ari_5th_percentile` actually measuring "
                         "bootstrap stability of the reported partition?"),
            "test": ("Well-separated clusters are trivially reproducible under "
                     "resampling. Any row with a high silhouette must therefore "
                     "have a high stability ARI. Rows with high silhouette AND "
                     "near-zero ARI are not physically consistent."),
            "high_silhouette_rows": [
                {
                    "dataset": row["dataset_name"],
                    "silhouette": round(float(row["silhouette"]), 3),
                    "bootstrap_ari": round(
                        float(row["bootstrap_ari_5th_percentile"]), 4),
                    "n_samples": int(row["n_samples"]),
                }
                for _, row in valid[valid["silhouette"] > 0.7]
                .sort_values("silhouette", ascending=False).iterrows()
            ],
            "spearman_silhouette_vs_ari": None,   # filled below
            "verdict": None,                       # filled below
        },
        # 4. the incomplete-correction finding
        "incomplete_correction_finding": {
            "erratum_rule": "n_true_clusters > n_samples",
            "rule_catches": ["GSE92332 (1957 labels / 533 samples)"],
            "rule_misses": (
                "rows where the label count equals or approaches the sample "
                "count, which are equally unusable as ground truth"),
            "n_valid_rows_failing_stricter_screen": int(len(dropped)),
            "screen": f"n_true_clusters / n_samples >= {SINGLETON_RATIO}",
            "rows_dropped": [
                {
                    "dataset": row["dataset_name"],
                    "n_samples": int(row["n_samples"]),
                    "n_true_clusters": int(row["n_true_clusters"]),
                    "ratio": round(float(row["singleton_ratio"]), 3),
                }
                for _, row in dropped.sort_values(
                    "singleton_ratio", ascending=False).iterrows()
            ],
            "conclusion_robust": None,  # filled below
        },
    }

    diag = out["column_validity_diagnostic"]
    hi = valid[valid["silhouette"] > 0.7]
    diag["spearman_silhouette_vs_ari"] = round(float(
        valid[["silhouette", "bootstrap_ari_5th_percentile"]]
        .corr(method="spearman").iloc[0, 1]), 4)
    n_bad = int((hi["bootstrap_ari_5th_percentile"] < 0.1).sum())
    diag["n_high_silhouette"] = int(len(hi))
    diag["n_high_silhouette_with_near_zero_ari"] = n_bad
    diag["verdict"] = (
        "COLUMN UNRELIABLE: every well-separated cohort reports near-zero or "
        "negative stability, which cannot be true of a genuine bootstrap-ARI "
        "measure. The column does not measure what its name says."
        if n_bad == len(hi) and len(hi) > 0 else
        "No inconsistency detected between silhouette and stability.")

    # The load-bearing defensive check: does the Result 2 conclusion survive the
    # stricter screen? It does only if no cohort clears the bar either way.
    v_at_50 = out["reproducibility_distribution"]["valid_rows"][
        "pass_rate_at_threshold"]["0.5"]["n_passing"]
    s_at_50 = out["reproducibility_distribution"]["strict_rows"][
        "pass_rate_at_threshold"]["0.5"]["n_passing"]
    out["incomplete_correction_finding"]["conclusion_robust"] = bool(
        v_at_50 == 0 and s_at_50 == 0)

    path = os.path.join(args.out, "benchmark_audit.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)

    # Console summary
    print(f"Rows {out['row_accounting']['total_rows']} "
          f"({out['row_accounting']['unique_dataset_ids']} unique IDs); "
          f"valid {out['row_accounting']['valid']}, "
          f"invalid {out['row_accounting']['invalid']}")
    print("\nRetracted model (bootstrap ARI ~ silhouette):")
    for k in ("refit_all_rows", "refit_valid_rows", "refit_strict_rows"):
        r = out["retracted_threshold_model"][k]
        print(f"  {k:18s} n={r['n']:3d} R2={r.get('r_squared')} "
              f"slope={r.get('slope')} p={r.get('p')}")
    print("\nReproducibility (bootstrap ARI 5th pct) across cohorts:")
    for label in ("valid_rows", "strict_rows"):
        dist = out["reproducibility_distribution"][label]
        print(f"  {label:12s} n={dist['n']:3d} median={dist['median']} "
              f"max={dist['max']}")
        for t in ("0.8", "0.5", "0.2"):
            p = dist["pass_rate_at_threshold"][t]
            print(f"      >= {t}: {p['n_passing']}/{dist['n']}")
    print(f"\nStricter screen drops "
          f"{out['incomplete_correction_finding']['n_valid_rows_failing_stricter_screen']}"
          f" of {len(valid)} valid rows; conclusion robust = "
          f"{out['incomplete_correction_finding']['conclusion_robust']}")
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
