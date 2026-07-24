#!/usr/bin/env python3
"""Separation sweep — the gate's operating characteristic (reviewer: sep was fixed
at 3.0, positives saturate the p-floor, no evidence of resolution).

The ablation used a single separation (well-separated clusters), so every positive
saturated the p-value floor and the gate's behaviour between "obvious structure" and
"no structure" was never probed. This sweep varies the between-component separation
delta from 0 (single Gaussian, a true negative) up to a clean 2-component mixture,
at several sample sizes, and records the gate's certification rate and its
single-Gaussian p-value at each step. A calibrated gate should transition from
reject/abstain to certify somewhere in the middle, not flip at the extremes only.

Deterministic (seed 42). No network. Slow (bootstrap per cell).

Usage: python separation_sweep.py [--reps N] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
try:
    from pathway_subtyping.discreteness import DiscretenessGateA
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../..")))
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../../src")))
    from pathway_subtyping.discreteness import DiscretenessGateA

import pandas as pd
from sklearn.mixture import GaussianMixture

SEPS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
NS = [120]        # single n keeps the sweep tractable; the curve is the point
P = 50            # features
N_REF = 100       # floor 1/101; modest so the sweep is tractable


def make(sep: float, n: int, p: int, rng) -> np.ndarray:
    """Two equal Gaussian blobs separated by `sep` SD along axis 0. sep=0 => one blob."""
    n1 = n // 2
    a = rng.normal(0, 1, size=(n1, p))
    b = rng.normal(0, 1, size=(n - n1, p))
    b[:, 0] += sep
    return np.vstack([a, b])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    floor = 1.0 / (N_REF + 1)

    rows = []
    for n in NS:
        for sep in SEPS:
            certs, ps, testable = 0, [], 0
            for rep in range(args.reps):
                rng = np.random.default_rng(1000 * n + int(sep * 10) + rep)
                X = make(sep, n, P, rng)
                Xdf = pd.DataFrame(X, columns=[f"f{i}" for i in range(P)])
                lab = GaussianMixture(2, covariance_type="full", n_init=1,
                                      random_state=42, reg_covar=1e-6).fit(X).predict(X)
                _ = lab
                res = DiscretenessGateA(seed=42, n_ref=N_REF, n_bootstrap=20).run(
                    "sweep", Xdf, 2, gmm_seed=42)
                if res.testable:
                    testable += 1
                    if res.passed:
                        certs += 1
                ps.append(float(res.sg_empirical_p))
            rows.append({
                "n": n, "sep": sep, "reps": args.reps,
                "certify_rate": round(certs / args.reps, 3),
                "testable_rate": round(testable / args.reps, 3),
                "median_sg_p": round(float(np.median(ps)), 4),
                "min_sg_p": round(float(np.min(ps)), 4),
                "p_floor": round(floor, 4),
                "median_p_above_floor": bool(np.median(ps) > floor + 1e-6),
            })
            print(f"n={n} sep={sep:.1f}: certify {rows[-1]['certify_rate']:.2f} "
                  f"testable {rows[-1]['testable_rate']:.2f} "
                  f"median_p {rows[-1]['median_sg_p']:.4f} "
                  f"(floor {floor:.4f})")

    out = {"design": {"seps": SEPS, "ns": NS, "p": P, "n_ref": N_REF,
                      "reps": args.reps, "p_floor": round(floor, 4)},
           "sweep": rows,
           "reads": ("A resolving gate should show certify_rate rising with sep and "
                     "median p descending through the floor region, not a flat "
                     "floor. Flat-at-floor across all sep would confirm the "
                     "reviewer's concern that the test has no resolution.")}
    path = os.path.join(args.out, "separation_sweep.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
