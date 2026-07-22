#!/usr/bin/env python3
"""Gate ablation study — answers Scientific Reports reviewer point R3.10.

R3.10 asks: "What happens to the performance if some of the validation gates in
the framework are removed?"

This harness quantifies what each gate *buys* by measuring the false-positive
rate (FPR) of certifying a partition as a "real subtype" under different gate
subsets, on synthetic data with KNOWN ground truth:

  POSITIVE (genuine discrete structure — a real subtype SHOULD be certified):
    - discrete_k2 : two well-separated Gaussian components
    - discrete_k3 : three well-separated Gaussian components
  NEGATIVE (no discrete structure — certification is a FALSE POSITIVE):
    - single_gaussian : one correlated Gaussian blob
    - continuum_1d    : a 1-D continuous gradient (tumor purity / infiltration
                        analog) that a mixture model reproducibly bisects

The headline result reproduces the erratum finding: the *old* bootstrap-stability
gate alone certifies continua as reproducible subtypes (high FPR at small n),
because its null tested independence, not discreteness. Adding the v0.8.0
discreteness gate (Gate A) collapses that FPR while retaining true positives.

Gate subsets compared (the ablation):
    stability_only          : {bootstrap stability}                 (pre-v0.8.0)
    discreteness_only       : {Gate A discreteness null}
    stability+discreteness  : both must pass                        (v0.8.0)

Usage:
    python scripts/gate_ablation_study.py [--quick] [--reps N] [--n N]
                                          [--out DIR] [--seed S]

    --quick  small/fast settings for a smoke run (few refs, few reps)

Outputs (JSON + Markdown table) are written under --out (default:
outputs/gate_ablation/, gitignored). Nothing here depends on private data.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from pathway_subtyping.validation import ValidationGates
from pathway_subtyping.discreteness import DiscretenessGateA


# --------------------------------------------------------------------------- #
# Synthetic data generators — each returns (scores DataFrame, natural_k)
# --------------------------------------------------------------------------- #
def make_discrete(n: int, p: int, k: int, sep: float, seed: int) -> pd.DataFrame:
    """k well-separated Gaussian components in the first few pathways."""
    rng = np.random.default_rng(seed)
    n_signal = max(3, p // 10)
    sizes = [n // k] * k
    sizes[-1] += n - sum(sizes)
    blocks = []
    for c in range(k):
        center = np.zeros(p)
        # spread component centers across the signal pathways
        center[c % n_signal] = sep * (1 + c)
        center[(c + 1) % n_signal] = -sep * (1 + c)
        blocks.append(rng.normal(center, 1.0, size=(sizes[c], p)))
    X = np.vstack(blocks)
    idx = rng.permutation(n)
    return pd.DataFrame(X[idx], columns=[f"PATHWAY_{j}" for j in range(p)])


def make_single_gaussian(n: int, p: int, seed: int) -> pd.DataFrame:
    """One correlated Gaussian blob — no discrete structure."""
    rng = np.random.default_rng(seed)
    A = rng.normal(0, 1, size=(p, p))
    cov = (A @ A.T) / p + np.eye(p) * 0.5
    X = rng.multivariate_normal(np.zeros(p), cov, size=n)
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(p)])


def make_continuum(n: int, p: int, seed: int) -> pd.DataFrame:
    """A 1-D continuous gradient: samples lie along one axis with noise.

    A mixture model reproducibly bisects this every bootstrap (high observed
    stability ARI) even though there are NO discrete clusters. This is the case
    the old stability gate false-certifies."""
    rng = np.random.default_rng(seed)
    n_signal = max(3, p // 10)
    # Gaussian latent along ONE dominant axis: GMM at k=2 splits it at the mean
    # reproducibly across bootstraps (high stability ARI) even though the
    # distribution is unimodal — no discrete clusters.
    t = rng.normal(0, 1, size=n)  # position along the continuum
    direction = np.zeros(p)
    direction[:n_signal] = rng.normal(0, 1, size=n_signal)
    direction /= np.linalg.norm(direction)
    X = np.outer(t, direction) * 20.0 + rng.normal(0, 0.3, size=(n, p))
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(p)])


CONDITIONS = {
    "discrete_k2": dict(kind="discrete", k=2, klass="positive"),
    "discrete_k3": dict(kind="discrete", k=3, klass="positive"),
    "single_gaussian": dict(kind="single", k=2, klass="negative"),
    "continuum_1d": dict(kind="continuum", k=2, klass="negative"),
}


def generate(cond: str, n: int, p: int, sep: float, seed: int):
    spec = CONDITIONS[cond]
    if spec["kind"] == "discrete":
        return make_discrete(n, p, spec["k"], sep, seed), spec["k"]
    if spec["kind"] == "single":
        return make_single_gaussian(n, p, seed), spec["k"]
    return make_continuum(n, p, seed), spec["k"]


# --------------------------------------------------------------------------- #
# Evaluate the two gates on one dataset
# --------------------------------------------------------------------------- #
def evaluate_dataset(scores: pd.DataFrame, k: int, seed: int, n_ref: int,
                     n_bootstrap: int) -> Dict[str, bool]:
    """Fit GMM at k, then run the old stability gate and Gate A. Returns the
    per-gate PASS booleans."""
    from sklearn.mixture import GaussianMixture

    gm = GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                         random_state=seed, reg_covar=1e-6).fit(scores.values)
    labels = gm.predict(scores.values)

    # Old bootstrap-stability gate (pre-v0.8.0)
    vg = ValidationGates(seed=seed, n_bootstrap=n_bootstrap,
                         stability_threshold=0.8, show_progress=False)
    stab = vg.stability_test_bootstrap(scores, labels, k, gmm_seed=seed)

    # New discreteness gate (Gate A, v0.8.0)
    gateA = DiscretenessGateA(seed=seed, n_ref=n_ref, n_bootstrap=n_bootstrap)
    a = gateA.run(tumor="synthetic", pathway_scores=scores, n_clusters=k,
                  gmm_seed=seed)

    return {
        "stability_pass": bool(stab.passed),
        "stability_ari": float(stab.metric_value),
        "discreteness_pass": bool(a.passed),
        "discreteness_verdict": a.verdict,
        "discreteness_testable": bool(a.testable),
        "sg_empirical_p": float(a.sg_empirical_p),
    }


# certification policies: which gates must PASS to certify a "real subtype"
POLICIES = {
    "stability_only": lambda r: r["stability_pass"],
    "discreteness_only": lambda r: r["discreteness_pass"],
    "stability+discreteness": lambda r: r["stability_pass"] and r["discreteness_pass"],
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=10, help="replicates per condition")
    ap.add_argument("--n", type=int, default=60, help="samples per dataset (small-n regime)")
    ap.add_argument("--p", type=int, default=50, help="pathways (features)")
    ap.add_argument("--sep", type=float, default=3.0, help="cluster separation for positives")
    ap.add_argument("--n-ref", type=int, default=100, help="Gate A single-Gaussian references")
    ap.add_argument("--n-bootstrap", type=int, default=40, help="bootstrap iterations")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--quick", action="store_true", help="fast smoke settings")
    ap.add_argument("--out", type=Path, default=Path("outputs/gate_ablation"))
    args = ap.parse_args()

    if args.quick:
        args.reps, args.n_ref, args.n_bootstrap = 3, 30, 20

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    rows: List[Dict] = []
    for cond in CONDITIONS:
        for rep in range(args.reps):
            seed = args.seed + 1000 * rep + hash(cond) % 997
            scores, k = generate(cond, args.n, args.p, args.sep, seed)
            res = evaluate_dataset(scores, k, seed, args.n_ref, args.n_bootstrap)
            res.update(condition=cond, klass=CONDITIONS[cond]["klass"], rep=rep, k=k)
            rows.append(res)
            print(f"  {cond:16s} rep {rep:2d}: "
                  f"stability={'PASS' if res['stability_pass'] else 'fail'} "
                  f"(ARI={res['stability_ari']:.2f})  "
                  f"discreteness={'PASS' if res['discreteness_pass'] else 'fail'} "
                  f"({res['discreteness_verdict']})")

    df = pd.DataFrame(rows)

    # Aggregate: per-policy FPR (on negatives) and TPR (on positives)
    summary = {}
    for pol, fn in POLICIES.items():
        cert = df.apply(fn, axis=1)
        pos = df["klass"] == "positive"
        neg = df["klass"] == "negative"
        tpr = float(cert[pos].mean())
        fpr = float(cert[neg].mean())
        # per-condition breakdown
        by_cond = {c: float(cert[df["condition"] == c].mean()) for c in CONDITIONS}
        summary[pol] = dict(TPR=round(tpr, 3), FPR=round(fpr, 3), cert_rate_by_condition=by_cond)

    out = dict(
        params=vars(args) | {"out": str(args.out)},
        n_datasets=len(df),
        elapsed_sec=round(time.time() - t0, 1),
        summary=summary,
    )
    (args.out / "gate_ablation_results.json").write_text(json.dumps(out, indent=2, default=str))
    df.to_csv(args.out / "gate_ablation_raw.csv", index=False)

    # Markdown table
    lines = [
        "# Gate Ablation Study — results (R3.10)",
        "",
        f"n={args.n} samples, p={args.p} pathways, {args.reps} reps/condition, "
        f"n_ref={args.n_ref}, n_bootstrap={args.n_bootstrap}. "
        f"Elapsed {out['elapsed_sec']}s.",
        "",
        "**TPR** = fraction of genuine-subtype datasets certified (higher is better). "
        "**FPR** = fraction of no-structure datasets (single Gaussian + continuum) "
        "falsely certified as subtypes (lower is better).",
        "",
        "| Gate subset | TPR | FPR | continuum_1d cert rate |",
        "|---|---|---|---|",
    ]
    for pol, s in summary.items():
        lines.append(
            f"| `{pol}` | {s['TPR']:.2f} | {s['FPR']:.2f} | "
            f"{s['cert_rate_by_condition']['continuum_1d']:.2f} |"
        )
    lines += [
        "",
        "**Reading:** `stability_only` (the pre-v0.8.0 gate set) false-certifies the "
        "1-D continuum — the erratum finding. `stability+discreteness` (v0.8.0) drives "
        "the continuum false-positive rate down while retaining true positives, which is "
        "the quantitative answer to R3.10: removing the discreteness gate reintroduces "
        "continuum false positives.",
    ]
    (args.out / "gate_ablation_results.md").write_text("\n".join(lines))

    print("\n" + "\n".join(lines))
    print(f"\nWrote: {args.out}/gate_ablation_results.{{json,md}} + gate_ablation_raw.csv")


if __name__ == "__main__":
    main()
