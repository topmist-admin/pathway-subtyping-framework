#!/usr/bin/env python3
"""JOB 1 — SigClust comparator (Paulus Point 3).

The paper positions its increment over Liu et al. as "apply the single-Gaussian null to a
bootstrap-ARI statistic instead of the 2-means cluster index" — but SigClust is never run,
so the comparison has never been made.

DESIGN. Hold the null FIXED and vary only the statistic. Both arms draw references from the
same single Gaussian N(mean, diag(var)) in the same PCA space, with the same n_ref. The only
difference is what is measured on each draw:

    arm "ari"      : bootstrap ARI at fixed k          (the paper's statistic)
    arm "sigclust" : 2-means cluster index             (Liu et al.'s statistic)

That isolates the choice of statistic, which is exactly the claim under test. A comparison
that also changed the null would answer a different question.

Cluster index (Liu et al. 2008): CI = within-cluster SS from 2-means / total SS about the
mean. LOWER = more clustered, so the p-value counts references at or BELOW the observed.

Outputs newline-delimited JSON so a killed run still leaves usable partial results.
Nothing here writes to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "scripts"))

from sklearn.cluster import KMeans
from pathway_subtyping.discreteness.gate_a_discreteness_null import (  # noqa: E402
    reduce_scores, reduced_dim, _single_gaussian_once, _stability_of,
)
import gate_ablation_study as gas  # noqa: E402


def cluster_index(X: np.ndarray, seed: int) -> float:
    """Liu et al. 2-means cluster index. Lower = more clustered."""
    km = KMeans(n_clusters=2, n_init=10, random_state=seed).fit(X)
    within = float(((X - km.cluster_centers_[km.labels_]) ** 2).sum())
    total = float(((X - X.mean(0)) ** 2).sum())
    return within / total if total > 0 else float("nan")


def run_dataset(scores, k, seed, n_ref, n_bootstrap, do_ari):
    """Both statistics against the SAME single-Gaussian reference."""
    Z, _ = reduce_scores(scores.values, reduced_dim(scores.shape[0]), seed)
    mean, var = Z.mean(0), Z.var(0)
    n = Z.shape[0]
    ss = np.random.SeedSequence(seed)
    seeds = [int(np.random.default_rng(s).integers(1, 2**31)) for s in ss.spawn(n_ref)]

    out = {}

    # --- SigClust arm (cheap: one 2-means per draw) ---
    obs_ci = cluster_index(Z, seed)
    ref_ci = []
    for i, sd in enumerate(seeds):
        rng = np.random.default_rng(sd)
        draw = rng.normal(loc=mean, scale=np.sqrt(var), size=(n, len(mean)))
        ref_ci.append(cluster_index(draw, sd))
    ref_ci = np.asarray([c for c in ref_ci if np.isfinite(c)])
    # lower CI = more clustered -> p counts references at or below observed
    out["sigclust"] = dict(
        observed=obs_ci,
        ref_mean=float(ref_ci.mean()), ref_p05=float(np.percentile(ref_ci, 5)),
        p=float((np.sum(ref_ci <= obs_ci) + 1) / (len(ref_ci) + 1)),
        n_ref=int(len(ref_ci)),
    )

    # --- bootstrap-ARI arm (expensive: n_bootstrap GMM fits per draw) ---
    if do_ari:
        obs_ari = _stability_of(Z, k, seed, n_bootstrap, seed)
        ref_ari = [
            _single_gaussian_once(mean, var, n, k, sd, n_bootstrap, seed + i + 1)
            for i, sd in enumerate(seeds)
        ]
        ref_ari = np.asarray([a for a in ref_ari if np.isfinite(a)])
        out["ari"] = dict(
            observed=float(obs_ari),
            ref_mean=float(ref_ari.mean()), ref_p95=float(np.percentile(ref_ari, 95)),
            p=float((np.sum(ref_ari >= obs_ari) + 1) / (len(ref_ari) + 1)),
            n_ref=int(len(ref_ari)),
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=15)
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--sep", type=float, default=3.0)
    ap.add_argument("--n-ref", type=int, default=100)
    ap.add_argument("--n-bootstrap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-ari", action="store_true",
                    help="SigClust arm only (fast); the ARI arm is already deposited")
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t0 = time.time()
    with open(a.out, "w") as fh:
        fh.write(json.dumps(dict(record="provenance", script=os.path.basename(__file__),
                                 argv=vars(a), seed=a.seed,
                                 design="same single-Gaussian null; statistic varied")) + "\n")
        fh.flush()
        for cond in gas.CONDITIONS:
            for rep in range(a.reps):
                seed = a.seed + 1000 * rep + hash(cond) % 997
                scores, k = gas.generate(cond, a.n, a.p, a.sep, seed)
                res = run_dataset(scores, k, seed, a.n_ref, a.n_bootstrap,
                                  do_ari=not a.skip_ari)
                rec = dict(record="dataset", condition=cond,
                           klass=gas.CONDITIONS[cond]["klass"], rep=rep, k=k,
                           n=a.n, p=a.p, **res)
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print(f"  {cond:16s} rep{rep:<2d} sigclust p={res['sigclust']['p']:.3f}"
                      + (f"  ari p={res['ari']['p']:.3f}" if "ari" in res else ""), flush=True)
        fh.write(json.dumps(dict(record="done", elapsed_sec=round(time.time() - t0, 1))) + "\n")
    print(f"JOB1 done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
