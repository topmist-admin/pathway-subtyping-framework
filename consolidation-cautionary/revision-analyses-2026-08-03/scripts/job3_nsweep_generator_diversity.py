#!/usr/bin/env python3
"""JOB 3 — n-sweep and negative-generator diversification (Paulus, "on the synthetic expansion").

He rejected a flat tenfold expansion: it buys precision, not unbiasedness. The two testable
negatives are not a random sample of the 30 — they are the ones the silhouette check happened
to pass — so more replicates at the same n reproduce that conditioning at higher n. And 600
draws from two generator families are pseudoreplicates with respect to a claim about continua
in general.

Two things are tested here instead.

(a) N-SWEEP. n in {60, 120, 300, 600} with replicates fixed. The prediction, if abstention is
    a power phenomenon rather than a property, is that the abstention rate FALLS as n rises
    while the false-certify rate stays at 0. That directly addresses the conditioning concern
    and feeds the testability x ground-truth two-by-two (Paulus P4).

(b) GENERATOR DIVERSITY. Five new negative families beyond the deposited two, so "continua in
    general" has some claim to generality:
      corr_weak / corr_strong   - correlated Gaussian at two correlation strengths
      nongaussian_gradient      - heavy-tailed (t_3) latent along one axis
      heteroscedastic_gradient  - variance grows along the gradient
      curved_manifold           - samples on a curved 1-D manifold (arc), no discrete groups
      three_component_gradient  - three overlapping shoulders, still continuous

All are NEGATIVES: none contains discrete groups, so any certification is a false certification.
`discrete_k2` is carried as the positive control so the sweep also shows sensitivity.

Newline-delimited JSON; nothing written to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from _provenance_safe import safe_argv  # noqa: E402
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from pathway_subtyping.discreteness import DiscretenessGateA  # noqa: E402



import zlib  # noqa: E402

# NOTE: was `hash(<str>) % 997`. Python salts str hashes per process (PYTHONHASHSEED),
# so that seeding produced DIFFERENT datasets on every invocation and nothing derived
# from it could be regenerated. zlib.crc32 is stable across processes and platforms.
def _stable_key(name: str) -> int:
    """Process-stable substitute for hash(name) % 997."""
    return zlib.crc32(name.encode()) % 997


def _df(X):
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(X.shape[1])])


def gen(kind: str, n: int, p: int, seed: int):
    rng = np.random.default_rng(seed)
    ns = max(3, p // 10)

    if kind == "discrete_k2":                      # POSITIVE control
        n1 = n // 2
        a = rng.normal(0, 1, size=(n1, p)); b = rng.normal(0, 1, size=(n - n1, p))
        b[:, 0] += 3.0
        return _df(np.vstack([a, b])), 2

    if kind in ("corr_weak", "corr_strong"):
        rho = 0.3 if kind == "corr_weak" else 0.85
        cov = np.full((p, p), rho) + np.eye(p) * (1 - rho)
        return _df(rng.multivariate_normal(np.zeros(p), cov, size=n)), 2

    if kind == "nongaussian_gradient":             # heavy-tailed latent
        t = rng.standard_t(df=3, size=n)
        d = np.zeros(p); d[:ns] = rng.normal(0, 1, ns); d /= np.linalg.norm(d)
        return _df(np.outer(t, d) + rng.normal(0, 1, size=(n, p))), 2

    if kind == "heteroscedastic_gradient":         # noise grows along the gradient
        t = rng.normal(0, 1, size=n)
        d = np.zeros(p); d[:ns] = rng.normal(0, 1, ns); d /= np.linalg.norm(d)
        scale = 0.5 + 1.5 * (t - t.min()) / (np.ptp(t) + 1e-9)   # np.ptp: ndarray.ptp() removed in NumPy 2
        return _df(np.outer(t, d) + rng.normal(0, 1, size=(n, p)) * scale[:, None]), 2

    if kind == "curved_manifold":                  # arc in the first two axes
        t = rng.uniform(-1.5, 1.5, size=n)
        X = rng.normal(0, 0.4, size=(n, p))
        X[:, 0] += 3.0 * np.sin(t); X[:, 1] += 3.0 * np.cos(t)
        return _df(X), 2

    if kind == "three_component_gradient":         # three overlapping shoulders, continuous
        t = np.concatenate([rng.normal(m, 1.0, size=n // 3) for m in (-1.2, 0.0, 1.2)])
        t = np.concatenate([t, rng.normal(0.0, 1.0, size=n - len(t))])
        d = np.zeros(p); d[:ns] = rng.normal(0, 1, ns); d /= np.linalg.norm(d)
        return _df(np.outer(t, d) + rng.normal(0, 1, size=(n, p))), 2

    raise ValueError(kind)


NEGATIVES = ["corr_weak", "corr_strong", "nongaussian_gradient",
             "heteroscedastic_gradient", "curved_manifold", "three_component_gradient"]
ALL = ["discrete_k2"] + NEGATIVES


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ns", default="60,120,300,600")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--n-ref", type=int, default=100)
    ap.add_argument("--n-bootstrap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--kinds", default=",".join(ALL))
    a = ap.parse_args()
    ns = [int(x) for x in a.ns.split(",")]
    kinds = a.kinds.split(",")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t0 = time.time()
    with open(a.out, "w") as fh:
        fh.write(json.dumps(dict(record="provenance", script=os.path.basename(__file__),
                                 argv=safe_argv(a),
                                 note="all non-discrete kinds are NEGATIVES; "
                                      "any certify is a false certification")) + "\n")
        fh.flush()
        for n in ns:
            for kind in kinds:
                tally = dict(certify=0, reject=0, abstain=0)
                for rep in range(a.reps):
                    seed = a.seed + 1000 * rep + _stable_key(kind) + n
                    scores, k = gen(kind, n, a.p, seed)
                    g = DiscretenessGateA(seed=seed, n_ref=a.n_ref,
                                          n_bootstrap=a.n_bootstrap)
                    r = g.run(tumor="synthetic", pathway_scores=scores,
                              n_clusters=k, gmm_seed=seed)
                    if r.testable and r.passed:
                        v = "certify"
                    elif str(r.verdict).startswith("not-testable"):
                        v = "abstain"
                    else:
                        v = "reject"
                    tally[v] += 1
                    fh.write(json.dumps(dict(record="run", n=n, kind=kind, rep=rep,
                                             verdict=v, raw_verdict=str(r.verdict))) + "\n")
                    fh.flush()
                fh.write(json.dumps(dict(record="cell_summary", n=n, kind=kind,
                                         reps=a.reps, **tally)) + "\n")
                fh.flush()
                print(f"  n={n:<4} {kind:26s} certify={tally['certify']:>2} "
                      f"reject={tally['reject']:>2} abstain={tally['abstain']:>2}", flush=True)
        fh.write(json.dumps(dict(record="done", elapsed_sec=round(time.time()-t0, 1))) + "\n")
    print(f"JOB3 done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
