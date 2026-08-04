#!/usr/bin/env python3
"""JOB 2 — is the gate's conservatism a designed property or an estimation artifact?
(Paulus Point 6.)

The gate builds its single-Gaussian reference from `var = Z.var(0)` — per-PC variances
estimated from the OBSERVED data. Under the alternative, those data contain the
between-component separation, so the leading PC's variance is inflated by exactly the
structure the null is meant to exclude. The reference cloud is then elongated along the
separation axis, which makes it split more reproducibly, which raises the reference ARI,
which raises the bar to certify — and it does so MORE as true separation grows. That is
the shape of the reported "resolves only at delta >= 2.5" finding.

DESIGN. Re-run the deposited separation sweep unchanged except for one knob: shrink the
per-PC reference variances toward their mean by lambda.

    var_shrunk = (1 - lam) * var + lam * mean(var)

    lam = 0.0  -> IDENTICAL to the deposited gate. This is the control: if lam=0 does not
                  reproduce the deposited curve, the harness is wrong and nothing else here
                  can be believed.
    lam = 1.0  -> fully isotropic reference at the average variance (no inflation possible)

If the delta floor moves materially as lam rises, the conservatism is an estimation
artifact and the Discussion's characterisation of the screen has to change. If it does not
move, the conservatism survives as a real property and the claim stands.

Reuses the deposited generator and the gate's own reference machinery, so the ONLY thing
varying is the covariance estimate. Writes newline-delimited JSON; nothing touches the
deposited tree.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))

from pathway_subtyping.discreteness.gate_a_discreteness_null import (  # noqa: E402
    reduce_scores, reduced_dim, _single_gaussian_once, _stability_of,
)

SEPS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]


def make(sep: float, n: int, p: int, rng) -> np.ndarray:
    """Two equal Gaussian blobs separated by `sep` SD along axis 0 (deposited generator)."""
    n1 = n // 2
    a = rng.normal(0, 1, size=(n1, p))
    b = rng.normal(0, 1, size=(n - n1, p))
    b[:, 0] += sep
    return np.vstack([a, b])


def cell(sep, n, p, seed, n_ref, n_bootstrap, lams, k=2):
    rng = np.random.default_rng(seed)
    X = make(sep, n, p, rng)
    Z, _ = reduce_scores(X, reduced_dim(n), seed)
    mean, var = Z.mean(0), Z.var(0)
    obs = _stability_of(Z, k, seed, n_bootstrap, seed)

    ss = np.random.SeedSequence(seed)
    seeds = [int(np.random.default_rng(s).integers(1, 2**31)) for s in ss.spawn(n_ref)]

    out = dict(sep=sep, n=n, observed=float(obs),
               var_ratio_pc1_to_mean=float(var[0] / var.mean()))
    for lam in lams:
        v = (1.0 - lam) * var + lam * var.mean()
        ref = np.asarray([
            _single_gaussian_once(mean, v, Z.shape[0], k, sd, n_bootstrap, seed + i + 1)
            for i, sd in enumerate(seeds)
        ])
        ref = ref[np.isfinite(ref)]
        p95 = float(np.percentile(ref, 95))
        pval = float((np.sum(ref >= obs) + 1) / (len(ref) + 1))
        out[f"lam{lam}"] = dict(ref_mean=float(ref.mean()), ref_p95=p95, p=pval,
                                certify=bool(obs > p95 and pval < 0.05))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--n-ref", type=int, default=100)
    ap.add_argument("--n-bootstrap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--lams", default="0.0,0.25,0.5,0.75,1.0")
    ap.add_argument("--seps", default=",".join(str(s) for s in SEPS))
    a = ap.parse_args()
    lams = [float(x) for x in a.lams.split(",")]
    seps = [float(x) for x in a.seps.split(",")]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t0 = time.time()
    with open(a.out, "w") as fh:
        fh.write(json.dumps(dict(record="provenance", script=os.path.basename(__file__),
                                 argv=vars(a),
                                 control="lam=0.0 must reproduce the deposited sweep")) + "\n")
        fh.flush()
        for sep in seps:
            cert = {lam: 0 for lam in lams}
            for rep in range(a.reps):
                seed = a.seed + 1000 * rep + int(sep * 10)
                r = cell(sep, a.n, a.p, seed, a.n_ref, a.n_bootstrap, lams)
                r.update(record="cell", rep=rep)
                fh.write(json.dumps(r) + "\n"); fh.flush()
                for lam in lams:
                    cert[lam] += int(r[f"lam{lam}"]["certify"])
            summary = dict(record="sep_summary", sep=sep, reps=a.reps,
                           certify_rate={str(l): cert[l] / a.reps for l in lams})
            fh.write(json.dumps(summary) + "\n"); fh.flush()
            print(f"  sep={sep:<4} " + "  ".join(
                f"lam{l}={cert[l]/a.reps:.2f}" for l in lams), flush=True)
        fh.write(json.dumps(dict(record="done", elapsed_sec=round(time.time()-t0, 1))) + "\n")
    print(f"JOB2 done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
