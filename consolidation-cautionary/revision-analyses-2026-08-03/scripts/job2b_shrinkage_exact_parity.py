#!/usr/bin/env python3
"""JOB 2b — covariance shrinkage, run through the DEPOSITED code path (Paulus Point 6).

WHY THIS SUPERSEDES JOB 2. Job 2 reimplemented the gate's internals (reduce_scores ->
_stability_of -> _single_gaussian_once) and varied `var` inside that reimplementation.
Its lam=0 control returned certify 0.30 at delta=3.0 where the deposited sweep reports
0.55, so by the standard job 2 set for itself the harness was not validated and its
numbers are not quotable. Diagnosed divergences: n_bootstrap 40 vs the deposited 20;
per-rep varying gate/gmm seeds vs the deposited fixed seed=42; the testability
precondition omitted; and — the substantive one — the deposited OBSERVED statistic comes
from `stability_test_bootstrap(Zdf, obs_labels, ...)` against labels from a GMM with
n_init=10, not from `_stability_of`, which fits its own.

DESIGN. Do not reimplement anything. Monkeypatch the module-level
`_single_gaussian_once` with a wrapper that shrinks the per-PC reference variances

    var_shrunk = (1 - lam) * var + lam * mean(var)

and delegates to the original. The gate calls that function by module global, so every
other line executed is deposited code. Then:

    lam = 0.0  -> the wrapper is the identity. This is not "an approximation of" the
                  deposited gate, it IS the deposited gate, and it must reproduce
                  certify 0.15 at delta=2.5 and 0.55 at delta=3.0 exactly. Parity is
                  asserted per-sep and recorded; if it fails, nothing else here holds.
    lam = 1.0  -> fully isotropic reference at the average variance, so the leading PC's
                  variance can no longer be inflated by the between-component separation.

Loop parameters are copied verbatim from the deposited
`consolidation-cautionary/cross-domain/gate_ablation/scripts/separation_sweep.py`:
data rng default_rng(1000*n + int(sep*10) + rep); DiscretenessGateA(seed=42, n_ref=100,
n_bootstrap=20); run("sweep", Xdf, 2, gmm_seed=42); certify iff res.testable and
res.passed.

Newline-delimited JSON; nothing is written to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))

from pathway_subtyping.discreteness import DiscretenessGateA          # noqa: E402
from pathway_subtyping.discreteness import gate_a_discreteness_null as G  # noqa: E402

_ORIG_SG_ONCE = G._single_gaussian_once

# Deposited constants (separation_sweep.py)
SEPS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
P, N_REF, N_BOOT, N = 50, 100, 20, 120
DEPOSITED = {0.0: 0.00, 0.5: 0.00, 1.0: 0.00, 1.5: 0.00,
             2.0: 0.00, 2.5: 0.15, 3.0: 0.55}


def make(sep: float, n: int, p: int, rng) -> np.ndarray:
    """Deposited generator: two equal Gaussian blobs separated by `sep` SD on axis 0."""
    n1 = n // 2
    a = rng.normal(0, 1, size=(n1, p))
    b = rng.normal(0, 1, size=(n - n1, p))
    b[:, 0] += sep
    return np.vstack([a, b])


def patch(lam: float):
    """Shrink ONLY the reference variance; delegate everything else to deposited code."""
    if lam == 0.0:
        G._single_gaussian_once = _ORIG_SG_ONCE          # exact identity, not a wrapper
        return
    def shrunk(mean, var, n, n_clusters, seed, n_bootstrap, gmm_seed):
        var = np.asarray(var, dtype=float)
        return _ORIG_SG_ONCE(mean, (1.0 - lam) * var + lam * var.mean(),
                             n, n_clusters, seed, n_bootstrap, gmm_seed)
    G._single_gaussian_once = shrunk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--lams", default="0.0,0.5,1.0")
    ap.add_argument("--seps", default=",".join(str(s) for s in SEPS))
    a = ap.parse_args()
    lams = [float(x) for x in a.lams.split(",")]
    seps = [float(x) for x in a.seps.split(",")]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t0 = time.time()
    parity_ok = True
    with open(a.out, "w") as fh:
        fh.write(json.dumps(dict(
            record="provenance", script=os.path.basename(__file__), argv=vars(a),
            design="monkeypatch reference variance only; all other lines are deposited code",
            control="lam=0.0 is the identity and must reproduce the deposited certify rates",
            deposited=DEPOSITED, supersedes="job2_covariance_shrinkage_sweep.py")) + "\n")
        fh.flush()
        for sep in seps:
            rates = {}
            for lam in lams:
                patch(lam)
                certs = 0
                for rep in range(a.reps):
                    rng = np.random.default_rng(1000 * N + int(sep * 10) + rep)
                    X = make(sep, N, P, rng)
                    Xdf = pd.DataFrame(X, columns=[f"f{i}" for i in range(P)])
                    res = DiscretenessGateA(seed=42, n_ref=N_REF,
                                            n_bootstrap=N_BOOT).run("sweep", Xdf, 2,
                                                                    gmm_seed=42)
                    ok = bool(res.testable and res.passed)
                    certs += int(ok)
                    fh.write(json.dumps(dict(record="run", sep=sep, lam=lam, rep=rep,
                                             certify=ok, testable=bool(res.testable),
                                             p=float(res.sg_empirical_p),
                                             p95=float(res.sg_ref_p95))) + "\n")
                    fh.flush()
                rates[lam] = certs / a.reps
            dep = DEPOSITED.get(sep)
            match = (dep is None) or abs(rates.get(0.0, -1) - dep) < 1e-9
            if 0.0 in rates and not match:
                parity_ok = False
            fh.write(json.dumps(dict(record="sep_summary", sep=sep, reps=a.reps,
                                     certify_rate={str(l): rates[l] for l in lams},
                                     deposited=dep, parity_ok=match)) + "\n")
            fh.flush()
            print(f"  sep={sep:<4} " + "  ".join(f"lam{l}={rates[l]:.2f}" for l in lams)
                  + f"   [deposited lam0={dep}  parity={'OK' if match else 'FAIL'}]",
                  flush=True)
        fh.write(json.dumps(dict(record="done", parity_ok=parity_ok,
                                 elapsed_sec=round(time.time() - t0, 1))) + "\n")
    print(f"JOB2b done in {time.time()-t0:.0f}s  parity_ok={parity_ok} -> {a.out}")
    if not parity_ok:
        print("PARITY FAILED — lam=0 did not reproduce the deposited curve; "
              "do not quote these numbers.", flush=True)


if __name__ == "__main__":
    main()
