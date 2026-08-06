#!/usr/bin/env python3
"""JOB 1d — canonical sigclust across the SEPARATION SWEEP (the decisive comparison).

Job 1 found the reimplemented cluster index at least as good as the gate on the ablation
grid; job 1b found the opposite across the separation sweep, where the cluster index never
certified at any δ including 3.0, while the gate certified 55%. Job 1b is the more
informative of the two, because the ablation grid saturates both statistics.

But job 1b used the reimplementation. If the canonical package behaves differently — and it
plausibly might, since Liu et al.'s soft-thresholded covariance is exactly the kind of
estimator that matters in the borderline regime — then the reversal is an artifact of the
reimplementation rather than a property of the statistic. **This job decides that.**

It runs the CRAN `sigclust` (icovest=2, full p=50 feature space) on the deposited separation
generator at each δ, with the deposited loop parameters (n=120, p=50, 20 reps, data rng
`default_rng(1000*n + int(sep*10) + rep)`), so the result is directly comparable to the
deposited gate curve:

    δ:                0.0   0.5   1.0   1.5   2.0   2.5   3.0
    gate (deposited)  0.00  0.00  0.00  0.00  0.00  0.15  0.55

Read alongside job 2b: if shrinking the gate's reference covariance raises its certify rate
(it does — 0.55 → 0.90 at δ=3.0), then the gate's own covariance estimator is costing it
sensitivity, and the fair comparison against a package that soft-thresholds its covariance
is the one run here.

Requires R with `sigclust`. Deterministic. Writes nothing to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, subprocess, tempfile
from _provenance_safe import safe_argv  # noqa: E402
import numpy as np
import pandas as pd

SEPS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
DEPOSITED_GATE = {0.0: 0.00, 0.5: 0.00, 1.0: 0.00, 1.5: 0.00,
                  2.0: 0.00, 2.5: 0.15, 3.0: 0.55}

R_SCRIPT = r"""
suppressMessages(library(sigclust))
args <- commandArgs(trailingOnly=TRUE)
csv <- args[1]; out <- args[2]; nsim <- as.integer(args[3]); seed <- as.integer(args[4])
nrep <- if (length(args) >= 5) as.integer(args[5]) else 1
x <- as.matrix(read.csv(csv, header=FALSE))
set.seed(seed)
r <- sigclust(x, nsim=nsim, nrep=nrep, labflag=0, icovest=2)
cat(sprintf("%.6f %.6f\n", r@pval, r@xcindex), file=out)
"""


def make(sep, n, p, rng):
    """Deposited generator: two equal Gaussian blobs separated by `sep` SD on axis 0."""
    n1 = n // 2
    a = rng.normal(0, 1, size=(n1, p))
    b = rng.normal(0, 1, size=(n - n1, p))
    b[:, 0] += sep
    return np.vstack([a, b])



def _require_sigclust(rscript_body):
    """Fail closed unless R AND the CRAN `sigclust` package actually work.

    `which Rscript` is not enough: an Rscript that exists but errors (wrong R, missing
    package, broken library path) passed that check, and every per-dataset call then failed
    into a `continue`. The job exited 0, printed its success banner, and wrote a
    deposited-shaped JSONL with ZERO certifications -- which is indistinguishable from this
    paper's actual Result 5 headline. A reviewer without a working sigclust would have
    "reproduced" our result by running nothing. Same class as the diptest fail-open, and
    worse, because the empty answer here looks exactly like the true one.
    """
    import subprocess as _sp, tempfile as _tf, os as _os, shutil as _sh
    if _sh.which("Rscript") is None:
        raise SystemExit("Rscript not found. This job needs R >= 4 with the CRAN `sigclust` "
                         "package. Refusing to run: a partial SigClust sweep is not a result.")
    with _tf.TemporaryDirectory() as _td:
        p = _os.path.join(_td, "probe.R")
        open(p, "w").write('suppressMessages(library(sigclust)); cat("ok\\n")\n')
        r = _sp.run(["Rscript", "--vanilla", p], capture_output=True, text=True)
    if r.returncode != 0 or "ok" not in (r.stdout or ""):
        raise SystemExit(
            "R is present but the CRAN `sigclust` package could not be loaded "
            "(install.packages(\"sigclust\")). Refusing to run: every SigClust call would "
            "fail into a skip and this job would emit ZERO certifications, which is "
            "indistinguishable from the paper's real Result 5 finding.\n"
            "  Rscript said: " + (r.stderr or "").strip()[:300])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--nsim", type=int, default=100)
    ap.add_argument("--nrep", type=int, default=1,
                    help="k-means restarts inside sigclust (CRAN default 1; use 100 for a converged comparison)")
    ap.add_argument("--seps", default=",".join(str(s) for s in SEPS))
    a = ap.parse_args()
    _require_sigclust(R_SCRIPT)
    seps = [float(x) for x in a.seps.split(",")]

    with tempfile.TemporaryDirectory() as td:
        rs = os.path.join(td, "sc.R"); open(rs, "w").write(R_SCRIPT)
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as fh:
            fh.write(json.dumps(dict(
                record="provenance", script=os.path.basename(__file__), argv=safe_argv(a),
                design="canonical CRAN sigclust, icovest=2, full feature space",
                deposited_gate=DEPOSITED_GATE,
                note="deposited loop params so the curve is directly comparable")) + "\n")
            fh.flush()
            for sep in seps:
                ps = []
                for rep in range(a.reps):
                    rng = np.random.default_rng(1000 * a.n + int(sep * 10) + rep)
                    X = make(sep, a.n, a.p, rng)
                    csv = os.path.join(td, "x.csv")
                    pd.DataFrame(X).to_csv(csv, header=False, index=False)
                    txt = os.path.join(td, "o.txt")
                    r = subprocess.run(["Rscript", "--vanilla", rs, csv, txt,
                                        str(a.nsim), str(42), str(a.nrep)],
                                       capture_output=True, text=True)
                    if r.returncode != 0 or not os.path.exists(txt):
                        raise SystemExit(
                            f"R/sigclust failed at sep={sep} rep{rep} -- refusing to continue.\n"
                            f"  Rscript said: {r.stderr.strip()[:300]}")
                    pval, cindex = [float(v) for v in open(txt).read().split()]
                    ps.append(pval)
                    fh.write(json.dumps(dict(record="run", sep=sep, rep=rep,
                                             p=pval, cindex=cindex)) + "\n")
                    fh.flush()
                    os.remove(txt)
                ps = np.asarray(ps)
                rate = float((ps < 0.05).mean()) if len(ps) else float("nan")
                fh.write(json.dumps(dict(record="sep_summary", sep=sep, reps=int(len(ps)),
                                         certify_rate=rate,
                                         median_p=float(np.median(ps)) if len(ps) else None,
                                         deposited_gate=DEPOSITED_GATE.get(sep))) + "\n")
                fh.flush()
                print(f"  sep={sep:<4} sigclust certify={rate:.2f} "
                      f"median_p={np.median(ps):.4f}   [deposited gate="
                      f"{DEPOSITED_GATE.get(sep)}]", flush=True)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
