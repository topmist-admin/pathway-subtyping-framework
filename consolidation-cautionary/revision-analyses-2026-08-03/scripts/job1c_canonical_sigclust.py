#!/usr/bin/env python3
"""JOB 1c — the CANONICAL sigclust R package (closes job 1's caveats 1 and 2).

Job 1 compared the gate against a Python reimplementation of Liu et al.'s 2-means cluster
index, deliberately run against the gate's own covariance and inside the gate's PCA-reduced
space so that the null was matched and only the statistic varied. That was the right design
for isolating the statistic, but it left two caveats that a reviewer would press:

  (1) it is not the `sigclust` package — Liu et al. use eigenvalue soft-thresholding for the
      covariance in the HDLSS regime, which the reimplementation did not;
  (2) it ran in a d = min(10, n/10) PCA space, not the full feature space where canonical
      SigClust operates, and the reduction may disadvantage the cluster index.

This job removes both. It exports the SAME 60 ablation datasets (deposited generator, same
seeds as job 1) and runs the actual CRAN `sigclust` on the full p=50 matrix with
`icovest=2`, the soft-thresholded estimator Liu et al. recommend for HDLSS.

Checklist item 20 allows "the `sigclust` R package or a reimplementation of its cluster-index
statistic". Having both, reported side by side, is strictly better than either: agreement
means the job 1 result is not an artifact of the reimplementation, and disagreement localises
the increment to the covariance estimator rather than to the statistic.

Requires R with `sigclust`. Deterministic seeds. Writes nothing to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, tempfile
from _provenance_safe import safe_argv  # noqa: E402
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import gate_ablation_study as gas  # noqa: E402


R_SCRIPT = r"""
suppressMessages(library(sigclust))
args <- commandArgs(trailingOnly=TRUE)
csv <- args[1]; out <- args[2]; nsim <- as.integer(args[3]); seed <- as.integer(args[4])
nrep <- if (length(args) >= 5) as.integer(args[5]) else 1
x <- as.matrix(read.csv(csv, header=FALSE))
set.seed(seed)
# icovest=2 : soft-thresholded covariance estimate, Liu et al.'s HDLSS recommendation
r <- sigclust(x, nsim=nsim, nrep=nrep, labflag=0, icovest=2)
cat(sprintf("%.6f %.6f\n", r@pval, r@xcindex), file=out)
"""


import zlib  # noqa: E402

# NOTE: was `hash(<str>) % 997`. Python salts str hashes per process (PYTHONHASHSEED),
# so that seeding produced DIFFERENT datasets on every invocation and nothing derived
# from it could be regenerated. zlib.crc32 is stable across processes and platforms.
def _stable_key(name: str) -> int:
    """Process-stable substitute for hash(name) % 997."""
    return zlib.crc32(name.encode()) % 997



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
    ap.add_argument("--reps", type=int, default=15)
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--sep", type=float, default=3.0)
    ap.add_argument("--nsim", type=int, default=100)
    ap.add_argument("--nrep", type=int, default=1,
                    help="k-means restarts inside sigclust (CRAN default 1; use 100 for a converged comparison)")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    _require_sigclust(R_SCRIPT)

    with tempfile.TemporaryDirectory() as td:
        rs = os.path.join(td, "sc.R")
        open(rs, "w").write(R_SCRIPT)
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as fh:
            fh.write(json.dumps(dict(
                record="provenance", script=os.path.basename(__file__), argv=safe_argv(a),
                design="canonical CRAN sigclust, full p-dim feature space, icovest=2",
                closes="job 1 caveats 1 (reimplementation) and 2 (PCA-reduced space)")) + "\n")
            fh.flush()
            tally = {}
            for cond in gas.CONDITIONS:
                klass = gas.CONDITIONS[cond]["klass"]
                ps = []
                for rep in range(a.reps):
                    seed = a.seed + 1000 * rep + _stable_key(cond)
                    scores, k = gas.generate(cond, a.n, a.p, a.sep, seed)
                    csv = os.path.join(td, "x.csv")
                    pd.DataFrame(scores.values).to_csv(csv, header=False, index=False)
                    txt = os.path.join(td, "o.txt")
                    r = subprocess.run(["Rscript", "--vanilla", rs, csv, txt,
                                        str(a.nsim), str(seed), str(a.nrep)],
                                       capture_output=True, text=True)
                    if r.returncode != 0 or not os.path.exists(txt):
                        raise SystemExit(
                            f"R/sigclust failed on {cond} rep{rep} -- refusing to continue.\n"
                            f"  A skipped dataset silently lowers the certification count,\n"
                            f"  which is the direction that flatters this paper.\n"
                            f"  Rscript said: {r.stderr.strip()[:300]}")
                    pval, cindex = [float(v) for v in open(txt).read().split()]
                    ps.append(pval)
                    fh.write(json.dumps(dict(record="dataset", condition=cond, klass=klass,
                                             rep=rep, p=pval, cindex=cindex)) + "\n")
                    fh.flush()
                    print(f"  {cond:16s} rep{rep:<2d} p={pval:.4f} cindex={cindex:.4f}",
                          flush=True)
                if ps:
                    ps = np.asarray(ps)
                    cert = int((ps < 0.05).sum())
                    tally[cond] = dict(klass=klass, n=len(ps), certify=cert,
                                       median_p=float(np.median(ps)))
                    fh.write(json.dumps(dict(record="condition_summary", condition=cond,
                                             **tally[cond])) + "\n")
                    fh.flush()
                    print(f"  -> {cond}: certify {cert}/{len(ps)} "
                          f"median p={np.median(ps):.4f}", flush=True)
            pos = [v for v in tally.values() if v["klass"] == "positive"]
            neg = [v for v in tally.values() if v["klass"] == "negative"]
            summary = dict(record="overall",
                           TPR=f"{sum(v['certify'] for v in pos)}/{sum(v['n'] for v in pos)}",
                           FPR=f"{sum(v['certify'] for v in neg)}/{sum(v['n'] for v in neg)}")
            fh.write(json.dumps(summary) + "\n")
            print(f"\nCANONICAL SIGCLUST — TPR {summary['TPR']}  FPR {summary['FPR']}")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
