#!/usr/bin/env python3
"""JOB 1e — canonical SigClust on the HETEROSCEDASTIC generators (closes a reporting gap).

WHY THIS EXISTS. Result 5's comparison table carried a "canonical SigClust" column for the
heteroscedastic continua, but no script in this deposit ever computed it: `grep -l icovest`
returned only job1c (the ablation grid) and job1d (the separation sweep), and
`job5_heteroscedastic_fpr.jsonl` contains no sigclust field. The column's "0/15" denominator
also disagreed with job5's 20 reps/cell, which is what exposed it. Those figures had been
computed ad hoc and never deposited — so a reviewer could not check the single concession the
paper leans on hardest ("on heteroscedastic continua SigClust is decisively better", and the
"run both" recommendation that follows from it).

DESIGN. Import job5's generator and reproduce its Phase-B seeding exactly
(`seed = a.seed + 1000*rep + n`), so the datasets scored here are **the same datasets** job5
scored, by construction rather than by re-implementation. Then run the canonical CRAN
`sigclust` on each — `icovest=2`, full p-dimensional feature space — which is the same
invocation job1c and job1d use.

Reporting both arms matters: the heteroscedastic arm is where the screen fails, and the matched
homoscedastic control is where it behaves. If SigClust is flat across both, that is the finding;
if it fails in the same place, the paper's concession is wrong and must be withdrawn.

Requires R >= 4 with the CRAN `sigclust` package. Deterministic. Writes nothing to the
deposited tree.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, tempfile
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from job5_heteroscedastic_fpr_confirmation import gen, KINDS, wilson  # noqa: E402

R_SCRIPT = r"""
suppressMessages(library(sigclust))
a <- commandArgs(trailingOnly=TRUE)
x <- as.matrix(read.csv(a[1], header=FALSE)); set.seed(as.integer(a[4]))
nrep <- if (length(a) >= 5) as.integer(a[5]) else 1
r <- sigclust(x, nsim=as.integer(a[3]), nrep=nrep, labflag=0, icovest=2)
cat(sprintf("%.6f %.6f\n", r@pval, r@xcindex), file=a[2])
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ns", default="120,300,600")
    ap.add_argument("--reps", type=int, default=20)   # matches job5
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--nsim", type=int, default=100)
    ap.add_argument("--nrep", type=int, default=1,
                    help="k-means restarts inside sigclust (CRAN default 1; use 100 for a converged comparison)")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    ns = [int(x) for x in a.ns.split(",")]

    if subprocess.run(["which", "Rscript"], capture_output=True).returncode != 0:
        sys.exit("Rscript not found. This job needs R >= 4 with the CRAN `sigclust` package.")

    with tempfile.TemporaryDirectory() as td:
        rs = os.path.join(td, "sc.R"); open(rs, "w").write(R_SCRIPT)
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as fh:
            w = lambda o: (fh.write(json.dumps(o) + "\n"), fh.flush())
            w(dict(record="provenance", script=os.path.basename(__file__),
                   argv={k: v for k, v in vars(a).items()},
                   design="canonical CRAN sigclust (icovest=2, full feature space) on job5's "
                          "generators, reproducing job5 Phase-B seeding exactly",
                   seeding="seed = args.seed + 1000*rep + n  (identical to job5)"))
            for n in ns:
                for kind in KINDS:
                    ps, cert = [], 0
                    for rep in range(a.reps):
                        seed = a.seed + 1000 * rep + n          # job5:134, verbatim
                        X = gen(kind, n, a.p, seed)
                        csv = os.path.join(td, "x.csv")
                        pd.DataFrame(X.values).to_csv(csv, header=False, index=False)
                        txt = os.path.join(td, "o.txt")
                        r = subprocess.run(["Rscript", "--vanilla", rs, csv, txt,
                                            str(a.nsim), str(seed), str(a.nrep)], capture_output=True, text=True)
                        if r.returncode != 0 or not os.path.exists(txt):
                            print(f"  R FAILED {kind} n={n} rep{rep}: {r.stderr.strip()[:160]}",
                                  flush=True)
                            continue
                        pval, ci = [float(v) for v in open(txt).read().split()]
                        os.remove(txt)
                        ps.append(pval); cert += int(pval < 0.05)
                        w(dict(record="run", n=n, kind=kind, rep=rep, p=pval, cindex=ci,
                               certify=bool(pval < 0.05)))
                    if not ps:
                        continue
                    lo, hi = wilson(cert, len(ps))
                    w(dict(record="cell_summary", n=n, kind=kind, reps=len(ps),
                           certify=cert, fpr=cert / len(ps), wilson95=[lo, hi],
                           median_p=float(np.median(ps))))
                    print(f"  n={n:<4} {kind:26s} SigClust FPR={cert:>2}/{len(ps)} "
                          f"= {cert/len(ps):.2f}  Wilson95 [{lo:.2f}, {hi:.2f}]  "
                          f"median p={np.median(ps):.4f}", flush=True)
            w(dict(record="done"))
    print(f"\nJOB1e done -> {a.out}")


if __name__ == "__main__":
    main()
