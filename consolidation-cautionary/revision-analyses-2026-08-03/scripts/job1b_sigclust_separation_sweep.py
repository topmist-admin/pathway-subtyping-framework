#!/usr/bin/env python3
"""JOB 1b — SigClust cluster index across the SEPARATION SWEEP.

Job 1 compared the two statistics only at the extremes (sep=3.0 positives, pure negatives),
where both saturate. The decisive question is the borderline regime: does bootstrap-ARI
resolve structure EARLIER or LATER than the 2-means cluster index, holding the null fixed?

If CI resolves earlier, the paper's statistic is worse than the one it replaced.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from sklearn.cluster import KMeans
from pathway_subtyping.discreteness.gate_a_discreteness_null import reduce_scores, reduced_dim

def ci(X, seed):
    km = KMeans(n_clusters=2, n_init=10, random_state=seed).fit(X)
    w = float(((X - km.cluster_centers_[km.labels_])**2).sum())
    t = float(((X - X.mean(0))**2).sum())
    return w/t if t > 0 else float("nan")

def make(sep, n, p, rng):
    n1 = n//2
    a = rng.normal(0,1,size=(n1,p)); b = rng.normal(0,1,size=(n-n1,p)); b[:,0]+=sep
    return np.vstack([a,b])

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--reps", type=int, default=20)
ap.add_argument("--n", type=int, default=120); ap.add_argument("--p", type=int, default=50)
ap.add_argument("--n-ref", type=int, default=100); ap.add_argument("--seed", type=int, default=42)
a = ap.parse_args()
SEPS=[0.0,0.5,1.0,1.5,2.0,2.5,3.0]
t0=time.time()
with open(a.out,"w") as fh:
    # Record basenames only: absolute paths are machine-specific and would leak into the
    # deposit. Content identity is carried by the deposited values, not by paths.
    _argv = {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v)
             for k, v in vars(a).items()}
    fh.write(json.dumps(dict(record="provenance", script="job1b", argv=_argv))+"\n")
    for sep in SEPS:
        cert=0; ps=[]
        for rep in range(a.reps):
            seed=a.seed+1000*rep+int(sep*10)
            rng=np.random.default_rng(seed)
            X=make(sep,a.n,a.p,rng)
            Z,_=reduce_scores(X, reduced_dim(a.n), seed)
            obs=ci(Z,seed); mean,var=Z.mean(0),Z.var(0)
            ss=np.random.SeedSequence(seed)
            sds=[int(np.random.default_rng(s).integers(1,2**31)) for s in ss.spawn(a.n_ref)]
            ref=np.array([ci(np.random.default_rng(sd).normal(mean,np.sqrt(var),size=Z.shape),sd) for sd in sds])
            ref=ref[np.isfinite(ref)]
            p=float((np.sum(ref<=obs)+1)/(len(ref)+1)); ps.append(p); cert+=int(p<0.05)
            fh.write(json.dumps(dict(record="cell",sep=sep,rep=rep,observed_ci=obs,p=p))+"\n")
        fh.write(json.dumps(dict(record="sep_summary",sep=sep,reps=a.reps,
                                 certify_rate=cert/a.reps,median_p=float(np.median(ps))))+"\n")
        fh.flush(); print(f"  sep={sep:<4} certify_rate={cert/a.reps:.2f}  median p={np.median(ps):.3f}",flush=True)
    fh.write(json.dumps(dict(record="done",elapsed_sec=round(time.time()-t0,1)))+"\n")
print(f"JOB1b done in {time.time()-t0:.0f}s")
