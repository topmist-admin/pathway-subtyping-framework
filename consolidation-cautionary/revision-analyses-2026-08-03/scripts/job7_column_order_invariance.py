#!/usr/bin/env python3
"""JOB 7 — is the Result 4 partition invariant to gene column order? (Paulus item 24.)

THE STANDING CAVEAT. Methods currently says: "Because ssGSEA ranks genes within a sample,
its output depends on gene column order wherever expression values are exactly tied;
log2(CPM + 1) maps every zero count to an exact tie. Deposited matrices must therefore be
read as shipped." That is an open caveat on the flagship Result, and it is testable.

WHAT THIS TESTS, AND WHAT IT DOES NOT. The deposited Result 4 scores were generated under
framework v0.3.0; only the stability figure was recomputed later. So comparing a fresh
v0.9.0 scoring against the deposited matrix would confound the scoring version with column
order and could not answer the question. Instead this holds EVERYTHING fixed — same raw
counts, same panel, same samples, same v0.9.0 scorer, same seed — and varies ONLY the order
of the gene columns:

    arm "shipped"   : genes in the order they appear in the GEO matrix
    arm "canonical" : the same genes, sorted lexicographically

If ssGSEA's tie handling is order-dependent in a way that matters, the two arms will differ.
The reported quantities are then the ones the caveat threatens: the pathway scores
themselves, the k=3 partition, and the bootstrap stability figure.

A third arm, "reversed", is included as a stronger perturbation: if shipped-vs-canonical
happens to be a benign permutation, reversal is the adversarial case.

Data: GSE80655 counts from GEO (public). Panel: the deposited 14-set schizophrenia GMT.
Samples: the 141 in the deposited Result 4 partition. Deterministic (seed 42).
Writes nothing to the deposited tree.
"""
from __future__ import annotations
import argparse, gzip, json, os, sys, time
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from sklearn.mixture import GaussianMixture                                    # noqa: E402
from sklearn.metrics import adjusted_rand_score                                # noqa: E402
from pathway_subtyping.expression import (                                     # noqa: E402
    score_pathways_from_expression, ExpressionScoringMethod,
)
from pathway_subtyping.validation import ValidationGates                       # noqa: E402

SEED = 42
PANEL = os.path.join(REPO, "consolidation-cautionary/panels/schizophrenia_pathways.gmt")
META = os.path.join(REPO, "consolidation-cautionary/data/partition/sample_metadata_with_subtypes.csv")
SYM2ENS = os.environ.get("SYM2ENS_JSON") or os.path.join(
    REPO, "consolidation-cautionary/revision-analyses-2026-08-03/inputs/"
          "symbol_to_ensembl_scz_panel.json")   # deposited; avoids a mygene.info call


def read_gmt(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                out[f[0]] = [g for g in f[2:] if g]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n-bootstrap", type=int, default=100)
    a = ap.parse_args()
    t0 = time.time()

    # ---------- expression: genes x samples -> samples x genes ---------- #
    with gzip.open(a.expr, "rt") as fh:
        expr = pd.read_csv(fh, sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated(keep="first")]
    print(f"  raw expression {expr.shape} (genes x samples)", flush=True)

    meta = pd.read_csv(META)
    sl = meta["title"].str.extract(r"_(SL\d+)$")[0]
    keep = [s for s in sl.dropna() if s in expr.columns]
    print(f"  matched {len(keep)} of {len(meta)} deposited samples to GEO columns", flush=True)

    X = expr[keep].T                                   # samples x genes, raw counts

    # ---------- panel: symbols -> ENSG ---------- #
    panel = read_gmt(PANEL)
    s2e = json.load(open(SYM2ENS)) if os.path.exists(SYM2ENS) else {}
    have = set(X.columns)
    pathways, cover = {}, {}
    for name, genes in panel.items():
        ids = []
        for g in genes:
            e = s2e.get(g)
            if isinstance(e, list):
                e = e[0] if e else None
            if e and e in have:
                ids.append(e)
            elif g in have:
                ids.append(g)
        ids = sorted(set(ids))
        cover[name] = (len(ids), len(genes))
        if len(ids) >= 2:
            pathways[name] = ids
    print(f"  panel: {len(pathways)}/{len(panel)} sets with >=2 mapped genes", flush=True)
    for n_, (a_, b_) in cover.items():
        print(f"     {n_:34s} {a_:>4}/{b_}", flush=True)
    if not pathways:
        print("  NO PATHWAY COVERAGE — cannot run; aborting rather than reporting a null result")
        sys.exit(2)

    used = sorted({g for ids in pathways.values() for g in ids})
    X = X[[c for c in X.columns if c in set(used)]]
    X = np.log2(X.div(X.sum(1), axis=0) * 1e6 + 1)     # log2(CPM+1): zeros -> exact ties
    ties = int((X == X.iloc[:, 0].min()).sum().sum())
    print(f"  scoring matrix {X.shape}; log2(CPM+1); tied-at-floor cells = {ties}", flush=True)

    # ---------- three column orders ---------- #
    orders = {"shipped": list(X.columns),
              "canonical": sorted(X.columns),
              "reversed": list(X.columns)[::-1]}
    scores, parts, stabs = {}, {}, {}
    for tag, cols in orders.items():
        Xo = X[cols]
        P = score_pathways_from_expression(
            Xo, pathways, method=ExpressionScoringMethod.SSGSEA,
            seed=SEED, show_progress=False).pathway_scores
        P = P.reindex(sorted(P.columns), axis=1)       # align for comparison only
        scores[tag] = P
        lab = GaussianMixture(a.k, covariance_type="full", n_init=10,
                              random_state=SEED, reg_covar=1e-6).fit(P.values).predict(P.values)
        parts[tag] = lab
        stabs[tag] = float(ValidationGates(seed=SEED, n_bootstrap=a.n_bootstrap,
                                           show_progress=False)
                           .stability_test_bootstrap(P, lab, a.k, gmm_seed=SEED).metric_value)
        print(f"  [{tag:9s}] scores {P.shape}  stability={stabs[tag]:.4f}", flush=True)

    base = "shipped"
    recs = [dict(record="provenance", script=os.path.basename(__file__), argv=vars(a),
                 n_samples=int(X.shape[0]), n_genes=int(X.shape[1]),
                 n_pathways=len(pathways), pathway_coverage=cover,
                 design="everything fixed except gene column order")]
    for tag in orders:
        if tag == base:
            continue
        d = float(np.abs(scores[tag].values - scores[base].values).max())
        ari = float(adjusted_rand_score(parts[base], parts[tag]))
        recs.append(dict(record="comparison", arm=tag, vs=base,
                         max_abs_score_diff=d,
                         partition_ari_vs_shipped=ari,
                         identical_partition=bool(ari == 1.0),
                         stability_shipped=stabs[base], stability_arm=stabs[tag],
                         stability_abs_diff=abs(stabs[tag] - stabs[base])))
        print(f"  {tag:9s} vs {base}: max|Δscore|={d:.3e}  partition ARI={ari:.4f}  "
              f"stability {stabs[base]:.4f} -> {stabs[tag]:.4f}", flush=True)
    recs.append(dict(record="done", elapsed_sec=round(time.time() - t0, 1)))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    print(f"\nJOB7 done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
