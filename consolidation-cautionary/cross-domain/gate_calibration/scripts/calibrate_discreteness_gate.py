#!/usr/bin/env python3
"""Real-data calibration of the discreteness gate (Gate A).

The cautionary-framework paper hinges on one claim a hostile reviewer will attack:
is the discreteness gate CALIBRATED, or just pessimistic? A synthetic ablation is
not enough — the gate must, on REAL public data, (1) CERTIFY genuinely discrete
structure and (2) REJECT a genuine continuous gradient.

Two real-data controls, both from public cBioPortal (no auth):

  DISCRETE positive control — pool three biologically distinct tumor types
  (colorectal, glioblastoma, lung adenocarcinoma). Ground truth = tumor type;
  this is unambiguously discrete. The gate MUST certify discrete structure and
  recovery of tumor-of-origin must be near-perfect. If it fails here, the gate
  is broken (blind to obvious structure).

  CONTINUUM negative control — within a single tumor type, isolate the immune-
  infiltration axis: a continuous cytolytic/T-cell signature score (canonical,
  Rooney et al. 2015-style). We (a) confirm the score is UNIMODAL (Hartigan dip),
  i.e. a gradient not two groups, (b) cluster k=2 in immune-pathway space, (c)
  show the k=2 split just bisects the continuous score, (d) ask the gate its
  verdict. A calibrated gate REJECTS this as a continuum.

The experiment REPORTS whatever the gate does — a mis-call is a real limitation to
disclose, not to hide.

Deterministic (seed 42). Network: www.cbioportal.org. Requires pathway-subtyping
(>=0.8 line), numpy, pandas, scikit-learn, scipy, requests.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import requests
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

_HERE = os.path.dirname(os.path.abspath(__file__))
try:
    from pathway_subtyping.discreteness import DiscretenessGateA
    from pathway_subtyping.discreteness.gate_a_discreteness_null import dip_of
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../..")))
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../../src")))
    from pathway_subtyping.discreteness import DiscretenessGateA
    from pathway_subtyping.discreteness.gate_a_discreteness_null import dip_of

API = "https://www.cbioportal.org/api"

# --- shared provenance recording (see consolidation-cautionary/scripts/_provenance.py) -
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "..", "..", "scripts"))
from _provenance import env_provenance, fetch_provenance  # noqa: E402

PANEL = os.path.join(_HERE, "../../../panels/hallmark_200genes.gmt")
# canonical cytolytic / T-cell infiltration signature (continuous axis)
IMMUNE_SIG = ["CD8A", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7", "CCL5",
              "CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "GBP1", "HLA-DRA"]

S = requests.Session()
S.headers.update({"Accept": "application/json"})


def post(url, payload, **p):
    import time
    last = None
    for attempt in range(5):
        try:
            r = S.post(f"{API}{url}", params=p, data=json.dumps(payload),
                       headers={"Content-Type": "application/json"}, timeout=300)
            if r.status_code in (502, 503, 504):        # transient server errors
                last = requests.HTTPError(f"{r.status_code}")
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            return r.json()
        except (requests.ConnectionError, requests.Timeout) as ex:
            last = ex
            time.sleep(2 ** attempt)
    raise last


def read_gmt(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = [g for g in f[2:] if g]
    return d


def fetch_expr(study, symbols, zscored=True):
    """Gene matrix (samples x genes). zscored=True -> within-study median Zscores
    (fine for a SINGLE cohort); zscored=False -> continuous RSEM (needed when
    POOLING cohorts, since within-study Zscores erase between-cohort signal)."""
    profile = (f"{study}_rna_seq_v2_mrna_median_all_sample_Zscores" if zscored
               else f"{study}_rna_seq_v2_mrna")
    genes = post("/genes/fetch", sorted(set(symbols)),
                 geneIdType="HUGO_GENE_SYMBOL", projection="SUMMARY")
    sym2ent = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in genes}
    ent2sym = {v: k for k, v in sym2ent.items()}
    md = post(f"/molecular-profiles/{profile}/molecular-data/fetch",
              {"entrezGeneIds": sorted(set(sym2ent.values())),
               "sampleListId": f"{study}_rna_seq_v2_mrna"}, projection="SUMMARY")
    return (pd.DataFrame([(m["sampleId"], ent2sym.get(m["entrezGeneId"]), m["value"])
                         for m in md], columns=["s", "g", "v"]).dropna()
            .pivot_table(index="s", columns="g", values="v", aggfunc="mean"))


def pathway_matrix(expr, pw):
    return pd.DataFrame({n: expr[[x for x in gs if x in expr.columns]].mean(axis=1)
                        for n, gs in pw.items()
                        if len([x for x in gs if x in expr.columns]) >= 2}).dropna()


def gateA_call(P, k):
    a = DiscretenessGateA(seed=42, n_ref=150).run("calib", P, k, gmm_seed=42)
    return {"certified": bool(a.testable and a.passed), "verdict": a.verdict,
            "testable": bool(a.testable), "sg_p": round(float(a.sg_empirical_p), 4),
            "dip_pc1_p": round(float(a.dip_pc1_p), 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    ap.add_argument("--continuum-study", default="luad_tcga_pan_can_atlas_2018")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    pw = read_gmt(PANEL)
    hallmark_genes = sorted({x for gs in pw.values() for x in gs})
    out = {}

    # ---- DISCRETE positive control: pool 3 distinct tumor types ----
    studies = {"COAD": "coadread_tcga_pan_can_atlas_2018",
               "GBM": "gbm_tcga_pan_can_atlas_2018",
               "LUAD": "luad_tcga_pan_can_atlas_2018"}
    print("DISCRETE control: pooling", list(studies))
    # fetch CONTINUOUS expression per study, pool, then log2 + pooled-gene z-score
    # (within-study Zscores would erase the between-tumor-type signal)
    raw, origin = [], []
    for name, st in studies.items():
        e = fetch_expr(st, hallmark_genes, zscored=False)
        raw.append(e)
        origin += [name] * len(e)
        print(f"  {name}: {e.shape}")
    E = pd.concat(raw).dropna(axis=1)
    origin = np.array(origin)
    Elog = np.log2(E.clip(lower=0) + 1.0)
    Ez = (Elog - Elog.mean(0)) / (Elog.std(0) + 1e-9)   # z-score each gene ACROSS pool
    Pd = pathway_matrix(Ez, pw)
    origin = origin[[E.index.get_loc(s) for s in Pd.index]]
    Xd = StandardScaler().fit_transform(Pd.values)
    lab = GaussianMixture(3, covariance_type="full", n_init=10, random_state=42,
                          reg_covar=1e-6).fit(Xd).predict(Xd)
    ari = float(adjusted_rand_score(pd.Series(origin).astype("category").cat.codes, lab))
    out["discrete_control"] = {"n": int(len(Pd)), "tumor_types": list(studies),
                               "recovery_of_origin_ari": round(ari, 3),
                               "gateA": gateA_call(Pd, 3),
                               "EXPECTED": "certified=True (real discrete structure)"}
    print(f"  origin recovery ARI={ari:.3f}; gateA={out['discrete_control']['gateA']}")

    # ---- CONTINUUM negative control: immune-infiltration gradient ----
    st = args.continuum_study
    print(f"\nCONTINUUM control: immune gradient in {st}")
    e = fetch_expr(st, hallmark_genes + IMMUNE_SIG)
    imm_present = [g for g in IMMUNE_SIG if g in e.columns]
    immune_score = e[imm_present].mean(axis=1)          # continuous infiltration axis
    dip = dip_of(immune_score.values)                    # unimodality of the axis
    # immune-pathway subspace so the isolated axis IS the immune gradient
    immune_pw = {n: gs for n, gs in pw.items()
                 if any(t in n.upper() for t in ["IMMUN", "INFLAMMAT", "INTERFERON",
                                                 "TNFA", "IL6", "IL2", "COMPLEMENT",
                                                 "ALLOGRAFT"])}
    Pc = pathway_matrix(e, immune_pw)
    Pc = Pc.loc[immune_score.index.intersection(Pc.index)]
    Xc = StandardScaler().fit_transform(Pc.values)
    labc = GaussianMixture(2, covariance_type="full", n_init=10, random_state=42,
                           reg_covar=1e-6).fit(Xc).predict(Xc)
    # does the k=2 split just bisect the continuous immune score?
    imm = immune_score.loc[Pc.index]
    sep = float(abs(imm[labc == 0].mean() - imm[labc == 1].mean()) / imm.std())
    out["continuum_control"] = {
        "study": st, "n": int(len(Pc)), "immune_pathways": len(immune_pw),
        "immune_score_unimodal_dip_p": round(float(dip["p"]), 4),
        "k2_split_tracks_immune_score_(std_gap)": round(sep, 2),
        "gateA": gateA_call(Pc, 2),
        "EXPECTED": "certified=False (continuous gradient, not discrete)"}
    print(f"  immune-score dip p={dip['p']:.3f} (unimodal if >0.05); "
          f"split std-gap={sep:.2f}; gateA={out['continuum_control']['gateA']}")

    # ---- verdict ----
    disc_ok = out["discrete_control"]["gateA"]["certified"] is True
    cont_ok = out["continuum_control"]["gateA"]["certified"] is False
    out["calibration_verdict"] = {
        "discrete_certified_correctly": disc_ok,
        "continuum_rejected_correctly": cont_ok,
        "both_correct": bool(disc_ok and cont_ok)}
    out["provenance"] = {"environment": env_provenance(),
                         "fetch": fetch_provenance(API, "pooled-3-study (see design_note)")}
    with open(os.path.join(args.out, "gate_calibration.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n=== CALIBRATION VERDICT ===")
    print(f"  discrete certified correctly: {disc_ok}")
    print(f"  continuum rejected correctly: {cont_ok}")
    print(f"  BOTH CORRECT: {out['calibration_verdict']['both_correct']}")
    print(f"\nWrote {args.out}/gate_calibration.json")


if __name__ == "__main__":
    main()
