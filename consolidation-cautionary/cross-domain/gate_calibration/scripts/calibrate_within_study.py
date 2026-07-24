#!/usr/bin/env python3
"""Real-data calibration of the discreteness gate — WITHIN-STUDY controls.

Supersedes `calibrate_discreteness_gate.py`. That version's discrete positive
control pooled three cBioPortal studies (COAD/GBM/LUAD), one tumour type each, so
study was perfectly confounded with the ground-truth label: within-study
normalisation collapsed recovery to ARI 0.05, and only across-study normalisation
(which reintroduces the batch axis) gave 0.92. A reviewer correctly flagged that
the "discrete" control was largely a batch effect — the very failure this paper
diagnoses elsewhere. This script removes the confound by drawing BOTH controls
from within a SINGLE study each.

  DISCRETE positive control (within one study) — IDH-mutant vs IDH-wildtype
  low-grade glioma (TCGA-LGG). IDH status is the single most discrete molecular
  dichotomy in adult glioma (a near-binary methylation/expression state), and both
  arms come from the SAME cBioPortal study, so study cannot stand in for the label.
  Within-study median z-scores are used — the normalisation the pooled version had
  to abandon. The gate MUST certify; recovery of IDH status must be strong.

  CONTINUUM negative control (within one study) — the immune-infiltration gradient
  in TCGA-LUAD, unchanged from the prior version (it was already within-study and
  is not affected by the confound). A continuous cytolytic signature, confirmed
  unimodal by Hartigan's dip; the k=2 split merely bisects it. The gate MUST reject.

Reports whatever the gate does. Deterministic (seed 42). Network: cbioportal.org.
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
PANEL = os.path.join(_HERE, "../../../panels/hallmark_200genes.gmt")
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
            if r.status_code in (502, 503, 504):
                last = requests.HTTPError(str(r.status_code)); time.sleep(2 ** attempt); continue
            r.raise_for_status(); return r.json()
        except (requests.ConnectionError, requests.Timeout) as ex:
            last = ex; time.sleep(2 ** attempt)
    raise last


def get(url, **p):
    r = S.get(f"{API}{url}", params=p, timeout=180); r.raise_for_status(); return r.json()


def read_gmt(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = [g for g in f[2:] if g]
    return d


def fetch_expr(study, symbols, zscored=True):
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


def fetch_subtype(study):
    """PATIENT-level SUBTYPE, mapped to each sample id in the study."""
    rows = get(f"/studies/{study}/clinical-data", clinicalDataType="PATIENT",
               attributeId="SUBTYPE", projection="SUMMARY")
    by_patient = {r["patientId"]: r["value"] for r in rows}
    # sample -> patient
    samples = post(f"/studies/{study}/samples/fetch", {}, projection="SUMMARY") \
        if False else get(f"/studies/{study}/samples", projection="SUMMARY")
    s2p = {s["sampleId"]: s["patientId"] for s in samples}
    return {s: by_patient.get(p) for s, p in s2p.items()}


def pathway_matrix(expr, pw):
    return pd.DataFrame({n: expr[[x for x in gs if x in expr.columns]].mean(axis=1)
                        for n, gs in pw.items()
                        if len([x for x in gs if x in expr.columns]) >= 2}).dropna()


def gateA_call(P, k):
    a = DiscretenessGateA(seed=42, n_ref=150).run("calib", P, k, gmm_seed=42)
    floor = 1.0 / (150 + 1)
    return {"certified": bool(a.testable and a.passed), "verdict": a.verdict,
            "testable": bool(a.testable),
            "sg_p": round(float(a.sg_empirical_p), 4),
            "sg_p_at_floor": bool(abs(float(a.sg_empirical_p) - floor) < 1e-4),
            "dip_pc1_p": round(float(a.dip_pc1_p), 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    ap.add_argument("--discrete-study", default="lgg_tcga_pan_can_atlas_2018")
    ap.add_argument("--continuum-study", default="luad_tcga_pan_can_atlas_2018")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    pw = read_gmt(PANEL)
    hallmark_genes = sorted({x for gs in pw.values() for x in gs})
    out = {"design_note": ("Both controls are WITHIN a single study each; study is "
                           "not confounded with the ground-truth label. Supersedes "
                           "the pooled-3-study version (batch-confounded).")}

    # ---- DISCRETE positive control: IDH-mut vs IDH-wt WITHIN TCGA-LGG ----
    st = args.discrete_study
    print(f"DISCRETE control (within-study): IDH status in {st}")
    e = fetch_expr(st, hallmark_genes, zscored=True)   # within-study z-scores: valid, one study
    sub = fetch_subtype(st)
    # collapse LGG SUBTYPE -> binary IDH status
    def idh(v):
        if not isinstance(v, str):
            return None
        if "IDHmut" in v:
            return "IDHmut"
        if "IDHwt" in v:
            return "IDHwt"
        return None
    label = pd.Series({s: idh(sub.get(s)) for s in e.index}).dropna()
    P = pathway_matrix(e, pw)
    common = P.index.intersection(label.index)
    P = P.loc[common]; label = label.loc[common]
    X = StandardScaler().fit_transform(P.values)
    lab = GaussianMixture(2, covariance_type="full", n_init=10, random_state=42,
                          reg_covar=1e-6).fit(X).predict(X)
    ari = float(adjusted_rand_score(label.astype("category").cat.codes, lab))
    out["discrete_control"] = {
        "study": st, "axis": "IDH-mutant vs IDH-wildtype (single study)",
        "n": int(len(P)),
        "n_idhmut": int((label == "IDHmut").sum()),
        "n_idhwt": int((label == "IDHwt").sum()),
        "normalization": "within-study median z-scores (no batch confound)",
        "recovery_of_idh_ari": round(ari, 3),
        "gateA": gateA_call(P, 2),
        "EXPECTED": "certified=True (IDH status is genuinely discrete)"}
    print(f"  n={len(P)} ({(label=='IDHmut').sum()} mut / {(label=='IDHwt').sum()} wt); "
          f"IDH recovery ARI={ari:.3f}; gateA={out['discrete_control']['gateA']}")

    # ---- CONTINUUM negative control: immune gradient (within-study, unchanged) ----
    st = args.continuum_study
    print(f"\nCONTINUUM control (within-study): immune gradient in {st}")
    e = fetch_expr(st, hallmark_genes + IMMUNE_SIG)
    imm_present = [g for g in IMMUNE_SIG if g in e.columns]
    immune_score = e[imm_present].mean(axis=1)
    dip = dip_of(immune_score.values)
    immune_pw = {n: gs for n, gs in pw.items()
                 if any(t in n.upper() for t in ["IMMUN", "INFLAMMAT", "INTERFERON",
                                                 "TNFA", "IL6", "IL2", "COMPLEMENT",
                                                 "ALLOGRAFT"])}
    Pc = pathway_matrix(e, immune_pw)
    Pc = Pc.loc[immune_score.index.intersection(Pc.index)]
    Xc = StandardScaler().fit_transform(Pc.values)
    labc = GaussianMixture(2, covariance_type="full", n_init=10, random_state=42,
                           reg_covar=1e-6).fit(Xc).predict(Xc)
    imm = immune_score.loc[Pc.index]
    sep = float(abs(imm[labc == 0].mean() - imm[labc == 1].mean()) / imm.std())
    out["continuum_control"] = {
        "study": st, "n": int(len(Pc)), "immune_pathways": len(immune_pw),
        "immune_score_unimodal_dip_p": round(float(dip["p"]), 4),
        "k2_split_tracks_immune_score_(std_gap)": round(sep, 2),
        "gateA": gateA_call(Pc, 2),
        "EXPECTED": "certified=False (continuous gradient, not discrete)"}
    print(f"  immune-score dip p={dip['p']:.3f}; split std-gap={sep:.2f}; "
          f"gateA={out['continuum_control']['gateA']}")

    disc_ok = out["discrete_control"]["gateA"]["certified"] is True
    cont_ok = out["continuum_control"]["gateA"]["certified"] is False
    out["calibration_verdict"] = {
        "discrete_certified_correctly": disc_ok,
        "continuum_rejected_correctly": cont_ok,
        "both_correct": bool(disc_ok and cont_ok),
        "caveat": ("Two anchors at opposite ends of the difficulty range; this "
                   "shows the gate is not degenerate, not that it resolves "
                   "borderline within-disease structure. See the separation sweep "
                   "for the operating characteristic.")}
    with open(os.path.join(args.out, "gate_calibration_within_study.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n=== VERDICT === discrete_ok={disc_ok} continuum_ok={cont_ok} "
          f"both={out['calibration_verdict']['both_correct']}")
    print(f"Wrote {args.out}/gate_calibration_within_study.json")


if __name__ == "__main__":
    main()
