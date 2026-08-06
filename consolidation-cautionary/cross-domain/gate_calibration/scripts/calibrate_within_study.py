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

# --- shared provenance recording (see consolidation-cautionary/scripts/_provenance.py) -
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "..", "..", "scripts"))
from _provenance import (env_provenance, fetch_provenance,  # noqa: E402
                         cache_read, cache_write)

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


def _require_diptest():
    """Fail closed if the optional `diptest` extra is missing.

    `dip_of()` returns NaN rather than raising when the extra is absent. Every dip value
    this script writes is reported as evidence, and a NaN would be serialised as a bare
    `NaN` token, which is not valid JSON (RFC 8259) and is accepted only by Python's own
    parser. Probe once, up front, so no partial output is ever written.
    """
    import numpy as _np
    if not _np.isfinite(dip_of(_np.linspace(0.0, 1.0, 32))["p"]):
        raise SystemExit(
            "Hartigan dip unavailable (install the `diptest` extra: pip install "
            "'pathway-subtyping[diptest]'). Refusing to run: this script reports dip "
            "p-values as evidence, and NaN would be written to the output as invalid JSON.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    ap.add_argument("--discrete-study", default="lgg_tcga_pan_can_atlas_2018")
    ap.add_argument("--continuum-study", default="luad_tcga_pan_can_atlas_2018")
    # Offline path. cBioPortal is an unversioned live service, so the fetched inputs are
    # deposited alongside the results and can be replayed with --cache-in. That turns this
    # package from "trust our fetch" into a checkable offline reproduction, the same
    # treatment already applied to the GSE28521 autism arm.
    ap.add_argument("--cache-dir", default=os.path.join(_HERE, "../results/cached_inputs"),
                    help="where deposited inputs live (written on a live run)")
    ap.add_argument("--cache-in", action="store_true",
                    help="replay from deposited inputs; no network access at all")
    args = ap.parse_args()
    _require_diptest()
    os.makedirs(args.out, exist_ok=True)
    pw = read_gmt(PANEL)
    hallmark_genes = sorted({x for gs in pw.values() for x in gs})
    out = {"design_note": ("Both controls are WITHIN a single study each; study is "
                           "not confounded with the ground-truth label. Supersedes "
                           "the pooled-3-study version (batch-confounded).")}

    # ---- DISCRETE positive control: IDH-mut vs IDH-wt WITHIN TCGA-LGG ----
    st = args.discrete_study
    print(f"DISCRETE control (within-study): IDH status in {st}")
    if args.cache_in:
        P_cached = cache_read(os.path.join(args.cache_dir, "discrete_pathways.csv.gz"))
        label = P_cached.pop("__label__")
        e = None
    else:
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

    if args.cache_in:
        P = P_cached
    else:
        label = pd.Series({s: idh(sub.get(s)) for s in e.index}).dropna()
        P = pathway_matrix(e, pw)
        common = P.index.intersection(label.index)
        P = P.loc[common]; label = label.loc[common]
        # Deposit the assembled inputs. Pathway scores + label are sufficient here:
        # nothing in this package resamples from the gene universe (unlike the autism
        # arm's random-gene-set control), so the 50-column matrix regenerates every
        # reported number and costs kilobytes rather than tens of megabytes.
        _c = P.copy(); _c["__label__"] = label
        cache_write(_c, os.path.join(args.cache_dir, "discrete_pathways.csv.gz"))
    X = StandardScaler().fit_transform(P.values)
    lab = GaussianMixture(2, covariance_type="full", n_init=10, random_state=42,
                          reg_covar=1e-6).fit(X).predict(X)
    ari = float(adjusted_rand_score(label.astype("category").cat.codes, lab))
    # Does the CERTIFIED partition actually align with IDH, or a different axis?
    # Cross-tab GMM cluster x IDH, and report the best cluster->IDH-wt mapping purity.
    ct = pd.crosstab(pd.Series(lab, index=P.index, name="cluster"), label)
    # which cluster is enriched for IDH-wt (the minority, transcriptionally distinct arm)
    wt_frac_by_cluster = (ct.get("IDHwt", pd.Series(0, index=ct.index)) /
                          ct.sum(axis=1)).to_dict()
    wt_cluster = max(wt_frac_by_cluster, key=wt_frac_by_cluster.get)
    # purity: of the samples the GMM puts in the IDH-wt-enriched cluster, what frac are IDH-wt;
    # and of all IDH-wt, what frac land in that cluster (recall)
    in_wt_cluster = (lab == wt_cluster)
    is_wt = (label.values == "IDHwt")
    precision = float((in_wt_cluster & is_wt).sum() / max(1, in_wt_cluster.sum()))
    recall = float((in_wt_cluster & is_wt).sum() / max(1, is_wt.sum()))
    # dip on PC1 of the pathway matrix (the diagnostic the reviewer wanted reported)
    from sklearn.decomposition import PCA
    pc1 = PCA(n_components=1, random_state=42).fit_transform(
        StandardScaler().fit_transform(P.values)).ravel()
    dip_pos = dip_of(pc1)
    out["discrete_control"] = {
        "study": st, "axis": "IDH-mutant vs IDH-wildtype (single study)",
        "n": int(len(P)),
        "n_idhmut": int((label == "IDHmut").sum()),
        "n_idhwt": int((label == "IDHwt").sum()),
        "normalization": "within-study median z-scores (no batch confound)",
        "recovery_of_idh_ari": round(ari, 3),
        "certified_partition_vs_idh": {
            "crosstab_cluster_x_idh": ct.to_dict(),
            "idhwt_enriched_cluster": int(wt_cluster),
            "precision_for_idhwt": round(precision, 3),
            "recall_of_idhwt": round(recall, 3),
            "note": ("does the CERTIFIED GMM partition track IDH? high precision/recall "
                     "=> the certified structure IS the IDH axis, and the modest ARI "
                     "reflects the 415-mut arm's internal codel/non-codel substructure "
                     "collapsed to binary, not a wrong axis")},
        # NOTE ON THE TWO DIP P-VALUES. This record carries two, and they legitimately
        # differ because they are computed on differently preprocessed matrices:
        #   * this key -- dip on PC1 of the pathway matrix after StandardScaler,
        #     i.e. pathway columns z-scored before PCA.
        #   * gateA.dip_pc1_p -- dip on PC1 of the gate's own internal reduction, which
        #     is CENTRED ONLY (reduce_scores() calls PCA directly on the raw scores).
        # Neither decides the gate; both are reported diagnostics. The key names now say
        # which is which so the pair cannot be read as an inconsistency.
        "pc1_dip_p_standardized_pathways": round(float(dip_pos["p"]), 4),
        "pc1_dip_p_note": ("dip on PC1 after z-scoring pathway columns; "
                           "gateA.dip_pc1_p is the same test on the gate's "
                           "centred-only reduction, so the two differ by design"),
        "gateA": gateA_call(P, 2),
        "EXPECTED": "certified=True (IDH status is genuinely discrete)"}
    print(f"  n={len(P)} ({(label=='IDHmut').sum()} mut / {(label=='IDHwt').sum()} wt); "
          f"IDH recovery ARI={ari:.3f}; gateA={out['discrete_control']['gateA']}")

    # ---- CONTINUUM negative control: immune gradient (within-study, unchanged) ----
    st = args.continuum_study
    print(f"\nCONTINUUM control (within-study): immune gradient in {st}")
    immune_pw = {n: gs for n, gs in pw.items()
                 if any(t in n.upper() for t in ["IMMUN", "INFLAMMAT", "INTERFERON",
                                                 "TNFA", "IL6", "IL2", "COMPLEMENT",
                                                 "ALLOGRAFT"])}
    if args.cache_in:
        Pc = cache_read(os.path.join(args.cache_dir, "continuum_pathways.csv.gz"))
        immune_score = Pc.pop("__immune_score__")
    else:
        e = fetch_expr(st, hallmark_genes + IMMUNE_SIG)
        imm_present = [g for g in IMMUNE_SIG if g in e.columns]
        immune_score = e[imm_present].mean(axis=1)
        Pc = pathway_matrix(e, immune_pw)
        Pc = Pc.loc[immune_score.index.intersection(Pc.index)]
        # The immune score rides along with the pathway matrix: it is the continuum this
        # control is built on, and re-deriving it needs the raw cytolytic genes.
        _c = Pc.copy(); _c["__immune_score__"] = immune_score.loc[Pc.index]
        cache_write(_c, os.path.join(args.cache_dir, "continuum_pathways.csv.gz"))
    dip = dip_of(immune_score.loc[Pc.index].values)
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
    out["provenance"] = {"environment": env_provenance(),
                         "fetch_discrete_control": fetch_provenance(API, args.discrete_study),
                         "fetch_continuum_control": fetch_provenance(API, args.continuum_study)}
    with open(os.path.join(args.out, "gate_calibration_within_study.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n=== VERDICT === discrete_ok={disc_ok} continuum_ok={cont_ok} "
          f"both={out['calibration_verdict']['both_correct']}")
    print(f"Wrote {args.out}/gate_calibration_within_study.json")


if __name__ == "__main__":
    main()
