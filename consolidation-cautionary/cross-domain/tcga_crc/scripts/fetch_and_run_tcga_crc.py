#!/usr/bin/env python3
"""
Real positive control for Gate 7 (somatic) — TCGA-CRC BRAF / KRAS / MSI.

Wires ``ValidationGates.somatic_anchoring_gate`` against real Colorectal
Adenocarcinoma data (TCGA PanCancer Atlas, via the public cBioPortal API). No CMS
labels are borrowed: the transcriptomic partition is built by the framework from
a DISEASE-AGNOSTIC Hallmark gene panel, so its alignment with the somatic strata
is a genuine, non-circular test.

Two results are produced (both on the same real cohort):

  A. THE WIRING (unbiased transcriptomic partition vs drivers).
     Score the Hallmark panel into 50 pathway means, cluster with a Gaussian
     mixture (k by BIC over 2..6 — not tuned to the somatic outcome), then run
     Gate 6 (confound: tissue) and Gate 7 (somatic: BRAF-V600E / KRAS / MSI).
     Expected/observed: the partition is NOT a tissue classifier (Gate 6 V~0),
     and its driver alignment is statistically significant but effect-size-modest
     (Cramer's V < 0.30) — so the somatic gate correctly returns NOT-anchored.
     This is the honest whole-transcriptome result (CMS-level recovery needs the
     dedicated CMS classifier, not generic Hallmark clustering) and it
     demonstrates the gate's effect-size discipline: it does not over-call a
     significant-but-weak association as an anchor.

  B. GATE-PASS POSITIVE CONTROL (real strong association).
     BRAF-V600E and MSI-high strongly co-occur in CRC (serrated / MLH1-methylation
     pathway) — a textbook strong association. Using MSI status as the partition
     and BRAF/KRAS as strata, the gate PASSES on BRAF-V600E (V ~ 0.55, q ~ 1e-36)
     and correctly rejects KRAS. This confirms the gate's PASS branch fires on
     real data, the somatic analog of the Voineagu control for feature-level
     anchoring.

Deterministic (seed 42). Requires: pathway-subtyping (>=0.8.0 line), numpy,
pandas, scikit-learn, requests. Network: www.cbioportal.org (public, no auth).
The deposited results/ snapshot lets you read the outcome without re-fetching;
cBioPortal is a live service, so re-fetched counts can drift slightly over time.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import requests
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

# Framework import: prefer an installed package; fall back to the repo src tree.
try:
    from pathway_subtyping.validation import ValidationGates
except ModuleNotFoundError:  # running from a source checkout
    _repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../.."))
    sys.path.insert(0, os.path.join(_repo, "src"))
    from pathway_subtyping.validation import ValidationGates

API = "https://www.cbioportal.org/api"
STUDY = "coadread_tcga_pan_can_atlas_2018"
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PANEL = os.path.join(_HERE, "../../../panels/hallmark_200genes.gmt")
S = requests.Session()
S.headers.update({"Accept": "application/json"})


def read_gmt(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = [g for g in f[2:] if g]
    return d


def post(url, payload, **params):
    r = S.post(f"{API}{url}", params=params, data=json.dumps(payload),
               headers={"Content-Type": "application/json"}, timeout=180)
    r.raise_for_status()
    return r.json()


def get(url, **params):
    r = S.get(f"{API}{url}", params=params, timeout=180)
    r.raise_for_status()
    return r.json()


def carrier(entrez_id, v600e=False):
    recs = post(f"/molecular-profiles/{STUDY}_mutations/mutations/fetch",
                {"entrezGeneIds": [entrez_id], "sampleListId": f"{STUDY}_sequenced"},
                projection="DETAILED")
    if v600e:
        return {m["sampleId"] for m in recs if "V600E" in (m.get("proteinChange") or "")}
    return {m["sampleId"] for m in recs}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=DEFAULT_PANEL)
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pw = read_gmt(args.panel)
    symbols = sorted({g for gs in pw.values() for g in gs})
    genes = post("/genes/fetch", symbols, geneIdType="HUGO_GENE_SYMBOL", projection="SUMMARY")
    sym2ent = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in genes}
    ent2sym = {v: k for k, v in sym2ent.items()}
    entrez = sorted(set(sym2ent.values()))

    prof = f"{STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores"
    md = post(f"/molecular-profiles/{prof}/molecular-data/fetch",
              {"entrezGeneIds": entrez, "sampleListId": f"{STUDY}_rna_seq_v2_mrna"},
              projection="SUMMARY")
    expr = (pd.DataFrame([(m["sampleId"], ent2sym.get(m["entrezGeneId"]), m["value"]) for m in md],
                         columns=["sample", "gene", "z"]).dropna()
            .pivot_table(index="sample", columns="gene", values="z", aggfunc="mean"))

    scores = {name: expr[[g for g in gs if g in expr.columns]].mean(axis=1)
              for name, gs in pw.items()
              if len([g for g in gs if g in expr.columns]) >= 2}
    P = pd.DataFrame(scores).dropna()
    X = StandardScaler().fit_transform(P.values)
    bic = {k: GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                              random_state=42, reg_covar=1e-6).fit(X).bic(X)
           for k in range(2, 7)}
    K = min(bic, key=bic.get)
    labels = pd.Series(
        GaussianMixture(n_components=K, covariance_type="full", n_init=10,
                        random_state=42, reg_covar=1e-6).fit(X).predict(X),
        index=P.index)

    braf = carrier(673, v600e=True)
    kras = carrier(3845)
    mantis = {c["sampleId"]: float(c["value"]) for c in
              get(f"/studies/{STUDY}/clinical-data", clinicalDataType="SAMPLE",
                  attributeId="MSI_SCORE_MANTIS", projection="SUMMARY")}
    ctype = {c["sampleId"]: c["value"] for c in
             get(f"/studies/{STUDY}/clinical-data", clinicalDataType="SAMPLE",
                 attributeId="CANCER_TYPE_DETAILED", projection="SUMMARY")}
    seq = set(get(f"/sample-lists/{STUDY}_sequenced")["sampleIds"])

    keep = [x for x in labels.index if x in seq]
    lab = labels.loc[keep].to_numpy()
    braf_s = np.array(["V600E" if x in braf else "wt" for x in keep])
    kras_s = np.array(["mut" if x in kras else "wt" for x in keep])
    msi_s = np.array(["MSI-H" if mantis.get(x, np.nan) > 0.4 else
                      ("MSS" if x in mantis else None) for x in keep], dtype=object)
    tissue = np.array(["COAD" if "Colon" in ctype.get(x, "") else
                       ("READ" if "Rectal" in ctype.get(x, "") else "other") for x in keep])

    g = ValidationGates(seed=42, show_progress=False)

    # A. the wiring
    g6 = g.confound_association_gate(lab, {"tissue": tissue}, cramers_v_max=0.30)
    g7 = g.somatic_anchoring_gate(lab, {"BRAF_V600E": braf_s, "KRAS": kras_s, "MSI": msi_s})

    # B. gate-PASS positive control (real strong BRAF<->MSI co-occurrence)
    known = [i for i, m in enumerate(msi_s) if m is not None]
    msi_part = np.array([0 if msi_s[i] == "MSS" else 1 for i in known])
    g7b = g.somatic_anchoring_gate(msi_part, {"BRAF_V600E": braf_s[known], "KRAS": kras_s[known]})

    print("=== TCGA-CRC somatic anchoring (real, cBioPortal) ===")
    print(f"n={len(keep)} sequenced tumors | k(BIC)={K} | "
          f"BRAF-V600E={int((braf_s=='V600E').sum())} KRAS={int((kras_s=='mut').sum())} "
          f"MSI-H={int((msi_s=='MSI-H').sum())}")
    print(f"A. Gate 6 (tissue confound): passed={g6.passed} V={g6.metric_value:.3f}")
    print(f"A. Gate 7 (somatic):         passed={g7.passed} anchored={g7.details['anchored_strata']} "
          f"max_V={g7.metric_value:.3f}")
    print(f"B. Gate 7 PASS control (MSI-partition vs BRAF): passed={g7b.passed} "
          f"anchored={g7b.details['anchored_strata']} max_V={g7b.metric_value:.3f}")

    out = {
        "study": STUDY, "source": "cBioPortal public API", "seed": 42,
        "n_sequenced_tumors": len(keep), "k_selected_by_bic": K,
        "bic_by_k": {k: round(v) for k, v in bic.items()},
        "strata_counts": {"BRAF_V600E": int((braf_s == "V600E").sum()),
                          "KRAS": int((kras_s == "mut").sum()),
                          "MSI_H": int((msi_s == "MSI-H").sum()),
                          "MSS": int((msi_s == "MSS").sum())},
        "A_wiring": {"gate6_confound_tissue": g6.to_dict(), "gate7_somatic": g7.to_dict()},
        "B_gate_pass_poscontrol_msi_vs_braf": g7b.to_dict(),
    }
    dest = os.path.join(args.out, "tcga_crc_somatic_result.json")
    json.dump(out, open(dest, "w"), indent=2, default=str)
    print(f"\nsaved -> {dest}")


if __name__ == "__main__":
    main()
