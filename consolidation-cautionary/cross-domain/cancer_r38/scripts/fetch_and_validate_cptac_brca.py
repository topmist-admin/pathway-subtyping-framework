#!/usr/bin/env python3
"""CPTAC-BRCA multi-omic subtyping — reviewer points R1.6 (multi-omic) + R3.8.

R1.6: "although VCF burden scoring is implemented, the validation is mainly based
on gene expression data." This closes the multi-omic gap on the PROTEIN side:
CPTAC breast (Krug et al., Cell 2020) provides matched mRNA AND mass-spec protein
quantification for the same tumors, both public via cBioPortal (no auth).

We score Hallmark pathways from each modality independently, subtype each, and:
  1. recover PAM50 from the expression-based and the protein-based partitions;
  2. measure expression<->protein subtype concordance (do the two modalities
     agree on the partition?);
  3. run the discreteness + stability gates on each modality.

This shows PSF's pathway subtyping operates on protein data, not only expression
(directly answering R1.6), and that the validation gates are modality-general.

Deterministic (seed 42). Network: www.cbioportal.org. Requires pathway-subtyping
(>=0.8 line), numpy, pandas, scikit-learn, requests.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import requests
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

_HERE = os.path.dirname(os.path.abspath(__file__))
try:
    from pathway_subtyping.validation import ValidationGates
    from pathway_subtyping.discreteness import DiscretenessGateA
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../..")))
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../../src")))
    from pathway_subtyping.validation import ValidationGates
    from pathway_subtyping.discreteness import DiscretenessGateA

API = "https://www.cbioportal.org/api"
STUDY = "brca_cptac_2020"
PROFILE = {"mrna": f"{STUDY}_mrna_median_Zscores",
           "protein": f"{STUDY}_protein_quantification"}
DEFAULT_PANEL = os.path.join(_HERE, "../../../panels/hallmark_200genes.gmt")

S = requests.Session()
S.headers.update({"Accept": "application/json"})


def post(url, payload, **params):
    r = S.post(f"{API}{url}", params=params, data=json.dumps(payload),
               headers={"Content-Type": "application/json"}, timeout=300)
    r.raise_for_status()
    return r.json()


def get(url, **params):
    r = S.get(f"{API}{url}", params=params, timeout=300)
    r.raise_for_status()
    return r.json()


def read_gmt(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = [g for g in f[2:] if g]
    return d


def pathway_scores(panel, profile_id):
    """Hallmark pathway-mean matrix for one modality (mRNA z or protein log2)."""
    pw = read_gmt(panel)
    symbols = sorted({x for gs in pw.values() for x in gs})
    genes = post("/genes/fetch", symbols, geneIdType="HUGO_GENE_SYMBOL",
                 projection="SUMMARY")
    sym2ent = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in genes}
    ent2sym = {v: k for k, v in sym2ent.items()}
    md = post(f"/molecular-profiles/{profile_id}/molecular-data/fetch",
              {"entrezGeneIds": sorted(set(sym2ent.values())),
               "sampleListId": f"{STUDY}_all"}, projection="SUMMARY")
    expr = (pd.DataFrame([(m["sampleId"], ent2sym.get(m["entrezGeneId"]), m["value"])
                          for m in md], columns=["sample", "gene", "v"]).dropna()
            .pivot_table(index="sample", columns="gene", values="v", aggfunc="mean"))
    P = pd.DataFrame({n: expr[[x for x in gs if x in expr.columns]].mean(axis=1)
                      for n, gs in pw.items()
                      if len([x for x in gs if x in expr.columns]) >= 2}).dropna()
    return P


def subtype_labels():
    rows = get(f"/studies/{STUDY}/clinical-data", clinicalDataType="SAMPLE",
               attributeId="PAM50", projection="SUMMARY")
    return pd.Series({r["sampleId"]: r["value"] for r in rows if r.get("value")})


def partition(P, k):
    X = StandardScaler().fit_transform(P.values)
    lab = GaussianMixture(k, covariance_type="full", n_init=10, random_state=42,
                          reg_covar=1e-6).fit(X).predict(X)
    return pd.Series(lab, index=P.index), X


def recovery(labels, truth):
    idx = truth.dropna().index.intersection(labels.index)
    if len(idx) < 10:
        return None
    t = truth.loc[idx].astype("category").cat.codes.to_numpy()
    return {"ari": float(adjusted_rand_score(t, labels.loc[idx])),
            "ami": float(adjusted_mutual_info_score(t, labels.loc[idx])),
            "n": int(len(idx))}


def gate_summary(P, lab, k, n_ref):
    g = ValidationGates(seed=42, show_progress=False)
    stab = g.stability_test_bootstrap(P, lab.to_numpy(), k)
    a = DiscretenessGateA(seed=42, n_ref=n_ref).run(STUDY, P, k, gmm_seed=42)
    return {"bootstrap_stability": {"passed": bool(stab.passed), "ari": float(stab.metric_value)},
            "discreteness_gateA": {"certified": bool(a.testable and a.passed), "verdict": a.verdict}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default=DEFAULT_PANEL)
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    ap.add_argument("--n-ref", type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pam50 = subtype_labels()
    k = max(2, int(pam50.dropna().nunique()))
    print(f"CPTAC-BRCA ({STUDY}); PAM50 classes={k}")

    out = {"study": STUDY, "k": k, "modalities": {}}
    parts = {}
    for mod, prof in PROFILE.items():
        print(f"  fetching {mod} pathway scores ({prof})...")
        P = pathway_scores(args.panel, prof)
        lab, _ = partition(P, k)
        parts[mod] = lab
        out["modalities"][mod] = {
            "n_samples": int(P.shape[0]),
            "recovery_vs_pam50": recovery(lab, pam50),
            "gates": gate_summary(P, lab, k, args.n_ref),
        }
        print(f"    {mod}: n={P.shape[0]}, recovery={out['modalities'][mod]['recovery_vs_pam50']}")

    # expression <-> protein subtype concordance (shared samples)
    shared = parts["mrna"].index.intersection(parts["protein"].index)
    conc = (float(adjusted_rand_score(parts["mrna"].loc[shared], parts["protein"].loc[shared]))
            if len(shared) >= 10 else None)
    out["expression_protein_concordance_ari"] = conc
    out["n_shared_samples"] = int(len(shared))

    with open(os.path.join(args.out, "cptac_brca_multiomic.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    print("\n=== CPTAC-BRCA multi-omic (R1.6 + R3.8) ===")
    for mod in PROFILE:
        m = out["modalities"][mod]
        print(f"  {mod:8s}: recovery {m['recovery_vs_pam50']}  gates {m['gates']}")
    print(f"  expression<->protein concordance ARI = {conc}  (n_shared={len(shared)})")
    print(f"\nWrote {args.out}/cptac_brca_multiomic.json")


if __name__ == "__main__":
    main()
