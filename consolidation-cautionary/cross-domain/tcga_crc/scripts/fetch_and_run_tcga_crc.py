#!/usr/bin/env python3
"""
Real positive control for Gate 7 (somatic) — TCGA-CRC BRAF / KRAS / MSI.

Wires ``ValidationGates.somatic_anchoring_gate`` against real Colorectal
Adenocarcinoma data (TCGA PanCancer Atlas, via the public cBioPortal API) with
actual BRAF-V600E / KRAS / MSI strata. The roadmap's Stage C, on real driver
calls. Two transcriptomic partitions are tested against the same real strata:

  A. CANONICAL CMS (the strong positive control).
     Uses the published Consensus Molecular Subtype labels (Guinney et al. 2015,
     CRC Subtyping Consortium public-tier file, vendored in ../inputs/) as the
     transcriptomic partition. CMS is expression-defined and independent of the
     somatic strata, so this is a genuine, non-circular test — and it reproduces
     textbook CRC biology: CMS1 (MSI-immune) is ~80% MSI-H / ~60% BRAF-V600E,
     CMS3 is ~71% KRAS-mutant. The somatic gate PASSES, anchoring on BRAF-V600E
     and MSI at Cramer's V ~ 0.66 (q ~ 1e-40). This is the somatic analog of the
     feature-level gate reproducing Voineagu's enrichment.

  B. UNBIASED PSF CLUSTERING (honest contrast).
     Builds a partition with the framework from a disease-agnostic Hallmark panel
     (pathway means -> Gaussian mixture, k by BIC). This generic whole-transcriptome
     partition couples to the drivers only modestly (significant but Cramer's V <
     0.30), so the gate correctly returns NOT-anchored. The contrast is the point:
     the gate faithfully reflects whether the partition captures driver-aligned
     biology (CMS does; a generic clustering does not) and does not over-call a
     significant-but-weak association. Recovering the driver-aligned structure
     de novo needs the dedicated CMS classifier, not generic clustering.

Deterministic (seed 42). Requires: pathway-subtyping (>=0.8.0 line), numpy,
pandas, scikit-learn, requests. Network: www.cbioportal.org (public, no auth).
The deposited results/ snapshot lets you read the outcome without re-fetching;
cBioPortal is a live service, so re-fetched counts can drift slightly over time.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
import pandas as pd
import requests
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

try:
    from pathway_subtyping.validation import ValidationGates
except ModuleNotFoundError:
    _repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../.."))
    sys.path.insert(0, os.path.join(_repo, "src"))
    from pathway_subtyping.validation import ValidationGates

API = "https://www.cbioportal.org/api"

# --- shared provenance recording (see consolidation-cautionary/scripts/_provenance.py) -
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "..", "..", "scripts"))
from _provenance import env_provenance, fetch_provenance  # noqa: E402

STUDY = "coadread_tcga_pan_can_atlas_2018"
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PANEL = os.path.join(_HERE, "../../../panels/hallmark_200genes.gmt")
DEFAULT_CMS = os.path.join(_HERE, "../inputs/cms_labels_public_all.txt")
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


def strata_for(samples, braf, kras, mantis):
    return (
        np.array(["V600E" if s in braf else "wt" for s in samples]),
        np.array(["mut" if s in kras else "wt" for s in samples]),
        np.array(["MSI-H" if mantis.get(s, np.nan) > 0.4 else
                  ("MSS" if s in mantis else None) for s in samples], dtype=object),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=DEFAULT_PANEL)
    ap.add_argument("--cms", default=DEFAULT_CMS)
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    g = ValidationGates(seed=42, show_progress=False)

    # --- shared: somatic strata + tissue from cBioPortal ---
    braf = carrier(673, v600e=True)
    kras = carrier(3845)
    mantis = {c["sampleId"]: float(c["value"]) for c in
              get(f"/studies/{STUDY}/clinical-data", clinicalDataType="SAMPLE",
                  attributeId="MSI_SCORE_MANTIS", projection="SUMMARY")}
    ctype = {c["sampleId"]: c["value"] for c in
             get(f"/studies/{STUDY}/clinical-data", clinicalDataType="SAMPLE",
                 attributeId="CANCER_TYPE_DETAILED", projection="SUMMARY")}
    seq = set(get(f"/sample-lists/{STUDY}_sequenced")["sampleIds"])
    tissue_of = lambda s: ("COAD" if "Colon" in ctype.get(s, "")
                           else "READ" if "Rectal" in ctype.get(s, "") else "other")

    # ================= A. canonical CMS partition (strong control) =================
    cms_rows = [r for r in csv.DictReader(open(args.cms), delimiter="\t")
                if r["dataset"] == "tcga"]
    cms = {r["sample"] + "-01": r["CMS_final_network_plus_RFclassifier_in_nonconsensus_samples"]
           for r in cms_rows}
    cms = {k: v for k, v in cms.items() if v != "NOLBL"}
    keepA = [s for s in cms if s in seq]
    labA = np.array([cms[s] for s in keepA])
    brafA, krasA, msiA = strata_for(keepA, braf, kras, mantis)
    tissueA = np.array([tissue_of(s) for s in keepA])
    g6A = g.confound_association_gate(labA, {"tissue": tissueA}, cramers_v_max=0.30)
    g7A = g.somatic_anchoring_gate(labA, {"BRAF_V600E": brafA, "KRAS": krasA, "MSI": msiA})

    # ================= B. unbiased PSF Hallmark clustering (contrast) =================
    pw = read_gmt(args.panel)
    symbols = sorted({x for gs in pw.values() for x in gs})
    genes = post("/genes/fetch", symbols, geneIdType="HUGO_GENE_SYMBOL", projection="SUMMARY")
    sym2ent = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in genes}
    ent2sym = {v: k for k, v in sym2ent.items()}
    md = post(f"/molecular-profiles/{STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores/molecular-data/fetch",
              {"entrezGeneIds": sorted(set(sym2ent.values())),
               "sampleListId": f"{STUDY}_rna_seq_v2_mrna"}, projection="SUMMARY")
    expr = (pd.DataFrame([(m["sampleId"], ent2sym.get(m["entrezGeneId"]), m["value"]) for m in md],
                         columns=["sample", "gene", "z"]).dropna()
            .pivot_table(index="sample", columns="gene", values="z", aggfunc="mean"))
    P = pd.DataFrame({n: expr[[x for x in gs if x in expr.columns]].mean(axis=1)
                      for n, gs in pw.items()
                      if len([x for x in gs if x in expr.columns]) >= 2}).dropna()
    X = StandardScaler().fit_transform(P.values)
    bic = {k: GaussianMixture(k, covariance_type="full", n_init=10, random_state=42,
                              reg_covar=1e-6).fit(X).bic(X) for k in range(2, 7)}
    K = min(bic, key=bic.get)
    lab = pd.Series(GaussianMixture(K, covariance_type="full", n_init=10, random_state=42,
                                    reg_covar=1e-6).fit(X).predict(X), index=P.index)
    keepB = [s for s in lab.index if s in seq]
    labB = lab.loc[keepB].to_numpy()
    brafB, krasB, msiB = strata_for(keepB, braf, kras, mantis)
    tissueB = np.array([tissue_of(s) for s in keepB])
    g6B = g.confound_association_gate(labB, {"tissue": tissueB}, cramers_v_max=0.30)
    g7B = g.somatic_anchoring_gate(labB, {"BRAF_V600E": brafB, "KRAS": krasB, "MSI": msiB})

    print("=== TCGA-CRC somatic anchoring (real, cBioPortal) ===")
    print(f"strata: BRAF-V600E={len(braf & seq)} KRAS={len(kras & seq)} "
          f"MSI-H={sum(1 for s in seq if mantis.get(s, 0) > 0.4)}")
    print(f"\nA. CANONICAL CMS partition (n={len(keepA)}): "
          f"CMS={dict(zip(*np.unique(labA, return_counts=True)))}")
    print(f"   Gate 6 (tissue): passed={g6A.passed} V={g6A.metric_value:.3f}")
    print(f"   Gate 7 (somatic): passed={g7A.passed} anchored={g7A.details['anchored_strata']} "
          f"maxV={g7A.metric_value:.3f}")
    for st, i in g7A.details["per_stratum"].items():
        print(f"      {st:12} V={i['cramers_v']:.3f} q={i['p_adjusted']:.2e} anchored={i['anchored']}")
    print(f"\nB. UNBIASED PSF clustering (n={len(keepB)}, k(BIC)={K}):")
    print(f"   Gate 6 (tissue): passed={g6B.passed} V={g6B.metric_value:.3f}")
    print(f"   Gate 7 (somatic): passed={g7B.passed} anchored={g7B.details['anchored_strata']} "
          f"maxV={g7B.metric_value:.3f}")

    def ct(lab_, strat_, order):
        return pd.crosstab(pd.Series(lab_, name="subtype"),
                           pd.Series(strat_, name="s"), normalize="index").round(3).to_dict("index")

    out = {
        "study": STUDY, "source": "cBioPortal public API", "seed": 42,
        "strata_counts": {"BRAF_V600E": len(braf & seq), "KRAS": len(kras & seq),
                          "MSI_H": sum(1 for s in seq if mantis.get(s, 0) > 0.4)},
        "A_canonical_cms": {
            "partition": "Consensus Molecular Subtype (Guinney 2015, consortium public labels)",
            "n": len(keepA), "cms_sizes": {k: int(v) for k, v in zip(*np.unique(labA, return_counts=True))},
            "gate6_confound_tissue": g6A.to_dict(), "gate7_somatic": g7A.to_dict(),
            "crosstabs_row_fraction": {"MSI": ct(labA, msiA, None),
                                        "BRAF_V600E": ct(labA, brafA, None),
                                        "KRAS": ct(labA, krasA, None)},
        },
        "B_unbiased_psf_clustering": {
            "partition": "PSF Hallmark pathway GMM", "n": len(keepB), "k_selected_by_bic": K,
            "bic_by_k": {k: round(v) for k, v in bic.items()},
            "gate6_confound_tissue": g6B.to_dict(), "gate7_somatic": g7B.to_dict(),
        },
    }
    dest = os.path.join(args.out, "tcga_crc_somatic_result.json")
    out["provenance"] = {"environment": env_provenance(),
                         "fetch": fetch_provenance(API, STUDY, matrix=P)}
    json.dump(out, open(dest, "w"), indent=2, default=str)
    print(f"\nsaved -> {dest}")


if __name__ == "__main__":
    main()
