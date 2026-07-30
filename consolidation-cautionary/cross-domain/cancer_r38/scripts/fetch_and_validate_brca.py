#!/usr/bin/env python3
"""Large-cohort validation on TCGA-BRCA / PAM50 (reviewer R3.8, cancer arm).

Independent large-cohort test of pathway-level subtyping on public data — no
controlled access required. Fetches Hallmark-panel mRNA and PAM50 subtype labels
for TCGA-BRCA from the public cBioPortal API (no auth), builds a PSF partition,
and reports:

  1. Recovery of PAM50 (ARI/AMI of PSF partition vs the published PAM50 calls),
     alongside classical (k-means) and deep (DEC, VAE-GMM) baselines — one table
     that answers both the R2.2 method comparison and the R3.8 recovery ask.
  2. The full validation-gate battery on the PSF partition: label shuffle,
     random gene set, bootstrap stability, and the v0.8.0 discreteness gate.
     BRCA (~1000) is large enough that Gate A is well-powered — the contrast with
     the small-n psychiatric caution.

Deterministic (seed 42). Network: www.cbioportal.org (public). Requires
pathway-subtyping (>=0.8 line), numpy, pandas, scikit-learn, requests; torch
optional (DL baselines skipped without it).

⚠️ DRAFT — study ID and the PAM50 clinical attribute id must be confirmed against
cBioPortal before trusting the numbers (see VERIFY notes). Run:
    python fetch_and_validate_brca.py --out ../results
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import requests
from sklearn.cluster import KMeans
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

# --- shared provenance recording (see consolidation-cautionary/scripts/_provenance.py) -
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "..", "..", "scripts"))
from _provenance import env_provenance, fetch_provenance  # noqa: E402

STUDY = "brca_tcga_pan_can_atlas_2018"          # VERIFY on cBioPortal
SUBTYPE_ATTR = "SUBTYPE"                          # VERIFY: PAM50 calls (e.g. "BRCA_LumA")
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


def fetch_pathway_scores(panel):
    """Hallmark-panel mRNA z-scores from cBioPortal -> pathway-mean matrix."""
    pw = read_gmt(panel)
    symbols = sorted({x for gs in pw.values() for x in gs})
    genes = post("/genes/fetch", symbols, geneIdType="HUGO_GENE_SYMBOL",
                 projection="SUMMARY")
    sym2ent = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in genes}
    ent2sym = {v: k for k, v in sym2ent.items()}
    md = post(f"/molecular-profiles/{STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores/molecular-data/fetch",
              {"entrezGeneIds": sorted(set(sym2ent.values())),
               "sampleListId": f"{STUDY}_rna_seq_v2_mrna"}, projection="SUMMARY")
    expr = (pd.DataFrame([(m["sampleId"], ent2sym.get(m["entrezGeneId"]), m["value"])
                          for m in md], columns=["sample", "gene", "z"]).dropna()
            .pivot_table(index="sample", columns="gene", values="z", aggfunc="mean"))
    P = pd.DataFrame({n: expr[[x for x in gs if x in expr.columns]].mean(axis=1)
                      for n, gs in pw.items()
                      if len([x for x in gs if x in expr.columns]) >= 2}).dropna()
    return P


def fetch_pam50():
    """PAM50 subtype ('SUBTYPE') is a PATIENT-level attribute in this study
    (values BRCA_LumA/LumB/Basal/Her2/Normal). Return {patientId: value}."""
    rows = get(f"/studies/{STUDY}/clinical-data", clinicalDataType="PATIENT",
               attributeId=SUBTYPE_ATTR, projection="SUMMARY")
    return {r["patientId"]: r["value"] for r in rows if r.get("value")}


def _patient_of(sample_id):
    """TCGA sample 'TCGA-3C-AAAU-01' -> patient 'TCGA-3C-AAAU'."""
    return "-".join(sample_id.split("-")[:3])


def recovery(labels, truth_series):
    """ARI/AMI of a partition vs the published subtype on shared samples."""
    idx = truth_series.dropna().index.intersection(labels.index)
    if len(idx) < 10:
        return None
    t = truth_series.loc[idx].astype("category").cat.codes.to_numpy()
    p = labels.loc[idx].to_numpy()
    return {"ari": float(adjusted_rand_score(t, p)),
            "ami": float(adjusted_mutual_info_score(t, p)),
            "n": int(len(idx))}


def single_subtype_enrichment(labels, truth_series):
    """Best-cluster enrichment per subtype, matching how the manuscript reported
    CMS4 (fraction of the most-enriched cluster that is the subtype + Fisher OR/p).
    This is the metric COMPARABLE to the paper's '75.9% CMS4' headline; the k-way
    ARI above is a harder, stricter test. Reporting BOTH is the honest fix."""
    from scipy.stats import fisher_exact
    idx = truth_series.dropna().index.intersection(labels.index)
    t = truth_series.loc[idx]
    lab = labels.loc[idx]
    per = {}
    for sub in sorted(t.unique()):
        is_sub = (t == sub).to_numpy()
        best = None
        for cl in sorted(lab.unique()):
            in_cl = (lab == cl).to_numpy()
            if in_cl.sum() == 0:
                continue
            frac = float((is_sub & in_cl).sum() / in_cl.sum())      # % of cluster that is sub
            a = int((is_sub & in_cl).sum()); b = int((~is_sub & in_cl).sum())
            c = int((is_sub & ~in_cl).sum()); d = int((~is_sub & ~in_cl).sum())
            try:
                orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            except Exception:
                orr, p = float("nan"), float("nan")
            cand = {"cluster": int(cl), "enrichment_frac": round(frac, 3),
                    "odds_ratio": round(float(orr), 2), "p": float(p)}
            if best is None or cand["enrichment_frac"] > best["enrichment_frac"]:
                best = cand
        per[str(sub)] = best
    best_sub = max(per, key=lambda s: per[s]["enrichment_frac"])
    return {"per_subtype": per, "best_subtype": best_sub, "best": per[best_sub]}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=DEFAULT_PANEL)
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    ap.add_argument("--n-ref", type=int, default=200)
    ap.add_argument("--no-dl", action="store_true", help="skip DEC/VAE baselines")
    ap.add_argument("--skip-gates", action="store_true",
                    help="skip the heavy Gate A (fast recovery-metric recompute)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Fetching TCGA-BRCA pathway scores from cBioPortal ({STUDY})...")
    P = fetch_pathway_scores(args.panel)
    X = StandardScaler().fit_transform(P.values)
    # map PAM50 (patient-level) onto the expression samples
    pam50_by_patient = fetch_pam50()
    pam50 = pd.Series({s: pam50_by_patient.get(_patient_of(s)) for s in P.index})
    n_pam = pam50.dropna().nunique()
    print(f"  pathway matrix {P.shape}; PAM50-labelled samples = "
          f"{int(pam50.notna().sum())}; classes = {n_pam}")

    # PSF partition: GMM at k = number of PAM50 classes (fair recovery test)
    k = max(2, n_pam) if n_pam else 5
    gmm = GaussianMixture(k, covariance_type="full", n_init=10, random_state=42,
                          reg_covar=1e-6).fit(X)
    lab_psf = pd.Series(gmm.predict(X), index=P.index)

    # framework's OWN k by BIC (so the stability result is not an artifact of
    # forcing k=5 to match PAM50)
    bic = {kk: GaussianMixture(kk, covariance_type="full", n_init=5, random_state=42,
                               reg_covar=1e-6).fit(X).bic(X) for kk in range(2, 8)}
    k_bic = min(bic, key=bic.get)
    lab_bic = pd.Series(GaussianMixture(k_bic, covariance_type="full", n_init=10,
                        random_state=42, reg_covar=1e-6).fit(X).predict(X), index=P.index)

    # --- recovery vs PAM50: PSF + baselines ---
    methods = {"PSF (pathway-GMM)": lab_psf,
               "k-means (pathway)": pd.Series(
                   KMeans(k, n_init=10, random_state=42).fit_predict(X), index=P.index)}
    if not args.no_dl:
        try:
            from pathway_subtyping.clustering_dl import run_dec, run_vae_gmm
            methods["DEC"] = pd.Series(run_dec(X, k, seed=42)[0], index=P.index)
            methods["VAE-GMM"] = pd.Series(run_vae_gmm(X, k, seed=42)[0], index=P.index)
        except Exception as e:
            print(f"  DL baselines skipped: {e}")

    rec = {name: recovery(lab, pam50) for name, lab in methods.items()}
    # single-subtype enrichment (metric-consistent with the manuscript's CMS4 report)
    enrich = {name: single_subtype_enrichment(lab, pam50) for name, lab in methods.items()}

    # --- validation gates: at forced k (=PAM50) AND at the framework's BIC k ---
    if args.skip_gates:
        gates = "skipped"
    else:
        g = ValidationGates(seed=42, show_progress=False)
        ns = g.negative_control_label_shuffle(P, lab_psf.to_numpy(), k)
        stab_k = g.stability_test_bootstrap(P, lab_psf.to_numpy(), k)
        stab_bic = g.stability_test_bootstrap(P, lab_bic.to_numpy(), k_bic)
        gateA = DiscretenessGateA(seed=42, n_ref=args.n_ref).run("BRCA", P, k, gmm_seed=42)
        gateA_bic = DiscretenessGateA(seed=42, n_ref=args.n_ref).run("BRCA", P, k_bic, gmm_seed=42)
        gates = {
            "label_shuffle": {"passed": bool(ns.passed), "ari": float(ns.metric_value)},
            "bootstrap_stability_at_k5": {"passed": bool(stab_k.passed), "ari": float(stab_k.metric_value)},
            "bootstrap_stability_at_bic_k": {"k": int(k_bic), "passed": bool(stab_bic.passed),
                                             "ari": float(stab_bic.metric_value)},
            "discreteness_gateA_at_k5": {"certified": bool(gateA.testable and gateA.passed),
                                         "verdict": gateA.verdict},
            "discreteness_gateA_at_bic_k": {"k": int(k_bic),
                                            "certified": bool(gateA_bic.testable and gateA_bic.passed),
                                            "verdict": gateA_bic.verdict},
        }

    out = {"study": STUDY, "n_samples": int(P.shape[0]), "k_pam50": int(k),
           "k_bic": int(k_bic), "pam50_classes": int(n_pam),
           "recovery_vs_pam50_ARI": rec,
           "recovery_vs_pam50_single_subtype_enrichment": enrich, "gates": gates}
    out["provenance"] = {"environment": env_provenance(),
                         "fetch": fetch_provenance(API, STUDY, n_samples=out.get("n_samples"))}
    with open(os.path.join(args.out, "brca_pam50_validation.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    print("\n=== TCGA-BRCA / PAM50 recovery (R3.8 + R2.2) ===")
    print("  k-way ARI (strict) | best-subtype enrichment (CMS4-comparable):")
    for name in rec:
        r = rec[name]; e = enrich[name]["best"]; bs = enrich[name]["best_subtype"]
        ari = "n/a" if r is None else f"ARI={r['ari']:.3f}"
        print(f"  {name:20s}: {ari:11s} | {bs} {e['enrichment_frac']*100:.1f}% "
              f"(OR={e['odds_ratio']}, p={e['p']:.1e})")
    print(f"\n=== Validation gates (forced k={k} vs BIC k={k_bic}) ===")
    if gates == "skipped":
        print("  (skipped)")
    else:
        for gname, gv in gates.items():
            print(f"  {gname}: {gv}")
    print(f"\nWrote {args.out}/brca_pam50_validation.json")


if __name__ == "__main__":
    main()
