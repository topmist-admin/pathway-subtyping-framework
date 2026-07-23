#!/usr/bin/env python3
"""Brain-region confound at scale — GTEx BRAIN (n~2931), the large-N answer to
the reviewers' small-N criticism (R1.5/R3.5/R3.6/R3.9).

Reads the recount3-exported pathway matrix (from fetch_gtex_brain_recount3.R) and
tests, on ~2931 public neurotypical brain samples across 13 subregions, whether
pathway-level structure is dominated by BRAIN REGION and whether that structure is
discreteness-gate-certified. If yes at this scale, it establishes that brain
region is a massive, real, discrete molecular axis — so any diagnosis-based
subtyping on mixed-region psychiatric tissue is at severe risk of recovering
region, not disease. Scales the "stable but confounded" finding from n=281 to
n~2931.

Deterministic (seed 42). Input: results/gtex_brain_pathway_scores.tsv.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
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


def region_enrichment(labels, region):
    """Best-cluster enrichment per region (metric-consistent with the cancer runs)."""
    from scipy.stats import fisher_exact
    best = {"enrichment_frac": -1}
    for reg in region.unique():
        iss = (region == reg).to_numpy()
        for cl in np.unique(labels):
            inc = (labels == cl)
            if inc.sum() == 0:
                continue
            frac = float((iss & inc).sum() / inc.sum())
            if frac > best["enrichment_frac"]:
                a, b = int((iss & inc).sum()), int((~iss & inc).sum())
                c, d = int((iss & ~inc).sum()), int((~iss & ~inc).sum())
                orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                best = {"region": str(reg), "enrichment_frac": round(frac, 3),
                        "odds_ratio": round(float(orr), 2), "p": float(p)}
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=os.path.join(_HERE, "../results/gtex_brain_pathway_scores.tsv"))
    ap.add_argument("--out", default=os.path.join(_HERE, "../results"))
    ap.add_argument("--n-ref", type=int, default=120)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(args.scores, sep="\t")
    region = df["region"].astype(str)
    P = df.drop(columns=["sample", "region"]).set_index(df["sample"])
    X = StandardScaler().fit_transform(P.values)
    n_reg = region.nunique()
    print(f"GTEx BRAIN: {P.shape[0]} samples x {P.shape[1]} pathways; {n_reg} regions")

    # cluster at k = #regions (recovery test) and at BIC-selected k
    k = n_reg
    lab_k = GaussianMixture(k, covariance_type="full", n_init=8, random_state=42,
                            reg_covar=1e-6).fit(X).predict(X)
    bic = {kk: GaussianMixture(kk, covariance_type="full", n_init=3, random_state=42,
                               reg_covar=1e-6).fit(X).bic(X) for kk in range(2, min(15, n_reg + 3))}
    k_bic = min(bic, key=bic.get)
    lab_bic = GaussianMixture(k_bic, covariance_type="full", n_init=8, random_state=42,
                              reg_covar=1e-6).fit(X).predict(X)

    reg_codes = region.astype("category").cat.codes.to_numpy()
    out = {
        "n_samples": int(P.shape[0]), "n_regions": int(n_reg),
        "k_region": int(k), "k_bic": int(k_bic),
        "region_recovery_at_k": {
            "ari": round(float(adjusted_rand_score(reg_codes, lab_k)), 3),
            "ami": round(float(adjusted_mutual_info_score(reg_codes, lab_k)), 3),
            "best_region_enrichment": region_enrichment(lab_k, region)},
        "region_recovery_at_bic_k": {
            "ari": round(float(adjusted_rand_score(reg_codes, lab_bic)), 3),
            "best_region_enrichment": region_enrichment(lab_bic, region)},
    }

    # discreteness gate at both k: is the region structure certified at scale?
    ga = DiscretenessGateA(seed=42, n_ref=args.n_ref).run("GTEx_brain", P, k, gmm_seed=42)
    ga_bic = DiscretenessGateA(seed=42, n_ref=args.n_ref).run("GTEx_brain", P, k_bic, gmm_seed=42)
    out["discreteness_gateA_at_k"] = {"certified": bool(ga.testable and ga.passed), "verdict": ga.verdict}
    out["discreteness_gateA_at_bic_k"] = {"k": int(k_bic),
                                          "certified": bool(ga_bic.testable and ga_bic.passed),
                                          "verdict": ga_bic.verdict}

    with open(os.path.join(args.out, "gtex_brain_region_confound.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"\n=== GTEx BRAIN region-confound at scale (n={P.shape[0]}) ===")
    print(f"  region recovery @k={k}: ARI={out['region_recovery_at_k']['ari']} "
          f"| best {out['region_recovery_at_k']['best_region_enrichment']}")
    print(f"  region recovery @BIC-k={k_bic}: ARI={out['region_recovery_at_bic_k']['ari']}")
    print(f"  discreteness gate @k: {out['discreteness_gateA_at_k']}")
    print(f"  discreteness gate @BIC-k: {out['discreteness_gateA_at_bic_k']}")
    print(f"\nWrote {args.out}/gtex_brain_region_confound.json")


if __name__ == "__main__":
    main()
