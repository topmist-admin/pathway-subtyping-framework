#!/usr/bin/env python3
"""
Cross-diagnosis region-gate reproduction (manuscript Table 2, section 3.3).

For each GSE80655 diagnosis-vs-Control contrast (SCZ / BD / MDD) and the pooled
4-diagnosis cohort, cluster ssGSEA pathway scores and report the region and
diagnosis association at BOTH a fixed k=3 (cross-contrast comparability) and the
BIC-selected k, to show the region artifact is invariant to cluster number while
stability is not.

Framework machinery (ssGSEA, GMM/BIC, confound gate) comes from the installed
`pathway_subtyping` package -- install the PUBLIC release: `pip install
pathway-subtyping==0.7.0`. Shared GEO/GMT parsers are imported from the sibling
`reproduce_consolidation_paper.py` in this same folder.

Seed 20260708. Deterministic given inputs + package versions.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling helpers
import reproduce_consolidation_paper as R  # noqa: E402
from pathway_subtyping.clustering import (  # noqa: E402
    ClusteringAlgorithm, run_clustering, select_n_clusters,
)
from pathway_subtyping.expression import (  # noqa: E402
    ExpressionScoringMethod, score_pathways_from_expression,
)
from pathway_subtyping.validation import ValidationGates  # noqa: E402

SEED = 20260708


def cv(labels, groups, bc=False):
    ct = pd.crosstab(pd.Series(labels), pd.Series(groups))
    chi2, p, _, _ = chi2_contingency(ct.values, correction=False)
    n = ct.values.sum(); r, k = ct.shape
    if bc:
        phi2 = chi2 / n
        phi2c = max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
        rc = r - (r - 1) ** 2 / (n - 1); kc = k - (k - 1) ** 2 / (n - 1)
        v = np.sqrt(phi2c / max(min(kc - 1, rc - 1), 1e-12))
    else:
        v = np.sqrt((chi2 / n) / (min(r, k) - 1))
    return chi2, p, v


def score(expr, cols, ens_pathways):
    lib = expr[cols].sum(axis=0)
    logcpm = np.log2(expr[cols].div(lib, axis=1) * 1e6 + 1.0)
    logcpm = logcpm.loc[logcpm.sum(axis=1) > 0].T
    scored = score_pathways_from_expression(
        gene_expression=logcpm, pathways=ens_pathways,
        method=ExpressionScoringMethod.SSGSEA, min_genes_per_pathway=2,
        alpha=0.25, seed=SEED, show_progress=False)
    return scored.pathway_scores.loc[cols]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gse80655-expr", required=True)
    ap.add_argument("--gse80655-matrix", required=True)
    ap.add_argument("--gmt", required=True, help="curated schizophrenia panel (.gmt)")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    pathways = R.parse_gmt(args.gmt)
    symbols = sorted({g for gs in pathways.values() for g in gs})
    sym2ens = R.map_symbols_to_ensembl(symbols)
    expr = pd.read_csv(args.gse80655_expr, sep="\t", index_col=0)
    expr.index = expr.index.str.replace(r"\.\d+$", "", regex=True)
    ens_pathways = {pw: [sym2ens[g] for g in genes
                         if g in sym2ens and sym2ens[g] in set(expr.index)]
                    for pw, genes in pathways.items()}
    meta = R.parse_series_matrix(args.gse80655_matrix)
    meta["sample_id"] = meta["title"].str.extract(r"(SL\d+)")

    contrasts = [("Schizophrenia", "SCZ+Control"), ("Bipolar Disorder", "BD+Control"),
                 ("Major Depression", "MDD+Control"), (None, "All-4-dx")]
    rows = []
    for dx, name in contrasts:
        keep = (["Schizophrenia", "Bipolar Disorder", "Major Depression", "Control"]
                if dx is None else [dx, "Control"])
        sub = meta[meta["clinical diagnosis"].isin(keep)]
        sub = sub[sub["sample_id"].isin(expr.columns)]
        cols = sub["sample_id"].tolist()
        ps = score(expr, cols, ens_pathways)
        sub = sub.set_index("sample_id").loc[cols]
        region = sub["brain region"].values
        dxv = sub["clinical diagnosis"].values
        bic_k = select_n_clusters(ps.values, list(range(2, 7)), method="bic", seed=SEED).optimal_k
        for kmode, kval in [("k3", 3), ("bicK", int(bic_k))]:
            lab = run_clustering(ps.values, kval, ClusteringAlgorithm.GMM, seed=SEED).labels
            g = ValidationGates(seed=SEED, n_bootstrap=100, show_progress=False)
            stab = g.stability_test_bootstrap(ps, lab, n_clusters=kval, gmm_seed=SEED)
            rc, rp, rv = cv(lab, region); _, _, rvc = cv(lab, region, True)
            _, dp, dv = cv(lab, dxv)
            rows.append(dict(contrast=name, n=len(cols), k_mode=kmode, k=kval, bic_k=int(bic_k),
                             stability_ari=round(float(stab.metric_value), 3),
                             stab_pass=bool(stab.passed),
                             region_chi2=round(rc, 1), region_V=round(rv, 3),
                             region_Vbc=round(rvc, 3), region_p=f"{rp:.2e}",
                             diagnosis_V=round(dv, 3), diagnosis_p=round(dp, 3)))
    out = pd.DataFrame(rows)
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print(out.to_string(index=False))
    out.to_csv(f"{args.out}/figure1_crossdiagnosis_k_robustness.csv", index=False)
    print(f"\nSEED {SEED} | panel {os.path.basename(args.gmt)} | mapped {len(sym2ens)}/{len(symbols)}")


if __name__ == "__main__":
    main()
