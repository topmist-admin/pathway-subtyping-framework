#!/usr/bin/env python3
"""Within-region SCHIZOPHRENIA marker-panel check (symmetric to the ASD §2.5 analysis).

Rationale (§3.2): the ASD analysis rescued a signal with a TARGETED canonical
cell-type marker panel (PVALB/SST/VIP/...), covariate-adjusted, BH-corrected, where
a genome-wide FDR null was empty. This applies the IDENTICAL panel + method to
schizophrenia, within each of the 3 true regions (AnCg/DLPFC/nAcc), so ASD and SCZ
get symmetric treatment.

REPRODUCIBILITY: this script now builds everything from RAW public GSE80655 (the
series matrix carries age / post-mortem interval / brain pH / sex in its
characteristics, and the expression file is Ensembl-indexed) -- no processed
intermediates. Reuses reproduce_consolidation_paper.parse_series_matrix. Public
inputs only; deterministic; framework not required (numpy/pandas/scipy/statsmodels).

Group: Schizophrenia vs Control, within each region. Model: OLS panel-score ~
diagnosis + age + PMI + brain_pH + sex; partial Cohen's d = 2t/sqrt(dof); BH within
region across the 9 panels. Canonical result (public GSE80655):
  AnCg : SST d=-1.38 q=2e-4 ; pan-GABA d=-1.20 q=6e-4 ; PVALB n.s.
  DLPFC: SST d=-1.54 q=4e-5 ; pan-GABA d=-0.67 q~0.06 ; PVALB n.s.
  nAcc : nothing survives (MSN-dominated, smaller n).
=> SST + pan-GABAergic depletion in both cortical regions, PVALB n.s., nAcc null.
"""
import argparse, os, sys, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import t as tdist
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reproduce_consolidation_paper as R  # noqa: E402  (parse_series_matrix)

# Ensembl IDs for the identical §2.5 canonical panel
ENS = {"PVALB": "ENSG00000100362", "SST": "ENSG00000157005", "VIP": "ENSG00000146469",
       "GAD1": "ENSG00000128683", "GAD2": "ENSG00000136750", "SLC32A1": "ENSG00000101438",
       "RORB": "ENSG00000198963", "PDGFRA": "ENSG00000134853", "MBP": "ENSG00000197971",
       "PLP1": "ENSG00000123560", "GFAP": "ENSG00000131095", "AQP4": "ENSG00000171885",
       "CX3CR1": "ENSG00000168329"}
PANELS = {"PV_interneuron": ["PVALB"], "SST_interneuron": ["SST"], "VIP_interneuron": ["VIP"],
          "GABAergic_pan": ["GAD1", "GAD2", "SLC32A1"], "L4_IT_excitatory": ["RORB"],
          "OPC": ["PDGFRA"], "Oligodendrocyte": ["MBP", "PLP1"],
          "Astrocyte_ctrl": ["GFAP", "AQP4"], "Microglia_ctrl": ["CX3CR1"]}
REGIONS = ["AnCg", "DLPFC", "nAcc"]


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def bh(pv):
    p = np.asarray(pv, float); o = np.argsort(p); m = len(p); q = np.empty(m); prev = 1.0
    for r, i in enumerate(o[::-1]):
        prev = min(prev, p[i] * m / (m - r)); q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gse80655-expr", required=True, help="GSE80655 raw gene expression (Ensembl x sample), .txt.gz")
    ap.add_argument("--gse80655-matrix", required=True, help="GSE80655_series_matrix.txt.gz")
    ap.add_argument("--out", default="outputs_scz_marker_panel")
    ap.add_argument("--seed", type=int, default=20260708)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

    expr = pd.read_csv(a.gse80655_expr, sep="\t", index_col=0)
    expr.index = expr.index.str.replace(r"\.\d+$", "", regex=True)
    m = R.parse_series_matrix(a.gse80655_matrix); m["sid"] = m["title"].str.extract(r"(SL\d+)")
    m = m[m["sid"].isin(expr.columns)].set_index("sid")
    cols = [c for c in expr.columns if c in m.index]; m = m.loc[cols]
    logcpm = np.log2(expr[cols].div(expr[cols].sum(axis=0), axis=1) * 1e6 + 1.0)  # Ensembl x samples
    cov = pd.DataFrame({"age": pd.to_numeric(m["age at death"], errors="coerce"),
                        "pmi": pd.to_numeric(m["post-mortem interval"], errors="coerce"),
                        "ph": pd.to_numeric(m["brain ph"], errors="coerce"),
                        "sex_M": m["gender"].astype(str).str.upper().str.startswith("M").astype(float)},
                       index=m.index)
    missing = sorted({g for gs in PANELS.values() for g in gs if ENS.get(g) not in logcpm.index})

    def panel_score(sids, genes):
        rows = [ENS[g] for g in genes if ENS.get(g) in logcpm.index]
        sub = logcpm.loc[rows, sids]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, 1e-9), axis=0)
        return z.mean(axis=0)

    def run_region(region):
        idx = [i for i in m.index[(m["brain region"] == region) &
               (m["clinical diagnosis"].isin(["Schizophrenia", "Control"]))] if cov.loc[i].notna().all()]
        dx = (m.loc[idx, "clinical diagnosis"] == "Schizophrenia").astype(int).values
        C = cov.loc[idx]; rows = []
        for name, genes in PANELS.items():
            s = panel_score(idx, genes).values
            X = np.column_stack([np.ones(len(idx)), dx, C["age"], C["pmi"], C["ph"], C["sex_M"]])
            beta, _, _, _ = np.linalg.lstsq(X, s, rcond=None)
            resid = s - X @ beta; dof = len(idx) - X.shape[1]
            se = np.sqrt((resid @ resid) / dof * np.linalg.pinv(X.T @ X)[1, 1])
            tval = beta[1] / se if se > 0 else np.nan
            rows.append(dict(region=region, n_SCZ=int((dx == 1).sum()), n_ctrl=int((dx == 0).sum()),
                             panel=name, partial_d_adj=round(2 * tval / np.sqrt(dof), 3),
                             p_adj=float(2 * tdist.sf(abs(tval), dof))))
        df = pd.DataFrame(rows); df["q_adj_BH"] = bh(df["p_adj"].values)
        return df

    res = pd.concat([run_region(r) for r in REGIONS], ignore_index=True)
    res.to_csv(f"{a.out}/scz_within_region_marker_panel_results.csv", index=False)
    json.dump({"generated_from": "raw public GSE80655", "expr_sha": sha(a.gse80655_expr),
               "matrix_sha": sha(a.gse80655_matrix), "missing_panel_genes": missing,
               "model": "panel ~ dx + age + pmi + brain_pH + sex", "seed": a.seed},
              open(f"{a.out}/provenance.json", "w"), indent=2)
    pd.set_option("display.width", 200)
    if missing: print("missing panel genes:", missing)
    for r in REGIONS:
        sub = res[res.region == r]
        print(f"\n### {r} (SCZ {sub.n_SCZ.iloc[0]} / ctrl {sub.n_ctrl.iloc[0]}) ###")
        print(sub[["panel", "partial_d_adj", "p_adj", "q_adj_BH"]].to_string(index=False))
    key = res[res.panel.isin(["SST_interneuron", "GABAergic_pan", "PV_interneuron"])]
    print("\nHeadline (SST / pan-GABA / PVALB by region):")
    print(key.pivot(index="region", columns="panel", values="partial_d_adj").to_string())
    print("Outputs ->", a.out)


if __name__ == "__main__":
    main()
