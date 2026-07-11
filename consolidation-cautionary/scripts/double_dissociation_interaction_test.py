#!/usr/bin/env python3
"""Double-dissociation test (§3.4): marker x diagnosis interaction.

Per-sample (PVALB - SST) log-expression contrast, tested case-vs-control within
each cohort by covariate-adjusted OLS. The diagnosis coefficient on the contrast
IS the marker x diagnosis interaction (which marker is depleted more), and cancels
platform/cohort baselines because PVALB and SST are measured in the same samples.
A double dissociation requires the interaction to be significant and OPPOSITE in
sign across ASD vs SCZ.

Result (2026-07-10):
  ASD GSE64018 (RNA-seq, adj)  d=-1.18  p=0.018   PVALB depleted > SST
  ASD GSE28521 (array,   adj)  d=-0.33  p=0.213   same direction, n.s.
  SCZ GSE80655 (cortex,  adj)  d=+1.01  p=9.5e-6  SST depleted > PVALB
  => opposite signs, both primary cohorts significant => double dissociation supported.

Inputs: same public GEO files as reproduce_consolidation_paper.py (framework:
`pip install pathway-subtyping==0.7.0`). Reuses that sibling module's
series-matrix / annotation parsers. Deterministic (no clustering; OLS).

Usage: python double_dissociation_interaction_test.py --data-dir <root>
where <root> contains GSE64018/, GSE28521/, GSE80655/ with the public GEO files.
"""
import argparse, os, sys, gzip, io, numpy as np, pandas as pd, statsmodels.api as sm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling helpers
import reproduce_consolidation_paper as R  # noqa: E402

PV, SST = "ENSG00000100362", "ENSG00000157005"
_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--data-dir", required=True,
                 help="root holding GSE64018/, GSE28521/, GSE80655/ public GEO files")
DD = _ap.parse_args().data_dir


def ols_delta(delta, design, label, n):
    fit = sm.OLS(delta, design, missing="drop").fit()
    b, p = fit.params[1], fit.pvalues[1]
    d = b / np.sqrt(fit.mse_resid)
    who = "PVALB more depleted" if b < 0 else "SST more depleted"
    print(f"{label:28s} n={n:3d}  beta={b:+.3f}  d={d:+.2f}  p={p:.2e}  {who}")
    return b, p


# GSE64018 ASD count model
m = R.parse_series_matrix(f"{DD}/GSE64018/GSE64018_series_matrix.txt.gz")
cnt = pd.read_csv(f"{DD}/GSE64018/GSE64018_countlevel_12asd_12ctl.txt.gz", sep="\t", index_col=0)
cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True)
m["col"] = [next((c for c in cnt.columns if b in c), None) for b in m["brainid"]]
m = m.set_index("col").loc[cnt.columns]
m["dx"] = (m["diagnosis"] == "ASD").astype(float); m["male"] = (m["Sex"] == "M").astype(float)
for c in ["age", "rin", "pmi"]:
    m[c] = pd.to_numeric(m[c], errors="coerce"); m[c] = m[c].fillna(m[c].median())
logcpm = np.log2(cnt.div(cnt.sum(axis=0), axis=1) * 1e6 + 1.0)
dlt = (logcpm.loc[PV] - logcpm.loc[SST]).values
des = np.column_stack([np.ones(len(m)), m["dx"], m["age"], m["male"], m["rin"], m["pmi"]])
b_asd1, p_asd1 = ols_delta(dlt, des, "ASD GSE64018 (count,adj)", len(m))

# GSE28521 ASD array
raw = gzip.open(f"{DD}/GSE28521/GSE28521_series_matrix.txt.gz", "rt").read()
tbl = raw.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
e28 = pd.read_csv(io.StringIO(tbl), sep="\t", index_col=0)
e28.index = e28.index.astype(str).str.replace('"', ""); e28.columns = [c.replace('"', '') for c in e28.columns]
m28 = R.parse_series_matrix(f"{DD}/GSE28521/GSE28521_series_matrix.txt.gz")
smp = pd.DataFrame({"geo": m28["geo"].values, "title": m28["title"].values})
smp["dx"] = (smp["title"].str[0] == "A").astype(float)
smp["region"] = smp["title"].str[-1].map({"C": "Cerebellum", "F": "Frontal", "T": "Temporal"})
smp = smp.set_index("geo"); smp = smp.loc[[c for c in e28.columns if c in smp.index]]
cx = smp[smp["region"].isin(["Frontal", "Temporal"])].copy(); cx["temporal"] = (cx["region"] == "Temporal").astype(float)
ann = gzip.open(f"{DD}/GSE28521/GPL6883.annot.gz", "rt", errors="ignore").read().split("\n")
hidx = next(i for i, l in enumerate(ann) if l.split("\t")[0] == "ID"); scol = ann[hidx].split("\t").index("Gene symbol")
sym2p = {}
for l in ann[hidx + 1:]:
    f = l.split("\t")
    if len(f) > scol and f[0].startswith("ILMN"):
        sym2p.setdefault(f[scol], []).append(f[0])
pm = lambda sym: e28.loc[[p for p in sym2p.get(sym, []) if p in e28.index], cx.index].astype(float).mean(axis=0).values
dlt28 = pm("PVALB") - pm("SST")
des28 = np.column_stack([np.ones(len(cx)), cx["dx"], cx["temporal"]])
b_asd2, p_asd2 = ols_delta(dlt28, des28, "ASD GSE28521 (array)", len(cx))

# GSE80655 SCZ cortex
e80 = pd.read_csv(f"{DD}/GSE80655/GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz", sep="\t", index_col=0)
e80.index = e80.index.str.replace(r"\.\d+$", "", regex=True)
m80 = R.parse_series_matrix(f"{DD}/GSE80655/GSE80655_series_matrix.txt.gz")
m80["sid"] = m80["title"].str.extract(r"(SL\d+)")
sc = m80[(m80["clinical diagnosis"].isin(["Schizophrenia", "Control"])) & (m80["brain region"].isin(["AnCg", "DLPFC"]))].copy()
sc = sc[sc["sid"].isin(e80.columns)]; cols = sc["sid"].tolist()
lc = np.log2(e80[cols].div(e80[cols].sum(axis=0), axis=1) * 1e6 + 1.0)
sc = sc.set_index("sid").loc[cols]
sc["dx"] = (sc["clinical diagnosis"] == "Schizophrenia").astype(float)
sc["dlpfc"] = (sc["brain region"] == "DLPFC").astype(float); sc["male"] = (sc["gender"] == "M").astype(float)
sc["age"] = pd.to_numeric(sc["age at death"], errors="coerce"); sc["age"] = sc["age"].fillna(sc["age"].median())
dlt80 = (lc.loc[PV] - lc.loc[SST]).values
des80 = np.column_stack([np.ones(len(sc)), sc["dx"], sc["dlpfc"], sc["male"], sc["age"]])
b_scz, p_scz = ols_delta(dlt80, des80, "SCZ GSE80655 (cortex,adj)", len(sc))

opp = (b_asd1 < 0) and (b_scz > 0); sig = (p_asd1 < 0.05) and (p_scz < 0.05)
print(f"\nDouble dissociation supported (opposite signs & both primary cohorts sig): {opp and sig}")
