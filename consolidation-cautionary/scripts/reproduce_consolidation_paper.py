#!/usr/bin/env python3
"""
Deposit-ready reproduction script for the consolidation / cautionary paper

    "Stable but confounded: brain region, not diagnosis, drives a
     validation-passing molecular subtype in postmortem psychiatric brain"

Regenerates every manuscript number for both figures from public inputs, using
the framework's own machinery (ssGSEA, GMM/BIC, and the first-class Confound
Association Gate added in pathway_subtyping.validation).

FIGURE 1 (the artifact).  For each diagnosis-vs-Control contrast in GSE80655, a
k=3 partition of ssGSEA pathway scores is scored by the bootstrap stability gate
and the confound association gate. In every contrast the confound gate FAILS on
brain REGION with diagnosis non-significant; the three single-diagnosis
contrasts (SCZ/BD/MDD vs Control) also PASS bootstrap stability (ARI ~0.85-0.91),
while the pooled all-4-diagnosis case falls below the stability threshold
(ARI ~0.71) because pooling fragments the region structure -- the region-confound
verdict holds regardless. The region check is run through
ValidationGates.confound_association_gate -- the gate whose absence let the
artifact through.

FIGURE 2 (the recovered biology).  A covariate-adjusted case-control marker test
-- run directly on the diagnosis labels, never on a cluster -- recovers real
disease signal: PVALB interneuron depletion in ASD (GSE64018, GSE28521) and SST
depletion in SCZ (GSE80655), a double dissociation (ASD spares SST; SCZ spares
PVALB).

PANEL PARAMETERISATION.  The ssGSEA pathway panel used for Figure 1 is supplied
via --gmt. It defaults to the public schizophrenia panel shipped in the repo
(data/pathways/schizophrenia_pathways.gmt). Supplying the manuscript's own
deposited panel regenerates the exact published Figure-1 numbers; the public
panel reproduces them to within reconstruction tolerance (same anatomy axis,
Cramer's V ~0.67). Figure 2 uses a fixed interneuron-marker gene list and does
not depend on --gmt.

All steps are seeded (20260708) and deterministic given the same inputs and
package versions. A provenance.json (input SHA-256s, seed, versions, results) is
written to --out.

Public inputs
-------------
  GSE80655 : GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz  (+ series matrix)
  GSE64018 : GSE64018_countlevel_12asd_12ctl.txt.gz,
             GSE64018_adjfpkm_12asd_12ctl.txt.gz                    (+ series matrix)
  GSE28521 : GSE28521_series_matrix.txt.gz  (+ GPL6883.annot.gz)
  panel    : schizophrenia_pathways.gmt (or the manuscript's deposited panel)
  Ensembl  : rest.ensembl.org (symbol -> Ensembl for the GSE80655 panel)

Research use only. Not for clinical decision-making.
"""
import argparse
import datetime
import gzip
import hashlib
import io
import json
import logging
import time

import numpy as np
import pandas as pd
import requests
import statsmodels.api as sm
from scipy.stats import chi2_contingency, ttest_ind

from pathway_subtyping.clustering import (
    ClusteringAlgorithm,
    run_clustering,
    select_n_clusters,
)
from pathway_subtyping.expression import (
    ExpressionScoringMethod,
    score_pathways_from_expression,
)
from pathway_subtyping.statistical_rigor import benjamini_hochberg
from pathway_subtyping.validation import ValidationGates, cramers_v

logging.getLogger("pathway_subtyping").setLevel(logging.ERROR)
SEED = 20260708

# Interneuron marker panel for Figure 2 (Ensembl, GRCh38) + specificity controls.
MARKERS = {
    "PVALB": "ENSG00000100362", "SST": "ENSG00000157005", "VIP": "ENSG00000146469",
    "GAD1": "ENSG00000128683", "GAD2": "ENSG00000136750", "CALB1": "ENSG00000104327",
    "CALB2": "ENSG00000172137", "RELN": "ENSG00000189056", "NPY": "ENSG00000122585",
    "SLC17A7": "ENSG00000104888", "GAPDH": "ENSG00000111640", "ACTB": "ENSG00000075624",
}
INTERNEURON = ["PVALB", "SST", "VIP", "GAD1", "GAD2", "CALB1", "CALB2", "RELN", "NPY"]


# ---------------------------------------------------------------- helpers
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_series_matrix(path):
    titles = geo = None
    char = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                char.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    meta = {}
    for cl in char:
        key = cl[0].split(":")[0].strip()
        meta[key] = [c.split(":", 1)[1].strip() if ":" in c else "" for c in cl]
    df = pd.DataFrame(meta)
    df.insert(0, "title", titles)
    df.insert(1, "geo", geo)
    return df


def parse_gmt(path):
    pw = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            pw[p[0]] = [g for g in p[2:] if g]
    return pw


def map_symbols_to_ensembl(symbols):
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    out = {}
    for i in range(0, len(symbols), 900):
        r = requests.post(url, headers=headers,
                          data=json.dumps({"symbols": symbols[i:i + 900]}))
        r.raise_for_status()
        for s, rec in r.json().items():
            if isinstance(rec, dict) and rec.get("id", "").startswith("ENSG"):
                out[s] = rec["id"]
        time.sleep(0.3)
    return out


def score_and_cluster(expr, cols, ens_pathways, k=3):
    """log2(CPM+1) -> ssGSEA -> GMM/BIC. Returns (pathway_scores, labels, bic_k)."""
    lib = expr[cols].sum(axis=0)
    logcpm = np.log2(expr[cols].div(lib, axis=1) * 1e6 + 1.0)
    logcpm = logcpm.loc[logcpm.sum(axis=1) > 0].T
    scored = score_pathways_from_expression(
        gene_expression=logcpm, pathways=ens_pathways,
        method=ExpressionScoringMethod.SSGSEA, min_genes_per_pathway=2,
        alpha=0.25, seed=SEED, show_progress=False)
    ps = scored.pathway_scores.loc[cols]
    ms = select_n_clusters(ps.values, list(range(2, 7)), method="bic", seed=SEED)
    cl = run_clustering(ps.values, k, ClusteringAlgorithm.GMM, seed=SEED)
    return ps, cl.labels, ms.optimal_k


# ---------------------------------------------------------------- Figure 1
def figure1_region_gate(args, ens_pathways):
    expr = pd.read_csv(args.gse80655_expr, sep="\t", index_col=0)
    expr.index = expr.index.str.replace(r"\.\d+$", "", regex=True)
    meta = parse_series_matrix(args.gse80655_matrix)
    meta["sample_id"] = meta["title"].str.extract(r"(SL\d+)")

    contrasts = [
        ("Schizophrenia", "SCZ+Control"),
        ("Bipolar Disorder", "BD+Control"),
        ("Major Depression", "MDD+Control"),
        (None, "All-4-dx"),
    ]
    rows = []
    for dx, name in contrasts:
        keep = (["Schizophrenia", "Bipolar Disorder", "Major Depression", "Control"]
                if dx is None else [dx, "Control"])
        sub = meta[meta["clinical diagnosis"].isin(keep)]
        sub = sub[sub["sample_id"].isin(expr.columns)]
        cols = sub["sample_id"].tolist()
        ps, labels, bic_k = score_and_cluster(expr, cols, ens_pathways, k=3)
        sub = sub.set_index("sample_id").loc[cols]
        regions = sub["brain region"].values
        diagnosis = sub["clinical diagnosis"].values

        gates = ValidationGates(seed=SEED, n_bootstrap=100, show_progress=False)
        stab = gates.stability_test_bootstrap(ps, labels, n_clusters=3, gmm_seed=SEED)
        gate = gates.confound_association_gate(
            labels, {"region": regions, "diagnosis": diagnosis})
        pc = gate.details["per_confound"]
        rows.append(dict(
            contrast=name, n=len(cols), bic_k=int(bic_k),
            stability_ari=round(float(stab.metric_value), 3),
            stability_pass=bool(stab.passed),
            region_chi2=pc["region"]["chi2"], region_cramers_v=pc["region"]["cramers_v"],
            region_p=pc["region"]["p_value"],
            diagnosis_cramers_v=pc["diagnosis"]["cramers_v"],
            diagnosis_p=pc["diagnosis"]["p_value"],
            confound_gate_pass=bool(gate.passed),
            gate_fails_on=";".join(gate.details["failing_confounds"]),
        ))
    return pd.DataFrame(rows).set_index("contrast")


# ---------------------------------------------------------------- Figure 2
def _adj_marker_test(logexpr, meta, formula_cols, dx_col, markers):
    """OLS log-expr ~ dx + covariates; adjusted Cohen's d = beta_dx / sqrt(mse)."""
    rows = []
    dxv = meta[dx_col].values.astype(float)
    cov = np.column_stack([meta[c].values.astype(float) for c in formula_cols]) \
        if formula_cols else np.empty((len(meta), 0))
    design = np.column_stack([np.ones(len(meta)), dxv, cov]) if cov.size \
        else np.column_stack([np.ones(len(meta)), dxv])
    for g, ens in markers.items():
        if ens not in logexpr.index:
            continue
        y = logexpr.loc[ens].reindex(meta.index).values.astype(float)
        fit = sm.OLS(y, design, missing="drop").fit()
        rows.append(dict(gene=g, adj_d=fit.params[1] / np.sqrt(fit.mse_resid),
                         p=fit.pvalues[1],
                         mean_case=float(y[dxv == 1].mean()),
                         mean_ctl=float(y[dxv == 0].mean())))
    res = pd.DataFrame(rows).set_index("gene")
    res["q"] = benjamini_hochberg(res["p"].values)
    return res


def figure2_markers(args):
    out = {}

    # --- GSE64018 (ASD, RNA-seq STG): count-level model ---
    m = parse_series_matrix(args.gse64018_matrix)
    cnt = pd.read_csv(args.gse64018_counts, sep="\t", index_col=0)
    cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True)
    m["col"] = [next((c for c in cnt.columns if b in c), None) for b in m["brainid"]]
    m = m.set_index("col").loc[cnt.columns]
    m["dx"] = (m["diagnosis"] == "ASD").astype(float)
    m["male"] = (m["Sex"] == "M").astype(float)
    for c in ["age", "rin", "pmi"]:
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m = m.assign(age=m["age"].fillna(m["age"].median()),
                 rin=m["rin"].fillna(m["rin"].median()),
                 pmi=m["pmi"].fillna(m["pmi"].median()))
    logcpm = np.log2(cnt.div(cnt.sum(axis=0), axis=1) * 1e6 + 1.0)
    out["GSE64018_asd_countmodel"] = _adj_marker_test(
        logcpm, m, ["age", "male", "rin", "pmi"], "dx", MARKERS)

    # --- GSE64018 published adjusted-FPKM (author covariate model) ---
    adj = pd.read_csv(args.gse64018_adjfpkm, sep="\t", index_col=0)
    adj.index = adj.index.str.replace(r"\.\d+$", "", regex=True)
    dx_by_bid = m.reset_index().set_index("brainid")["diagnosis"].to_dict()
    acols = {c: dx_by_bid.get(c.split("_")[0]) for c in adj.columns}
    asd = [c for c, d in acols.items() if d == "ASD"]
    ctl = [c for c, d in acols.items() if d == "CTL"]
    frows = []
    for g, ens in MARKERS.items():
        if ens not in adj.index:
            continue
        a = adj.loc[ens, asd].astype(float).values
        c = adj.loc[ens, ctl].astype(float).values
        sp = np.sqrt(((len(a) - 1) * a.std(ddof=1) ** 2 +
                      (len(c) - 1) * c.std(ddof=1) ** 2) / (len(a) + len(c) - 2))
        frows.append(dict(gene=g, d=(a.mean() - c.mean()) / sp,
                          p=ttest_ind(a, c)[1]))
    fpk = pd.DataFrame(frows).set_index("gene")
    fpk["q"] = benjamini_hochberg(fpk["p"].values)
    out["GSE64018_asd_pubadjfpkm"] = fpk

    # --- GSE28521 (ASD, array cortex) ---
    raw = gzip.open(args.gse28521_matrix, "rt").read()
    tbl = raw.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    expr28 = pd.read_csv(io.StringIO(tbl), sep="\t", index_col=0)
    expr28.index = expr28.index.astype(str).str.replace('"', "")
    expr28.columns = [c.replace('"', "") for c in expr28.columns]
    m28 = parse_series_matrix(args.gse28521_matrix)
    smp = pd.DataFrame({"geo": m28["geo"].values, "title": m28["title"].values})
    smp["dx"] = (smp["title"].str[0] == "A").astype(float)
    smp["region"] = smp["title"].str[-1].map({"C": "Cerebellum", "F": "Frontal", "T": "Temporal"})
    smp = smp.set_index("geo")
    smp = smp.loc[[c for c in expr28.columns if c in smp.index]]
    cortex = smp[smp["region"].isin(["Frontal", "Temporal"])].copy()
    cortex["temporal"] = (cortex["region"] == "Temporal").astype(float)
    # probe -> symbol from GPL6883 annotation
    ann = gzip.open(args.gse28521_annot, "rt", errors="ignore").read().split("\n")
    hidx = next(i for i, l in enumerate(ann) if l.split("\t")[0] == "ID")
    hdr = ann[hidx].split("\t")
    scol = hdr.index("Gene symbol")
    sym2probes = {}
    for l in ann[hidx + 1:]:
        f = l.split("\t")
        if len(f) > scol and f[0].startswith("ILMN"):
            sym2probes.setdefault(f[scol], []).append(f[0])
    rows = []
    dxv = cortex["dx"].values
    # UNADJUSTED case-vs-control (matches §2.5: "GSE28521 was tested unadjusted").
    # Pooled frontal+temporal cortex; no region covariate. `unadj_d` = Cohen's d.
    design = np.column_stack([np.ones(len(cortex)), dxv])
    for g in list(MARKERS):
        probes = [p for p in sym2probes.get(g, []) if p in expr28.index]
        if not probes:
            continue
        y = expr28.loc[probes, cortex.index].astype(float).mean(axis=0).reindex(cortex.index).values
        if np.isnan(y).all():
            continue
        fit = sm.OLS(y, design, missing="drop").fit()
        rows.append(dict(gene=g, unadj_d=fit.params[1] / np.sqrt(fit.mse_resid), p=fit.pvalues[1]))
    r28 = pd.DataFrame(rows).set_index("gene")
    r28["q"] = benjamini_hochberg(r28["p"].values)
    out["GSE28521_asd_cortex"] = r28

    # --- GSE80655 (SCZ pole, cortex) ---
    expr80 = pd.read_csv(args.gse80655_expr, sep="\t", index_col=0)
    expr80.index = expr80.index.str.replace(r"\.\d+$", "", regex=True)
    m80 = parse_series_matrix(args.gse80655_matrix)
    m80["sample_id"] = m80["title"].str.extract(r"(SL\d+)")
    sc = m80[(m80["clinical diagnosis"].isin(["Schizophrenia", "Control"])) &
             (m80["brain region"].isin(["AnCg", "DLPFC"]))].copy()
    sc = sc[sc["sample_id"].isin(expr80.columns)]
    cols = sc["sample_id"].tolist()
    logcpm80 = np.log2(expr80[cols].div(expr80[cols].sum(axis=0), axis=1) * 1e6 + 1.0)
    sc = sc.set_index("sample_id").loc[cols]
    sc["dx"] = (sc["clinical diagnosis"] == "Schizophrenia").astype(float)
    sc["dlpfc"] = (sc["brain region"] == "DLPFC").astype(float)
    sc["male"] = (sc["gender"] == "M").astype(float)
    sc["age"] = pd.to_numeric(sc["age at death"], errors="coerce")
    sc["age"] = sc["age"].fillna(sc["age"].median())
    out["GSE80655_scz_cortex"] = _adj_marker_test(
        logcpm80, sc, ["dlpfc", "male", "age"], "dx",
        {g: MARKERS[g] for g in ["PVALB", "SST", "VIP", "GAD1", "ACTB", "GAPDH"]})
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gse80655-expr", required=True)
    ap.add_argument("--gse80655-matrix", required=True)
    ap.add_argument("--gse64018-counts", required=True)
    ap.add_argument("--gse64018-adjfpkm", required=True)
    ap.add_argument("--gse64018-matrix", required=True)
    ap.add_argument("--gse28521-matrix", required=True)
    ap.add_argument("--gse28521-annot", required=True)
    ap.add_argument("--gmt", required=True,
                    help="ssGSEA pathway panel for Figure 1. Default = public "
                         "schizophrenia panel; supply the manuscript's deposited "
                         "panel to regenerate the exact published numbers.")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    # Figure 1: translate the panel to Ensembl (GSE80655 is Ensembl-keyed).
    pathways = parse_gmt(args.gmt)
    symbols = sorted({g for gs in pathways.values() for g in gs})
    sym2ens = map_symbols_to_ensembl(symbols)
    expr80_idx = pd.read_csv(args.gse80655_expr, sep="\t", index_col=0, usecols=[0]).index
    expr80_idx = expr80_idx.str.replace(r"\.\d+$", "", regex=True)
    ens_pathways = {
        pw: [sym2ens[g] for g in genes if g in sym2ens and sym2ens[g] in set(expr80_idx)]
        for pw, genes in pathways.items()
    }

    print("=== FIGURE 1: region artifact across diagnoses (confound gate) ===")
    fig1 = figure1_region_gate(args, ens_pathways)
    print(fig1.to_string())
    fig1.to_csv(f"{args.out}/figure1_region_gate_all_diagnoses.csv")

    print("\n=== FIGURE 2: recovered biology (covariate-adjusted marker test) ===")
    fig2 = figure2_markers(args)
    for name, tbl in fig2.items():
        print(f"\n--- {name} ---")
        print(tbl.round(4).to_string())
        tbl.to_csv(f"{args.out}/figure2_{name}.csv")

    # Provenance
    prov = {
        "script": "reproduce_consolidation_paper.py",
        "date": datetime.date.today().isoformat(),
        "seed": SEED,
        "panel_gmt": args.gmt,
        "panel_symbols_mapped": f"{len(sym2ens)}/{len(symbols)}",
        "inputs_sha256": {
            "gse80655_expr": sha256(args.gse80655_expr),
            "gse80655_matrix": sha256(args.gse80655_matrix),
            "gse64018_counts": sha256(args.gse64018_counts),
            "gse64018_adjfpkm": sha256(args.gse64018_adjfpkm),
            "gse64018_matrix": sha256(args.gse64018_matrix),
            "gse28521_matrix": sha256(args.gse28521_matrix),
            "gse28521_annot": sha256(args.gse28521_annot),
            "gmt": sha256(args.gmt),
        },
        "packages": {"statsmodels": sm.__version__, "numpy": np.__version__,
                     "pandas": pd.__version__},
        "figure1": json.loads(fig1.reset_index().to_json(orient="records")),
        "figure2": {k: json.loads(v.reset_index().to_json(orient="records"))
                    for k, v in fig2.items()},
    }
    with open(f"{args.out}/provenance.json", "w") as f:
        json.dump(prov, f, indent=2)
    print(f"\nProvenance written to {args.out}/provenance.json")


if __name__ == "__main__":
    main()
