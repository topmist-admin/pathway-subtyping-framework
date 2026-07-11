#!/usr/bin/env python3
"""Figure: the PVALB/SST double dissociation (§3.4).

Two panels of the per-sample *PVALB* − *SST* log-expression contrast, by group:
  A) ASD (GSE64018): cases fall BELOW controls -> PVALB depleted more than SST.
  B) SCZ (GSE80655 cortex): cases rise ABOVE controls -> SST depleted more than PVALB.
The within-sample contrast cancels platform/cohort baselines; the case-vs-control
shift in it IS the marker x diagnosis interaction. Opposite signs = double dissociation.
Annotated d / p are the covariate-adjusted interaction estimates from
double_dissociation_interaction_test.py. Public GEO inputs only.
"""
import argparse, os, sys, gzip, io, numpy as np, pandas as pd, statsmodels.api as sm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reproduce_consolidation_paper as R  # noqa: E402
PV, SST = "ENSG00000100362", "ENSG00000157005"


def adj_interaction(delta, design):
    fit = sm.OLS(delta, design, missing="drop").fit()
    return fit.params[1] / np.sqrt(fit.mse_resid), fit.pvalues[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True, help="root with GSE64018/, GSE80655/")
    ap.add_argument("--out", default="figures/figure_double_dissociation")
    a = ap.parse_args(); os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    DD = a.data_dir

    # ASD GSE64018
    m = R.parse_series_matrix(f"{DD}/GSE64018/GSE64018_series_matrix.txt.gz")
    cnt = pd.read_csv(f"{DD}/GSE64018/GSE64018_countlevel_12asd_12ctl.txt.gz", sep="\t", index_col=0)
    cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True)
    m["col"] = [next((c for c in cnt.columns if b in c), None) for b in m["brainid"]]
    m = m.set_index("col").loc[cnt.columns]
    dx = (m["diagnosis"] == "ASD").astype(float).values
    male = (m["Sex"] == "M").astype(float).values
    covs = {c: pd.to_numeric(m[c], errors="coerce").fillna(pd.to_numeric(m[c], errors="coerce").median()).values
            for c in ["age", "rin", "pmi"]}
    logcpm = np.log2(cnt.div(cnt.sum(axis=0), axis=1) * 1e6 + 1.0)
    dA = (logcpm.loc[PV] - logcpm.loc[SST]).values
    desA = np.column_stack([np.ones(len(m)), dx, covs["age"], male, covs["rin"], covs["pmi"]])
    dA_d, dA_p = adj_interaction(dA, desA)
    asd = pd.DataFrame({"contrast": dA, "grp": np.where(dx == 1, "ASD", "Control")})

    # SCZ GSE80655 cortex
    e80 = pd.read_csv(f"{DD}/GSE80655/GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz", sep="\t", index_col=0)
    e80.index = e80.index.str.replace(r"\.\d+$", "", regex=True)
    m80 = R.parse_series_matrix(f"{DD}/GSE80655/GSE80655_series_matrix.txt.gz")
    m80["sid"] = m80["title"].str.extract(r"(SL\d+)")
    sc = m80[(m80["clinical diagnosis"].isin(["Schizophrenia", "Control"])) &
             (m80["brain region"].isin(["AnCg", "DLPFC"]))].copy()
    sc = sc[sc["sid"].isin(e80.columns)]; cols = sc["sid"].tolist()
    lc = np.log2(e80[cols].div(e80[cols].sum(axis=0), axis=1) * 1e6 + 1.0)
    sc = sc.set_index("sid").loc[cols]
    dxS = (sc["clinical diagnosis"] == "Schizophrenia").astype(float).values
    dlpfc = (sc["brain region"] == "DLPFC").astype(float).values
    maleS = (sc["gender"] == "M").astype(float).values
    ageS = pd.to_numeric(sc["age at death"], errors="coerce").fillna(
        pd.to_numeric(sc["age at death"], errors="coerce").median()).values
    dS = (lc.loc[PV] - lc.loc[SST]).values
    desS = np.column_stack([np.ones(len(sc)), dxS, dlpfc, maleS, ageS])
    dS_d, dS_p = adj_interaction(dS, desS)
    scz = pd.DataFrame({"contrast": dS, "grp": np.where(dxS == 1, "SCZ", "Control")})

    # ---- plot ----
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.0), sharey=False)
    panels = [("A  Autism (GSE64018)", asd, ["Control", "ASD"], dA_d, dA_p, "#4C72B0"),
              ("B  Schizophrenia (GSE80655)", scz, ["Control", "SCZ"], dS_d, dS_p, "#C44E52")]
    for ax, (ttl, df, order, dval, pval, ccase) in zip(axes, panels):
        for i, g in enumerate(order):
            v = df.loc[df.grp == g, "contrast"].values
            ax.scatter(np.full(len(v), i) + rng.uniform(-.08, .08, len(v)), v, s=26,
                       color=("#999999" if g == "Control" else ccase), alpha=.75,
                       edgecolor="white", linewidth=.5, zorder=3)
            ax.plot([i - .22, i + .22], [v.mean()] * 2, color="black", lw=2, zorder=4)
            ax.errorbar(i, v.mean(), yerr=v.std(ddof=1) / np.sqrt(len(v)), color="black", capsize=4, zorder=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels(order)
        ax.set_title(ttl, fontsize=10, loc="left", fontweight="bold")
        star = "****" if pval < 1e-4 else ("*" if pval < .05 else "n.s.")
        ax.text(.5, .96, f"interaction $d$ = {dval:+.2f}, $p$ = {pval:.1e}  {star}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8.5)
        ax.axhline(df.contrast.mean(), color="#cccccc", ls=":", lw=.8, zorder=1)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("per-sample  $PVALB - SST$  (log$_2$ CPM)")
    fig.suptitle("Double dissociation: PVALB depleted in autism, SST in schizophrenia",
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=200, bbox_inches="tight")
    print(f"ASD interaction d={dA_d:+.2f} p={dA_p:.2e} | SCZ interaction d={dS_d:+.2f} p={dS_p:.2e}")
    print("saved ->", f"{a.out}.png / .pdf")


if __name__ == "__main__":
    main()
