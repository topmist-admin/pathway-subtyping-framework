#!/usr/bin/env python3
"""JOB 4 — donor-level inference on a CONTINUOUS statistic (Paulus Points 2, items 14-16).

THE PROBLEM. The deposited donor-level test collapses each donor to a single majority
subtype and computes Cramer's V on the resulting 3x2 table at n=48. Two consequences:

  (a) DEGENERATE NULL. Under donor permutation the table can only take a small number of
      configurations, so the permutation distribution of V is a handful of atoms. A p-value
      read off it is coarse and its resolution is an artifact of the discretisation, not of
      the evidence. The manuscript's p=0.234 inherits that.

  (b) INFORMATION DISCARDED. 44 of 48 donors (92%) have samples in more than one cluster,
      so the majority rule throws away the composition of almost every donor. One donor
      (X4296, 2 samples) has no majority at all and is assigned by whichever label pandas
      `value_counts()` happens to order first — an undocumented tie-break deciding a data
      point in an n=48 analysis.

THE REPLACEMENT. Keep the unit of analysis at the donor, but carry a continuous quantity:
each donor's SUBTYPE COMPOSITION, the fraction of that donor's samples in each cluster.
Nothing is discarded and no tie-break is needed.

Three statistics, all under the same donor-permutation null (diagnosis permuted across
donors, 10,000 draws, seed 42):

  T1  AUC of donor cluster-1 fraction against diagnosis        (univariate, continuous)
  T2  LOO-CV AUC of logistic regression on the full composition (multivariate, continuous)
  T3  GEE logistic on SAMPLES, subtype ~ diagnosis, exchangeable working correlation
      clustered by donor, robust SE                            (uses all 141, respects nesting)

An earlier version of T3 used `BinomialBayesMixedGLM.fit_vb` and formed a likelihood-ratio
statistic from the fitted objects. That is invalid: `fit_vb` returns variational-Bayes
pseudo-likelihoods, and differences between them are not chi-square-distributed. It is
replaced by GEE, which is a proper estimator here and needs no likelihood. The two happened
to agree to three decimals (p=0.0616 vs 0.0617), but the agreement was luck, not validation.
Because n=48 donors is small for the GEE asymptotics, T3 also carries a donor-permutation
empirical p on |z| and that is the one to quote.

T3 is the estimand Paulus named. T1/T2 are reported alongside because a mixed model at
n=48 donors is itself fragile, and three statistics agreeing is worth more than one.

Also emitted, to document (a) directly: the SIZE OF THE SUPPORT of the deposited V
statistic under the same permutations — how many distinct values it can take at all.

Deterministic. No network. Reads the deposited partition; writes nothing to it.
"""
from __future__ import annotations
import argparse, json, os, warnings
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, chi2 as chi2_dist
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut

SEED, N_PERM = 42, 10000

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
LABELS = os.path.join(REPO, "consolidation-cautionary/data/partition/"
                            "sample_metadata_with_subtypes.csv")


def cramers_v(ct):
    if min(ct.shape) < 2 or ct.sum() == 0:
        return 0.0
    chi2 = chi2_contingency(ct, correction=False)[0]
    return float(np.sqrt((chi2 / ct.sum()) / (min(ct.shape) - 1)))


def loo_auc(X, y, seed=SEED):
    if len(np.unique(y)) < 2:
        return np.nan
    pred = np.empty(len(y), float)
    for tr, te in LeaveOneOut().split(X):
        if len(np.unique(y[tr])) < 2:
            pred[te] = y[tr].mean()
            continue
        m = LogisticRegression(max_iter=2000, random_state=seed).fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
    return float(roc_auc_score(y, pred))


def gee_dx(samples, dx_by_donor, seed=SEED):
    """GEE logistic, exchangeable by donor, robust SE. Returns (beta, z, asymptotic p)."""
    import statsmodels.api as sm
    d = samples.copy()
    d["dx"] = d["donor"].map(dx_by_donor).astype(float)
    d["y"] = (d["subtype"] == 1).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = sm.GEE.from_formula("y ~ dx", groups="donor", data=d,
                                family=sm.families.Binomial(),
                                cov_struct=sm.cov_struct.Exchangeable()).fit()
    return float(m.params["dx"]), float(m.tvalues["dx"]), float(m.pvalues["dx"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default=LABELS)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    a = ap.parse_args()

    df = pd.read_csv(a.labels)
    df["donor"] = df["title"].str.extract(r"^(X?\d+)_")
    df = df.dropna(subset=["donor", "subtype", "brain region", "diagnosis"])

    # ---------- item 15: tie-breaking, quantified ---------- #
    g = df.groupby("donor")["subtype"]
    tied = g.apply(lambda s: (lambda vc: len(vc) > 1 and vc.iloc[0] == vc.iloc[1])
                   (s.value_counts()))
    split = g.apply(lambda s: s.nunique() > 1)
    tie_info = dict(n_donors=int(df["donor"].nunique()), n_samples=int(len(df)),
                    donors_split_across_clusters=int(split.sum()),
                    donors_with_no_majority=int(tied.sum()),
                    tied_donor_ids=[str(x) for x in tied[tied].index],
                    rule="pandas value_counts().index[0] — ordering-dependent, undocumented")

    # ---------- donor table ---------- #
    comp = (pd.crosstab(df["donor"], df["subtype"])
              .pipe(lambda t: t.div(t.sum(1), axis=0)))          # continuous composition
    dxmap = df.groupby("donor")["diagnosis"].first()
    donors = comp.index
    y = (dxmap.loc[donors].values == "SCZ").astype(int)
    maj = g.apply(lambda s: s.value_counts().index[0]).loc[donors].values   # deposited rule

    # ---------- observed statistics ---------- #
    f1 = comp.iloc[:, list(comp.columns).index(1)].values if 1 in comp.columns \
        else comp.iloc[:, 0].values
    X = comp.values
    obs = dict(
        T1_auc_cluster1_fraction=float(roc_auc_score(y, f1)),
        T2_loocv_auc_full_composition=loo_auc(X, y),
        V_majority_deposited=cramers_v(pd.crosstab(maj, y).values),
    )
    t3_b, t3_z, t3_p = gee_dx(df, pd.Series(y, index=donors))
    obs["T3_gee_beta"], obs["T3_gee_z"], obs["T3_gee_asymptotic_p"] = t3_b, t3_z, t3_p

    # ---------- donor-permutation null ---------- #
    rng = np.random.default_rng(SEED)
    n1 = np.empty(a.n_perm); n2 = np.empty(a.n_perm); nv = np.empty(a.n_perm)
    for i in range(a.n_perm):
        yp = rng.permutation(y)
        n1[i] = roc_auc_score(yp, f1)
        nv[i] = cramers_v(pd.crosstab(maj, yp).values)
        n2[i] = np.nan
    # T2 permutation is expensive (LOO-CV per draw); use a reduced count
    n_small = min(1000, a.n_perm)
    for i in range(n_small):
        n2[i] = loo_auc(X, rng.permutation(y))

    n_t3 = min(2000, a.n_perm)
    n3 = np.empty(n_t3)
    for i in range(n_t3):
        n3[i] = abs(gee_dx(df, pd.Series(rng.permutation(y), index=donors))[1])
    t3_perm_p = float((np.sum(n3 >= abs(t3_z)) + 1) / (len(n3) + 1))

    def emp_p(o, null, two_sided_center=None):
        null = null[np.isfinite(null)]
        if two_sided_center is not None:
            return float((np.sum(np.abs(null - two_sided_center)
                                 >= abs(o - two_sided_center)) + 1) / (len(null) + 1))
        return float((np.sum(null >= o) + 1) / (len(null) + 1))

    res = dict(
        tie_breaking=tie_info,
        observed=obs,
        p_values=dict(
            T1_auc=emp_p(obs["T1_auc_cluster1_fraction"], n1, two_sided_center=0.5),
            T2_loocv_auc=emp_p(obs["T2_loocv_auc_full_composition"],
                               n2[:n_small], two_sided_center=0.5),
            T3_gee_asymptotic=obs["T3_gee_asymptotic_p"],
            T3_gee_donor_permutation=t3_perm_p,
            V_majority_deposited=emp_p(obs["V_majority_deposited"], nv),
        ),
        degeneracy=dict(
            distinct_values_of_deposited_V_under_permutation=int(
                len(np.unique(np.round(nv[np.isfinite(nv)], 10)))),
            distinct_values_of_T1_AUC=int(len(np.unique(np.round(n1, 10)))),
            n_perm=a.n_perm, n_perm_T2=n_small,
            reads=("A permutation null supported on a handful of atoms cannot resolve a "
                   "p-value; the continuous statistics can."),
        ),
    )
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(json.dumps(res, indent=2))
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
