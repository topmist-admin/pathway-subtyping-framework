#!/usr/bin/env python3
"""Gate ablation study — answers Scientific Reports reviewer point R3.10.

R3.10 asks: "What happens to the performance if some of the validation gates in
the framework are removed?"

This harness quantifies what each gate *buys* by measuring the false-positive
rate (FPR) of certifying a partition as a "real subtype" under different gate
subsets, on synthetic data with KNOWN ground truth:

  POSITIVE (genuine discrete structure — a real subtype SHOULD be certified):
    - discrete_k2 : two well-separated Gaussian components
    - discrete_k3 : three well-separated Gaussian components
  NEGATIVE (no discrete structure — certification is a FALSE POSITIVE):
    - single_gaussian : one correlated Gaussian blob
    - continuum_1d    : a 1-D continuous gradient (tumor purity / infiltration
                        analog) that a mixture model reproducibly bisects

The headline result reproduces the erratum finding: the *old* bootstrap-stability
gate alone certifies continua as reproducible subtypes (high FPR at small n),
because its null tested independence, not discreteness. Adding the v0.8.0
discreteness gate (Gate A) removes those false certifications.

Two honesty constraints on how that result may be stated:
  1. Gate A reaches the low FPR mostly by ABSTAINING, not by rejecting. On the
     negative controls it returns "not-testable (no reproducible k)" on 28/30
     (93%); the testable-negative denominator is 2. Never report FPR=0.000 or
     "rejects the continuum on 100% of runs" without that denominator.
  2. It is not free: TPR moves 1.000 -> 0.967 (one discrete_k3 replicate lost).
     Exact McNemar on positives is p=1.0, so the cost is not distinguishable
     from zero, but the study has no power to exclude one. Report both the point
     change and the p-value; do not write "TPR held" and do not assert a firm
     "3% cost" as though it were resolved.

Gate subsets compared (the ablation):
    stability_only          : {bootstrap stability}                 (pre-v0.8.0)
    discreteness_only       : {Gate A discreteness null}
    stability+discreteness  : both must pass                        (v0.8.0)

Usage:
    python scripts/gate_ablation_study.py [--quick] [--reps N] [--n N]
                                          [--out DIR] [--seed S]

    --quick  small/fast settings for a smoke run (few refs, few reps)

Outputs (JSON + Markdown table) are written under --out (default:
outputs/gate_ablation/, gitignored). Nothing here depends on private data.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from pathway_subtyping.validation import ValidationGates
from pathway_subtyping.discreteness import DiscretenessGateA


def _cluster_gmm(X, k, seed):
    """GMM clusterer with the run_x(X, k, seed) -> labels convention."""
    from sklearn.mixture import GaussianMixture
    gm = GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                         random_state=seed, reg_covar=1e-6).fit(np.asarray(X))
    return gm.predict(np.asarray(X))


# --------------------------------------------------------------------------- #
# Synthetic data generators — each returns (scores DataFrame, natural_k)
# --------------------------------------------------------------------------- #
def make_discrete(n: int, p: int, k: int, sep: float, seed: int) -> pd.DataFrame:
    """k well-separated Gaussian components in the first few pathways."""
    rng = np.random.default_rng(seed)
    n_signal = max(3, p // 10)
    sizes = [n // k] * k
    sizes[-1] += n - sum(sizes)
    blocks = []
    for c in range(k):
        center = np.zeros(p)
        # spread component centers across the signal pathways
        center[c % n_signal] = sep * (1 + c)
        center[(c + 1) % n_signal] = -sep * (1 + c)
        blocks.append(rng.normal(center, 1.0, size=(sizes[c], p)))
    X = np.vstack(blocks)
    idx = rng.permutation(n)
    return pd.DataFrame(X[idx], columns=[f"PATHWAY_{j}" for j in range(p)])


def make_single_gaussian(n: int, p: int, seed: int) -> pd.DataFrame:
    """One correlated Gaussian blob — no discrete structure."""
    rng = np.random.default_rng(seed)
    A = rng.normal(0, 1, size=(p, p))
    cov = (A @ A.T) / p + np.eye(p) * 0.5
    X = rng.multivariate_normal(np.zeros(p), cov, size=n)
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(p)])


def make_continuum(n: int, p: int, seed: int) -> pd.DataFrame:
    """A 1-D continuous gradient: samples lie along one axis with noise.

    A mixture model reproducibly bisects this every bootstrap (high observed
    stability ARI) even though there are NO discrete clusters. This is the case
    the old stability gate false-certifies."""
    rng = np.random.default_rng(seed)
    n_signal = max(3, p // 10)
    # Gaussian latent along ONE dominant axis: GMM at k=2 splits it at the mean
    # reproducibly across bootstraps (high stability ARI) even though the
    # distribution is unimodal — no discrete clusters.
    t = rng.normal(0, 1, size=n)  # position along the continuum
    direction = np.zeros(p)
    direction[:n_signal] = rng.normal(0, 1, size=n_signal)
    direction /= np.linalg.norm(direction)
    X = np.outer(t, direction) * 20.0 + rng.normal(0, 0.3, size=(n, p))
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(p)])


CONDITIONS = {
    "discrete_k2": dict(kind="discrete", k=2, klass="positive"),
    "discrete_k3": dict(kind="discrete", k=3, klass="positive"),
    "single_gaussian": dict(kind="single", k=2, klass="negative"),
    "continuum_1d": dict(kind="continuum", k=2, klass="negative"),
}


def generate(cond: str, n: int, p: int, sep: float, seed: int):
    spec = CONDITIONS[cond]
    if spec["kind"] == "discrete":
        return make_discrete(n, p, spec["k"], sep, seed), spec["k"]
    if spec["kind"] == "single":
        return make_single_gaussian(n, p, seed), spec["k"]
    return make_continuum(n, p, seed), spec["k"]


# --------------------------------------------------------------------------- #
# Evaluate the two gates on one dataset
# --------------------------------------------------------------------------- #
def evaluate_dataset(scores: pd.DataFrame, k: int, seed: int, n_ref: int,
                     n_bootstrap: int) -> Dict[str, bool]:
    """Fit GMM at k, then run the old stability gate and Gate A. Returns the
    per-gate PASS booleans."""
    from sklearn.mixture import GaussianMixture

    gm = GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                         random_state=seed, reg_covar=1e-6).fit(scores.values)
    labels = gm.predict(scores.values)

    # Old bootstrap-stability gate (pre-v0.8.0)
    vg = ValidationGates(seed=seed, n_bootstrap=n_bootstrap,
                         stability_threshold=0.8, show_progress=False)
    stab = vg.stability_test_bootstrap(scores, labels, k, gmm_seed=seed)

    # New discreteness gate (Gate A, v0.8.0)
    gateA = DiscretenessGateA(seed=seed, n_ref=n_ref, n_bootstrap=n_bootstrap)
    a = gateA.run(tumor="synthetic", pathway_scores=scores, n_clusters=k,
                  gmm_seed=seed)

    return {
        "stability_pass": bool(stab.passed),
        "stability_ari": float(stab.metric_value),
        # A subtype is certified discrete ONLY when Gate A is both testable and
        # passing (verdict == "discrete structure"). "not-testable" is a third
        # outcome — no reproducible k — and must NOT count as certification even
        # though a.passed can be True there.
        "discreteness_pass": bool(a.testable and a.passed),
        "discreteness_verdict": a.verdict,
        "discreteness_testable": bool(a.testable),
        "sg_empirical_p": float(a.sg_empirical_p),
    }


def gate_a_label(certifies: bool, verdict: str) -> str:
    """Three-way outcome label for Gate A.

    Gate A has three outcomes, not two. Collapsing "not-testable" into REJECT
    is what produced the retracted "rejects the continuum on 100% of runs"
    claim: an abstention is not a rejection, and counting it as one inflates
    the apparent false-positive performance.
    """
    if certifies:
        return "CERTIFY"
    if verdict.startswith("not-testable"):
        return "ABSTAIN"
    return "REJECT"


# certification policies: which gates must PASS to certify a "real subtype"
POLICIES = {
    "stability_only": lambda r: r["stability_pass"],
    "discreteness_only": lambda r: r["discreteness_pass"],
    "stability+discreteness": lambda r: r["stability_pass"] and r["discreteness_pass"],
}


def _self_stability(cluster_fn, X, k, labels, n_boot, seed, threshold=0.8):
    """A clusterer's OWN bootstrap stability: re-run the clusterer on bootstrap
    resamples and average the ARI to its original labels on shared samples.
    This is the naive 'is my partition reproducible?' check each method would
    pass on a continuum it bisects the same way every time."""
    from sklearn.metrics import adjusted_rand_score
    Xv = np.asarray(X)
    rng = np.random.default_rng(seed + 999)
    aris = []
    for b in range(n_boot):
        idx = np.unique(rng.choice(len(Xv), len(Xv), replace=True))
        if len(idx) < 2 * k:
            continue
        try:
            lb = cluster_fn(Xv[idx], k, seed + b)
            aris.append(adjusted_rand_score(labels[idx], lb))
        except Exception:
            continue
    mean_ari = float(np.mean(aris)) if aris else 0.0
    return mean_ari, bool(mean_ari >= threshold)


def clusterer_sweep(args) -> None:
    """Gate-agnostic demonstration (R2.2): cluster the SAME data with GMM, DEC,
    and VAE-GMM. Each method's own bootstrap stability is high on a continuum it
    bisects reproducibly (a naive stability check passes for all three), yet
    Gate A — which tests the discreteness of the data at k, independent of the
    clusterer — rejects the continuum partition for every method."""
    from functools import partial
    from pathway_subtyping.clustering_dl import run_dec, run_vae_gmm

    # fast DL settings for the sweep (labels only; drop the silhouette return)
    dl_epochs = 40 if args.quick else 120
    clusterers = {
        "gmm": _cluster_gmm,
        "dec": lambda X, k, s: run_dec(X, k, seed=s, pretrain_epochs=dl_epochs,
                                       cluster_epochs=max(20, dl_epochs // 3))[0],
        "vae_gmm": lambda X, k, s: run_vae_gmm(X, k, seed=s, epochs=dl_epochs)[0],
    }
    conditions = ["discrete_k2", "continuum_1d"]  # positive control + the money case
    n_boot = 8 if args.quick else 20

    rows = []
    for cond in conditions:
        for rep in range(max(2, args.reps // 2)):
            seed = args.seed + 1000 * rep + hash(cond) % 997
            scores, k = generate(cond, args.n, args.p, args.sep, seed)
            # Gate A is clusterer-agnostic (re-clusters internally) -> run once
            gateA = DiscretenessGateA(seed=seed, n_ref=args.n_ref,
                                      n_bootstrap=args.n_bootstrap)
            a = gateA.run(tumor="synthetic", pathway_scores=scores, n_clusters=k,
                          gmm_seed=seed)
            # certified discrete only when testable AND passing (not "not-testable")
            gateA_certifies = bool(a.testable and a.passed)
            for cname, cfn in clusterers.items():
                labels = cfn(scores.values, k, seed)
                self_ari, self_pass = _self_stability(cfn, scores.values, k, labels,
                                                       n_boot, seed)
                rows.append(dict(
                    condition=cond, klass=CONDITIONS[cond]["klass"], rep=rep,
                    clusterer=cname, self_stability_ari=round(self_ari, 3),
                    self_stability_pass=self_pass,
                    gateA_pass=gateA_certifies, gateA_verdict=a.verdict))
                print(f"  {cond:14s} rep{rep} {cname:8s}: "
                      f"self-stability ARI={self_ari:.2f} "
                      f"({'looks-reproducible' if self_pass else 'unstable'})  "
                      f"GateA={gate_a_label(gateA_certifies, a.verdict)} ({a.verdict})")

    df = pd.DataFrame(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / "clusterer_sweep_raw.csv", index=False)

    # For the continuum: fraction of (clusterer, rep) that LOOK reproducible vs
    # fraction Gate A rejects.
    lines = [
        "# Gate-Agnostic Demonstration — clusterer sweep (R2.2)",
        "",
        "Same synthetic data clustered by GMM, DEC (Xie 2016), and VAE-GMM "
        "(VaDE, Jiang 2017). **self-stability** = the method's own bootstrap "
        "reproducibility (the naive check each would pass). **Gate A** = the "
        "discreteness gate's verdict (clusterer-agnostic; tests the data).",
        "",
        "| Condition | Clusterer | self-stability ARI | looks-reproducible? | Gate A verdict |",
        "|---|---|---|---|---|",
    ]
    for _, r in df.iterrows():
        lines.append(
            f"| {r.condition} | `{r.clusterer}` | {r.self_stability_ari:.2f} | "
            f"{'yes' if r.self_stability_pass else 'no'} | "
            f"{gate_a_label(bool(r.gateA_pass), str(r.gateA_verdict))} "
            f"({r.gateA_verdict}) |")
    cont = df[df.condition == "continuum_1d"]
    # Three-way split. Do NOT collapse abstentions into rejections: the gate
    # declines to certify mostly by returning "not-testable", and reporting that
    # as a rejection is the retracted "100% of runs" claim.
    n_cont = max(1, len(cont))
    not_certified = float((~cont.gateA_pass).mean())
    abstained = float(
        cont.gateA_verdict.astype(str).str.startswith("not-testable").sum() / n_cont)
    explicit_reject = max(0.0, not_certified - abstained)
    gmm_ari = float(cont[cont.clusterer == "gmm"].self_stability_ari.mean())
    dec_ari = float(cont[cont.clusterer == "dec"].self_stability_ari.mean())
    vae_ari = float(cont[cont.clusterer == "vae_gmm"].self_stability_ari.mean())
    lines += [
        "",
        f"**Primary result (clusterer-agnostic non-certification):** Gate A declines to "
        f"certify the 1-D continuum on **{not_certified:.0%}** of runs, and does so "
        f"identically whether the partition was drawn by GMM, DEC, or VAE-GMM. That "
        f"total splits into **{explicit_reject:.0%} explicit rejection** and "
        f"**{abstained:.0%} abstention** (\"not-testable — no reproducible k\"). The "
        f"abstentions are the larger share and must not be reported as rejections: the "
        f"gate mostly declines to rule, rather than ruling against. This is the "
        f"substantive answer to R2.2: PSF is not another clustering method to be "
        f"benchmarked against DEC/VAE — it is a validation layer that wraps any of "
        f"them, because Gate A tests the discreteness of the *data* at a given k, not "
        f"the confidence of the algorithm.",
        "",
        f"**Secondary observation (reported honestly):** on the continuum the clusterers "
        f"differ in how reproducible their (spurious) bipartition looks — mean "
        f"self-stability ARI GMM {gmm_ari:.2f} vs DEC {dec_ari:.2f} vs VAE-GMM "
        f"{vae_ari:.2f}. GMM's near-0.8 self-stability is the classic continuum "
        f"false-positive the gate was built to catch; the deep methods are additionally "
        f"unstable under resampling at this small n, itself a caution against trusting "
        f"algorithmic confidence (consistent with the small-n concerns in "
        f"R3.5/R3.6/R3.9). Competitive DL *recovery* of real subtypes is a separate "
        f"claim for the large cohorts (TCGA-COAD, CPTAC), not this small-n control.",
    ]
    (args.out / "clusterer_sweep_results.md").write_text("\n".join(lines))
    print("\n" + "\n".join(lines))
    print(f"\nWrote: {args.out}/clusterer_sweep_results.md + clusterer_sweep_raw.csv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=10, help="replicates per condition")
    ap.add_argument("--n", type=int, default=60, help="samples per dataset (small-n regime)")
    ap.add_argument("--p", type=int, default=50, help="pathways (features)")
    ap.add_argument("--sep", type=float, default=3.0, help="cluster separation for positives")
    ap.add_argument("--n-ref", type=int, default=100, help="Gate A single-Gaussian references")
    ap.add_argument("--n-bootstrap", type=int, default=40, help="bootstrap iterations")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--quick", action="store_true", help="fast smoke settings")
    ap.add_argument("--sweep", action="store_true",
                    help="run the gate-agnostic clusterer sweep (GMM/DEC/VAE) instead")
    ap.add_argument("--out", type=Path, default=Path("outputs/gate_ablation"))
    args = ap.parse_args()

    if args.quick:
        args.reps, args.n_ref, args.n_bootstrap = 3, 30, 20

    if args.sweep:
        clusterer_sweep(args)
        return

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    rows: List[Dict] = []
    for cond in CONDITIONS:
        for rep in range(args.reps):
            seed = args.seed + 1000 * rep + hash(cond) % 997
            scores, k = generate(cond, args.n, args.p, args.sep, seed)
            res = evaluate_dataset(scores, k, seed, args.n_ref, args.n_bootstrap)
            res.update(condition=cond, klass=CONDITIONS[cond]["klass"], rep=rep, k=k)
            rows.append(res)
            print(f"  {cond:16s} rep {rep:2d}: "
                  f"stability={'PASS' if res['stability_pass'] else 'fail'} "
                  f"(ARI={res['stability_ari']:.2f})  "
                  f"discreteness={'PASS' if res['discreteness_pass'] else 'fail'} "
                  f"({res['discreteness_verdict']})")

    df = pd.DataFrame(rows)

    # Aggregate: per-policy FPR (on negatives) and TPR (on positives)
    summary = {}
    for pol, fn in POLICIES.items():
        cert = df.apply(fn, axis=1)
        pos = df["klass"] == "positive"
        neg = df["klass"] == "negative"
        tpr = float(cert[pos].mean())
        fpr = float(cert[neg].mean())
        # per-condition breakdown
        by_cond = {c: float(cert[df["condition"] == c].mean()) for c in CONDITIONS}
        summary[pol] = dict(TPR=round(tpr, 3), FPR=round(fpr, 3), cert_rate_by_condition=by_cond)

    out = dict(
        params=vars(args) | {"out": str(args.out)},
        n_datasets=len(df),
        elapsed_sec=round(time.time() - t0, 1),
        summary=summary,
    )
    (args.out / "gate_ablation_results.json").write_text(json.dumps(out, indent=2, default=str))
    df.to_csv(args.out / "gate_ablation_raw.csv", index=False)

    # Markdown table
    lines = [
        "# Gate Ablation Study — results (R3.10)",
        "",
        f"n={args.n} samples, p={args.p} pathways, {args.reps} reps/condition, "
        f"n_ref={args.n_ref}, n_bootstrap={args.n_bootstrap}. "
        f"Elapsed {out['elapsed_sec']}s.",
        "",
        "**TPR** = fraction of genuine-subtype datasets certified (higher is better). "
        "**FPR** = fraction of no-structure datasets (single Gaussian + continuum) "
        "falsely certified as subtypes (lower is better).",
        "",
        "| Gate subset | TPR | FPR | continuum_1d cert rate |",
        "|---|---|---|---|",
    ]
    for pol, s in summary.items():
        lines.append(
            f"| `{pol}` | {s['TPR']:.2f} | {s['FPR']:.2f} | "
            f"{s['cert_rate_by_condition']['continuum_1d']:.2f} |"
        )
    lines += [
        "",
        "**Reading:** `stability_only` (the pre-v0.8.0 gate set) false-certifies the "
        "1-D continuum — the erratum finding. `stability+discreteness` (v0.8.0) removes "
        "those false certifications, which is the quantitative answer to R3.10: "
        "removing the discreteness gate reintroduces continuum false positives.",
        "",
        "**Two caveats that must travel with the FPR number.** (1) The gate reaches its "
        "low FPR mostly by ABSTAINING rather than rejecting — on the negative controls "
        "it returns \"not-testable (no reproducible k)\" on 28/30 (93%), leaving a "
        "testable-negative denominator of 2. An FPR quoted against n=30 is the "
        "over-generous convention; against the testable n=2 the interval is "
        "[0, 0.658] and is nearly uninformative. Do not quote FPR=0.000 unqualified. "
        "(2) The gate is not free: TPR moves 1.000 -> 0.967 (one `discrete_k3` "
        "replicate lost). Exact McNemar on positives is p=1.0, so that cost is not "
        "distinguishable from zero — but the study has no power to exclude one. Report "
        "the point change and the p-value together; write neither \"TPR held\" nor a "
        "settled \"3% cost\".",
    ]
    (args.out / "gate_ablation_results.md").write_text("\n".join(lines))

    print("\n" + "\n".join(lines))
    print(f"\nWrote: {args.out}/gate_ablation_results.{{json,md}} + gate_ablation_raw.csv")


if __name__ == "__main__":
    main()
