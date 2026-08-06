#!/usr/bin/env python3
"""JOB 6 — the TCGA-BRCA certification statistics (Paulus item 17).

TCGA-BRCA is the one real-data case in the paper where the discreteness screen is MORE
permissive than bootstrap stability: the screen certifies at both k=5 and the BIC-selected
k=2, while stability fails the 0.80 bar at both (ARI ~0.39 / ~0.43). That is the single
place a reader is most likely to suspect over-certification — and it is reported as two bare
verdicts. The deposited record (`brca_pam50_validation_with_DL.json`) contains only:

    discreteness_gateA_at_k5    : {certified: true, verdict: "discrete structure"}
    discreteness_gateA_at_bic_k : {k: 2, certified: true, verdict: "discrete structure"}

No observed statistic, no reference percentile, no p-value — while every other certification
in the paper reports one. The numbers were computed and discarded, not never computed: the
script builds full `GateAv2Result` objects and then serialises two fields.

This job re-runs the same two gate calls and records the whole result, plus the two PC1
diagnostics reported for the LUAD continuum control so the three real-data cases can be read
side by side:

    observed bootstrap ARI, single-Gaussian 95th percentile, empirical p, p-at-floor flag,
    testability, Hartigan dip on PC1, standardised gap on PC1

Fetch and scoring are delegated to the deposited script's own functions — this file does not
re-implement them, so the pathway matrix is identical to the deposited one by construction.

Network: www.cbioportal.org (public, read-only). Deterministic (seed 42).
Writes nothing to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
BRCA_SCRIPTS = os.path.join(REPO, "consolidation-cautionary/cross-domain/cancer_r38/scripts")
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, BRCA_SCRIPTS)

from sklearn.mixture import GaussianMixture                                    # noqa: E402
from sklearn.preprocessing import StandardScaler                               # noqa: E402
from pathway_subtyping.discreteness import DiscretenessGateA                   # noqa: E402
from pathway_subtyping.discreteness.gate_a_discreteness_null import (          # noqa: E402
    reduce_scores, reduced_dim, dip_of,
)
import fetch_and_validate_brca as B                                            # noqa: E402


def std_gap_pc1(Z: np.ndarray, labels: np.ndarray) -> float:
    """Standardised gap between the two largest groups' PC1 means, as for LUAD."""
    u, c = np.unique(labels, return_counts=True)
    a, b = u[np.argsort(-c)[:2]]
    xa, xb = Z[labels == a, 0], Z[labels == b, 0]
    sd = np.sqrt(((len(xa) - 1) * xa.var(ddof=1) + (len(xb) - 1) * xb.var(ddof=1))
                 / max(len(xa) + len(xb) - 2, 1))
    return float(abs(xa.mean() - xb.mean()) / sd) if sd > 0 else float("nan")


def _require_diptest():
    """Fail closed if the optional `diptest` extra is missing.

    `dip_of()` returns NaN rather than raising when the extra is absent. Every dip value
    this script writes is reported as evidence, and a NaN would be serialised as a bare
    `NaN` token, which is not valid JSON (RFC 8259) and is accepted only by Python's own
    parser. Probe once, up front, so no partial output is ever written.
    """
    import numpy as _np
    if not _np.isfinite(dip_of(_np.linspace(0.0, 1.0, 32))["p"]):
        raise SystemExit(
            "Hartigan dip unavailable (install the `diptest` extra: pip install "
            "'pathway-subtyping[discreteness]'). Refusing to run: this script reports dip "
            "p-values as evidence, and NaN would be written to the output as invalid JSON.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--panel", default=B.DEFAULT_PANEL)
    ap.add_argument("--n-ref", type=int, default=200)     # deposited default
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    _require_diptest()

    t0 = time.time()
    print(f"Fetching TCGA-BRCA pathway scores from cBioPortal ({B.STUDY})...", flush=True)
    P = B.fetch_pathway_scores(a.panel)
    X = StandardScaler().fit_transform(P.values)
    print(f"  pathway matrix {P.shape}", flush=True)

    pam50_by_patient = B.fetch_pam50()
    pam50 = pd.Series({s: pam50_by_patient.get(B._patient_of(s)) for s in P.index})
    n_pam = pam50.dropna().nunique()
    k5 = max(2, n_pam) if n_pam else 5

    bic = {kk: GaussianMixture(kk, covariance_type="full", n_init=5, random_state=a.seed,
                               reg_covar=1e-6).fit(X).bic(X) for kk in range(2, 8)}
    k_bic = min(bic, key=bic.get)
    print(f"  k(PAM50)={k5}  k(BIC)={k_bic}  PAM50-labelled={int(pam50.notna().sum())}",
          flush=True)

    Z, _ = reduce_scores(P.values, reduced_dim(P.shape[0]), a.seed)
    # Record basenames only: absolute paths are machine-specific and would leak into the
    # deposit. Content identity is carried by the deposited values, not by paths.
    _argv = {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v)
             for k, v in vars(a).items()}
    out = dict(record="provenance", script=os.path.basename(__file__), argv=_argv,
               study=B.STUDY, matrix_shape=list(P.shape),
               n_pam50_labelled=int(pam50.notna().sum()),
               k_pam50=int(k5), k_bic=int(k_bic),
               reduced_d=int(Z.shape[1]),
               note="fetch+scoring delegated to the deposited script; gate re-run in full")
    recs = [out]

    for tag, k in (("k5", k5), ("bic_k", k_bic)):
        print(f"  running gate at k={k} (n_ref={a.n_ref}) ...", flush=True)
        r = DiscretenessGateA(seed=a.seed, n_ref=a.n_ref).run("BRCA", P, k, gmm_seed=a.seed)
        lab = GaussianMixture(k, covariance_type="full", n_init=10, random_state=a.seed,
                              reg_covar=1e-6).fit(Z).predict(Z)
        rec = dict(record="gate", arm=tag, k=int(k),
                   certified=bool(r.testable and r.passed), verdict=str(r.verdict),
                   testable=bool(r.testable),
                   observed_bootstrap_ari=float(r.observed_stability),
                   observed_ci95=[float(x) for x in r.observed_ci95],
                   sg_ref_p95=float(r.sg_ref_p95),
                   sg_empirical_p=float(r.sg_empirical_p),
                   # `sg_empirical_p` is stored rounded to 4dp, so comparing it against the UNROUNDED
                   # floor reported False for values that are exactly at the floor
                   # (1/201 = 0.0049751 -> stored 0.005 > floor). Compare at the same precision.
                   p_at_floor=bool(round(float(r.sg_empirical_p), 4)
                                   <= round(1.0 / (a.n_ref + 1), 4) + 1e-12),
                   p_floor=1.0 / (a.n_ref + 1),
                   dip_pc1_p=float(dip_of(Z[:, 0])["p"]),
                   std_gap_pc1=std_gap_pc1(Z, lab))
        recs.append(rec)
        print(f"    k={k}: obs={rec['observed_bootstrap_ari']:.4f} "
              f"p95={rec['sg_ref_p95']:.4f} p={rec['sg_empirical_p']:.4f} "
              f"(floor {rec['p_floor']:.4f}) certified={rec['certified']} "
              f"dip_pc1_p={rec['dip_pc1_p']:.4f} std_gap={rec['std_gap_pc1']:.3f}",
              flush=True)

    recs.append(dict(record="done", elapsed_sec=round(time.time() - t0, 1)))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    print(f"\nJOB6 done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
