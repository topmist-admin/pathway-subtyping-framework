#!/usr/bin/env python3
"""JOB 6b — PC1 diagnostics for the three real-data cohorts, on a stated basis.

WHY THIS EXISTS. Result 3 compares the certified BRCA partition against the two within-study
controls using two PC1 diagnostics (Hartigan's dip and a standardised gap). The LGG gap value
quoted there had **no deposited source**: `std_gap_pc1` is defined only in job6, and job6 is
BRCA-only. The number was computed ad hoc, which is exactly the situation the manuscript
criticises elsewhere.

It also matters that these be comparable. An earlier version of that table placed the LUAD
control alongside LGG and BRCA as if "computed identically"; it is not — the deposited LUAD
control is an **8-pathway immune panel plus a derived `__immune_score__` aggregate** (this
script selects all numeric columns and so reports 9) while LGG and BRCA are scored on
50-pathway matrices,
so its PC1 is a different object. This job therefore records the feature count for every row,
so a reader can see which comparisons are like-for-like rather than having to trust a caption.

Both diagnostics are computed in the gate's own reduced space (`reduce_scores` at
`reduced_dim(n)`), which is where the gate's `dip_pc1_p` is measured — not on z-scored pathway
columns, which is a different statistic that the gate-calibration package also reports and
which differs by construction (0.86 vs 0.9963 for LGG).

Reads the deposited cached control inputs; no network. Deterministic (seed 42).
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from sklearn.mixture import GaussianMixture                                     # noqa: E402
from pathway_subtyping.discreteness.gate_a_discreteness_null import (           # noqa: E402
    reduce_scores, reduced_dim, dip_of,
)

CACHED = os.path.join(REPO, "consolidation-cautionary/cross-domain/gate_calibration/"
                            "results/cached_inputs")
COHORTS = [
    ("TCGA-LGG",  "discrete_pathways.csv.gz",  2, "discrete control (IDH axis)"),
    ("TCGA-LUAD", "continuum_pathways.csv.gz", 2, "continuum control (immune gradient)"),
]


def std_gap_pc1(Z, labels):
    """Standardised gap between the two largest groups' PC1 means (pooled SD)."""
    u, c = np.unique(labels, return_counts=True)
    a, b = u[np.argsort(-c)[:2]]
    xa, xb = Z[labels == a, 0], Z[labels == b, 0]
    sd = np.sqrt(((len(xa) - 1) * xa.var(ddof=1) + (len(xb) - 1) * xb.var(ddof=1))
                 / max(len(xa) + len(xb) - 2, 1))
    return float(abs(xa.mean() - xb.mean()) / sd) if sd > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    # Fail closed on a missing `diptest` extra. `dip_of` returns NaN when the optional
    # dependency is absent; every downstream use here is a reported diagnostic, so a NaN
    # would be written straight into the deposit as a bare `NaN` token — which is not valid
    # JSON (RFC 8259) and is silently accepted only by Python's own parser. Probe once, up
    # front, so we never emit a partial file.
    _probe = dip_of(np.linspace(0.0, 1.0, 32))["p"]
    if not np.isfinite(_probe):
        raise SystemExit(
            "Hartigan dip unavailable (install the `diptest` extra: pip install "
            "'pathway-subtyping[discreteness]'). Refusing to run: this job reports dip p-values "
            "as evidence, and NaN would be written to the output as invalid JSON.")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        w = lambda o: (fh.write(json.dumps(o) + "\n"), fh.flush())
        # Record basenames only: absolute paths are machine-specific and would leak into the
        # deposit. Content identity is carried by the deposited values, not by paths.
        _argv = {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v)
                 for k, v in vars(a).items()}
        w(dict(record="provenance", script=os.path.basename(__file__), argv=_argv,
               note="dip and standardised gap on PC1 of the gate's reduced space; "
                    "n_features recorded per cohort because the cohorts are NOT on a common "
                    "feature basis and only same-basis rows are comparable"))
        print(f"{'cohort':26s}{'n':>6}{'features':>10}{'d':>4}{'dip p (PC1)':>13}{'std gap':>10}",
              flush=True)
        for name, fn, k, role in COHORTS:
            p = os.path.join(CACHED, fn)
            if not os.path.exists(p):
                print(f"  {name}: cached input missing ({fn}) — skipped", flush=True)
                w(dict(record="skipped", cohort=name, reason=f"missing {fn}")); continue
            P = pd.read_csv(p, index_col=0).select_dtypes(include=[np.number])
            Z, _ = reduce_scores(P.values, reduced_dim(P.shape[0]), a.seed)
            lab = GaussianMixture(k, covariance_type="full", n_init=10,
                                  random_state=a.seed, reg_covar=1e-6).fit(Z).predict(Z)
            dip = float(dip_of(Z[:, 0])["p"]); gap = std_gap_pc1(Z, lab)
            w(dict(record="cohort", cohort=name, role=role, k=k,
                   n=int(P.shape[0]), n_features=int(P.shape[1]), reduced_d=int(Z.shape[1]),
                   dip_pc1_p=dip, std_gap_pc1=gap, source=fn))
            print(f"{name:26s}{P.shape[0]:>6}{P.shape[1]:>10}{Z.shape[1]:>4}"
                  f"{dip:>13.4f}{gap:>10.3f}", flush=True)
        w(dict(record="comparability_note",
               text="TCGA-LGG (50 features) and TCGA-BRCA (50 features, see job6) are on a "
                    "common basis and may be compared. The deposited TCGA-LUAD control is a "
                    "8-pathway immune panel plus a derived __immune_score__ aggregate (this script "
                    "selects all numeric columns and so reports 9); its PC1 is a "
                    "different object and its gap must "
                    "not be placed in the same column."))
        w(dict(record="done"))
    print(f"\nJOB6b done -> {a.out}")


if __name__ == "__main__":
    main()
