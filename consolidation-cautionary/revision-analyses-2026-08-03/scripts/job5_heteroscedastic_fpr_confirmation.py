#!/usr/bin/env python3
"""JOB 5 — confirm the heteroscedastic false-certification, and isolate its cause.

WHAT JOB 3b FOUND. On dip-validated negatives, `heteroscedastic_gradient` — a unimodal 1-D
continuum whose noise scale grows along the gradient — was certified as discrete
1/5 at n=60, 1/5 at n=120, and **4/5 at n=300**. The error rate RISES with sample size. That
is the signature of a systematic false positive, not of limited power: more data makes the
departure from the null easier to detect, and the gate reports every detectable departure as
"discrete structure".

WHY THIS WOULD MATTER. The gate's reference is a single Gaussian with diagonal covariance.
Its null hypothesis is therefore "these data are one homoscedastic Gaussian blob", and
rejecting it licenses the conclusion "discrete structure". But a continuum with non-constant
variance also rejects that null while containing no groups whatsoever. If this holds, the
screen tests **departure from a single Gaussian**, which is necessary but NOT sufficient for
discreteness — and variance heterogeneity along a gradient is entirely ordinary in real
expression data.

DESIGN — two things job 3b could not establish with 5 replicates per cell.

  1. POWER. 20 replicates per cell at n = 120, 300, 600, so the trend in n has an interval
     around it rather than being read off 5 draws.

  2. CAUSAL ISOLATION. A matched `homoscedastic_gradient` control: identical latent, identical
     signal direction, identical seeds — the ONLY difference is that the noise scale is held
     constant at the heteroscedastic version's mean. If the heteroscedastic arm false-certifies
     and the homoscedastic arm does not, variance heterogeneity is the cause, not the gradient,
     not the dimensionality, not the seeds.

Both arms are dip-screened at every n before scoring, on the same one-sided rule as job 3b
(>30% rejection disqualifies), so neither can be read as a false certification unless it has
first been shown to be a continuum.

Newline-delimited JSON; nothing written to the deposited tree.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
# Resolve the repo root from this file so the deposit reproduces anywhere.
# Layout: <root>/consolidation-cautionary/revision-analyses-2026-08-03/scripts/
REPO = os.environ.get("PSF_REPO") or os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from pathway_subtyping.discreteness import DiscretenessGateA                  # noqa: E402
from pathway_subtyping.discreteness.gate_a_discreteness_null import (         # noqa: E402
    reduce_scores, reduced_dim, dip_of,
)

DIP_MAX_REJECT_FRAC = 0.30


def gen(kind: str, n: int, p: int, seed: int):
    """Matched pair. Same latent t, same direction d, same seed — noise scale is the only
    difference, so any verdict difference is attributable to variance heterogeneity."""
    rng = np.random.default_rng(seed)
    ns = max(3, p // 10)
    d = np.zeros(p); d[:ns] = rng.normal(0, 1, ns); d /= np.linalg.norm(d)
    t = rng.normal(0, 1, size=n)
    noise = rng.normal(0, 1, size=(n, p))
    scale = 0.5 + 1.5 * (t - t.min()) / (np.ptp(t) + 1e-9)     # ranges 0.5 -> 2.0
    if kind == "heteroscedastic_gradient":
        X = np.outer(t, d) + noise * scale[:, None]
    elif kind == "homoscedastic_gradient":
        X = np.outer(t, d) + noise * float(scale.mean())        # same mean noise SCALE (not power: E[s^2] > (E[s])^2, so the control carries ~4% less total variance — the manuscript states scale, which is correct)
    else:
        raise ValueError(kind)
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(p)])


KINDS = ["heteroscedastic_gradient", "homoscedastic_gradient"]


def validate(kind, n, p, reps, seed):
    ps = [float(dip_of(reduce_scores(gen(kind, n, p, seed + 7919 * s).values,
                                     reduced_dim(n), 42)[0][:, 0])["p"])
          for s in range(reps)]
    ps = np.asarray(ps)
    # FAIL CLOSED — see job3b. dip_of() returns NaN when the optional `diptest` extra is
    # missing, and `NaN < 0.05` is False, so frac would be 0.0 and every generator would
    # silently "qualify". This flag gates which cells are scored below, so a screen that
    # returns PASS when it did not run would invert the result rather than blur it.
    if not np.isfinite(ps).all():
        raise SystemExit(
            "Hartigan dip unavailable (install the `diptest` extra: "
            "pip install 'pathway-subtyping[discreteness]'). Refusing to qualify negative "
            "generators without the validation screen.")
    frac = float((ps < 0.05).mean())
    return dict(kind=kind, n=n, dip_reject_frac=frac,
                dip_median_p=float(np.median(ps)),
                qualified=bool(frac <= DIP_MAX_REJECT_FRAC))


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    ph = k / n; den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (round(max(0.0, c - h), 4), round(min(1.0, c + h), 4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ns", default="120,300,600")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--dip-reps", type=int, default=20)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--n-ref", type=int, default=100)
    ap.add_argument("--n-bootstrap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    ns = [int(x) for x in a.ns.split(",")]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t0 = time.time()
    with open(a.out, "w") as fh:
        w = lambda o: (fh.write(json.dumps(o) + "\n"), fh.flush())
        # Record basenames only: absolute paths are machine-specific and would leak into the
        # deposit. Content identity is carried by the deposited values, not by paths.
        _argv = {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v)
                 for k, v in vars(a).items()}
        w(dict(record="provenance", script=os.path.basename(__file__), argv=_argv,
               design="matched hetero/homoscedastic pair; identical latent, direction, seeds",
               claim_under_test="gate tests departure-from-single-Gaussian, not discreteness"))

        print("PHASE A — both arms must pass the dip screen to be scored", flush=True)
        ok = {}
        for n in ns:
            for kind in KINDS:
                v = validate(kind, n, a.p, a.dip_reps, a.seed)
                ok[(n, kind)] = v["qualified"]; w(dict(record="negative_validation", **v))
                print(f"  n={n:<4} {kind:26s} dip_reject_frac={v['dip_reject_frac']:.2f} "
                      f"median_p={v['dip_median_p']:.4f}  "
                      f"{'qualified' if v['qualified'] else '** DISQUALIFIED **'}", flush=True)

        print("\nPHASE B — false-certification rate, 20 reps/cell", flush=True)
        for n in ns:
            for kind in KINDS:
                if not ok[(n, kind)]:
                    w(dict(record="excluded", n=n, kind=kind, reason="failed dip screen"))
                    print(f"  n={n:<4} {kind:26s} EXCLUDED", flush=True)
                    continue
                tally = dict(certify=0, reject=0, abstain=0)
                for rep in range(a.reps):
                    seed = a.seed + 1000 * rep + n
                    r = DiscretenessGateA(seed=seed, n_ref=a.n_ref,
                                          n_bootstrap=a.n_bootstrap).run(
                        tumor="synthetic", pathway_scores=gen(kind, n, a.p, seed),
                        n_clusters=2, gmm_seed=seed)
                    v = ("certify" if (r.testable and r.passed)
                         else "abstain" if not r.testable else "reject")
                    tally[v] += 1
                    w(dict(record="run", n=n, kind=kind, rep=rep, verdict=v))
                lo, hi = wilson(tally["certify"], a.reps)
                w(dict(record="cell_summary", n=n, kind=kind, reps=a.reps,
                       fpr=tally["certify"] / a.reps, wilson95=[lo, hi], **tally))
                print(f"  n={n:<4} {kind:26s} FPR={tally['certify']:>2}/{a.reps} "
                      f"= {tally['certify']/a.reps:.2f}  Wilson95 [{lo:.2f}, {hi:.2f}]  "
                      f"(abstain {tally['abstain']})", flush=True)
        w(dict(record="done", elapsed_sec=round(time.time() - t0, 1)))
    print(f"\nJOB5 done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
