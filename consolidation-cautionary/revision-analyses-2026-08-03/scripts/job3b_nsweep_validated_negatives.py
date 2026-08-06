#!/usr/bin/env python3
"""JOB 3b — n-sweep with VALIDATED negative generators (Paulus, "on the synthetic expansion").

WHY THIS SUPERSEDES JOB 3. Job 3 introduced a `curved_manifold` negative — samples on an
arc, t ~ U(-1.5, 1.5), X0 += 3 sin t, X1 += 3 cos t — and the gate certified it 4 of 5 at
n=60. Read naively that is an 80% false-certification rate on a continuum, which would
have been the headline finding. It is not one. Because sin flattens as |t| -> 1.5
(cos 1.5 = 0.07), a uniform t piles density at both ends of the arc, so the PC1 marginal
is genuinely bimodal: Hartigan's dip rejects unimodality in 50% of replicates at n=120 and
70% at n=300 (median p 0.055 and 0.015). The generator was mislabelled. Certifying it was
arguably correct behaviour, and reporting it as a false certification would have repeated,
inside the revision analysis, precisely the error this manuscript exists to document —
asserting a ground-truth label without testing whether the data support it.

THE FIX IS PROCEDURAL, NOT COSMETIC. Every generator labelled NEGATIVE is now validated
before any of its verdicts are read, and the check is recorded in the output alongside the
results rather than performed by hand and forgotten:

    Phase A  for each (generator, n): 20 replicates, Hartigan dip on PC1 of the gate's own
             reduced space. A generator is DISQUALIFIED at that n if the dip rejects
             unimodality in more than 30% of replicates.
    Phase B  only qualified generators are scored. Disqualified ones are written out with
             their statistics and excluded, so the exclusion is auditable.

The dip screen is deliberately ONE-SIDED. Failing it disqualifies a negative; passing it
is necessary, not sufficient. Hartigan's dip on PC1 is underpowered here — it does not
reject unimodality for the `discrete_k2` POSITIVE at n=60 either — so it is used only to
catch generators that are visibly not continua, never to certify that one is.

`curved_manifold` is retained in fixed form: Gaussian t (density peaks mid-arc rather than
at the ends), a narrower angular range, and larger isotropic noise. It passes the screen at
both n (frac 0.00, median dip p 0.83/0.88) while keeping the curvature that made it worth
testing. The original is kept as `curved_manifold_INVALID` so the disqualification is
visible in the record rather than silently dropped.

Newline-delimited JSON; nothing is written to the deposited tree.
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

DIP_MAX_REJECT_FRAC = 0.30   # >30% of replicates rejecting unimodality disqualifies


import zlib  # noqa: E402

# NOTE: was `hash(<str>) % 997`. Python salts str hashes per process (PYTHONHASHSEED),
# so that seeding produced DIFFERENT datasets on every invocation and nothing derived
# from it could be regenerated. zlib.crc32 is stable across processes and platforms.
def _stable_key(name: str) -> int:
    """Process-stable substitute for hash(name) % 997."""
    return zlib.crc32(name.encode()) % 997


def _df(X):
    return pd.DataFrame(X, columns=[f"PATHWAY_{j}" for j in range(X.shape[1])])


def gen(kind: str, n: int, p: int, seed: int):
    rng = np.random.default_rng(seed)
    ns = max(3, p // 10)
    d = np.zeros(p); d[:ns] = rng.normal(0, 1, ns); d /= np.linalg.norm(d)

    if kind == "discrete_k2":                       # POSITIVE control
        n1 = n // 2
        a = rng.normal(0, 1, size=(n1, p)); b = rng.normal(0, 1, size=(n - n1, p))
        b[:, 0] += 3.0
        return _df(np.vstack([a, b])), 2

    if kind in ("corr_weak", "corr_strong"):
        rho = 0.3 if kind == "corr_weak" else 0.85
        cov = np.full((p, p), rho) + np.eye(p) * (1 - rho)
        return _df(rng.multivariate_normal(np.zeros(p), cov, size=n)), 2

    if kind == "nongaussian_gradient":              # heavy-tailed (t_3) latent
        return _df(np.outer(rng.standard_t(df=3, size=n), d)
                   + rng.normal(0, 1, size=(n, p))), 2

    if kind == "heteroscedastic_gradient":          # noise grows along the gradient
        t = rng.normal(0, 1, size=n)
        scale = 0.5 + 1.5 * (t - t.min()) / (np.ptp(t) + 1e-9)
        return _df(np.outer(t, d) + rng.normal(0, 1, size=(n, p)) * scale[:, None]), 2

    if kind == "curved_manifold":                   # FIXED: Gaussian t, narrow arc
        t = np.clip(rng.normal(0, 0.40, size=n), -0.9, 0.9)
        X = rng.normal(0, 0.7, size=(n, p))
        X[:, 0] += 3.0 * np.sin(t); X[:, 1] += 3.0 * np.cos(t)
        return _df(X), 2

    if kind == "curved_manifold_INVALID":           # ORIGINAL: uniform t, endpoint piling
        t = rng.uniform(-1.5, 1.5, size=n)
        X = rng.normal(0, 0.4, size=(n, p))
        X[:, 0] += 3.0 * np.sin(t); X[:, 1] += 3.0 * np.cos(t)
        return _df(X), 2

    if kind == "three_component_gradient":          # overlapping shoulders, still continuous
        t = np.concatenate([rng.normal(m, 1.0, size=n // 3) for m in (-1.2, 0.0, 1.2)])
        t = np.concatenate([t, rng.normal(0.0, 1.0, size=n - len(t))])
        return _df(np.outer(t, d) + rng.normal(0, 1, size=(n, p))), 2

    raise ValueError(kind)


NEGATIVES = ["corr_weak", "corr_strong", "nongaussian_gradient", "heteroscedastic_gradient",
             "curved_manifold", "curved_manifold_INVALID", "three_component_gradient"]
ALL = ["discrete_k2"] + NEGATIVES


def validate_negative(kind, n, p, reps, seed):
    """One-sided screen: reject generators that are visibly not continua."""
    ps = []
    for s in range(reps):
        Z, _ = reduce_scores(gen(kind, n, p, seed + 7919 * s)[0].values, reduced_dim(n), 42)
        ps.append(float(dip_of(Z[:, 0])["p"]))
    ps = np.asarray(ps)
    # FAIL CLOSED. dip_of() returns NaN when the optional `diptest` extra is missing, and
    # `NaN < 0.05` is False -> frac would be 0.0 -> every generator "qualifies", including
    # curved_manifold_INVALID, whose exclusion is the entire premise of this job. A screen
    # that silently returns PASS when it did not run is worse than no screen.
    if not np.isfinite(ps).all():
        raise SystemExit(
            "Hartigan dip unavailable (install the `diptest` extra: "
            "pip install 'pathway-subtyping[discreteness]'). Refusing to score negatives "
            "without the validation screen — it would silently qualify every generator.")
    frac = float((ps < 0.05).mean())
    return dict(kind=kind, n=n, dip_reject_frac=frac, dip_median_p=float(np.median(ps)),
                qualified=bool(frac <= DIP_MAX_REJECT_FRAC))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ns", default="60,120,300")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--dip-reps", type=int, default=20)
    ap.add_argument("--p", type=int, default=50)
    ap.add_argument("--n-ref", type=int, default=100)
    ap.add_argument("--n-bootstrap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--kinds", default=",".join(ALL))
    a = ap.parse_args()
    ns = [int(x) for x in a.ns.split(",")]
    kinds = a.kinds.split(",")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t0 = time.time()
    with open(a.out, "w") as fh:
        w = lambda o: (fh.write(json.dumps(o) + "\n"), fh.flush())
        # Record basenames only: absolute paths are machine-specific and would leak into the
        # deposit. Content identity is carried by the deposited values, not by paths.
        _argv = {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v)
                 for k, v in vars(a).items()}
        w(dict(record="provenance", script=os.path.basename(__file__), argv=_argv,
               supersedes="job3_nsweep_generator_diversity.py",
               dip_max_reject_frac=DIP_MAX_REJECT_FRAC,
               note="negatives are dip-screened before scoring; screen is one-sided"))

        # ---- Phase A: validate the negatives ---- #
        print("PHASE A — negative-generator validation", flush=True)
        qualified = {}
        for n in ns:
            for kind in kinds:
                if kind == "discrete_k2":
                    qualified[(n, kind)] = True
                    continue
                v = validate_negative(kind, n, a.p, a.dip_reps, a.seed)
                qualified[(n, kind)] = v["qualified"]
                w(dict(record="negative_validation", **v))
                flag = "qualified" if v["qualified"] else "** DISQUALIFIED **"
                print(f"  n={n:<4} {kind:26s} dip_reject_frac={v['dip_reject_frac']:.2f} "
                      f"median_p={v['dip_median_p']:.4f}  {flag}", flush=True)

        # ---- Phase B: score only qualified generators ---- #
        print("\nPHASE B — gate verdicts on qualified generators", flush=True)
        for n in ns:
            for kind in kinds:
                if not qualified[(n, kind)]:
                    w(dict(record="excluded", n=n, kind=kind,
                           reason="failed dip screen; not a valid negative"))
                    print(f"  n={n:<4} {kind:26s} EXCLUDED (failed dip screen)", flush=True)
                    continue
                tally = dict(certify=0, reject=0, abstain=0)
                for rep in range(a.reps):
                    seed = a.seed + 1000 * rep + _stable_key(kind) + n
                    scores, k = gen(kind, n, a.p, seed)
                    r = DiscretenessGateA(seed=seed, n_ref=a.n_ref,
                                          n_bootstrap=a.n_bootstrap).run(
                        tumor="synthetic", pathway_scores=scores, n_clusters=k,
                        gmm_seed=seed)
                    v = ("certify" if (r.testable and r.passed)
                         else "abstain" if not r.testable else "reject")
                    tally[v] += 1
                    w(dict(record="run", n=n, kind=kind, rep=rep, verdict=v,
                           raw_verdict=str(r.verdict)))
                w(dict(record="cell_summary", n=n, kind=kind, reps=a.reps, **tally))
                note = ""
                if kind != "discrete_k2" and tally["certify"]:
                    note = f"  <-- {tally['certify']}/{a.reps} FALSE CERTIFICATIONS"
                print(f"  n={n:<4} {kind:26s} certify={tally['certify']:>2} "
                      f"reject={tally['reject']:>2} abstain={tally['abstain']:>2}{note}",
                      flush=True)
        w(dict(record="done", elapsed_sec=round(time.time() - t0, 1)))
    print(f"\nJOB3b done in {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
