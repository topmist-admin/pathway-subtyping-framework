#!/usr/bin/env python3
"""Honest re-analysis of the gate ablation — three-way accounting, CIs, McNemar,
and a head-to-head of the gate against the methods it is built from.

Consumes the deposited raw CSVs from the ablation run (no recompute):
  ../results/gate_ablation_raw.csv     (60 datasets x gate subsets)
  ../results/clusterer_sweep_raw.csv   (clusterer-generality sweep)

Answers three reviewer points that the first write-up got wrong:

  1. THREE-WAY ACCOUNTING. The gate returns certify / reject / not-testable. The
     original ablation folded `not-testable` (an abstention) into "correct
     rejection", manufacturing FPR=0.000 on negatives it never actually judged.
     Here abstentions are reported as abstentions and excluded from FPR.

  2. CONFIDENCE INTERVALS + PAIRED TEST. Every rate gets a Wilson 95% CI. The
     stability-only vs stability+discreteness contrast is a PAIRED design (same
     datasets), so significance is by exact McNemar, not an unpaired comparison.

  3. HEAD-TO-HEAD. The gate's decision rule is `obs > sg_p95 and sg_p < alpha`
     (single-Gaussian null on a bootstrap-ARI statistic). Its cited components —
     the gap statistic and Hartigan's dip — are computed but never enter the rule.
     We show what the SigClust-style p-value ALONE achieves versus the full gate,
     which quantifies the actual increment.

Deterministic; no network.
"""
from __future__ import annotations

import argparse
import json
import math
import os

import pandas as pd

ALPHA = 0.05
POS = {"discrete_k2", "discrete_k3"}
NEG = {"single_gaussian", "continuum_1d"}


def wilson(k: int, n: int):
    """Wilson score 95% CI for a binomial proportion."""
    if n == 0:
        return [None, None]
    z = 1.959963985
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [round((c - h) / d, 4), round((c + h) / d, 4)]


def mcnemar_exact(b: int, c: int):
    """Exact two-sided McNemar p from discordant counts b, c."""
    n = b + c
    if n == 0:
        return 1.0
    from math import comb
    m = min(b, c)
    tail = sum(comb(n, i) for i in range(0, m + 1)) / (2 ** n)
    return round(min(1.0, 2 * tail), 4)


def three_way(df: pd.DataFrame) -> dict:
    """certify / explicit-reject / not-testable per condition for the gate."""
    out = {}
    for cond, s in df.groupby("condition"):
        cert = int((s["discreteness_pass"] == True).sum())  # noqa: E712
        testable = s["discreteness_testable"] == True        # noqa: E712
        rej = int(((s["discreteness_pass"] == False) & testable).sum())  # noqa: E712
        abst = int((~testable).sum())
        out[cond] = {"n": int(len(s)), "certify": cert,
                     "explicit_reject": rej, "not_testable": abst}
    return out


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(here, "../results/gate_ablation_raw.csv"))
    ap.add_argument("--sweep", default=os.path.join(here, "../results/clusterer_sweep_raw.csv"))
    ap.add_argument("--out", default=os.path.join(here, "../results"))
    args = ap.parse_args()

    d = pd.read_csv(args.raw)
    d["is_pos"] = d["condition"].isin(POS)
    d["testable"] = d["discreteness_testable"] == True   # noqa: E712
    d["gate_certify"] = d["discreteness_pass"] == True    # noqa: E712
    d["sig_alone"] = d["sg_empirical_p"] < ALPHA          # SigClust null alone
    d["stab_certify"] = d["stability_pass"] == True       # noqa: E712

    pos, neg = d[d.is_pos], d[~d.is_pos]

    # --- 1. three-way accounting ---
    tw = three_way(d)

    # --- gate metrics, abstentions excluded from FPR ---
    neg_testable = neg[neg.testable]
    gate = {
        "TPR": {"k": int(pos.gate_certify.sum()), "n": int(len(pos)),
                "rate": round(pos.gate_certify.mean(), 4),
                "wilson95": wilson(int(pos.gate_certify.sum()), len(pos))},
        "negatives_total": int(len(neg)),
        "negatives_abstained": int((~neg.testable).sum()),
        "negatives_abstain_rate": round((~neg.testable).mean(), 4),
        "negatives_testable": int(len(neg_testable)),
        "FPR_excluding_abstentions": {
            "k": int(neg_testable.gate_certify.sum()),
            "n": int(len(neg_testable)),
            "rate": (round(neg_testable.gate_certify.mean(), 4)
                     if len(neg_testable) else None),
            "wilson95": wilson(int(neg_testable.gate_certify.sum()), len(neg_testable)),
            "note": ("denominator is tiny because the gate abstains on almost all "
                     "negatives; FPR here is nearly uninformative — that is the "
                     "honest finding, not FPR=0")},
        "FPR_abstention_as_reject": {
            "k": int(neg.gate_certify.sum()), "n": int(len(neg)),
            "rate": round(neg.gate_certify.mean(), 4),
            "wilson95": wilson(int(neg.gate_certify.sum()), len(neg)),
            "note": "the original (over-generous) convention, shown for contrast"},
    }

    # --- stability-only baseline ---
    stab = {
        "TPR": {"k": int(pos.stab_certify.sum()), "n": int(len(pos)),
                "rate": round(pos.stab_certify.mean(), 4),
                "wilson95": wilson(int(pos.stab_certify.sum()), len(pos))},
        "FPR": {"k": int(neg.stab_certify.sum()), "n": int(len(neg)),
                "rate": round(neg.stab_certify.mean(), 4),
                "wilson95": wilson(int(neg.stab_certify.sum()), len(neg))},
    }

    # --- 2. paired McNemar: stability-only vs full gate, on the NEGATIVES ---
    # (abstention counted as non-certify, i.e. a pass on a negative = good)
    neg_stab = neg.stab_certify.values
    neg_gate = neg.gate_certify.values
    b = int(((neg_stab) & (~neg_gate)).sum())   # stability certifies, gate doesn't
    c = int(((~neg_stab) & (neg_gate)).sum())   # gate certifies, stability doesn't
    mcnemar_fpr = {"discordant_stab_only": b, "discordant_gate_only": c,
                   "exact_p": mcnemar_exact(b, c),
                   "interpretation": ("large b, ~0 c, small p => the gate removes "
                                      "false certifications the stability gate made")}
    # on POSITIVES: did the gate lose true positives?
    pos_stab = pos.stab_certify.values
    pos_gate = pos.gate_certify.values
    b2 = int(((pos_stab) & (~pos_gate)).sum())
    c2 = int(((~pos_stab) & (pos_gate)).sum())
    mcnemar_tpr = {"discordant_stab_only": b2, "discordant_gate_only": c2,
                   "exact_p": mcnemar_exact(b2, c2),
                   "interpretation": ("p ~ 1 with 1 discordant pair => the ~3% TPR "
                                      "difference is not statistically distinguishable "
                                      "from zero; do not claim a 'cost'")}

    # --- 3. head-to-head: SigClust p ALONE vs the full composite gate ---
    head = {
        "sigclust_p_alone": {
            "TPR": {"k": int(pos.sig_alone.sum()), "n": int(len(pos)),
                    "rate": round(pos.sig_alone.mean(), 4)},
            "FPR": {"k": int(neg.sig_alone.sum()), "n": int(len(neg)),
                    "rate": round(neg.sig_alone.mean(), 4)}},
        "full_composite_gate": {
            "TPR": {"k": int(pos.gate_certify.sum()), "n": int(len(pos)),
                    "rate": round(pos.gate_certify.mean(), 4)},
            "FPR_abstention_as_reject": {"k": int(neg.gate_certify.sum()),
                                         "n": int(len(neg)),
                                         "rate": round(neg.gate_certify.mean(), 4)}},
        "increment": ("The single-Gaussian p-value alone reproduces the composite "
                      "gate's TPR and FPR exactly. The gap statistic and dip test "
                      "are computed but not in the decision rule; the k-stability "
                      "routing only adds abstentions. The honest contribution is the "
                      "null recalibration (feature-permutation -> single-Gaussian), "
                      "not a new multi-test instrument."),
    }

    # p-value floor disclosure
    ref_n = 120
    floor = 1.0 / (ref_n + 1)
    cert_ps = sorted(d[d.gate_certify]["sg_empirical_p"].unique().tolist())
    floor_note = {
        "n_ref": ref_n, "p_floor": round(floor, 4),
        "unique_p_among_certified": cert_ps,
        "all_certified_at_floor": bool(len(cert_ps) == 1 and abs(cert_ps[0] - floor) < 1e-4),
        "note": ("every certified positive sits at the p-floor; report as "
                 "p <= 1/(n_ref+1), and see the separation sweep for whether the "
                 "gate has resolution between borderline and obvious structure")}

    out = {
        "three_way_accounting": tw,
        "gate": gate,
        "stability_only_baseline": stab,
        "paired_mcnemar": {"on_negatives_FPR": mcnemar_fpr,
                           "on_positives_TPR": mcnemar_tpr},
        "head_to_head": head,
        "p_value_floor": floor_note,
    }
    path = os.path.join(args.out, "ablation_honest.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)

    # console
    print("THREE-WAY (certify / reject / not-testable):")
    for cond, v in tw.items():
        print(f"  {cond:16s} certify {v['certify']:2d}  reject {v['explicit_reject']:2d}"
              f"  not-testable {v['not_testable']:2d}  (n={v['n']})")
    print(f"\nGate: TPR {gate['TPR']['rate']} {gate['TPR']['wilson95']}")
    print(f"  negatives: {gate['negatives_abstained']}/{gate['negatives_total']} "
          f"abstained ({gate['negatives_abstain_rate']*100:.0f}%)")
    print(f"  FPR excl. abstentions: {gate['FPR_excluding_abstentions']['rate']} "
          f"on n={gate['FPR_excluding_abstentions']['n']} "
          f"{gate['FPR_excluding_abstentions']['wilson95']}")
    print(f"Stability-only FPR: {stab['FPR']['rate']} {stab['FPR']['wilson95']}")
    print(f"\nMcNemar negatives (FPR): b={mcnemar_fpr['discordant_stab_only']} "
          f"c={mcnemar_fpr['discordant_gate_only']} p={mcnemar_fpr['exact_p']}")
    print(f"McNemar positives (TPR): b={mcnemar_tpr['discordant_stab_only']} "
          f"c={mcnemar_tpr['discordant_gate_only']} p={mcnemar_tpr['exact_p']}")
    print(f"\nHead-to-head: SigClust-alone TPR {head['sigclust_p_alone']['TPR']['rate']} "
          f"FPR {head['sigclust_p_alone']['FPR']['rate']}  ==  "
          f"full gate TPR {head['full_composite_gate']['TPR']['rate']} "
          f"FPR {head['full_composite_gate']['FPR_abstention_as_reject']['rate']}")
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
