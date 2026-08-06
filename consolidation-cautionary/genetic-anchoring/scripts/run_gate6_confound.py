#!/usr/bin/env python3
"""
Gate 6 — Confound Association Gate on the deposited partition.

Thin driver over the shipped ``ValidationGates.confound_association_gate``. Applies
the gate to the deposited SCZ+Control partition (39/92/10; framework v0.3.0, seed
42) and shows it FAILS on brain region while exempting diagnosis. This is the gate
whose absence let the anatomy artifact pass the original stability battery — it is now a
first-class framework gate, not a bolt-on. (An earlier revision quoted "bootstrap ARI
~0.914" here. No deposited artifact contains that number, and the deposit elsewhere
enumerates the values on record for this quantity; an unsourced fourth one does not belong
in a docstring. For the canonical list see cross-domain/flagship_stats/README.md.)

Reproduces (deposited reference: ../results/gate6_confound_association.json):
  passed=False, max nuisance Cramer's V=0.66 (brain region, chi2=125.1, p~4e-26),
  diagnosis exempt (biology-of-interest, V=0.0, p=0.408).

INPUT:
  --labels  deposited partition CSV with per-sample `subtype`, 3-level
            `brain region`, `diagnosis`, `gender`, `ethnicity`. On main this is
            ../data/partition/sample_metadata_with_subtypes.csv (parent harness).
  --out     output dir (default: ../results).

Requires: pathway-subtyping, numpy, pandas.
"""
import argparse
import json
import os

import pandas as pd

from pathway_subtyping.validation import ValidationGates


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--labels",
                    default=os.path.join(root, "../data/partition/sample_metadata_with_subtypes.csv"))
    ap.add_argument("--out", default=os.path.join(root, "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(args.labels, index_col=0)
    labels = df["subtype"].to_numpy()
    confounds = {
        "brain region": df["brain region"].to_numpy(),  # nuisance (anatomy)
        "diagnosis": df["diagnosis"].to_numpy(),         # biology-of-interest (exempt)
        "gender": df["gender"].to_numpy(),
        "ethnicity": df["ethnicity"].to_numpy(),
    }
    res = ValidationGates(seed=42, show_progress=False).confound_association_gate(
        labels, confounds, cramers_v_max=0.30, alpha=0.05
    )

    print(f"GATE: {res.name}")
    print(f"passed           : {res.passed}")
    print(f"{res.metric_name}  : {res.metric_value:.4f} "
          f"(threshold {res.comparison} {res.threshold})")
    print(f"worst_confound   : {res.details['worst_confound']}")
    print(f"failing_confounds: {res.details['failing_confounds']}")
    print(f"interpretation   : {res.details['interpretation']}")

    dest = os.path.join(args.out, "gate6_confound_association.json")
    json.dump(res.to_dict(), open(dest, "w"), indent=2, default=str)
    print(f"\nsaved -> {dest}")


if __name__ == "__main__":
    main()
