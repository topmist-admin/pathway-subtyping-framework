#!/usr/bin/env python3
"""
Confound-remapping demonstration: the SAME confounded-subtype detector, three domains.

The autism/psychiatry reproduction showed the framework's Gate 6 (Confound
Association Gate) catching a partition that was really a *brain-region*
classifier. The core insight for cross-disorder / cancer expansion is that
Gate 6 does not change — only its CONFOUND REGISTRY changes. The gate's job is
always "is this partition explained by a nuisance axis rather than the biology
of interest?"; the disease sets what the nuisance axis IS.

This script demonstrates the remap on SYNTHETIC partitions with a planted
nuisance structure, so a developer can see the gate behave identically across
domains before wiring real cohorts. It uses the installed framework gate
(pathway-subtyping==0.7.0) exactly as the autism Gate-6 demo did.

DOMAIN -> dominant nuisance confound (what the naive partition recovers):
  psychiatry (multi-region brain)  -> brain region        (shown in autism work)
  cancer (pan-cancer pooled)       -> tissue of origin      + tumor purity
  cancer (single tumor type)       -> tumor purity / stromal-immune infiltration

BIOLOGY-OF-INTEREST (exempt from the gate, per v0.7.0 keys diagnosis/dx/
disease/condition/phenotype): the disease label in psychiatry; the cancer
type or the known molecular subtype (e.g. CMS) in oncology. NOTE the exempt
key must be present under one of those names for the gate to spare it — see
`--help` and the remap table in docs/01_cross_disorder_cancer_roadmap.md.

Deterministic (seeded). Requires pathway-subtyping==0.7.0, numpy.
"""
import argparse
import json
import os

import numpy as np
from pathway_subtyping.validation import ValidationGates


def planted(n, k_part, nuisance_levels, align, seed):
    """Make a partition that is `align`-correlated with a nuisance variable."""
    rng = np.random.default_rng(seed)
    nuis = rng.integers(0, nuisance_levels, n)
    # partition follows nuisance with prob `align`, else random
    part = np.where(rng.random(n) < align, nuis % k_part, rng.integers(0, k_part, n))
    return part, nuis


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=141)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="out_confound_remap")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    g = ValidationGates(seed=args.seed, show_progress=False)

    scenarios = {
        # domain label : (n_partition_clusters, nuisance name, nuisance levels, alignment)
        "psychiatry_region":      (3, "brain region", 3, 0.92),   # ~ the autism/SCZ case
        "cancer_pan_tissue":      (4, "tissue of origin", 5, 0.90),
        "cancer_single_purity":   (2, "tumor purity bin", 4, 0.55),  # weaker, realistic
    }
    results = {}
    rng = np.random.default_rng(args.seed)
    for name, (kk, nuis_name, lv, align) in scenarios.items():
        part, nuis = planted(args.n, kk, lv, align, args.seed)
        # biology-of-interest: a diagnosis/subtype label that is (correctly) NOT the driver
        dx = rng.integers(0, 2, args.n)
        confounds = {nuis_name: nuis, "diagnosis": dx}
        res = g.confound_association_gate(part, confounds, cramers_v_max=0.30, alpha=0.05)
        results[name] = {
            "nuisance": nuis_name, "planted_alignment": align,
            "passed": bool(res.passed),
            "max_nuisance_cramers_v": float(res.metric_value),
            "worst_confound": res.details.get("worst_confound"),
            "failing_confounds": res.details.get("failing_confounds"),
        }
        print(f"{name:24} nuisance={nuis_name:18} V={res.metric_value:.3f} "
              f"passed={res.passed} worst={res.details.get('worst_confound')}")

    json.dump(results, open(f"{args.out}/confound_remap_results.json", "w"), indent=2)
    print(f"\nsaved -> {args.out}/confound_remap_results.json")
    print("\nInterpretation: the gate flags the domain's dominant nuisance axis in each "
          "case\nwhile exempting 'diagnosis' — identical logic, remapped registry.")


if __name__ == "__main__":
    main()
