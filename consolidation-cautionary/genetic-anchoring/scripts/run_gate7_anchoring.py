#!/usr/bin/env python3
"""
Gate 7 — feature-level genetic-anchoring positive control (Voineagu calibration).

REFACTORED DRIVER. The enrichment logic no longer lives in this script: it is now
the shipped framework gate ``ValidationGates.genetic_anchoring_gate`` backed by
``pathway_subtyping.genetics`` (added on the v0.8.0 line). This driver only does
I/O — build the risk set, the two gene sets, and the background-matched null — and
then calls the framework. The earlier standalone prototype
(``psf-repro-bundle .../gate7_feature_level_anchoring.py``) reimplemented the
hypergeometric test inline; this replaces it, exercising the same code a third
party gets from the package.

Ground-truth calibration (design audit s0.5): a NEURONAL/SYNAPTIC gene set is
enriched for ASD genetic risk while a GLIAL/IMMUNE set is not, under a
brain-expressed-matched null. If the pipeline cannot recover this known answer
(Voineagu 2011 F5, from public data), its verdict on any novel subtype is
worthless. The null must be background-matched, not genome-wide: brain-expressed
genes are already enriched for ASD risk, so a genome-wide null over-states
enrichment for a region-identity axis. Both nulls are reported; the matched null
decides the gate.

INPUTS (all deposited under ../inputs, so this runs offline; re-fetch with the
fetch_*.py utilities if you want to rebuild them):
  --risk-json      ../inputs/asd_risk_ensembl.json      (Open Targets, MONDO_0005258)
  --pc-json        ../inputs/protein_coding_ensembl.json (genome-wide reference null)
  --symbol-map     ../inputs/testset_symbol2ens.json     (test-set symbol -> Ensembl)
  --gse80655-expr  GSE80655 expression .txt.gz           (brain-expressed matched null)
  --gmt            ../../panels/schizophrenia_pathways.gmt (cell-type/pathway panels)
  --out            output dir (default: ../results)

Requires: pathway-subtyping (>=0.8.0 line, for the genetics subpackage), numpy,
pandas, scipy.
"""
import argparse
import gzip
import io
import json
import os

import pandas as pd

from pathway_subtyping.genetics import feature_level_anchoring
from pathway_subtyping.validation import ValidationGates


def read_gmt(path):
    d = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = [g for g in f[2:] if g]
    return d


def brain_expressed_universe(expr_path, min_frac=0.5):
    """Genes detected in >= min_frac of GSE80655 samples (Ensembl, version-stripped)."""
    raw = gzip.open(expr_path, "rt").read()
    e = pd.read_csv(io.StringIO(raw), sep="\t", index_col=0)
    e.index = e.index.astype(str).str.split(".").str[0]
    return set(e[(e > 0).mean(axis=1) >= min_frac].index)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)  # consolidation-cautionary/genetic-anchoring
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--risk-json", default=os.path.join(root, "inputs/asd_risk_ensembl.json"))
    ap.add_argument("--pc-json", default=os.path.join(root, "inputs/protein_coding_ensembl.json"))
    ap.add_argument("--symbol-map", default=os.path.join(root, "inputs/testset_symbol2ens.json"))
    ap.add_argument("--gse80655-expr", required=True,
                    help="GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz")
    ap.add_argument("--gmt", default=os.path.join(root, "../panels/schizophrenia_pathways.gmt"))
    ap.add_argument("--out", default=os.path.join(root, "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    risk = set(json.load(open(args.risk_json)))
    pc = set(json.load(open(args.pc_json)))
    brain = brain_expressed_universe(args.gse80655_expr)
    m = json.load(open(args.symbol_map))
    scz = read_gmt(args.gmt)

    neuronal_syms = sorted(set(scz["SYNAPTIC_TRANSMISSION"]) |
                           set(scz["GLUTAMATE_SIGNALING"]) |
                           set(scz["GABA_SIGNALING"]))
    glial_syms = sorted(set(scz["IMMUNE_COMPLEMENT"]) |
                        {"GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CX3CR1",
                         "AIF1", "CSF1R", "MOG", "MOBP", "SLC1A3", "S100B"})
    gene_sets = {
        "NEURONAL/SYNAPTIC": {m[g] for g in neuronal_syms if m.get(g)},
        "GLIAL/IMMUNE": {m[g] for g in glial_syms if m.get(g)},
    }

    # --- framework call: raw enrichment (both nulls, for the "null matters" table) ---
    enrich = feature_level_anchoring(gene_sets, risk, brain, reference_universe=pc)
    print("=== Gate 7 feature-level anchoring (framework: pathway_subtyping.genetics) ===")
    print("    neuronal SHOULD enrich for ASD risk; glial should NOT\n")
    for label in ("NEURONAL/SYNAPTIC", "GLIAL/IMMUNE"):
        for null in ("background-matched", "genome-wide"):
            r = enrich[f"{label}|{null}"]
            star = ("***" if r.p_value < 1e-3 else "**" if r.p_value < 1e-2
                    else "*" if r.p_value < 5e-2 else "n.s.")
            print(f"{label:20} vs {null:20}: {r.risk_hits}/{r.testset_n} hits, "
                  f"exp {r.expected:.2f}, fold={r.fold:.2f}, p={r.p_value:.2e}  {star}")
        print()

    # --- framework call: the gate verdict (matched null decides) ---
    gate = ValidationGates(show_progress=False).genetic_anchoring_gate(
        subtype_gene_sets=gene_sets,
        risk_genes=risk,
        background_universe=brain,
        reference_universe=pc,
    )
    print(f"GATE: {gate.name}")
    print(f"passed          : {gate.passed}")
    print(f"anchored subtypes: {gate.details['anchored_subtypes']}")
    print(f"max enrichment fold: {gate.metric_value:.2f} (threshold {gate.comparison} "
          f"{gate.threshold})")

    out = {
        "enrichment": {k: v.to_dict() for k, v in enrich.items()},
        "gate": gate.to_dict(),
        "discrimination_confirmed": gate.passed and "NEURONAL/SYNAPTIC" in gate.details["anchored_subtypes"],
    }
    dest = os.path.join(args.out, "gate7_genetic_anchoring_poscontrol.json")
    json.dump(out, open(dest, "w"), indent=2)
    print(f"\nsaved -> {dest}")


if __name__ == "__main__":
    main()
