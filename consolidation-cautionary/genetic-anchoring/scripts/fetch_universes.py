#!/usr/bin/env python3
"""
Fetch the two background universes + the test-set symbol->Ensembl map for
gate7_feature_level_anchoring.py, from Ensembl BioMart (one bulk call each).

Writes to --out:
  protein_coding_ensembl.json   list of protein-coding Ensembl gene IDs (genome-wide null)
  testset_symbol2ens.json       {symbol: ensembl_id} for the neuronal+glial test-set genes

The brain-expressed universe is NOT fetched here — it is derived at run time
from GSE80655 inside gate7_feature_level_anchoring.py (genes detected in >=50%
of samples), so it stays tied to the actual expression input.

Requires: requests. Network: www.ensembl.org (BioMart).
"""
import argparse
import json
import os

import requests

BIOMART = "https://www.ensembl.org/biomart/martservice"
PC_XML = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
          '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" '
          'count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" '
          'interface="default"><Filter name="biotype" value="protein_coding"/>'
          '<Attribute name="ensembl_gene_id"/></Dataset></Query>')
MAP_XML = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
           '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" '
           'count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" '
           'interface="default"><Attribute name="ensembl_gene_id"/>'
           '<Attribute name="hgnc_symbol"/></Dataset></Query>')


def read_gmt(path):
    d = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            d[f[0]] = [g for g in f[2:] if g]
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gmt", required=True, help="schizophrenia_pathways.gmt")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pc = [x.strip() for x in requests.get(BIOMART, params={"query": PC_XML}, timeout=180)
          .text.splitlines() if x.startswith("ENSG")]
    json.dump(pc, open(f"{args.out}/protein_coding_ensembl.json", "w"))
    print(f"protein-coding universe: {len(pc)}")

    sym2ens = {}
    for line in requests.get(BIOMART, params={"query": MAP_XML}, timeout=180).text.splitlines():
        p = line.split("\t")
        if len(p) == 2 and p[0].startswith("ENSG") and p[1]:
            sym2ens[p[1]] = p[0]

    scz = read_gmt(args.gmt)
    syms = sorted(set(scz["SYNAPTIC_TRANSMISSION"]) | set(scz["GLUTAMATE_SIGNALING"]) |
                  set(scz["GABA_SIGNALING"]) | set(scz["IMMUNE_COMPLEMENT"]) |
                  {"GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CX3CR1",
                   "AIF1", "CSF1R", "MOG", "MOBP", "SLC1A3", "S100B"})
    m = {g: sym2ens.get(g) for g in syms}
    json.dump(m, open(f"{args.out}/testset_symbol2ens.json", "w"), indent=0)
    print(f"test-set symbols mapped: {sum(1 for v in m.values() if v)}/{len(syms)}")
    print(f"saved -> {args.out}/protein_coding_ensembl.json, testset_symbol2ens.json")


if __name__ == "__main__":
    main()
