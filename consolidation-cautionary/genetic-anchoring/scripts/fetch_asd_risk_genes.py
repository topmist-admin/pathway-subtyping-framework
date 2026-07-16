#!/usr/bin/env python3
"""
Fetch the ASD genetic-risk reference gene set from the Open Targets Platform.

Open Targets is used as an independent, genetics-grounded ASD risk-gene
catalogue (SFARI Gene is the canonical source but its download endpoint was
allowlist-blocked in the environment this was built in; Open Targets'
`genetic_association` datatype is derived from GWAS + rare-variant evidence and
is causally-prior in exactly the same sense).

Writes two files to --out:
  asd_risk_ensembl.json        list of Ensembl gene IDs with genetic_association>0
  opentargets_asd_targets.json full pull (symbol, overall score, genetic_assoc)

Requires: requests. Network: api.platform.opentargets.org.
"""
import argparse
import json
import os

import requests

API = "https://api.platform.opentargets.org/api/v4/graphql"
EFO = "MONDO_0005258"  # autism spectrum disorder

QUERY = """
query ($efo:String!, $idx:Int!, $size:Int!){
  disease(efoId:$efo){
    id name
    associatedTargets(page:{index:$idx,size:$size}, orderByScore:"score"){
      count
      rows{ target{id approvedSymbol} score datatypeScores{id score} }
    }
  }
}"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=".")
    ap.add_argument("--efo", default=EFO)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rows, idx, size, total = [], 0, 500, None
    while True:
        d = requests.post(API, json={"query": QUERY,
                                     "variables": {"efo": args.efo, "idx": idx, "size": size}},
                          timeout=60).json()
        at = d["data"]["disease"]["associatedTargets"]
        total = at["count"]
        rows += at["rows"]
        idx += 1
        if idx * size >= total:
            break

    def ga(r):
        return next((x["score"] for x in r["datatypeScores"]
                     if x["id"] == "genetic_association"), 0.0)

    risk_ens = sorted({r["target"]["id"] for r in rows if ga(r) > 0})
    full = [(r["target"]["approvedSymbol"], round(r["score"], 4), round(ga(r), 4))
            for r in rows]
    json.dump(risk_ens, open(f"{args.out}/asd_risk_ensembl.json", "w"))
    json.dump({"efo": args.efo, "total": total, "genes": full},
              open(f"{args.out}/opentargets_asd_targets.json", "w"), indent=0)
    print(f"disease total targets: {total}")
    print(f"genetic_association>0 (risk set): {len(risk_ens)}")
    print(f"saved -> {args.out}/asd_risk_ensembl.json, opentargets_asd_targets.json")


if __name__ == "__main__":
    main()
