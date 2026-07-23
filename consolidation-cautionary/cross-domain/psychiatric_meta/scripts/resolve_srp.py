#!/usr/bin/env python3
"""Reliably resolve GSE -> SRP via direct SRA search (the GDS extrelations field
is unreliable and undercounts). Updates candidate_studies.tsv keepers in place
with a resolved SRP where missing.
"""
import json
import os
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
KEY = os.environ.get("NCBI_API_KEY", "")
EMAIL = "rohit.chauhan@topmist.com"
_HERE = os.path.dirname(os.path.abspath(__file__))
TSV = os.path.join(_HERE, "../results/candidate_studies.tsv")


def _get(endpoint, **params):
    params.update(retmode="json", email=EMAIL)
    if KEY:
        params["api_key"] = KEY
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception:
            if a == 3:
                return None
            time.sleep(2 ** a)


def srp_for(gse):
    """GSE -> SRP via SRA: esearch sra for the GSE, esummary first hit -> study acc."""
    es = _get("esearch.fcgi", db="sra", term=f"{gse}[All Fields]", retmax=1)
    ids = (es or {}).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return ""
    summ = _get("esummary.fcgi", db="sra", id=ids[0])
    if not summ:
        return ""
    raw = json.dumps(summ)
    import re
    m = re.search(r"(SRP\d+)", raw)
    return m.group(1) if m else ""


def main():
    rows = [l.rstrip("\n").split("\t") for l in open(TSV)]
    header, data = rows[0], rows[1:]
    gi, si, ki = header.index("gse"), header.index("srp"), header.index("keep")
    resolved = 0
    for r in data:
        if r[ki] == "True" and not r[si]:
            srp = srp_for(r[gi])
            if srp:
                r[si] = srp
                resolved += 1
                print(f"  {r[gi]} -> {srp}")
            time.sleep(0.15)
    with open(TSV, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in data:
            fh.write("\t".join(r) + "\n")
    print(f"resolved {resolved} additional SRPs; updated {TSV}")


if __name__ == "__main__":
    main()
