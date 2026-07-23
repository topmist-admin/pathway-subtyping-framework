#!/usr/bin/env python3
"""Enumerate PUBLIC psychiatric postmortem-brain bulk RNA-seq GEO series and map
each to its SRA study (SRP), for the Track-A recount3 meta-cohort.

Reproducible NCBI E-utilities search + curation. Writes a candidate table; a
companion R step (check_recount3.R) then flags which SRPs are in recount3.

Curation rules (bulk postmortem RNA-seq only):
  KEEP  gdsType contains "Expression profiling by high throughput sequencing"
  DROP  titles/summaries indicating iPSC / hiPSC / organoid / single-cell /
        single-nucleus / snRNA / scRNA (not bulk postmortem tissue)

Env: reads NCBI_API_KEY if set (raises the rate ceiling); works without it.
Output: results/candidate_studies.tsv  (gse, n_samples, srp, gdstype, keep, reason)
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
KEY = os.environ.get("NCBI_API_KEY", "")
EMAIL = "rohit.chauhan@topmist.com"
QUERY = ('(schizophrenia[Title] OR bipolar[Title] OR "major depress"[Title] OR '
         'psychiatric[Title]) AND brain[All Fields] AND '
         '"Expression profiling by high throughput sequencing"[DataSet Type] AND '
         '"Homo sapiens"[Organism] AND gse[Filter]')
DROP_TERMS = ("ipsc", "hipsc", "organoid", "single-cell", "single cell",
              "single-nucleus", "single nucleus", "snrna", "scrna", "sc-rna")
_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "../results/candidate_studies.tsv")


def _get(endpoint, **params):
    params.update(retmode="json", email=EMAIL)
    if KEY:
        params["api_key"] = KEY
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    return None


def main():
    es = _get("esearch.fcgi", db="gds", term=QUERY, retmax=200)
    uids = es["esearchresult"]["idlist"]
    print(f"{len(uids)} candidate GEO series")
    rows = []
    # batch esummary in chunks
    for i in range(0, len(uids), 50):
        chunk = uids[i:i + 50]
        summ = _get("esummary.fcgi", db="gds", id=",".join(chunk))["result"]
        for uid in chunk:
            d = summ.get(uid, {})
            gse = "GSE" + str(d.get("accession", "")).replace("GSE", "")
            gse = d.get("accession", "")
            n = int(d.get("n_samples", 0) or 0)
            gtype = d.get("gdstype", "")
            title = (d.get("title", "") + " " + d.get("summary", "")).lower()
            srp = ""
            for rel in d.get("extrelations", []):
                tgt = rel.get("targetobject", "")
                if tgt.startswith("SRP"):
                    srp = tgt
            keep, reason = True, ""
            if "expression profiling by high throughput sequencing" not in gtype.lower():
                keep, reason = False, "not-RNAseq"
            elif any(t in title for t in DROP_TERMS):
                keep = False
                reason = "iPSC/single-cell"
            rows.append((gse, n, srp, gtype, keep, reason))
        time.sleep(0.2)

    rows.sort(key=lambda r: (-int(r[4]), -r[1]))  # keepers first, then by n
    with open(OUT, "w") as fh:
        fh.write("gse\tn_samples\tsrp\tgdstype\tkeep\treason\n")
        for gse, n, srp, gt, keep, reason in rows:
            fh.write(f"{gse}\t{n}\t{srp}\t{gt}\t{keep}\t{reason}\n")

    keepers = [r for r in rows if r[4]]
    with_srp = [r for r in keepers if r[2]]
    print(f"  bulk-RNA-seq keepers: {len(keepers)}  (with SRP link: {len(with_srp)})")
    print(f"  total keeper samples (GEO-reported): {sum(r[1] for r in keepers)}")
    print("  top keepers:")
    for gse, n, srp, gt, keep, reason in keepers[:12]:
        print(f"    {gse:12s} n={n:4d} {srp or '(no SRP)'}")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
