#!/usr/bin/env python3
"""Build the Zenodo reproduction-deposit archive for the cautionary-framework paper.

This is the PERMANENT, CITABLE archive a journal reviewer reaches from the paper's
Data & Code Availability section. Unlike the pre-release stress-test bundle
(build_full_reviewer_bundle.py), it:

  - is self-contained and public-safe — it does NOT include the withdrawn prior
    manuscript or the raw journal reviewer comments;
  - mirrors the repo layout so every reproduction script runs unchanged;
  - ships the v0.8.0 framework source (belt-and-braces permanence: the primary
    install path is `pip install pathway-subtyping==0.8.0`, but the source is here
    too so the deposit reproduces even if PyPI/GitHub change);
  - includes the concern-handling matrix (each Scientific Reports comment ->
    resolution -> artifact) and a deposit README + Zenodo metadata.

Usage: python build_zenodo_deposit.py --out DIR [--stamp YYYY-MM-DD]
       [--concern-matrix PATH]   (the REVIEWER-CONCERN-HANDLING-MATRIX.md)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tarfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

SKIP_DIR = {"__pycache__", ".ipynb_checkpoints", ".git", ".pytest_cache",
            ".mypy_cache", "node_modules", ".venv", "pathwayenv"}
SKIP_EXT = {".pyc", ".pyo", ".so", ".o"}
SKIP_NAME = {".DS_Store"}

# repo-relative layout so scripts' ../../../../ paths resolve; pip install -e . at root
ITEMS = [
    ("src", "src"),
    ("tests", "tests"),
    ("pyproject.toml", "pyproject.toml"),
    ("CITATION.cff", "CITATION.cff"),
    ("CHANGELOG.md", "CHANGELOG.md"),
    ("README.md", "README-framework.md"),
    ("consolidation-cautionary", "consolidation-cautionary"),
    ("research-results/GSE28521", "research-results/GSE28521"),
    ("research-results/GSE80655", "research-results/GSE80655"),
    ("CORRECTION_2026-07", "CORRECTION_2026-07"),
    ("scripts/gate_ablation_study.py", "scripts/gate_ablation_study.py"),
    ("scripts/plot_gate_ablation.py", "scripts/plot_gate_ablation.py"),
]

VERSION = "0.8.0"

DEPOSIT_README = """# Reproduction package — cautionary molecular-subtype validation framework

Permanent, self-contained reproduction package for the Scientific Reports manuscript
on validation gates for molecular patient stratification (the discreteness / conservative
screen and the cross-disease audit). Everything here uses **public data only**; no
controlled-access data is included or used.

## Install and reproduce (verified from a clean environment)

```
python -m venv .venv && . .venv/bin/activate
pip install pathway-subtyping==%(version)s        # the framework (also on PyPI)
pip install -r consolidation-cautionary/requirements.txt
pip install torch                                  # only for the DL baselines
python -c "import pathway_subtyping as p; print(p.__version__)"   # %(version)s
```

If PyPI is unavailable, the framework source is bundled here — install it instead:
`pip install -e .` at the archive root.

Then follow **`consolidation-cautionary/RUNME.md`**, which indexes every result to a
runnable script and its deposited reference output. The no-network results reproduce in
minutes; the rest fetch public cBioPortal / GEO / recount3 data with no authentication.

## What is here
- `src/`, `tests/`, `pyproject.toml`  the v0.8.0 framework source (the discreteness gate + DL baselines)
- `consolidation-cautionary/`         reproduction packages (scripts + deposited results + READMEs); `RUNME.md` is the index
- `research-results/`, `CORRECTION_2026-07/`  cached public data + the corrected benchmark
- `REVIEWER-CONCERN-HANDLING-MATRIX.md`  each reviewer comment -> resolution -> the artifact that demonstrates it
- `MANIFEST.txt`                      every file with its SHA-256

## Cite
Framework: `pip install pathway-subtyping==%(version)s` · RRID:SCR_028051 ·
GitHub/Codeberg tag `v%(version)s`. This deposit: the Zenodo DOI shown on the record.
""" % {"version": VERSION}


def zenodo_metadata(stamp):
    return {
        "metadata": {
            "title": ("Reproduction package: a discreteness-aware validation "
                      "framework for molecular patient stratification (v%s)" % VERSION),
            "upload_type": "software",
            "description": (
                "Self-contained, public-data-only reproduction package for the "
                "Scientific Reports manuscript on mandatory validation gates for "
                "molecular subtyping. Contains the v%s framework source (installable "
                "from PyPI as pathway-subtyping==%s), the cross-disease reproduction "
                "packages (synthetic gate ablation, within-study real-data calibration, "
                "corrected 47-dataset benchmark audit, TCGA-BRCA/CPTAC cancer worked "
                "example, GTEx large-N, and the postmortem-psychiatric flagship), the "
                "cached public inputs, and a concern-handling matrix mapping each "
                "reviewer comment to its resolution and artifact. Reproduction verified "
                "from a clean PyPI install; see RUNME.md."
                % (VERSION, VERSION)),
            "version": VERSION,
            "language": "eng",
            "keywords": ["molecular subtyping", "cluster validation",
                         "reproducibility", "bioinformatics", "precision medicine",
                         "SigClust", "discreteness"],
            "access_right": "open",
            "license": "MIT",
            "creators": [
                {"name": "Chauhan, Rohit", "affiliation": "Topmist LLC",
                 "orcid": "0009-0003-9895-4629"},
                {"name": "Chauhan, Mohit",
                 "affiliation": "Mayo Clinic Jacksonville",
                 "orcid": "0000-0001-8848-4385"},
            ],
            "related_identifiers": [
                {"identifier": "https://pypi.org/project/pathway-subtyping/%s/" % VERSION,
                 "relation": "isSupplementedBy", "scheme": "url"},
                {"identifier": "https://github.com/topmist-admin/pathway-subtyping-framework",
                 "relation": "isDerivedFrom", "scheme": "url"},
                # Concept DOI, not the v2.0 version DOI: this relation means "the
                # software references that dataset", and it should keep pointing at
                # the current version after the benchmark record is re-issued.
                {"identifier": "10.5281/zenodo.19323753",
                 "relation": "references", "scheme": "doi"},
                {"identifier": "RRID:SCR_028051",
                 "relation": "isIdenticalTo", "scheme": "url"},
            ],
            "notes": ("Public data only; no controlled-access data included. Built %s. "
                      "Primary install path is pip install pathway-subtyping==%s; the "
                      "bundled source is a permanence fallback." % (stamp, VERSION)),
        }
    }


def keep(path):
    parts = set(path.split(os.sep))
    if parts & SKIP_DIR:
        return False
    if os.path.basename(path) in SKIP_NAME:
        return False
    return os.path.splitext(path)[1] not in SKIP_EXT


def collect(items):
    out = []
    for src_rel, dest_rel in items:
        src = os.path.join(REPO, src_rel)
        if not os.path.exists(src):
            print("  ! missing:", src_rel, file=sys.stderr)
            continue
        if os.path.isfile(src):
            if keep(src):
                out.append((src, dest_rel))
            continue
        for r, dirs, files in os.walk(src):
            dirs[:] = [d for d in dirs if d not in SKIP_DIR]
            for f in files:
                full = os.path.join(r, f)
                if keep(full):
                    out.append((full, os.path.join(dest_rel, os.path.relpath(full, src))))
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stamp", default="2026-07-25")
    ap.add_argument("--concern-matrix", default="")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pairs = collect(ITEMS)

    # generated files (temp -> tar -> cleanup)
    temp = []
    readme_p = os.path.join(args.out, "_DEPOSIT_README.md")
    with open(readme_p, "w") as fh:
        fh.write(DEPOSIT_README)
    pairs.append((readme_p, "README.md")); temp.append(readme_p)

    meta_p = os.path.join(args.out, "zenodo-metadata.json")   # kept OUTSIDE the tar too
    with open(meta_p, "w") as fh:
        json.dump(zenodo_metadata(args.stamp), fh, indent=2)
    pairs.append((meta_p, "zenodo-metadata.json"))

    if args.concern_matrix and os.path.exists(args.concern_matrix):
        pairs.append((args.concern_matrix, "REVIEWER-CONCERN-HANDLING-MATRIX.md"))
    else:
        print("  ! concern matrix not supplied/found — deposit will lack it",
              file=sys.stderr)

    pairs = sorted(set(pairs), key=lambda p: p[1])
    root = "psf-cautionary-REPRODUCTION-v%s-%s" % (VERSION, args.stamp)

    lines = [f"{root} — Zenodo reproduction deposit", "",
             f"{'sha256':<64}  size  path", "-" * 100]
    total = 0
    for src, dest in pairs:
        sz = os.path.getsize(src); total += sz
        lines.append(f"{sha256(src)}  {sz:>9}  {dest}")
    lines += ["", f"{len(pairs)} files, {total/1e6:.1f} MB uncompressed"]
    man_p = os.path.join(args.out, f"{root}.MANIFEST.txt")
    with open(man_p, "w") as fh:
        fh.write("\n".join(lines))

    out_path = os.path.join(args.out, f"{root}.tar.gz")
    with tarfile.open(out_path, "w:gz") as tar:
        for src, dest in pairs:
            tar.add(src, arcname=os.path.join(root, dest))
        tar.add(man_p, arcname=os.path.join(root, "MANIFEST.txt"))

    for p in temp:
        os.remove(p)
    print(f"Wrote {out_path}  ({os.path.getsize(out_path)/1e6:.1f} MB, {len(pairs)} files)")
    print(f"  manifest:        {man_p}")
    print(f"  zenodo metadata: {meta_p}  (upload the .tar.gz, paste these fields)")


if __name__ == "__main__":
    main()
