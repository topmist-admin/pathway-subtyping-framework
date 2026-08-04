#!/usr/bin/env python3
"""Build ONE self-contained archive for third-party stress testing of the
cautionary-framework manuscript, BEFORE v0.9.0 is publicly released.

Unlike the split CODE / RESULTS archives (build_review_archives.py), this bundles
everything a reviewer needs to run and audit the claims from scratch, with no PyPI
dependency on the unreleased gate:

  1. The full v0.9 framework SOURCE (src/ + pyproject + tests) — so the gate is
     importable via `pip install -e .` without the PyPI release.
  2. All reproduction packages (scripts + deposited results + READMEs).
  3. Redistributable CACHED public data (GSE28521 / GSE80655 pathway scores,
     partition labels, pathway panels, the corrected benchmark CSV) so the
     no-network analyses run immediately; cBioPortal/GEO fetches are documented.
  4. The manuscript, the verified reference list, the concern-handling matrix,
     the internal hostile-review logs, and the erratum.
  5. A top-level STRESS-TEST guide + a MANIFEST with per-file SHA-256.

Controlled-access data is never included; none is used.

Usage: python build_full_reviewer_bundle.py --out DIR [--stamp YYYY-MM-DD]
       --manuscript-dir PATH
         (or set the PSF_MANUSCRIPT_DIR environment variable; the manuscript lives
          outside this repository, so there is no baked-in default path)
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tarfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# The manuscript is maintained outside this public repo; supply its location via
# --manuscript-dir or the PSF_MANUSCRIPT_DIR env var. No private path is hard-coded.
DEFAULT_MS = os.environ.get("PSF_MANUSCRIPT_DIR", "")

# *.egg-info is build metadata, and `pip install -e .` -- step 1 of the reviewer
# guide -- REWRITES it in place. Shipping it meant every reviewer who followed the
# instructions then saw two MANIFEST.txt hash mismatches and had to work out that
# the bundle was fine and the instructions had corrupted it. Excluded.
SKIP_DIR = {"__pycache__", ".ipynb_checkpoints", ".git", ".pytest_cache",
            ".mypy_cache", "node_modules", ".venv", "pathwayenv"}
SKIP_DIR_SUFFIX = (".egg-info",)
SKIP_EXT = {".pyc", ".pyo", ".so", ".o"}
SKIP_NAME = {".DS_Store"}

# Self-audit documents, withheld from the reviewer bundle on purpose. These state what WE
# concluded was weak and how we fixed it. Handing them over before an independent review
# anchors the reviewer to our findings and turns agreement into an echo rather than
# corroboration. Instructions (RUNME.md, DATA-ACQUISITION.md, CLEAN-ROOM-VALIDATION.md)
# still ship -- the line is "how to check" versus "what we already concluded".
SKIP_SELF_AUDIT = {"FULL-RERUN-REVIEWER-PLAN.md", "VERIFICATION-MATRIX.md"}

# The bundle MIRRORS THE REPO ROOT LAYOUT so every reproduction script's
# repo-relative paths (e.g. ../../../../CORRECTION_2026-07/...) resolve unchanged.
# `pip install -e .` runs at the bundle root (pyproject.toml is there).

# --- framework source (the unreleased v0.9 line) — at bundle root ---
SOURCE_ITEMS = [
    ("src", "src"),
    ("tests", "tests"),
    ("pyproject.toml", "pyproject.toml"),
    ("CITATION.cff", "CITATION.cff"),
    ("CHANGELOG.md", "CHANGELOG.md"),
    ("README.md", "README-framework.md"),
]

# --- reproduction bundle (scripts + results + READMEs), keep its real name ---
BUNDLE_ITEMS = [
    ("consolidation-cautionary", "consolidation-cautionary"),
]

# --- redistributable cached public data — at the repo-relative paths scripts expect ---
DATA_ITEMS = [
    ("research-results/GSE28521", "research-results/GSE28521"),
    ("research-results/GSE80655", "research-results/GSE80655"),
    ("CORRECTION_2026-07", "CORRECTION_2026-07"),
]

# repo-root framework tooling the bundle references — at repo-relative paths
TOOLING_ITEMS = [
    ("scripts/gate_ablation_study.py", "scripts/gate_ablation_study.py"),
    ("scripts/plot_gate_ablation.py", "scripts/plot_gate_ablation.py"),
]


def manuscript_items(ms_dir):
    # (source relative to ms_dir, destination inside the bundle) — explicit so the
    # PROPOSED new version and the WITHDRAWN prior version land in clearly-labelled
    # subfolders and can never be confused for one another.
    files = [
        # --- shared audit material at manuscript/ top level ---
        ("REVIEWER-CONCERN-HANDLING-MATRIX.md",
         "manuscript/REVIEWER-CONCERN-HANDLING-MATRIX.md"),
        ("REVIEW-BLOCKERS-2026-07-23.md",
         "manuscript/internal-hostile-review-log.md"),
        ("REFERENCES-VERIFIED-2026-07-23.md",
         "manuscript/REFERENCES-VERIFIED.md"),
        ("correspondence/2026-07-22_REBUILD_architecture_cautionary_framework.md",
         "manuscript/rebuild-architecture-rationale.md"),
        # --- PROPOSED new version (not yet submitted) ---
        ("REBUILD-DRAFT-v2-2026-07-24-FULL.md",
         "manuscript/proposed-new-version/PROPOSED-manuscript-v2.2-2026-07-24.md"),
        # Point at the CURRENT typeset PDF. The old path
        # (rebuild-deliverables-2026-07-24/manuscript-v2.2-2026-07-24.pdf) no longer
        # exists, so the builder silently skipped it and every bundle shipped the
        # markdown with no PDF beside it. This file is built by
        # build_submission_pdf.py, which strips the internal draft banner.
        ("SUBMISSION-PACKAGE-2026-07-25/01_manuscript/Manuscript-Chauhan-Paulus-SciReports-2026-08-02-CURRENT.pdf",
         "manuscript/proposed-new-version/PROPOSED-manuscript-v2.2-2026-07-24.pdf"),
        # --- WITHDRAWN prior version (do not cite) ---
        ("Chauhan-Mandatory-Validation-Gates-Manuscript-RevisionR1-2026-05-25.docx",
         "manuscript/withdrawn-prior-version/WITHDRAWN-manuscript-as-reviewed-2026-05-25.docx"),
        ("correspondence/Comments.docx",
         "manuscript/withdrawn-prior-version/Scientific-Reports-reviewer-comments-ORIGINAL.docx"),
        ("correspondence/2026-07-21_reviewer_comments_MAJOR_REVISION_post-withdrawal.md",
         "manuscript/withdrawn-prior-version/Scientific-Reports-reviewer-comments-TRANSCRIBED.md"),
    ]
    out = []
    for src_rel, dest in files:
        src = os.path.join(ms_dir, src_rel)
        if os.path.exists(src):
            out.append((src, dest))
        else:
            print(f"  ! manuscript file missing: {src_rel}", file=sys.stderr)
    # the erratum ships inside the withdrawn folder too (it documents the retractions)
    erratum = os.path.join(REPO, "CORRECTION_2026-07/ERRATUM_2026-07-08.md")
    if os.path.exists(erratum):
        out.append((erratum,
                    "manuscript/withdrawn-prior-version/ERRATUM_2026-07-08.md"))
    return out


def _keep(path):
    parts = path.split(os.sep)
    if set(parts) & SKIP_DIR:
        return False
    if any(p.endswith(SKIP_DIR_SUFFIX) for p in parts):
        return False
    if os.path.basename(path) in SKIP_NAME:
        return False
    if os.path.basename(path) in SKIP_SELF_AUDIT:
        return False
    return os.path.splitext(path)[1] not in SKIP_EXT


def collect(items, root=REPO):
    out = []
    for src_rel, dest_rel in items:
        src = src_rel if os.path.isabs(src_rel) else os.path.join(root, src_rel)
        if not os.path.exists(src):
            print(f"  ! missing, skipped: {src_rel}", file=sys.stderr)
            continue
        if os.path.isfile(src):
            if _keep(src):
                out.append((src, dest_rel))
            continue
        for r, dirs, files in os.walk(src):
            dirs[:] = [d for d in dirs
                       if d not in SKIP_DIR and not d.endswith(SKIP_DIR_SUFFIX)]
            for f in files:
                full = os.path.join(r, f)
                if _keep(full):
                    out.append((full, os.path.join(dest_rel, os.path.relpath(full, src))))
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


PROPOSED_README = """# PROPOSED new version — NOT YET SUBMITTED

This folder holds the **proposed** rebuild of the manuscript. It has NOT been
submitted to any journal. It supersedes the withdrawn prior version
(`../withdrawn-prior-version/`).

- `PROPOSED-manuscript-v2.2-2026-07-24.md` — the draft (markdown source of record)
- `PROPOSED-manuscript-v2.2-2026-07-24.pdf` — typeset PDF of the same

Status: draft v2.2, after two internal hostile-review passes and one independent
third-party stress-test (MAJOR REVISION; its analytical findings applied). Known
pre-submission gaps (v0.9 release + code DOI, Zenodo supersession, funding/author
statements) are listed in the manuscript's Data/Code section and in Part C of
`../REVIEWER-CONCERN-HANDLING-MATRIX.md`. Every number in this draft is
reproducible from this bundle — see `../../STRESS-TEST.md`.
"""

WITHDRAWN_README = """# WITHDRAWN prior version — DO NOT CITE

⚠️ **The manuscript in this folder was withdrawn before publication and its
headline figures are RETRACTED. Do not cite any number from it.** It is included
for one purpose only: so a reviewer can verify, by direct comparison, that the
retracted numbers do **not** leak into the proposed new version — a de-leakage
check. It is not a current claim of any kind.

Contents:
- `WITHDRAWN-manuscript-as-reviewed-2026-05-25.docx` — the withdrawn manuscript,
  the exact version the Scientific Reports reviewers saw.
- `Scientific-Reports-reviewer-comments-ORIGINAL.docx` — the reviewer comments as
  received from the journal.
- `Scientific-Reports-reviewer-comments-TRANSCRIBED.md` — the same comments
  transcribed to text with point IDs (R1.1…R3.11) for cross-reference.
- `ERRATUM_2026-07-08.md` — the public erratum documenting exactly which figures
  are retracted and why the benchmark did not reproduce.

**Retracted figures that must appear NOWHERE in the proposed version** (use this
list to run the de-leakage check): CMS4 recovery 75.9% (OR 16.71, p 1.4×10⁻²⁵),
cross-platform ARI 0.870, adaptive-threshold R²=0.889 (RMSE 0.051). NOTE: the
GSE80655 psychiatric-partition bootstrap stability (≈0.92) is a *valid* figure and
is legitimately used in the proposed version's Result 4 — it is NOT retracted; see
the erratum and the concern-handling matrix.

The withdrawal was driven by the 47-dataset benchmark failing to reproduce and by
an unsupported adaptive-threshold model, not by misconduct; the authors self-issued
the erratum after an external reproducibility review.
"""

STRESS_GUIDE = """# STRESS-TEST THIS BUNDLE — third-party reviewer guide

Self-contained bundle to independently reproduce and stress-test every claim in the
cautionary-framework manuscript. Nothing here needs controlled-access data.

## Layout (mirrors the source repo so every script's relative paths resolve)
- `src/`, `tests/`, `pyproject.toml`  the v0.9.0 framework source (also on PyPI; bundled for permanence)
- `consolidation-cautionary/`         reproduction packages: scripts + deposited results + per-package READMEs (`RUNME.md` indexes them)
- `research-results/`, `CORRECTION_2026-07/`  redistributable cached public data (pathway scores, partition labels, corrected benchmark)
- `scripts/`                          repo-root framework tooling (gate ablation study + plotter)
- `manuscript/`                       audit material: concern-handling matrix, verified references, internal hostile-review log, rebuild rationale
- `manuscript/proposed-new-version/`  the PROPOSED rebuild (md + PDF) — NOT yet submitted
- `manuscript/withdrawn-prior-version/`  the WITHDRAWN prior manuscript + the journal's reviewer comments (original .docx + transcription) + the erratum — for the de-leakage check; DO NOT CITE
- `MANIFEST.txt`                      every file with its SHA-256

## Setup (5 minutes) — run at the bundle root
v0.9.0 IS on PyPI (`pip install pathway-subtyping==0.9.0`), but this bundle also
CONTAINS the v0.9.0 source under `src/`, and `pip install -e .` builds it from the files
in this folder (the `.` = "the project in this directory"; pip reads the bundled
`pyproject.toml`). Installing from the bundle is preferred for a clean-room check: it
proves the archive is self-sufficient rather than trusting the index.
```
python -m venv .venv && . .venv/bin/activate
pip install -e .                    # builds v0.9 from the local src/ HERE — not PyPI, not git
python -c "import pathway_subtyping as p; print(p.__version__)"   # 0.9.0  (if it says 0.8.0 you installed the PyPI release by mistake)
pip install numpy pandas scikit-learn scipy statsmodels requests torch   # analysis deps (torch only for DL baselines)
```

## Fastest checks (no network, seconds to minutes) — from the bundle root
- **Benchmark audit + column-validity diagnostic** (Result 2):
  `python consolidation-cautionary/cross-domain/benchmark_audit/scripts/audit_corrected_benchmark.py`
  → confirms the retracted R²=0.889 does not refit, and that the benchmark's stability
  column is unsound (5 high-silhouette cohorts with ≤0 stability).
- **Honest ablation re-analysis** (Result 1):
  `python consolidation-cautionary/cross-domain/gate_ablation/scripts/reanalyze_ablation_honest.py`
  → three-way accounting, Wilson CIs, McNemar, and the SigClust-p-alone head-to-head,
  all from the deposited raw CSV.
- **Donor-level flagship** (Result 4):
  `python consolidation-cautionary/cross-domain/flagship_stats/scripts/flagship_donor_level.py`
  → region V 0.66, diagnosis permutation p 0.234 at n=48 donors.
- **Flagship stability** (Result 4):
  `python consolidation-cautionary/cross-domain/flagship_stats/scripts/flagship_stability.py`
  → recomputes 0.921 from cached scores (confirms the 0.923 provenance).

## Checks that fetch public data (cBioPortal/GEO, no auth)
- **Within-study calibration** (Result 1):
  `python consolidation-cautionary/cross-domain/gate_calibration/scripts/calibrate_within_study.py`
  → IDH-glioma certified (87% wildtype-arm recall), immune continuum rejected.
- **Cancer worked example** (Result 3):
  `python consolidation-cautionary/cross-domain/cancer_r38/scripts/fetch_and_validate_brca.py`
  → ARI 0.218 vs enrichment 87.6%; DL baselines.
- **Separation sweep** (Result 1, slow — nested bootstrap; single-thread BLAS recommended):
  `OMP_NUM_THREADS=1 python consolidation-cautionary/cross-domain/gate_ablation/scripts/separation_sweep.py`

## De-leakage check (verify the reframe is real, not cosmetic)
The withdrawn prior manuscript is included so you can confirm the retracted numbers do
NOT survive into the proposed version:
- `manuscript/withdrawn-prior-version/` — the withdrawn manuscript (as reviewed), the
  Scientific Reports reviewer comments (original .docx + a transcription with point
  IDs), and the erratum. **Read its `README` first — DO NOT CITE anything in it.**
- `manuscript/proposed-new-version/` — the PROPOSED rebuild (not yet submitted).
- The retracted figures that must appear nowhere in the proposed version: **CMS4 75.9%
  (OR 16.71), cross-platform ARI 0.870, adaptive-threshold R²=0.889.** Grep the proposed
  draft for them — they appear only in its Prior-Work Disclosure, as retractions. (One
  subtlety flagged in the withdrawn folder's README: the GSE80655 stability ≈0.92 is a
  *valid* figure, legitimately used in Result 4, and is NOT among the retractions.)

## Where to start auditing
1. `manuscript/REVIEWER-CONCERN-HANDLING-MATRIX.md` — every SR reviewer point (R1.1–R3.11)
   and every internal hostile-review finding, mapped to the artifact that proves the
   resolution. Its Part C lists what is still unresolved.
2. `manuscript/internal-hostile-review-log.md` — the internal adversarial findings in full.
3. `consolidation-cautionary/RUNME.md` — headline numbers traced to deposited files.

## Known limits (disclosed)
Large-N psychiatric validation is access-gated and remains future work; the real-data
control is two anchors, not a per-domain calibration. (Resolved since an earlier draft of
this guide: v0.9.0 is released on PyPI with a Zenodo concept DOI, the Zenodo benchmark
record has been superseded by v2.1, and the funding and author-contribution statements
are supplied.) See `manuscript/REVIEWER-CONCERN-HANDLING-MATRIX.md` Part C.

## Determinism
Everything is seed 42. No-network scripts reproduce exactly. Network scripts draw from
live public endpoints whose contents can drift; the deposited JSONs are the reference.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stamp", default="2026-07-24")
    ap.add_argument("--manuscript-dir", default=DEFAULT_MS,
                    help="Directory holding the manuscript files (maintained outside "
                         "this repo). Or set PSF_MANUSCRIPT_DIR.")
    args = ap.parse_args()
    if not args.manuscript_dir:
        sys.exit("error: no manuscript directory given. Pass --manuscript-dir PATH "
                 "or set the PSF_MANUSCRIPT_DIR environment variable "
                 "(the manuscript is maintained outside this repository).")
    os.makedirs(args.out, exist_ok=True)

    pairs = []
    pairs += collect(SOURCE_ITEMS)
    pairs += collect(TOOLING_ITEMS)
    pairs += collect(BUNDLE_ITEMS)
    pairs += collect(DATA_ITEMS)
    pairs += manuscript_items(args.manuscript_dir)

    # The two reviewer-facing guides are hand-maintained next to the tarball. Ship
    # them INSIDE it as well: CLEAN-ROOM-VALIDATION.md tells the reviewer the
    # bundle is self-contained, which was not true of the entry-point documents
    # themselves -- anyone handed only the .tar.gz got the archive without the
    # instructions for it.
    # VERIFICATION-MATRIX.md is deliberately NOT bundled. It is our own self-audit: it
    # names every gap we found and how we resolved it. Shipping it means the reviewer
    # reads our answers before forming their own, which anchors them and weakens exactly
    # the independence the exercise is for. Withhold it, let them audit cold, then
    # compare -- agreement is then corroboration rather than an echo.
    #
    # CLEAN-ROOM-VALIDATION.md IS bundled: it is instructions (how to install, what to
    # run, what the acceptance bar is per result), not conclusions about correctness.
    for guide in ("CLEAN-ROOM-VALIDATION.md",):
        g = os.path.join(args.out, guide)
        if os.path.isfile(g):
            pairs.append((g, guide))
        else:
            print(f"  ! guide not found, not bundled: {guide}", file=sys.stderr)

    pairs = sorted(set(pairs), key=lambda p: p[1])

    root = f"psf-cautionary-REVIEWER-BUNDLE-{args.stamp}"
    # generated labelling files (written to temp, added to tar, cleaned up after)
    temp_files = []
    for content, dest, tmpname in [
        (STRESS_GUIDE, "STRESS-TEST.md", "_STRESS_TEST.md"),
        (PROPOSED_README, "manuscript/proposed-new-version/README.md", "_PROPOSED.md"),
        (WITHDRAWN_README, "manuscript/withdrawn-prior-version/README.md", "_WITHDRAWN.md"),
    ]:
        p = os.path.join(args.out, tmpname)
        with open(p, "w") as fh:
            fh.write(content)
        pairs.append((p, dest))
        temp_files.append(p)

    # manifest
    lines = [f"{root} — full third-party reviewer bundle ({args.stamp})", "",
             f"{'sha256':<64}  size  path", "-" * 100]
    total = 0
    for src, dest in pairs:
        sz = os.path.getsize(src)
        total += sz
        lines.append(f"{sha256(src)}  {sz:>9}  {dest}")
    lines += ["", f"{len(pairs)} files, {total/1e6:.1f} MB uncompressed"]
    man = os.path.join(args.out, f"{root}.MANIFEST.txt")
    with open(man, "w") as fh:
        fh.write("\n".join(lines))

    out_path = os.path.join(args.out, f"{root}.tar.gz")
    with tarfile.open(out_path, "w:gz") as tar:
        for src, dest in pairs:
            tar.add(src, arcname=os.path.join(root, dest))
        tar.add(man, arcname=os.path.join(root, "MANIFEST.txt"))

    for p in temp_files:    # they now live inside the tar
        os.remove(p)
    print(f"Wrote {out_path}")
    print(f"  {len(pairs)} files, {os.path.getsize(out_path)/1e6:.1f} MB compressed "
          f"({total/1e6:.1f} MB uncompressed)")
    print(f"  manifest: {man}")


if __name__ == "__main__":
    main()
