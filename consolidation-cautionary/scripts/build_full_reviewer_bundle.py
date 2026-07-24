#!/usr/bin/env python3
"""Build ONE self-contained archive for third-party stress testing of the
cautionary-framework manuscript, BEFORE v0.8.0 is publicly released.

Unlike the split CODE / RESULTS archives (build_review_archives.py), this bundles
everything a reviewer needs to run and audit the claims from scratch, with no PyPI
dependency on the unreleased gate:

  1. The full v0.8 framework SOURCE (src/ + pyproject + tests) — so the gate is
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
       [--manuscript-dir PATH]   (defaults to the AIForAutismOutReach SR folder)
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tarfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_MS = ("/Users/rohitchauhan/Downloads/AIForAutismOutReach/manuscripts/"
              "SCIENTIFIC-REPORTS-SUBMISSION-2026")

SKIP_DIR = {"__pycache__", ".ipynb_checkpoints", ".git", ".pytest_cache",
            ".mypy_cache", "node_modules", ".venv", "pathwayenv"}
SKIP_EXT = {".pyc", ".pyo", ".so", ".o"}
SKIP_NAME = {".DS_Store"}

# The bundle MIRRORS THE REPO ROOT LAYOUT so every reproduction script's
# repo-relative paths (e.g. ../../../../CORRECTION_2026-07/...) resolve unchanged.
# `pip install -e .` runs at the bundle root (pyproject.toml is there).

# --- framework source (the unreleased v0.8 line) — at bundle root ---
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
    files = [
        "REBUILD-DRAFT-v2-2026-07-24-FULL.md",
        "REFERENCES-VERIFIED-2026-07-23.md",
        "REVIEWER-CONCERN-HANDLING-MATRIX.md",
        "REVIEW-BLOCKERS-2026-07-23.md",
        "correspondence/2026-07-21_reviewer_comments_MAJOR_REVISION_post-withdrawal.md",
        "correspondence/2026-07-22_REBUILD_architecture_cautionary_framework.md",
        "rebuild-deliverables-2026-07-24/manuscript-v2.1-2026-07-24.pdf",
    ]
    out = []
    for f in files:
        src = os.path.join(ms_dir, f)
        if os.path.exists(src):
            out.append((src, os.path.join("manuscript", os.path.basename(f))))
        else:
            print(f"  ! manuscript file missing: {f}", file=sys.stderr)
    return out


def _keep(path):
    parts = set(path.split(os.sep))
    if parts & SKIP_DIR:
        return False
    if os.path.basename(path) in SKIP_NAME:
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
            dirs[:] = [d for d in dirs if d not in SKIP_DIR]
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


STRESS_GUIDE = """# STRESS-TEST THIS BUNDLE — third-party reviewer guide

Self-contained bundle to independently reproduce and stress-test every claim in the
cautionary-framework manuscript, BEFORE the framework's v0.8 line is published to
PyPI. Nothing here needs controlled-access data or the unreleased release.

## Layout (mirrors the source repo so every script's relative paths resolve)
- `src/`, `tests/`, `pyproject.toml`  the v0.8 framework source (the gate lives here; not yet on PyPI)
- `consolidation-cautionary/`         reproduction packages: scripts + deposited results + per-package READMEs (`RUNME.md` indexes them)
- `research-results/`, `CORRECTION_2026-07/`  redistributable cached public data (pathway scores, partition labels, corrected benchmark)
- `scripts/`                          repo-root framework tooling (gate ablation study + plotter)
- `manuscript/`                       the draft, verified references, the concern-handling matrix, the reviewer comments, the internal hostile-review log, and the built PDF
- `MANIFEST.txt`                      every file with its SHA-256

## Setup (5 minutes) — run at the bundle root
```
python -m venv .venv && . .venv/bin/activate
pip install -e .                    # installs the v0.8 line from source (pyproject.toml is here)
python -c "import pathway_subtyping as p; print(p.__version__)"   # 0.8.0
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

## Where to start auditing
1. `manuscript/REVIEWER-CONCERN-HANDLING-MATRIX.md` — every SR reviewer point and every
   internal hostile-review finding, mapped to the artifact that proves the resolution.
2. `manuscript/REVIEW-BLOCKERS-2026-07-23.md` — the internal adversarial findings in full.
3. `consolidation-cautionary/RUNME.md` — headline numbers traced to deposited files.

## Known limits (disclosed)
v0.8 is unreleased (that is why this bundle ships the source); the Zenodo benchmark
record needs superseding; large-N psychiatric validation is access-gated; the
real-data control is two anchors; funding/contributions statements are pending. See
`manuscript/REVIEWER-CONCERN-HANDLING-MATRIX.md` Part C.

## Determinism
Everything is seed 42. No-network scripts reproduce exactly. Network scripts draw from
live public endpoints whose contents can drift; the deposited JSONs are the reference.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stamp", default="2026-07-24")
    ap.add_argument("--manuscript-dir", default=DEFAULT_MS)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pairs = []
    pairs += collect(SOURCE_ITEMS)
    pairs += collect(TOOLING_ITEMS)
    pairs += collect(BUNDLE_ITEMS)
    pairs += collect(DATA_ITEMS)
    pairs += manuscript_items(args.manuscript_dir)
    pairs = sorted(set(pairs), key=lambda p: p[1])

    root = f"psf-cautionary-REVIEWER-BUNDLE-{args.stamp}"
    guide_path = os.path.join(args.out, "_STRESS_TEST.md")
    with open(guide_path, "w") as fh:
        fh.write(STRESS_GUIDE)
    pairs.append((guide_path, "STRESS-TEST.md"))

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

    print(f"Wrote {out_path}")
    print(f"  {len(pairs)} files, {os.path.getsize(out_path)/1e6:.1f} MB compressed "
          f"({total/1e6:.1f} MB uncompressed)")
    print(f"  manifest: {man}")


if __name__ == "__main__":
    main()
