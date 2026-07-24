#!/usr/bin/env python3
"""Build the two third-party-reviewer archives for the cautionary-framework paper.

Produces:

  psf-cautionary-CODE-<stamp>.tar.gz
      Everything needed to RUN the analyses: the framework source, the
      reproduction-bundle scripts, environment pins, and the top-level RUNME.
      No result files, no data.

  psf-cautionary-RESULTS-<stamp>.tar.gz
      Everything needed to CHECK the analyses without running them: deposited
      result JSON/TSV/CSV, the figures, the corrected benchmark input, and every
      package README stating what each number does and does not support.

Both archives get a MANIFEST.txt listing every file with its SHA-256, so a
reviewer can verify nothing shifted between the paper and the archive.

Controlled-access data is never included; none is used in this paper.

Usage:
    python build_review_archives.py --out DIR [--stamp YYYY-MM-DD]
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tarfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# (source path relative to repo root, path inside the archive)
CODE_ITEMS = [
    ("src/pathway_subtyping", "code/framework/src/pathway_subtyping"),
    ("pyproject.toml", "code/framework/pyproject.toml"),
    ("CITATION.cff", "code/framework/CITATION.cff"),
    ("CHANGELOG.md", "code/framework/CHANGELOG.md"),
    ("consolidation-cautionary/RUNME.md", "code/RUNME.md"),
    ("consolidation-cautionary/README.md", "code/README-stable-but-confounded.md"),
    ("consolidation-cautionary/requirements.txt", "code/requirements-v070-package.txt"),
    ("consolidation-cautionary/scripts", "code/scripts"),
    ("consolidation-cautionary/panels", "code/panels"),
    ("scripts/gate_ablation_study.py", "code/scripts_framework/gate_ablation_study.py"),
    ("scripts/plot_gate_ablation.py", "code/scripts_framework/plot_gate_ablation.py"),
]
# every cross-domain package contributes its scripts to CODE and its results to RESULTS
PACKAGES = [
    "gate_ablation", "gate_calibration", "benchmark_audit",
    "cancer_r38", "gtex_brain", "psychiatric_meta", "tcga_crc",
]

RESULTS_ITEMS = [
    ("consolidation-cautionary/RUNME.md", "results/RUNME.md"),
    ("CORRECTION_2026-07/corrected_benchmark_47datasets_v2.csv",
     "results/inputs/corrected_benchmark_47datasets_v2.csv"),
    ("CORRECTION_2026-07/ERRATUM_2026-07-08.md",
     "results/inputs/ERRATUM_2026-07-08.md"),
    ("CORRECTION_2026-07/_honest_numbers.json",
     "results/inputs/_honest_numbers.json"),
    ("consolidation-cautionary/data/partition",
     "results/inputs/partition_labels"),
]

SKIP_DIR = {"__pycache__", ".ipynb_checkpoints", ".git"}
SKIP_EXT = {".pyc", ".pyo"}


def _keep(path: str) -> bool:
    parts = set(path.split(os.sep))
    if parts & SKIP_DIR:
        return False
    return os.path.splitext(path)[1] not in SKIP_EXT


def collect(items):
    """Expand (src, dest) pairs into concrete file pairs."""
    out = []
    for src_rel, dest_rel in items:
        src = os.path.join(REPO, src_rel)
        if not os.path.exists(src):
            print(f"  ! missing, skipped: {src_rel}", file=sys.stderr)
            continue
        if os.path.isfile(src):
            if _keep(src):
                out.append((src, dest_rel))
            continue
        for root, dirs, files in os.walk(src):
            dirs[:] = [d for d in dirs if d not in SKIP_DIR]
            for f in files:
                full = os.path.join(root, f)
                if not _keep(full):
                    continue
                out.append((full, os.path.join(
                    dest_rel, os.path.relpath(full, src))))
    return out


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_archive(pairs, out_path: str, header: str) -> None:
    pairs = sorted(set(pairs), key=lambda p: p[1])
    lines = [header, "", f"{'sha256':<64}  size  path", "-" * 100]
    for src, dest in pairs:
        lines.append(f"{sha256(src)}  {os.path.getsize(src):>8}  {dest}")
    lines.append("")
    lines.append(f"{len(pairs)} files")
    manifest = "\n".join(lines)

    man_path = out_path + ".MANIFEST.txt"
    with open(man_path, "w") as fh:
        fh.write(manifest)

    root = os.path.basename(out_path).replace(".tar.gz", "")
    with tarfile.open(out_path, "w:gz") as tar:
        for src, dest in pairs:
            tar.add(src, arcname=os.path.join(root, dest))
        tar.add(man_path, arcname=os.path.join(root, "MANIFEST.txt"))
    print(f"  {out_path}  ({os.path.getsize(out_path)/1e6:.1f} MB, "
          f"{len(pairs)} files)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stamp", default="2026-07-23")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    code_items = list(CODE_ITEMS)
    results_items = list(RESULTS_ITEMS)
    for p in PACKAGES:
        base = f"consolidation-cautionary/cross-domain/{p}"
        code_items.append((f"{base}/scripts", f"code/packages/{p}/scripts"))
        code_items.append((f"{base}/README.md", f"code/packages/{p}/README.md"))
        results_items.append((f"{base}/results", f"results/packages/{p}/results"))
        results_items.append((f"{base}/README.md", f"results/packages/{p}/README.md"))
    # bundle-level cross-domain readme + gate-6 remap results
    code_items.append(("consolidation-cautionary/cross-domain/README.md",
                       "code/packages/README-cross-domain.md"))
    results_items.append(("consolidation-cautionary/cross-domain/results",
                          "results/packages/confound_remap/results"))
    results_items.append(("consolidation-cautionary/cross-domain/README.md",
                          "results/packages/README-cross-domain.md"))

    print("Building archives...")
    write_archive(
        collect(code_items),
        os.path.join(args.out, f"psf-cautionary-CODE-{args.stamp}.tar.gz"),
        "PSF cautionary-framework paper — CODE archive\n"
        "Framework source + reproduction scripts + environment pins.\n"
        "Run order and per-result mapping: RUNME.md.\n"
        "NOTE: requires the v0.8 line of the framework "
        "(pathway_subtyping.discreteness, clustering_dl); the source here IS "
        "that version. Do not substitute an earlier PyPI release.")
    write_archive(
        collect(results_items),
        os.path.join(args.out, f"psf-cautionary-RESULTS-{args.stamp}.tar.gz"),
        "PSF cautionary-framework paper — RESULTS archive\n"
        "Deposited outputs, figures, and the corrected benchmark input, with "
        "each package's README stating what its numbers do and do not support.\n"
        "Every figure in the manuscript traces to a file here; see RUNME.md.\n"
        "All inputs are public. No controlled-access data is included or used.")


if __name__ == "__main__":
    main()
