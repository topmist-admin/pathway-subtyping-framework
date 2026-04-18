# TCGA-ACC Data Pipeline Extension

**Date:** April 17, 2026
**Scope:** Extend the existing TCGA ingestion infrastructure to include TCGA-ACC (Adrenocortical Carcinoma) as a small-cohort rare-cancer indication for PSF methodology validation.

---

## Rationale

TCGA-ACC is a small-cohort (~92 samples) rare-cancer dataset with established C1A/C1B transcriptional subtype definitions (Assié et al. NEJM 2014; TCGA 2016). Its small N and existing subtype literature make it a useful target for demonstrating PSF's pathway-level subtyping methodology in a regime where conventional gene-level clustering is underpowered. Adding TCGA-ACC to the supported-projects whitelist is a single-line change and a natural extension of the TCGA coverage PSF already provides.

## Existing infrastructure

The repository has two independent paths for downloading TCGA data:

1. **`download_tcga_datasets.py` (API-direct path)** — queries the GDC API for each listed project and downloads TSV files directly. Does not require `gdc-client`. Uses a hardcoded `TCGA_PROJECTS` list at lines 34-47 as the whitelist of supported projects.
2. **`generate_gdc_manifest.py` + `gdc-client` (manifest path)** — generates a `gdc_manifest_TCGA-<project>.txt` file by querying the GDC API, then the `gdc-client` binary downloads by manifest. Pre-existing manifest files for other TCGA projects are present in the repo root.

Both paths target the same underlying GDC resources. The API-direct path is the lighter-weight option for cohorts under approximately 200 samples.

## Required changes

### A. Add TCGA-ACC to the API-direct whitelist

Edit `download_tcga_datasets.py` lines 34-47. Add TCGA-ACC at the top of the list:

```python
TCGA_PROJECTS = [
    'TCGA-ACC',   # Adrenocortical Carcinoma
    'TCGA-BRCA',  # Breast Cancer
    'TCGA-LUAD',  # Lung Adenocarcinoma
    'TCGA-OV',    # Ovarian Cancer
    'TCGA-GBM',   # Glioblastoma
    'TCGA-PAAD',  # Pancreatic Cancer
    'TCGA-UCEC',  # Uterine Corpus Endometrial Cancer
    'TCGA-HNSC',  # Head and Neck SCC
    'TCGA-COAD',  # Colon Adenocarcinoma
    'TCGA-STAD',  # Stomach Adenocarcinoma
    'TCGA-KIRC',  # Kidney Clear Cell Carcinoma
    'TCGA-KIRP',  # Kidney Papillary Carcinoma
    'TCGA-LIHC',  # Liver Hepatocellular Carcinoma
]
```

This is the entire change required for the API-direct path. The existing `fetch_file_ids()`, `download_project()`, and `download_file()` functions handle TCGA-ACC identically to other projects.

### B. Generate the manifest file (optional — only if `gdc-client` path is preferred)

Run the existing generator:

```bash
python3 generate_gdc_manifest.py TCGA-ACC
```

This produces `gdc_manifest_TCGA-ACC.txt` at the repo root.

**Known pre-existing issue with the manifest generator** (not introduced by this plan): `generate_gdc_manifest.py` lines 44 and 48 write newline characters between the `size` and `state` columns instead of tab characters. The produced manifest file has a broken column alignment that `gdc-client` will reject.

If the `gdc-client` path is required, fix the generator first:

```python
# Line 44 — change from \n to \t between 'size' and 'state'
f.write("id\tfilename\tmd5\tsize\tstate\n")

# Line 48 — same correction
f.write(f"{file_id}\t{filename}\t\t\treleased\n")
```

For a 92-sample cohort like TCGA-ACC, the API-direct path is adequate and this fix is not required.

### C. Convenience wrapper (optional)

A thin wrapper that invokes the API-direct download and performs a sanity check on the resulting file count. Save as `fetch_tcga_acc.py` at the repo root:

```python
#!/usr/bin/env python3
"""
Fetch TCGA-ACC expression data.

Thin wrapper around download_tcga_datasets.py that adds a sanity check on
expected sample count and a ready-state log.

Usage:
    python3 fetch_tcga_acc.py [--output tcga_data]

Expected outcome:
    ~92 .tsv files in <output>/TCGA-ACC/, each ~500 KB.
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
)
log = logging.getLogger("fetch_tcga_acc")

EXPECTED_SAMPLES_MIN = 80  # TCGA-ACC has ~92; allow headroom for GDC variability


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="tcga_data",
        help="Output directory root (default: tcga_data)",
    )
    args = parser.parse_args()

    output_root = Path(args.output)
    acc_dir = output_root / "TCGA-ACC"

    log.info("Invoking download_tcga_datasets.py for TCGA-ACC...")
    result = subprocess.run(
        [
            sys.executable,
            "download_tcga_datasets.py",
            "--project",
            "TCGA-ACC",
            "--output",
            str(output_root),
        ],
        check=False,
    )
    if result.returncode != 0:
        log.error("Download failed with exit code %d", result.returncode)
        return result.returncode

    if not acc_dir.exists():
        log.error("Output directory %s not created — download likely failed", acc_dir)
        return 1

    tsv_files = list(acc_dir.glob("*.tsv"))
    log.info("Downloaded %d .tsv files to %s", len(tsv_files), acc_dir)

    if len(tsv_files) < EXPECTED_SAMPLES_MIN:
        log.warning(
            "Sample count %d is below expected minimum %d — verify completeness "
            "before running analysis",
            len(tsv_files),
            EXPECTED_SAMPLES_MIN,
        )
        return 2

    log.info("TCGA-ACC ready. Next steps:")
    log.info("  1. Load via existing scripts/benchmark_loaders/loaders.py pattern")
    log.info("  2. Run PSF 12-pathway scoring on TCGA-ACC")
    log.info("  3. Validate against C1A/C1B classes (Assié 2014, TCGA 2016)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

Return codes:
- `0` — success, full cohort downloaded
- `1` — download failed outright
- `2` — download succeeded but sample count suspicious (manual check required)

## Execution sequence

```bash
cd /Users/rohitchauhan/Downloads/AI-Genetic-Research/pathway-subtyping-framework
source pathwayenv/bin/activate

# 1. Apply Patch A (one-line edit to download_tcga_datasets.py)

# 2. Run the download
python3 download_tcga_datasets.py --project TCGA-ACC --output tcga_data
# OR via the convenience wrapper (if Patch C applied):
python3 fetch_tcga_acc.py

# 3. Verify sample count
ls tcga_data/TCGA-ACC/*.tsv | wc -l   # expect ~92
```

Expected runtime: 5-15 minutes depending on GDC API response. Total download size approximately 50 MB.

## Integration with PSF analysis pipeline

No changes to PSF core are required. The existing `scripts/benchmark_loaders/loaders.py` pattern handles any TCGA project code once the data is on disk. TCGA-ACC loads identically to TCGA-KIRC, TCGA-COAD, or any other project in the whitelist.

After download, the standard analysis sequence is:

1. Load TCGA-ACC expression matrices via the benchmark loader
2. Run PSF 12-pathway scoring
3. Derive pathway-level subtypes via Gaussian mixture clustering with stability gates
4. Validate against published C1A/C1B labels and clinical outcomes

## Validation target

Minimum acceptance criteria for the downloaded cohort:

- ≥80 sample-level TSV files in `tcga_data/TCGA-ACC/`
- Each file parses as a valid gene-expression quantification TSV
- Clinical metadata (stage, survival, hormonal status, Weiss score) is retrievable from the accompanying clinical.tsv via GDC API (separate query)

## Scope boundaries

This plan covers data ingestion infrastructure only. It does not modify the PSF core framework, commit to a specific analytical outcome, or add new dependencies beyond those already used by the repository. Once the download completes, subsequent analysis is conducted using existing PSF primitives on the new cohort.

## Status tracking

- [ ] Patch A applied to `download_tcga_datasets.py`
- [ ] `python3 download_tcga_datasets.py --project TCGA-ACC --output tcga_data` executed
- [ ] Sample count ≥80 verified
- [ ] (Optional) Patch B manifest generator fix applied
- [ ] (Optional) Patch C `fetch_tcga_acc.py` convenience wrapper added
- [ ] PSF pathway scoring attempted on downloaded cohort (sanity-check run)
