"""
Fetch GSE123753 — Rett syndrome iPSC-derived cells with MECP2 deletion.

Cohort: Boxer et al. 2020 (pubmed 30923256) — WT vs isogenic MECP2-deletion
iPSCs, NPCs, and Neurons profiled by TRAP-seq (ribosome engagement) and
Input RNA-seq. We use the Input (bulk RNA-seq) NPC and Neuron samples to
build a WT vs MECP2-KO pair for F5 real-data acceptance.

Downloads:
    data/f5_gse123753/GSE123753_counts_gene_replicates.csv.gz  (raw counts)
    data/f5_gse123753/series_matrix.txt.gz                       (metadata)

Writes a parsed pandas-friendly samples-by-genes CSV:
    data/f5_gse123753/input_rnaseq_samples_x_genes.csv
with companion sample metadata at:
    data/f5_gse123753/input_rnaseq_sample_metadata.csv
"""

from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger("fetch_gse123753")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "f5_gse123753"

COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123753/suppl/"
    "GSE123753_counts_gene_replicates.csv.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123753/matrix/"
    "GSE123753_series_matrix.txt.gz"
)


def _curl(url: str, dest: Path) -> None:
    if dest.exists():
        logger.info("[cache] %s", dest.name)
        return
    logger.info("[fetch] %s -> %s", url, dest)
    subprocess.check_call(["curl", "-s", "-L", "-o", str(dest), url])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    _curl(COUNTS_URL, DATA_DIR / "GSE123753_counts_gene_replicates.csv.gz")
    _curl(MATRIX_URL, DATA_DIR / "series_matrix.txt.gz")

    logger.info("[parse] counts matrix")
    df = pd.read_csv(
        DATA_DIR / "GSE123753_counts_gene_replicates.csv.gz", index_col=0,
    )
    # "gene" column holds the symbol; the other columns are either
    # NumReads.* (raw counts) or TPM.* (normalised).
    gene_col = df["gene"]
    count_cols = [c for c in df.columns if c.startswith("NumReads.")]
    counts = df[count_cols].copy()
    counts.index = gene_col.values
    counts.index.name = "gene"
    counts = counts.dropna()
    # Drop any duplicate symbols by max-aggregation (keeps the best-covered tx)
    counts = counts.groupby(level=0).max()

    # Keep only Input.* (RNA-seq, not TRAP) NPC + Neu samples — these are
    # the bulk RNA-seq replicates where WT and RTT are directly comparable.
    input_cols = [c for c in counts.columns if c.startswith("NumReads.Input_")]
    keep = [c for c in input_cols if "_NPC_" in c or "_Neu_" in c]
    filtered = counts[keep].copy()

    def _describe(col: str) -> dict:
        # e.g. "NumReads.Input_Neu_WT_2"
        tag = col.replace("NumReads.Input_", "")
        parts = tag.split("_")
        cell_type, genotype, rep = parts[0], parts[1], parts[2]
        return {
            "sample_id": tag,
            "cell_type": cell_type,
            "genotype": "WT" if genotype == "WT" else "MECP2_KO",
            "replicate": int(rep),
        }

    meta = pd.DataFrame([_describe(c) for c in keep]).set_index("sample_id")

    # Rename columns to the clean sample_id tag for downstream use.
    filtered.columns = [c.replace("NumReads.Input_", "") for c in filtered.columns]

    # Samples × genes transpose (consistent with other cohort layouts)
    samples_x_genes = filtered.T.astype(float)
    samples_x_genes.index.name = "sample"

    out_counts = DATA_DIR / "input_rnaseq_samples_x_genes.csv"
    out_meta = DATA_DIR / "input_rnaseq_sample_metadata.csv"
    samples_x_genes.to_csv(out_counts)
    meta.to_csv(out_meta)

    logger.info("[wrote] %s shape=%s", out_counts, samples_x_genes.shape)
    logger.info("[wrote] %s", out_meta)
    logger.info("\n%s", meta.value_counts(["cell_type", "genotype"]).to_string())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
