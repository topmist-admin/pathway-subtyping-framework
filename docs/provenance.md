# Dataset Provenance

Tracked public datasets consumed by validation scripts and notebooks, with
exact accessions and file hashes where feasible. Aimed at reproducibility
of real-data acceptance results.

---

## F1 — Conformal coverage (Phase 1)

- **TCGA-COAD** (STAR gene counts, GDC). Consumed by
  `scripts/validate_f1_real_data.py`. 57 matched STAR TSVs under
  `tcga_data/TCGA-COAD/` — access via the GDC portal at
  https://portal.gdc.cancer.gov/projects/TCGA-COAD (bulk RNA-seq,
  STAR gene counts).
- **GSE28521** — autism post-mortem frontal cortex microarray. 79
  samples, pre-computed pathway scores consumed from
  `research-results/GSE28521/pathway_scores.csv`. Canonical accession:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28521.

## F2 — Cross-platform harmonization (Phase 1)

- **GSE28521** (Platform A) — Affymetrix U133 Plus 2.0 microarray,
  post-mortem frontal cortex, n=79. File:
  `research-results/GSE28521/gene_expression_processed.csv`
  (pre-normalised log2 intensity, 9,552 genes). Accession:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28521.
- **GSE80655** (Platform B) — Illumina HiSeq 2000 bulk RNA-seq,
  DLPFC, n=281. File:
  `research-results/GSE80655/gene_expression_processed.csv`
  (log2 abundance, 40,461 genes). Accession:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655.

## F5 — In-silico perturbation directional test (Phase 2)

- **TCGA-COAD** (as above, 57 samples) — serves as the WT cohort for
  `scripts/validate_f5_real_data.py`. The directional edges are
  literature-anchored (MYC / TP53 / E2F1 / CCNE1 / CDK1 → hallmark
  pathways where each gene is a direct member and a known driver).

## F10 — Multi-omics fusion (Phase 3)

- **10x Genomics `pbmc_1k_protein_v3`** — CITE-seq PBMC reference,
  713 cells × 33,538 genes × 17 antibodies. Public 10x tutorial
  dataset. File:
  `data/f10_citeseq/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5`
  (~5.2 MB h5).
  - URL: https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5
  - SHA-256: `b7f8d7c26b132c972edfc8624ee8387aa9d8a771258746a7285968df4192f094`
  - License: 10x Genomics sample data — free for research and
    educational use.
- **Antibody → pathway mapping** at `data/omics/cite_adt_to_pathway.yaml`
  — hand-curated from CellMarker 2.0 for the 17-antibody panel; seven
  immune-function pathways matched across RNA (gene symbols) and
  protein (TotalSeqB barcode labels).

---

## Notes

- Fetchable URL-backed datasets are not redistributed in this repo; fetch
  scripts (or copy-pasteable `curl` commands in each validation script's
  error path) reconstruct them on demand.
- Datasets pre-packaged under `research-results/` were produced by prior
  research runs using the same upstream GEO accessions; the processed CSVs
  are checked in for reproducibility of the validation scripts that
  consume them.
