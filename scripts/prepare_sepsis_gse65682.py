#!/usr/bin/env python3
"""
Prepare GSE65682 (MARS Sepsis Cohort) for PSF analysis.

Downloads expression data + clinical metadata from GEO,
maps probe IDs to gene symbols, and creates input files
for the pathway-subtyping pipeline.
"""

import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "sepsis"
)


def download_geo_data():
    """Download GSE65682 from GEO using GEOparse."""
    import GEOparse

    print("=" * 60)
    print("Phase 1: Downloading GSE65682 from GEO")
    print("=" * 60)

    os.environ["GEOPARSE_USE_HTTP_FOR_FTP"] = "yes"
    gse = GEOparse.get_GEO(geo="GSE65682", destdir=DATA_DIR, silent=True)
    print(f"Downloaded GSE65682: {len(gse.gsms)} samples")
    return gse


def extract_expression_matrix(gse):
    """Extract expression matrix from GEO samples, map probes to gene symbols."""
    print("\nExtracting expression matrix...")

    # Get platform annotation for probe-to-gene mapping
    platform_id = list(gse.gpls.keys())[0]
    gpl = gse.gpls[platform_id]
    print(f"Platform: {platform_id}")

    # Build probe-to-gene mapping from platform table
    gpl_table = gpl.table
    print(f"GPL table columns: {list(gpl_table.columns)[:10]}")

    # HG-U219 uses 'Gene Symbol' column
    gene_col = None
    for col in ["Gene Symbol", "gene_assignment", "Symbol", "GENE_SYMBOL"]:
        if col in gpl_table.columns:
            gene_col = col
            break

    if gene_col is None:
        print("WARNING: No gene symbol column found. Available columns:")
        print(gpl_table.columns.tolist())
        # Try to find any column containing gene info
        for col in gpl_table.columns:
            sample_vals = gpl_table[col].dropna().head(5).tolist()
            if any(
                isinstance(v, str) and re.match(r"^[A-Z][A-Z0-9]+$", str(v))
                for v in sample_vals
            ):
                gene_col = col
                print(f"Using '{col}' as gene column")
                break

    print(f"Gene symbol column: {gene_col}")

    # Build probe-to-gene map
    probe_to_gene = {}
    for _, row in gpl_table.iterrows():
        probe_id = str(row["ID"])
        gene = str(row.get(gene_col, ""))
        if gene and gene != "nan" and gene != "---":
            # Handle multi-gene annotations: "GENE1 /// GENE2"
            genes = [g.strip() for g in gene.split("///")]
            if genes:
                probe_to_gene[probe_id] = genes[0]  # Take first gene

    print(f"Mapped {len(probe_to_gene)} probes to gene symbols")

    # Extract expression values from each sample
    sample_ids = list(gse.gsms.keys())
    print(f"Processing {len(sample_ids)} samples...")

    # Get probe IDs from first sample
    first_gsm = gse.gsms[sample_ids[0]]
    probe_ids = first_gsm.table["ID_REF"].values

    # Build expression matrix (probes x samples)
    expr_data = {}
    for i, gsm_id in enumerate(sample_ids):
        gsm = gse.gsms[gsm_id]
        table = gsm.table
        if "VALUE" in table.columns:
            expr_data[gsm_id] = table.set_index("ID_REF")["VALUE"].values
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{len(sample_ids)} samples")

    expr_df = pd.DataFrame(expr_data, index=probe_ids)
    print(f"Raw expression matrix: {expr_df.shape} (probes x samples)")

    # Convert to numeric
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    # Map probes to genes (average duplicate genes)
    expr_df["gene"] = expr_df.index.map(lambda x: probe_to_gene.get(str(x), None))
    expr_df = expr_df.dropna(subset=["gene"])
    expr_df = expr_df.groupby("gene").mean()

    # Remove genes with low variance
    gene_var = expr_df.var(axis=1)
    expr_df = expr_df[gene_var > 0.01]

    print(f"Gene expression matrix: {expr_df.shape} (genes x samples)")

    # Check if data is log-transformed (microarray RMA is typically log2)
    max_val = expr_df.max().max()
    min_val = expr_df.min().min()
    print(f"Value range: [{min_val:.2f}, {max_val:.2f}]")
    if max_val < 20:
        print("Data appears to be log2-transformed (RMA normalized)")
        input_type = "log2"
    else:
        print("Data appears to be raw intensities")
        input_type = "raw"

    # Transpose: PSF wants (samples x genes) or auto-detects
    # Save as genes-as-rows (PSF auto-detects)
    out_path = os.path.join(DATA_DIR, "GSE65682_expression_matrix.csv")
    expr_df.to_csv(out_path)
    print(f"Saved expression matrix to {out_path}")

    return expr_df, input_type, sample_ids


def extract_clinical_metadata(gse, sample_ids):
    """Extract clinical phenotype data from GEO sample metadata."""
    print("\nExtracting clinical metadata...")

    records = []
    for gsm_id in sample_ids:
        gsm = gse.gsms[gsm_id]
        meta = gsm.metadata

        record = {"sample_id": gsm_id}

        # Extract all characteristics
        chars = meta.get("characteristics_ch1", [])
        for ch in chars:
            if ":" in ch:
                key, val = ch.split(":", 1)
                key = key.strip().lower().replace(" ", "_")
                val = val.strip()
                record[key] = val

        # Also grab title and description
        record["title"] = meta.get("title", [""])[0]
        record["source"] = meta.get("source_name_ch1", [""])[0]

        records.append(record)

    pheno_df = pd.DataFrame(records)
    print(f"Clinical metadata: {pheno_df.shape}")
    print(f"Columns: {list(pheno_df.columns)}")

    # Print sample of key clinical variables
    for col in pheno_df.columns:
        n_unique = pheno_df[col].nunique()
        if n_unique <= 10:
            print(f"  {col}: {pheno_df[col].value_counts().to_dict()}")
        else:
            print(f"  {col}: {n_unique} unique values")

    out_path = os.path.join(DATA_DIR, "GSE65682_phenotypes.csv")
    pheno_df.to_csv(out_path, index=False)
    print(f"Saved clinical metadata to {out_path}")

    return pheno_df


def verify_data(expr_df, pheno_df, input_type):
    """Verify data integrity."""
    print("\n" + "=" * 60)
    print("Data Verification")
    print("=" * 60)

    n_samples = expr_df.shape[1]
    n_genes = expr_df.shape[0]
    n_pheno = len(pheno_df)

    print(f"Expression matrix: {n_genes} genes x {n_samples} samples")
    print(f"Phenotype records: {n_pheno}")
    print(f"Expression type: {input_type}")

    # Check overlap
    expr_samples = set(expr_df.columns)
    pheno_samples = set(pheno_df["sample_id"])
    overlap = expr_samples & pheno_samples
    print(f"Sample overlap: {len(overlap)} / {n_samples}")

    # Missing values
    missing_pct = (expr_df.isna().sum().sum()) / (n_genes * n_samples) * 100
    print(f"Missing values in expression: {missing_pct:.2f}%")

    # Check for mortality data
    mortality_cols = [
        c
        for c in pheno_df.columns
        if any(
            kw in c.lower()
            for kw in ["mortal", "death", "survival", "alive", "outcome", "icu"]
        )
    ]
    print(f"Mortality-related columns: {mortality_cols}")

    return True


def main():
    # Step 1: Download
    gse = download_geo_data()

    # Step 2: Extract expression
    expr_df, input_type, sample_ids = extract_expression_matrix(gse)

    # Step 3: Extract clinical metadata
    pheno_df = extract_clinical_metadata(gse, sample_ids)

    # Step 4: Verify
    verify_data(expr_df, pheno_df, input_type)

    print("\n" + "=" * 60)
    print("Phase 1 COMPLETE")
    print("=" * 60)
    print(f"Expression matrix: data/sepsis/GSE65682_expression_matrix.csv")
    print(f"Phenotype data:    data/sepsis/GSE65682_phenotypes.csv")
    print(f"Expression type:   {input_type}")


if __name__ == "__main__":
    main()
