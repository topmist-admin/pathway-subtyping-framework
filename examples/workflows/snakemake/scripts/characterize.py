"""
Snakemake script for subtype characterization.

Called by the 'characterize' rule in the Snakefile.
Reads pipeline outputs and generates characterization results.
"""

import json

import numpy as np
import pandas as pd

from pathway_subtyping import (
    characterize_subtypes,
    export_characterization,
    generate_subtype_heatmap,
)

# Snakemake provides input/output/params via the `snakemake` object
pathway_scores = pd.read_csv(snakemake.input.scores, index_col=0)
assignments = pd.read_csv(snakemake.input.assignments)

cluster_labels = assignments["cluster_id"].values

confidence = None
if "confidence" in assignments.columns:
    confidence = assignments["confidence"].values

# Build cluster name mapping
names = {}
for _, row in assignments.drop_duplicates("cluster_id").iterrows():
    names[int(row["cluster_id"])] = row["cluster_label"]

# Run characterization
result = characterize_subtypes(
    pathway_scores=pathway_scores,
    cluster_labels=cluster_labels,
    cluster_names=names,
    confidence_scores=confidence,
)

# Export CSV files
export_characterization(result, snakemake.params.outdir, formats=["csv"])

# Generate heatmap
generate_subtype_heatmap(result, output_path=snakemake.output.heatmap)

# Save JSON summary
with open(snakemake.output.report_json, "w") as f:
    json.dump(result.to_dict(), f, indent=2)

print(result.format_report())
