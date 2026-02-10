#!/usr/bin/env nextflow

/*
 * Pathway Subtyping Framework — Nextflow Workflow
 *
 * A minimal Nextflow DSL2 pipeline that runs pathway-based molecular
 * subtype discovery and characterization.
 *
 * Usage:
 *   nextflow run main.nf \
 *     --vcf data/cohort.vcf.gz \
 *     --phenotypes data/phenotypes.csv \
 *     --pathways data/pathways/autism_pathways.gmt \
 *     --outdir results/
 *
 * Requirements:
 *   - Nextflow >= 22.04
 *   - Python 3.9+ with pathway-subtyping installed
 *   - Or use the provided conda environment (environment.yaml)
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------

params.vcf        = null
params.phenotypes = null
params.pathways   = null
params.outdir     = 'results'
params.n_clusters = null   // auto-select if not set
params.seed       = 42

// ---------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------

process SUBTYPE_DISCOVERY {
    tag "subtyping"
    publishDir "${params.outdir}/pipeline", mode: 'copy'

    input:
    path vcf
    path phenotypes
    path pathways

    output:
    path "pipeline_output/**", emit: results
    path "pipeline_output/report.json", emit: report_json

    script:
    def clusters_arg = params.n_clusters ? "--n-clusters ${params.n_clusters}" : ""
    """
    python3 -c "
import json, yaml

config = {
    'pipeline': {
        'name': 'nextflow_run',
        'output_dir': 'pipeline_output',
        'seed': ${params.seed},
    },
    'data': {
        'vcf_path': '${vcf}',
        'phenotype_path': '${phenotypes}',
        'pathway_db': '${pathways}',
    },
    'clustering': {
        'n_clusters': ${params.n_clusters ?: 'None'},
    },
    'output': {
        'disclaimer': 'Research use only. Not for clinical decision-making.',
    },
}

with open('config.yaml', 'w') as f:
    yaml.dump(config, f)
"
    python3 -m pathway_subtyping.cli run --config config.yaml
    """
}

process CHARACTERIZE {
    tag "characterization"
    publishDir "${params.outdir}/characterization", mode: 'copy'

    input:
    path pipeline_results

    output:
    path "characterization_output/*", emit: results

    script:
    """
    python3 -c "
import json
import pandas as pd
import numpy as np
from pathway_subtyping import characterize_subtypes, export_characterization, generate_subtype_heatmap

# Load pipeline outputs
pathway_scores = pd.read_csv('pipeline_output/pathway_scores.csv', index_col=0)
assignments = pd.read_csv('pipeline_output/cluster_assignments.csv')

cluster_labels = assignments['cluster_id'].values
confidence = assignments['confidence'].values if 'confidence' in assignments.columns else None

# Build cluster name mapping
names = {}
for _, row in assignments.drop_duplicates('cluster_id').iterrows():
    names[int(row['cluster_id'])] = row['cluster_label']

# Characterize
result = characterize_subtypes(
    pathway_scores=pathway_scores,
    cluster_labels=cluster_labels,
    cluster_names=names,
    confidence_scores=confidence,
    seed=${params.seed},
)

# Export
import os
os.makedirs('characterization_output', exist_ok=True)
export_characterization(result, 'characterization_output', formats=['csv'])
generate_subtype_heatmap(result, output_path='characterization_output/subtype_heatmap.png')

# Save JSON summary
with open('characterization_output/characterization.json', 'w') as f:
    json.dump(result.to_dict(), f, indent=2)

print(result.format_report())
"
    """
}

// ---------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------

workflow {
    // Validate inputs
    if (!params.vcf || !params.phenotypes || !params.pathways) {
        error "Required parameters: --vcf, --phenotypes, --pathways"
    }

    vcf_ch        = Channel.fromPath(params.vcf, checkIfExists: true)
    phenotypes_ch = Channel.fromPath(params.phenotypes, checkIfExists: true)
    pathways_ch   = Channel.fromPath(params.pathways, checkIfExists: true)

    // Run pipeline
    SUBTYPE_DISCOVERY(vcf_ch, phenotypes_ch, pathways_ch)

    // Characterize subtypes
    CHARACTERIZE(SUBTYPE_DISCOVERY.out.results.collect())
}
