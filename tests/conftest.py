"""
Pytest configuration and shared fixtures for pathway_subtyping tests.
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def project_root():
    """Return the project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def sample_data_dir(project_root):
    """Return the sample data directory."""
    return project_root / "data" / "sample"


@pytest.fixture
def config_dir(project_root):
    """Return the configs directory."""
    return project_root / "configs"


@pytest.fixture
def pathways_dir(project_root):
    """Return the pathways directory."""
    return project_root / "data" / "pathways"


@pytest.fixture
def synthetic_vcf_path(sample_data_dir):
    """Return path to synthetic VCF file."""
    return sample_data_dir / "synthetic_cohort.vcf"


@pytest.fixture
def synthetic_phenotypes_path(sample_data_dir):
    """Return path to synthetic phenotypes file."""
    return sample_data_dir / "synthetic_phenotypes.csv"


@pytest.fixture
def autism_pathways_path(pathways_dir):
    """Return path to autism pathways GMT file."""
    return pathways_dir / "autism_pathways.gmt"


@pytest.fixture
def test_config_path(config_dir):
    """Return path to test synthetic config."""
    return config_dir / "test_synthetic.yaml"


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_pathway_scores():
    """Generate sample pathway scores for testing."""
    np.random.seed(42)
    n_samples = 60
    n_pathways = 5

    # Create pathway scores with 4 distinct clusters
    scores = np.zeros((n_samples, n_pathways))

    # Cluster 1: High in pathway 0
    scores[:15, 0] = np.random.normal(2, 0.3, 15)
    scores[:15, 1:] = np.random.normal(0, 0.3, (15, 4))

    # Cluster 2: High in pathway 1
    scores[15:30, 1] = np.random.normal(2, 0.3, 15)
    scores[15:30, [0, 2, 3, 4]] = np.random.normal(0, 0.3, (15, 4))

    # Cluster 3: High in pathway 2
    scores[30:45, 2] = np.random.normal(2, 0.3, 15)
    scores[30:45, [0, 1, 3, 4]] = np.random.normal(0, 0.3, (15, 4))

    # Cluster 4: High in pathway 3
    scores[45:60, 3] = np.random.normal(2, 0.3, 15)
    scores[45:60, [0, 1, 2, 4]] = np.random.normal(0, 0.3, (15, 4))

    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    pathway_names = ["PATHWAY_A", "PATHWAY_B", "PATHWAY_C", "PATHWAY_D", "PATHWAY_E"]

    return pd.DataFrame(scores, index=sample_ids, columns=pathway_names)


@pytest.fixture
def sample_cluster_labels():
    """Generate sample cluster labels for testing."""
    labels = np.array([0] * 15 + [1] * 15 + [2] * 15 + [3] * 15)
    return labels


@pytest.fixture
def sample_pathways():
    """Generate sample pathway definitions for testing."""
    return {
        "PATHWAY_A": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
        "PATHWAY_B": ["GENE6", "GENE7", "GENE8", "GENE9", "GENE10"],
        "PATHWAY_C": ["GENE11", "GENE12", "GENE13", "GENE14", "GENE15"],
        "PATHWAY_D": ["GENE16", "GENE17", "GENE18", "GENE19", "GENE20"],
        "PATHWAY_E": ["GENE21", "GENE22", "GENE23", "GENE24", "GENE25"],
    }


@pytest.fixture
def sample_gene_burdens():
    """Generate sample gene burden matrix for testing."""
    np.random.seed(42)
    n_samples = 60
    n_genes = 25

    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    gene_names = [f"GENE{i}" for i in range(1, n_genes + 1)]

    burdens = np.random.exponential(0.5, (n_samples, n_genes))

    return pd.DataFrame(burdens, index=sample_ids, columns=gene_names)


@pytest.fixture
def minimal_vcf_content():
    """Return minimal VCF content for testing."""
    return """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_001\tSAMPLE_002
chr1\t100\tvar1\tA\tG\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=missense_variant;CADD=28.5\tGT\t0/1\t0/0
chr1\t200\tvar2\tC\tT\t99\tPASS\tGENE=CHD8;CONSEQUENCE=frameshift_variant;CADD=35.2\tGT\t0/0\t0/1
"""


@pytest.fixture
def minimal_phenotypes_content():
    """Return minimal phenotypes CSV content for testing."""
    return """sample_id,sex,age,planted_subtype
SAMPLE_001,M,5.2,synaptic
SAMPLE_002,F,4.8,chromatin
"""


@pytest.fixture
def minimal_gmt_content():
    """Return minimal GMT content for testing."""
    return """SYNAPTIC\thttp://example.com\tSHANK3\tSHANK2\tNRXN1
CHROMATIN\thttp://example.com\tCHD8\tCHD2\tARID1B
"""
