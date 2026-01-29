"""
Integration tests for the complete pathway subtyping pipeline.

These tests verify end-to-end functionality with real sample data.
"""

import json
import os
from pathlib import Path

import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig


@pytest.mark.integration
class TestEndToEndPipeline:
    """End-to-end integration tests."""

    def test_synthetic_data_pipeline(
        self, synthetic_vcf_path, synthetic_phenotypes_path, autism_pathways_path, tmp_path
    ):
        """Test complete pipeline with synthetic data."""
        config = PipelineConfig(
            name="e2e_test",
            output_dir=str(tmp_path / "e2e_output"),
            vcf_path=str(synthetic_vcf_path),
            phenotype_path=str(synthetic_phenotypes_path),
            pathway_db=str(autism_pathways_path),
            seed=42,
            n_clusters_range=[2, 6],
        )

        pipeline = DemoPipeline(config)
        pipeline.run()

        # Verify outputs
        output_dir = tmp_path / "e2e_output"

        # Check files exist
        assert (output_dir / "pathway_scores.csv").exists()
        assert (output_dir / "subtype_assignments.csv").exists()
        assert (output_dir / "report.json").exists()
        assert (output_dir / "report.md").exists()
        assert (output_dir / "run_metadata.yaml").exists()
        assert (output_dir / "pipeline.log").exists()
        assert (output_dir / "figures" / "summary.png").exists()

        # Check pathway scores
        scores = pd.read_csv(output_dir / "pathway_scores.csv", index_col=0)
        assert len(scores) == 60
        assert scores.shape[1] >= 1  # At least one pathway

        # Check assignments
        assignments = pd.read_csv(output_dir / "subtype_assignments.csv")
        assert len(assignments) == 60
        assert "sample_id" in assignments.columns
        assert "cluster_id" in assignments.columns
        assert "cluster_label" in assignments.columns
        assert "confidence" in assignments.columns

        # Check report
        with open(output_dir / "report.json") as f:
            report = json.load(f)

        assert report["pipeline_name"] == "e2e_test"
        assert report["summary"]["n_samples"] == 60
        assert "validation_gates" in report

    def test_planted_subtype_recovery(
        self, synthetic_vcf_path, synthetic_phenotypes_path, autism_pathways_path, tmp_path
    ):
        """Test that pipeline recovers planted subtypes."""
        config = PipelineConfig(
            name="recovery_test",
            output_dir=str(tmp_path / "recovery_output"),
            vcf_path=str(synthetic_vcf_path),
            phenotype_path=str(synthetic_phenotypes_path),
            pathway_db=str(autism_pathways_path),
            seed=42,
        )

        pipeline = DemoPipeline(config)
        pipeline.run()

        # Load assignments and phenotypes
        assignments = pd.read_csv(tmp_path / "recovery_output" / "subtype_assignments.csv")

        # The synthetic data has planted_subtype column
        assert "planted_subtype" in assignments.columns

        # Calculate ARI between discovered and planted
        ari = adjusted_rand_score(assignments["planted_subtype"], assignments["cluster_label"])

        # Should have high agreement (ARI > 0.7)
        # Note: With perfect synthetic data, should be very high
        assert ari > 0.7, f"ARI {ari} is too low - subtypes not recovered well"

    def test_validation_gates_run(
        self, synthetic_vcf_path, synthetic_phenotypes_path, autism_pathways_path, tmp_path
    ):
        """Test that validation gates execute and report results."""
        config = PipelineConfig(
            name="validation_test",
            output_dir=str(tmp_path / "validation_output"),
            vcf_path=str(synthetic_vcf_path),
            phenotype_path=str(synthetic_phenotypes_path),
            pathway_db=str(autism_pathways_path),
            seed=42,
        )

        pipeline = DemoPipeline(config)
        pipeline.run()

        # Load report
        with open(tmp_path / "validation_output" / "report.json") as f:
            report = json.load(f)

        # Check validation gates results
        gates = report["validation_gates"]
        assert "all_passed" in gates
        assert "tests" in gates
        assert len(gates["tests"]) == 3

        # Check individual tests
        test_names = [t["name"] for t in gates["tests"]]
        assert "Negative Control 1: Label Shuffle" in test_names
        assert "Negative Control 2: Random Gene Sets" in test_names
        assert "Stability Test: Bootstrap" in test_names


@pytest.mark.integration
class TestCLIIntegration:
    """Tests for CLI integration."""

    def test_cli_help(self):
        """Test CLI help command."""
        import subprocess

        result = subprocess.run(
            ["python", "-m", "pathway_subtyping", "--help"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "Pathway Subtyping Framework" in result.stdout
        assert "--config" in result.stdout

    def test_cli_version(self):
        """Test CLI version command."""
        import subprocess

        result = subprocess.run(
            ["python", "-m", "pathway_subtyping", "--version"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "0.1.0" in result.stdout

    def test_cli_run_with_config(
        self, synthetic_vcf_path, synthetic_phenotypes_path, autism_pathways_path, tmp_path
    ):
        """Test CLI run with config file."""
        import subprocess

        import yaml

        # Create config file
        config = {
            "pipeline": {
                "name": "cli_test",
                "output_dir": str(tmp_path / "cli_output"),
                "seed": 42,
            },
            "data": {
                "vcf_path": str(synthetic_vcf_path),
                "phenotype_path": str(synthetic_phenotypes_path),
                "pathway_db": str(autism_pathways_path),
            },
        }

        config_file = tmp_path / "test_config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config, f)

        result = subprocess.run(
            ["python", "-m", "pathway_subtyping", "--config", str(config_file)],
            capture_output=True,
            text=True,
            cwd=str(Path(__file__).parent.parent),
        )

        assert result.returncode == 0
        assert (tmp_path / "cli_output" / "report.json").exists()


@pytest.mark.integration
class TestDataFormats:
    """Tests for different data format handling."""

    def test_vcf_parsing(self, synthetic_vcf_path, tmp_path):
        """Test VCF file is parsed correctly."""
        from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

        config = PipelineConfig(
            name="vcf_test",
            output_dir=str(tmp_path),
            vcf_path=str(synthetic_vcf_path),
        )

        pipeline = DemoPipeline(config)
        pipeline.setup()
        pipeline._load_vcf()

        # Check variants
        assert pipeline.variants_df is not None
        assert "chrom" in pipeline.variants_df.columns
        assert "pos" in pipeline.variants_df.columns
        assert "gene" in pipeline.variants_df.columns
        assert "consequence" in pipeline.variants_df.columns
        assert "cadd" in pipeline.variants_df.columns

        # Check genotypes
        assert pipeline.genotypes_df is not None
        assert len(pipeline.genotypes_df.columns) == 60  # 60 samples

        # Check samples
        assert len(pipeline.samples) == 60
        assert all(s.startswith("SAMPLE_") for s in pipeline.samples)

    def test_gmt_parsing(self, autism_pathways_path, tmp_path):
        """Test GMT file is parsed correctly."""
        from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

        config = PipelineConfig(
            name="gmt_test",
            output_dir=str(tmp_path),
            pathway_db=str(autism_pathways_path),
        )

        pipeline = DemoPipeline(config)
        pipeline.setup()
        pipeline._load_pathways()

        assert len(pipeline.pathways) == 15

        # Check pathway structure
        for name, genes in pipeline.pathways.items():
            assert isinstance(name, str)
            assert isinstance(genes, list)
            assert len(genes) > 0
            assert all(isinstance(g, str) for g in genes)

        # Check specific pathways
        assert "SYNAPTIC_TRANSMISSION" in pipeline.pathways
        assert "SHANK3" in pipeline.pathways["SYNAPTIC_TRANSMISSION"]

    def test_phenotype_parsing(self, synthetic_phenotypes_path, tmp_path):
        """Test phenotype file is parsed correctly."""
        from pathway_subtyping.pipeline import DemoPipeline, PipelineConfig

        config = PipelineConfig(
            name="pheno_test",
            output_dir=str(tmp_path),
            phenotype_path=str(synthetic_phenotypes_path),
        )

        pipeline = DemoPipeline(config)
        pipeline.setup()
        pipeline._load_phenotypes()

        assert pipeline.phenotypes_df is not None
        assert len(pipeline.phenotypes_df) == 60

        # Check index is sample_id
        assert pipeline.phenotypes_df.index.name == "sample_id"

        # Check columns
        assert "sex" in pipeline.phenotypes_df.columns
        assert "planted_subtype" in pipeline.phenotypes_df.columns


@pytest.mark.integration
class TestOutputFormats:
    """Tests for output format correctness."""

    def test_json_report_structure(
        self, synthetic_vcf_path, synthetic_phenotypes_path, autism_pathways_path, tmp_path
    ):
        """Test JSON report has correct structure."""
        config = PipelineConfig(
            name="json_test",
            output_dir=str(tmp_path / "json_output"),
            vcf_path=str(synthetic_vcf_path),
            phenotype_path=str(synthetic_phenotypes_path),
            pathway_db=str(autism_pathways_path),
            seed=42,
        )

        pipeline = DemoPipeline(config)
        pipeline.run()

        with open(tmp_path / "json_output" / "report.json") as f:
            report = json.load(f)

        # Required fields
        assert "pipeline_name" in report
        assert "timestamp" in report
        assert "seed" in report
        assert "input_files" in report
        assert "summary" in report
        assert "clusters" in report
        assert "disclaimer" in report

        # Summary fields
        assert "n_variants" in report["summary"]
        assert "n_samples" in report["summary"]
        assert "n_genes" in report["summary"]
        assert "n_pathways" in report["summary"]
        assert "n_clusters" in report["summary"]

    def test_markdown_report_content(
        self, synthetic_vcf_path, synthetic_phenotypes_path, autism_pathways_path, tmp_path
    ):
        """Test Markdown report has expected sections."""
        config = PipelineConfig(
            name="md_test",
            output_dir=str(tmp_path / "md_output"),
            vcf_path=str(synthetic_vcf_path),
            phenotype_path=str(synthetic_phenotypes_path),
            pathway_db=str(autism_pathways_path),
            seed=42,
        )

        pipeline = DemoPipeline(config)
        pipeline.run()

        with open(tmp_path / "md_output" / "report.md") as f:
            content = f.read()

        # Check sections
        assert "# Pathway Subtyping Framework" in content
        assert "## Input Summary" in content
        assert "## Clustering Results" in content
        assert "## Validation Gates" in content
        assert "## Disclaimer" in content
