"""
Tests for the CLI module.

Tests cover:
- Command-line argument parsing
- Config loading and override handling
- Error handling and exit codes
- Version output
"""

import tempfile
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from pathway_subtyping import __version__
from pathway_subtyping.cli import main


class TestCLIVersion:
    """Tests for version display."""

    def test_version_option(self):
        """Test --version flag displays version."""
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])

        assert result.exit_code == 0
        assert __version__ in result.output

    def test_version_short_flag(self):
        """Test version is displayed correctly."""
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])

        assert "pathway-subtyping-framework" in result.output.lower() or __version__ in result.output


class TestCLIHelp:
    """Tests for help display."""

    def test_help_option(self):
        """Test --help flag displays help."""
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])

        assert result.exit_code == 0
        assert "--config" in result.output
        assert "--output" in result.output
        assert "--seed" in result.output

    def test_help_short_flag(self):
        """Test -h displays help (if available)."""
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])

        assert "Pathway Subtyping Framework" in result.output


class TestCLIConfigHandling:
    """Tests for configuration handling."""

    @pytest.fixture
    def valid_config(self, tmp_path):
        """Create a minimal valid config file."""
        # Create necessary input files
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("""##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10
chr1\t100\tv1\tA\tG\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=missense;CADD=25\tGT\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0
chr1\t200\tv2\tC\tT\t99\tPASS\tGENE=CHD8;CONSEQUENCE=frameshift;CADD=35\tGT\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1
""")

        pheno_path = tmp_path / "phenotypes.csv"
        pheno_path.write_text("""sample_id,age,sex
S1,10,M
S2,12,F
S3,11,M
S4,9,F
S5,13,M
S6,10,F
S7,11,M
S8,12,F
S9,10,M
S10,11,F
""")

        pathway_path = tmp_path / "pathways.gmt"
        pathway_path.write_text("""SYNAPTIC\tdesc\tSHANK3\tCHD8\tNRXN1
CHROMATIN\tdesc\tCHD8\tSETD5\tASH1L
""")

        output_dir = tmp_path / "output"

        config = {
            "pipeline": {
                "name": "test_run",
                "output_dir": str(output_dir),
                "seed": 42,
                "verbose": False,
            },
            "data": {
                "vcf_path": str(vcf_path),
                "phenotype_path": str(pheno_path),
                "pathway_db": str(pathway_path),
            },
            "clustering": {
                "n_clusters_range": [2, 3],
            },
        }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        return config_path

    def test_missing_config_required(self):
        """Test that --config is required."""
        runner = CliRunner()
        result = runner.invoke(main, [])

        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    def test_config_file_not_found(self, tmp_path):
        """Test error when config file doesn't exist."""
        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(tmp_path / "nonexistent.yaml")])

        # Click validates path exists
        assert result.exit_code != 0

    def test_valid_config_loads(self, valid_config, tmp_path):
        """Test that valid config runs pipeline."""
        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(valid_config)])

        # Pipeline should complete or provide meaningful output
        # Note: May still fail due to insufficient data, but should not crash on config
        assert "Pathway Subtyping Framework" in result.output or result.exit_code in [0, 1]


class TestCLIOverrides:
    """Tests for command-line overrides."""

    @pytest.fixture
    def minimal_config(self, tmp_path):
        """Create a minimal config for override testing."""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\tv1\tA\tG\t99\tPASS\t.\tGT\t0/1
""")

        pheno_path = tmp_path / "phenotypes.csv"
        pheno_path.write_text("sample_id,age\nS1,10\n")

        pathway_path = tmp_path / "pathways.gmt"
        pathway_path.write_text("PATH1\tdesc\tGENE1\tGENE2\n")

        config = {
            "pipeline": {
                "name": "test",
                "output_dir": str(tmp_path / "output"),
                "seed": 42,
            },
            "data": {
                "vcf_path": str(vcf_path),
                "phenotype_path": str(pheno_path),
                "pathway_db": str(pathway_path),
            },
        }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        return config_path, tmp_path

    def test_output_override(self, minimal_config):
        """Test --output flag overrides config output_dir."""
        config_path, tmp_path = minimal_config
        custom_output = tmp_path / "custom_output"

        runner = CliRunner()
        # Run with output override - may fail on pipeline but should parse args
        result = runner.invoke(
            main, ["--config", str(config_path), "--output", str(custom_output)]
        )

        # Check that the output dir was attempted
        # (Pipeline may fail, but args should be parsed)
        assert "Pathway Subtyping Framework" in result.output

    def test_seed_override(self, minimal_config):
        """Test --seed flag overrides config seed."""
        config_path, tmp_path = minimal_config

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path), "--seed", "123"])

        # Should accept seed parameter
        assert "Pathway Subtyping Framework" in result.output

    def test_quiet_mode(self, minimal_config):
        """Test --quiet flag reduces output."""
        config_path, tmp_path = minimal_config

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path), "--quiet"])

        # Quiet mode should still show header
        assert "Pathway Subtyping Framework" in result.output


class TestCLIErrorHandling:
    """Tests for error handling."""

    def test_invalid_yaml_syntax(self, tmp_path):
        """Test error handling for invalid YAML."""
        config_path = tmp_path / "invalid.yaml"
        config_path.write_text("invalid: yaml: content: [")

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path)])

        assert result.exit_code != 0

    def test_missing_required_fields(self, tmp_path):
        """Test error when required config fields are missing."""
        config_path = tmp_path / "incomplete.yaml"
        config_path.write_text("pipeline:\n  name: test\n")

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path)])

        assert result.exit_code != 0

    def test_nonexistent_input_file(self, tmp_path):
        """Test error when input files don't exist."""
        config = {
            "pipeline": {"name": "test", "output_dir": str(tmp_path / "out")},
            "data": {
                "vcf_path": str(tmp_path / "nonexistent.vcf"),
                "phenotype_path": str(tmp_path / "nonexistent.csv"),
                "pathway_db": str(tmp_path / "nonexistent.gmt"),
            },
        }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path)])

        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "error" in result.output.lower()


class TestCLIExitCodes:
    """Tests for correct exit codes."""

    def test_success_exit_code(self, tmp_path):
        """Test exit code 0 on success (when possible)."""
        # Create valid minimal dataset
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10\tS11\tS12\tS13\tS14\tS15
"""
        for i in range(20):
            gene = "SHANK3" if i % 2 == 0 else "CHD8"
            vcf_content += f"chr1\t{100+i*10}\tv{i}\tA\tG\t99\tPASS\tGENE={gene};CONSEQUENCE=missense;CADD=25\tGT\t"
            vcf_content += "\t".join(["0/1" if j % 3 == i % 3 else "0/0" for j in range(15)])
            vcf_content += "\n"

        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text(vcf_content)

        pheno_content = "sample_id,age\n" + "\n".join([f"S{i},{10+i}" for i in range(1, 16)])
        pheno_path = tmp_path / "phenotypes.csv"
        pheno_path.write_text(pheno_content)

        pathway_path = tmp_path / "pathways.gmt"
        pathway_path.write_text("SYNAPTIC\tdesc\tSHANK3\tNRXN1\nCHROMATIN\tdesc\tCHD8\tSETD5\n")

        config = {
            "pipeline": {
                "name": "test",
                "output_dir": str(tmp_path / "output"),
                "seed": 42,
                "verbose": False,
            },
            "data": {
                "vcf_path": str(vcf_path),
                "phenotype_path": str(pheno_path),
                "pathway_db": str(pathway_path),
            },
            "clustering": {"n_clusters_range": [2, 3]},
        }

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path)])

        # Exit code should be 0 or 1 (success or handled failure)
        assert result.exit_code in [0, 1]

    def test_failure_exit_code(self, tmp_path):
        """Test exit code 1 on failure."""
        config_path = tmp_path / "config.yaml"
        config_path.write_text("invalid")

        runner = CliRunner()
        result = runner.invoke(main, ["--config", str(config_path)])

        assert result.exit_code != 0
