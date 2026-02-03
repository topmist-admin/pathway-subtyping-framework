"""
Tests for configuration loading and validation.
"""

import tempfile
from pathlib import Path

import pytest
import yaml

from pathway_subtyping.config import (
    ConfigValidationError,
    load_config,
    validate_config,
    validate_gmt_file,
)


class TestLoadConfig:
    """Tests for load_config function."""

    def test_load_config_file_not_found(self):
        """Test that missing config file raises error."""
        with pytest.raises(FileNotFoundError):
            load_config("nonexistent.yaml")

    def test_load_config_valid_yaml(self, tmp_path):
        """Test loading a valid YAML config."""
        config_content = {
            "pipeline": {"name": "test", "seed": 42},
            "data": {"vcf_path": "test.vcf"},
        }
        config_file = tmp_path / "config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_content, f)

        config = load_config(str(config_file))

        assert config["pipeline"]["name"] == "test"
        assert config["pipeline"]["seed"] == 42
        assert config["data"]["vcf_path"] == "test.vcf"

    def test_load_config_empty_file(self, tmp_path):
        """Test loading an empty YAML file raises error."""
        config_file = tmp_path / "empty.yaml"
        config_file.touch()

        with pytest.raises(ConfigValidationError, match="empty or invalid"):
            load_config(str(config_file))

    def test_load_config_invalid_yaml(self, tmp_path):
        """Test loading invalid YAML raises error."""
        config_file = tmp_path / "invalid.yaml"
        with open(config_file, "w") as f:
            f.write("invalid: yaml: content: [")

        with pytest.raises(yaml.YAMLError):
            load_config(str(config_file))

    def test_load_config_with_all_sections(self, tmp_path):
        """Test loading config with all expected sections."""
        config_content = {
            "pipeline": {
                "name": "full_test",
                "output_dir": "outputs/test",
                "seed": 123,
                "verbose": True,
            },
            "data": {
                "vcf_path": "data/test.vcf",
                "phenotype_path": "data/phenotypes.csv",
                "pathway_db": "data/pathways.gmt",
            },
            "clustering": {
                "method": "gmm",
                "n_clusters": None,
                "n_clusters_range": [2, 8],
            },
            "validation": {
                "run_gates": True,
                "bootstrap_iterations": 50,
            },
        }
        config_file = tmp_path / "full_config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_content, f)

        config = load_config(str(config_file))

        assert config["pipeline"]["name"] == "full_test"
        assert config["clustering"]["n_clusters_range"] == [2, 8]
        assert config["validation"]["run_gates"] is True


class TestValidateConfig:
    """Tests for validate_config function."""

    def test_validate_config_valid(self):
        """Test validation passes for valid config (without file checks)."""
        config = {
            "pipeline": {"name": "test"},
            "data": {
                "vcf_path": "test.vcf",
                "phenotype_path": "test.csv",
                "pathway_db": "test.gmt",
            },
        }
        # Should not raise when file checks disabled
        result = validate_config(config, check_files=False)
        assert result is True

    def test_validate_config_missing_pipeline(self):
        """Test validation fails without pipeline section."""
        config = {"data": {"vcf_path": "test.vcf"}}

        with pytest.raises(ValueError, match="pipeline"):
            validate_config(config)

    def test_validate_config_missing_data(self):
        """Test validation fails without data section."""
        config = {"pipeline": {"name": "test"}}

        with pytest.raises(ValueError, match="data"):
            validate_config(config)

    def test_validate_config_empty(self):
        """Test validation fails for empty config."""
        with pytest.raises(ValueError):
            validate_config({})

    def test_validate_config_none(self):
        """Test validation fails for None config."""
        with pytest.raises(TypeError):
            validate_config(None)

    def test_validate_config_missing_vcf_path(self):
        """Test validation fails without vcf_path."""
        config = {
            "pipeline": {"name": "test"},
            "data": {"phenotype_path": "test.csv", "pathway_db": "test.gmt"},
        }
        with pytest.raises(ValueError, match="vcf_path"):
            validate_config(config)

    def test_validate_config_missing_phenotype_path(self):
        """Test validation fails without phenotype_path."""
        config = {
            "pipeline": {"name": "test"},
            "data": {"vcf_path": "test.vcf", "pathway_db": "test.gmt"},
        }
        with pytest.raises((ValueError, ConfigValidationError), match="phenotype_path"):
            validate_config(config, check_files=False)

    def test_validate_config_missing_pathway_db(self):
        """Test validation fails without pathway_db."""
        config = {
            "pipeline": {"name": "test"},
            "data": {"vcf_path": "test.vcf", "phenotype_path": "test.csv"},
        }
        with pytest.raises((ValueError, ConfigValidationError), match="pathway_db"):
            validate_config(config, check_files=False)


class TestPipelineConfig:
    """Tests for PipelineConfig dataclass."""

    def test_from_yaml_basic(self, test_config_path):
        """Test loading PipelineConfig from YAML."""
        from pathway_subtyping.pipeline import PipelineConfig

        config = PipelineConfig.from_yaml(str(test_config_path))

        assert config.name == "synthetic_test"
        assert config.seed == 42
        assert config.verbose is True

    def test_from_yaml_defaults(self, tmp_path):
        """Test PipelineConfig uses defaults for missing values."""
        from pathway_subtyping.pipeline import PipelineConfig

        config_content = {
            "pipeline": {"name": "minimal"},
            "data": {},
        }
        config_file = tmp_path / "minimal.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_content, f)

        config = PipelineConfig.from_yaml(str(config_file))

        assert config.name == "minimal"
        assert config.seed == 42  # Default
        assert config.n_clusters_range == [2, 8]  # Default

    def test_from_yaml_cluster_range(self, tmp_path):
        """Test PipelineConfig parses cluster range correctly."""
        from pathway_subtyping.pipeline import PipelineConfig

        config_content = {
            "pipeline": {"name": "test"},
            "data": {},
            "clustering": {"n_clusters_range": [3, 10]},
        }
        config_file = tmp_path / "cluster_config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_content, f)

        config = PipelineConfig.from_yaml(str(config_file))

        assert config.n_clusters_range == [3, 10]

    def test_from_yaml_fixed_clusters(self, tmp_path):
        """Test PipelineConfig with fixed number of clusters."""
        from pathway_subtyping.pipeline import PipelineConfig

        config_content = {
            "pipeline": {"name": "test"},
            "data": {},
            "clustering": {"n_clusters": 4},
        }
        config_file = tmp_path / "fixed_clusters.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_content, f)

        config = PipelineConfig.from_yaml(str(config_file))

        assert config.n_clusters == 4
