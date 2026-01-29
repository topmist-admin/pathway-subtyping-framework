"""
Configuration loading and validation.

Provides utility functions for loading pipeline configuration from YAML files.
"""

from pathlib import Path
from typing import Any, Dict

import yaml


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load pipeline configuration from YAML file.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Configuration dictionary
    """
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(path) as f:
        config = yaml.safe_load(f)

    return config


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate pipeline configuration.

    Args:
        config: Configuration dictionary

    Returns:
        True if valid

    Raises:
        ValueError: If configuration is invalid
    """
    required_sections = ["pipeline", "data"]
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required config section: {section}")

    data = config.get("data", {})
    if not data.get("vcf_path"):
        raise ValueError("data.vcf_path is required")
    if not data.get("phenotype_path"):
        raise ValueError("data.phenotype_path is required")
    if not data.get("pathway_db"):
        raise ValueError("data.pathway_db is required")

    return True
