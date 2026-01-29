"""
Configuration loading and validation.

This is a placeholder - copy/adapt from autism-pathway-framework.
"""

from pathlib import Path
from typing import Dict, Any
import yaml


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load and validate pipeline configuration.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Validated configuration dictionary
    """
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(path) as f:
        config = yaml.safe_load(f)

    # TODO: Add validation logic
    return config
