"""
Tests for configuration loading.
"""

import pytest
from pathway_subtyping.config import load_config


def test_load_config_file_not_found():
    """Test that missing config file raises error."""
    with pytest.raises(FileNotFoundError):
        load_config("nonexistent.yaml")


# TODO: Add more tests after implementing config validation
