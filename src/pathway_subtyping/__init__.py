"""
Pathway Subtyping Framework

A disease-agnostic tool for pathway-based molecular subtype discovery.
"""

__version__ = "0.1.0"
__author__ = "Rohit Chauhan"
__email__ = "info@topmist.com"

from .pipeline import run_pipeline
from .config import load_config

__all__ = ["run_pipeline", "load_config", "__version__"]
