"""
Reproducibility utilities for deterministic execution.

This module provides functions to set random seeds across all
stochastic components for reproducible pipeline runs.
"""

import random
from typing import Optional

import numpy as np


def set_global_seed(seed: int) -> None:
    """Set seed for all random number generators.

    Args:
        seed: Integer seed value for reproducibility
    """
    random.seed(seed)
    np.random.seed(seed)

    # Try to set PyTorch seed if available
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass  # PyTorch not installed


def get_rng(seed: Optional[int], module_name: str = "") -> np.random.RandomState:
    """Get a module-specific random state for isolated reproducibility.

    Creates a deterministic offset from the module name to ensure different
    modules get different but reproducible random sequences.

    Args:
        seed: Base seed (None for random)
        module_name: Module identifier for offset calculation

    Returns:
        NumPy RandomState instance
    """
    if seed is None:
        return np.random.RandomState()

    # Create deterministic offset from module name
    offset = sum(ord(c) for c in module_name) % 1000
    return np.random.RandomState(seed + offset)
