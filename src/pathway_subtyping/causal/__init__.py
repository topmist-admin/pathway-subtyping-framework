"""
Causal inference layer for the Pathway Subtyping Framework.

v0.6 Phase 3 F11: PSF's cascade scoring (F1) is correlational by design
— "upstream and downstream are co-activated." A causal statement — "up
causes down" — requires either intervention or invariant prediction
across environments.

This subpackage implements **invariant causal prediction (ICP)** for
pathway-level relationships when users can partition samples into
environments (e.g., cell line, batch, donor, stimulus condition). The
core intuition: if X causes Y, the conditional distribution P(Y | X)
is stable across environments even when P(X) changes; non-causal
correlates lose their predictive link across environments.

Submodules:
    - invariant: InvariantPathwayPredictor (subset-search ICP)

Research use only. Not for clinical decision-making.
"""

from .invariant import (
    CausalParentReport,
    InvariantPathwayPredictor,
    invariance_pvalue,
)

__all__ = [
    "CausalParentReport",
    "InvariantPathwayPredictor",
    "invariance_pvalue",
]
