"""
Active-learning sample selection for the Pathway Subtyping Framework.

v0.6 Phase 3 F12: given a pool of unlabelled samples and a fixed label
budget, select the samples whose labels would most reduce PSF subtype
uncertainty. Built on the Phase 1 F1 uncertainty layer, which provides
the calibrated per-sample confidence that drives the uncertainty
strategy.

Submodules:
    - selector: ``ActiveSampleSelector`` with three strategies —
      ``uncertainty``, ``diversity``, ``hybrid``.

Research use only. Not for clinical decision-making.
"""

from .selector import (
    ActiveSampleSelector,
    SelectionResult,
    SelectionStrategy,
)

__all__ = [
    "ActiveSampleSelector",
    "SelectionResult",
    "SelectionStrategy",
]
