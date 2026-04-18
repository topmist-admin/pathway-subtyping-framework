"""
Gene-set curation utilities for the Pathway Subtyping Framework.

v0.6 Phase 2 F7: regulatory-evidence-based gene-set expansion. Users
curating custom gene sets can consult this module for candidate
additions informed by a DNA-sequence-predicted regulatory model
(Borzoi / Enformer in production, co-expression fallback for testing
and offline runs).

Submodules:
    - regulatory: ``RegulatoryGeneSetExpander``

Install: pip install pathway-subtyping[genesets]

Research use only. Not for clinical decision-making.
"""

from .regulatory import (
    BorzoiBackend,
    CoexpressionBackend,
    ExpansionCandidate,
    ExpansionResult,
    RegulatoryBackend,
    RegulatoryGeneSetExpander,
)

__all__ = [
    "BorzoiBackend",
    "CoexpressionBackend",
    "ExpansionCandidate",
    "ExpansionResult",
    "RegulatoryBackend",
    "RegulatoryGeneSetExpander",
]
