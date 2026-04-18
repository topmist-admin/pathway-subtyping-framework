"""
Cell-embedding wrappers for the Pathway Subtyping Framework.

v0.6 Phase 2 F6: thin, backend-agnostic wrappers around single-cell
foundation models whose outputs complement PSF pathway scores. The
first supported backend is scGPT (Cui et al. 2024); additional models
plug in behind the same ``Embedder`` interface.

Submodules:
    - base: Abstract ``Embedder`` interface shared by all backends.
    - scgpt: ``scGPTEmbedder`` — thin scGPT wrapper with opt-in
      model loading and a deterministic PCA-based fallback for tests.
    - cache: Content-hashed embedding cache; safe across package
      upgrades because the hash covers the model checkpoint version,
      tokenizer version, and the input expression bytes.

Install: pip install pathway-subtyping[embed]

Research use only. Not for clinical decision-making.
"""

from .base import Embedder, EmbeddingResult
from .cache import EmbeddingCache, cache_key_for
from .scgpt import FallbackSCGPTEmbedder, OfficialSCGPTBackend, scGPTEmbedder

__all__ = [
    "Embedder",
    "EmbeddingCache",
    "EmbeddingResult",
    "FallbackSCGPTEmbedder",
    "OfficialSCGPTBackend",
    "cache_key_for",
    "scGPTEmbedder",
]
