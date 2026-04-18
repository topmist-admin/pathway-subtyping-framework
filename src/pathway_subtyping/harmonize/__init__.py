"""
Cross-Platform Harmonization Layer for the Pathway Subtyping Framework.

Addresses the v0.6 roadmap Phase 1 F2 limitation: PSF currently assumes
pathway scores from 10x Chromium, Smart-seq2, bulk RNA-seq, and spatial
transcriptomics are comparable after standard preprocessing. Empirically
this holds within a platform but cross-platform correlations are lower
than within-platform, leading to score drift when users combine public
datasets.

Submodules:
    - uce: UCE (Universal Cell Embeddings) embedding wrapper
    - align: Platform-specific -> reference alignment of pathway scores
    - confidence: Per-cell harmonization confidence scoring
    - benchmark: Cross-platform evaluation suite

The UCE backend is opt-in via ``pip install pathway-subtyping[harmonize]``.
When UCE is not installed, a deterministic PCA-based fallback embedder is
available for testing and development; it is not a substitute for the real
UCE model but preserves the public API so downstream code runs uniformly.

Research use only. Not for clinical decision-making.
"""

from .align import AlignmentResult, CrossPlatformAligner
from .benchmark import CrossPlatformBenchmark, simulate_platform_distortion
from .confidence import HarmonizationReport
from .uce import FallbackEmbedder, UCEEmbedder

__all__ = [
    "AlignmentResult",
    "CrossPlatformAligner",
    "CrossPlatformBenchmark",
    "FallbackEmbedder",
    "HarmonizationReport",
    "UCEEmbedder",
    "simulate_platform_distortion",
]
