"""
Multi-omics fusion for the Pathway Subtyping Framework.

v0.6 Phase 3 F10: PSF v0.5 is RNA-centric — pathway activity inferred
from transcript abundance alone misses post-transcriptional and
post-translational regulation. This subpackage adds ATAC-seq and
proteomics pathway scoring that can be fused with the standard
RNA-derived scores, plus an explicit discordance flag for pathways
where RNA and protein disagree.

Submodules:
    - atac: ATAC-seq chromatin-accessibility pathway scoring.
    - proteomics: Proteomics pathway scoring.
    - fusion: Weighted fusion of RNA + ATAC + protein pathway scores.
    - discordance: RNA-vs-protein disagreement detection.

All four layers assume pathway-level inputs with matching sample indices;
raw-data ingestion and gene/protein/peak mapping are the caller's
responsibility (bring your own scanpy / muon pipeline).

Research use only. Not for clinical decision-making.
"""

from .atac import ATACScorer, score_atac_pathways
from .discordance import DiscordanceReport, flag_discordant_pathways
from .fusion import FusionResult, FusionWeights, MultiOmicsFusion
from .proteomics import ProteomicsScorer, score_proteomics_pathways

__all__ = [
    "ATACScorer",
    "DiscordanceReport",
    "FusionResult",
    "FusionWeights",
    "MultiOmicsFusion",
    "ProteomicsScorer",
    "flag_discordant_pathways",
    "score_atac_pathways",
    "score_proteomics_pathways",
]
