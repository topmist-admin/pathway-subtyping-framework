"""
In-silico gene perturbation for the Pathway Subtyping Framework.

v0.6 Phase 2 F5: given an expression matrix and a gene of interest
(knockout or over-expression), return a perturbed molecular-state vector
(MSV) with calibrated conformal prediction intervals.

Pipeline:
    1. `GeneformerPerturber.perturb(expr, gene, mode)` produces a
       perturbed per-cell embedding. Real backend wraps the published
       Geneformer foundation model (Theodoris et al. 2023, Apache-2.0)
       via its official loader; a `FallbackPerturber` provides a
       deterministic substitute that preserves the API for testing.
    2. `MSVFromEmbedding.transform(embedding)` maps Geneformer's 256-dim
       embeddings to PSF's per-pathway MSV scores via a small learned
       head fit on reference data.
    3. `ConformalPathwayPredictor` (F1) wraps the chain for calibrated
       prediction intervals.
    4. `PerturbationScreen.run(expr, gene_list)` iterates over a panel
       of genes and ranks them by delta-MSV impact.

Install: pip install pathway-subtyping[perturb]

Research use only. Not for clinical decision-making.
"""

from .geneformer_wrapper import (
    FallbackPerturber,
    GeneformerBackend,
    GeneformerPerturber,
    OfficialBackend,
    PerturbationMode,
    PerturbationResult,
)
from .msv_from_embedding import MSVFromEmbedding
from .report import PerturbationReport
from .screen import PerturbationScreen, ScreenResult

__all__ = [
    "FallbackPerturber",
    "GeneformerBackend",
    "GeneformerPerturber",
    "MSVFromEmbedding",
    "OfficialBackend",
    "PerturbationMode",
    "PerturbationResult",
    "PerturbationReport",
    "PerturbationScreen",
    "ScreenResult",
]
