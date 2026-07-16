"""
Genetics-anchoring subpackage (feature-level; added on the v0.8.0 line).

Supplies the enrichment core behind validation **Gate 7 — Genetic Anchoring**
(``ValidationGates.genetic_anchoring_gate``). Unlike every other validation gate,
Gate 7 is a *positive-evidence* test: a pass is evidence a subtype axis is
genetically implicated in disease, because a germline-variant enrichment cannot
be manufactured by any postmortem or technical confound (PMI, RIN, dissection,
batch). It is the external, causally-prior validator the internal transcriptomic
gates (stability, nulls, confound gate) structurally cannot provide.

Scope on this line is **feature-level** anchoring only — "are the genes defining a
subtype enriched for disease genetic risk?" — which needs only public risk-gene
catalogues and a background-matched gene universe (zero data-use agreement).
*Subject-level* anchoring ("do subtype donors carry more rare-variant burden /
higher polygenic score?") needs transcriptome+genotype on the same donors behind
dbGaP/SFARI/Synapse data-use agreements and is a later, access-gated addition.

Import path::

    from pathway_subtyping.genetics import (
        feature_level_anchoring,
        hypergeometric_enrichment,
        EnrichmentResult,
    )

Research use only. Not for clinical decision-making.
"""
from .gwas_enrichment import (
    EnrichmentResult,
    feature_level_anchoring,
    hypergeometric_enrichment,
)

__all__ = [
    "EnrichmentResult",
    "feature_level_anchoring",
    "hypergeometric_enrichment",
]
