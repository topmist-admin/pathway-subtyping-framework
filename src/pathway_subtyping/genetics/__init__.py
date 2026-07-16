"""
Genetics-anchoring subpackage (added on the v0.8.0 line).

Supplies the anchoring cores behind validation **Gate 7 — Genetic Anchoring**.
Unlike every other validation gate, Gate 7 is a *positive-evidence* test: a pass
is evidence a subtype axis is genetically implicated in disease, because a
germline- or somatic-variant association cannot be manufactured by any postmortem
or technical confound (PMI, RIN, dissection, batch). It is the external,
causally-prior validator the internal transcriptomic gates (stability, nulls,
confound gate) structurally cannot provide.

Two modes, matched to the disease's genetic architecture:

- **Feature-level** (``ValidationGates.genetic_anchoring_gate``, ``gwas_enrichment``):
  "are the genes *defining* a subtype enriched for germline disease risk?" — a
  gene-set over-representation test against a background-matched universe. Needs
  only public risk-gene catalogues (zero data-use agreement). The
  psychiatry-appropriate mode; low power because germline risk is diffuse.

- **Somatic subject-level** (``ValidationGates.somatic_anchoring_gate``,
  ``somatic_anchoring``): "do the *tumors* in a subtype carry a somatic stratum
  (driver mutation / CNA / mutational signature) more than the others?" — a
  subject-level categorical association (the confound gate's statistic with
  inverted polarity). The cancer-appropriate mode; high power because somatic
  signal is strong per-tumor, and somatic calls are often available at a lower
  access tier than germline genotypes.

Still out of scope on this line: **germline** subject-level anchoring ("do subtype
donors carry more rare-variant burden / higher polygenic score?"), which needs
transcriptome+genotype on the same donors behind dbGaP/SFARI/Synapse data-use
agreements — a later, access-gated addition.

Import path::

    from pathway_subtyping.genetics import (
        feature_level_anchoring, hypergeometric_enrichment, EnrichmentResult,
        somatic_alignment, StratumAlignment,
    )

Research use only. Not for clinical decision-making.
"""
from .gwas_enrichment import (
    EnrichmentResult,
    feature_level_anchoring,
    hypergeometric_enrichment,
)
from .somatic_anchoring import (
    StratumAlignment,
    somatic_alignment,
)

__all__ = [
    "EnrichmentResult",
    "feature_level_anchoring",
    "hypergeometric_enrichment",
    "StratumAlignment",
    "somatic_alignment",
]
