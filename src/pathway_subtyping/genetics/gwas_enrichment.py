"""
Feature-level genetic anchoring — enrichment core (Gate 7 support).

Tests whether a gene set (e.g. a subtype's defining genes) is over-represented
for disease genetic-risk genes, by a hypergeometric over-representation test
against an explicit gene universe. This is the reusable primitive behind
``ValidationGates.genetic_anchoring_gate`` (validation.py), factored out so it can
also be called directly for calibration / positive-control work.

Why this belongs in the validation battery
-------------------------------------------
Every other gate is an *internal* transcriptomic check — stability, nulls, the
confound gate — and can only certify that a partition is *reproducible*, never
that it is *real* or *causal*. A germline-variant enrichment is the exception: it
cannot be manufactured by any postmortem or technical confound (PMI, RIN,
dissection, batch, sequencing platform), because the variants are fixed at
conception, upstream of every one of those. So a subtype whose defining genes
carry disease genetic risk is genetically implicated — necessary, not sufficient,
but confound-immune, and computable from public risk-gene catalogues with zero
data-use agreement.

The null must be background-matched
-----------------------------------
The single make-or-break methodological point: the universe used as the
background must be matched to the tissue, not genome-wide. Brain-expressed genes
are already enriched for brain-disease genetic risk, so a region-identity /
brain-identity subtype could show spurious enrichment against a genome-wide
protein-coding background. The discriminating test uses a background-matched
universe (e.g. brain-expressed genes). ``feature_level_anchoring`` reports both a
matched null and an optional genome-wide reference so the difference is visible;
only the matched null decides anchoring.

Research use only. Not for clinical decision-making.
"""

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional

import numpy as np
from scipy.stats import hypergeom


@dataclass
class EnrichmentResult:
    """Result of one hypergeometric over-representation test.

    Attributes:
        label: Name of the gene set tested (e.g. a subtype label).
        null: Name of the background universe used ("background-matched" or
            "genome-wide").
        universe_n: Size of the background universe after intersection (N).
        risk_in_universe: Risk genes present in the universe (K).
        testset_n: Test-set genes present in the universe (n).
        risk_hits: Test-set genes that are also risk genes (k).
        expected: Risk hits expected under the null (n * K / N).
        fold: Observed / expected fold enrichment ((k / n) / (K / N)).
        p_value: Hypergeometric survival-function p (P[X >= k]); NaN if the
            test is undefined (empty test set or no risk genes in the universe).
    """

    label: str
    null: str
    universe_n: int
    risk_in_universe: int
    testset_n: int
    risk_hits: int
    expected: float
    fold: float
    p_value: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "label": self.label,
            "null": self.null,
            "universe_n": self.universe_n,
            "risk_in_universe": self.risk_in_universe,
            "testset_n": self.testset_n,
            "risk_hits": self.risk_hits,
            "expected": round(self.expected, 2) if np.isfinite(self.expected) else None,
            "fold": round(self.fold, 4) if np.isfinite(self.fold) else None,
            "p_value": self.p_value,
        }


def hypergeometric_enrichment(
    test_set: Iterable[str],
    risk_genes: Iterable[str],
    universe: Iterable[str],
    label: str = "",
    null: str = "",
) -> EnrichmentResult:
    """Hypergeometric over-representation of ``risk_genes`` in ``test_set``.

    Both the test set and the risk-gene set are intersected with ``universe``
    before testing, so the universe fully defines the background. Gene
    identifiers must be in the same namespace (e.g. all Ensembl gene IDs, with
    version suffixes already stripped).

    Args:
        test_set: Genes defining the axis under test (e.g. a subtype's genes).
        risk_genes: Disease genetic-risk reference genes.
        universe: Background universe of testable genes (N).
        label: Optional name for the test set (carried into the result).
        null: Optional name for the background (carried into the result).

    Returns:
        EnrichmentResult. ``p_value`` is P[X >= k] under the hypergeometric null;
        NaN when the test is undefined (n == 0 or K == 0).
    """
    universe_set = set(universe)
    test = set(test_set) & universe_set
    risk = set(risk_genes) & universe_set

    n = len(test)
    k = len(test & risk)
    big_k = len(risk)
    big_n = len(universe_set)

    if n == 0 or big_k == 0 or big_n == 0:
        return EnrichmentResult(
            label=label,
            null=null,
            universe_n=big_n,
            risk_in_universe=big_k,
            testset_n=n,
            risk_hits=k,
            expected=float("nan"),
            fold=float("nan"),
            p_value=float("nan"),
        )

    expected = n * big_k / big_n
    fold = (k / n) / (big_k / big_n)
    # Survival function is P[X > k-1] = P[X >= k].
    p_value = float(hypergeom.sf(k - 1, big_n, big_k, n))

    return EnrichmentResult(
        label=label,
        null=null,
        universe_n=big_n,
        risk_in_universe=big_k,
        testset_n=n,
        risk_hits=k,
        expected=float(expected),
        fold=float(fold),
        p_value=p_value,
    )


def feature_level_anchoring(
    gene_sets: Dict[Any, Iterable[str]],
    risk_genes: Iterable[str],
    universe: Iterable[str],
    reference_universe: Optional[Iterable[str]] = None,
) -> Dict[str, EnrichmentResult]:
    """Test every gene set for risk-gene over-representation against a universe.

    Convenience wrapper over :func:`hypergeometric_enrichment` that runs each set
    against the (background-matched) ``universe`` and, optionally, against a
    ``reference_universe`` (e.g. genome-wide protein-coding) reported for
    contrast. This is the primitive used both by the Gate-7 wrapper and by the
    Voineagu-style calibration (neuronal set should enrich for autism risk, glial
    set should not, under the brain-expressed-matched null).

    Args:
        gene_sets: Mapping of set label -> gene identifiers.
        risk_genes: Disease genetic-risk reference genes.
        universe: Background-matched universe (the discriminating null).
        reference_universe: Optional genome-wide universe reported alongside.

    Returns:
        Dict keyed ``"{label}|{null}"`` -> EnrichmentResult. When
        ``reference_universe`` is given, each set contributes two entries
        (``background-matched`` and ``genome-wide``).
    """
    results: Dict[str, EnrichmentResult] = {}
    for label, genes in gene_sets.items():
        matched = hypergeometric_enrichment(
            genes, risk_genes, universe, label=str(label), null="background-matched"
        )
        results[f"{label}|background-matched"] = matched
        if reference_universe is not None:
            ref = hypergeometric_enrichment(
                genes, risk_genes, reference_universe, label=str(label), null="genome-wide"
            )
            results[f"{label}|genome-wide"] = ref
    return results
