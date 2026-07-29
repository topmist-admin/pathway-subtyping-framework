"""
Pathway Crosstalk & Interference Detection (F9).

Intended purpose: detect resource competition at shared pathway nodes. When
multiple pathways share intermediate transducers, one pathway can starve another.

⚠️ EXPERIMENTAL / UNSTABLE — THE SCORE DOES NOT YET MEASURE WHAT IT CLAIMS.
=========================================================================
``competition_score`` is **not a function of the shared nodes**, despite naming
and documentation that say otherwise. Do not interpret current values, and do
not use ``dominant`` or the PASS/FAIL summary to make a release decision.

What the code actually does. ``_compute_competition`` builds ``genes_a`` and
``genes_b`` by *subtracting* the shared set from each pathway, then returns
``a_mean - b_mean``. That compares the two pathways' **exclusive** machinery —
roughly "which pathway is more active overall" — and the shared genes never
enter the arithmetic. They are used only as a gate: if no shared gene is present
the function short-circuits to 0.0.

Demonstrated consequences (2026-07-28):
  * Varying shared-node expression across four orders of magnitude (0.1 -> 500)
    leaves ``competition_score`` bit-identical and ``dominant`` unchanged.
  * Holding shared nodes perfectly balanced and varying only the exclusive genes
    swings the score from +8.0 to -8.0 and flips ``dominant``.
  * ``dominant`` is assigned by ``score > 0`` with no dead zone, so sampling
    noise between two identical distributions still names a deterministic winner.
  * ``n_significant`` thresholds ``abs(score)`` at 0.3 by default; on log2-scale
    expression a >0.3 difference in mean pathway activity is routine, so the
    PASS/FAIL verdict would fail most real batches for an unrelated reason.

This is not a regression. The module shipped this way in the v0.5 QC layer; a
``shared_mean`` local was computed and never used from the first commit, and was
removed in the 2026-07-28 lint cleanup as dead code.

The replacement is specified in ``docs/roadmap-f9-competition-model.md``
(**True Competition / Starvation Interaction Model**, targeted at v0.9.0): a
partial correlation ``rho_AB.S`` as a screen, confirmed by a moderation model
``E_B ~ E_A + S + E_A:S`` in which starvation predicts a negative A->B slope that
attenuates as the shared transducer becomes abundant. That is an interaction
test and a structurally different computation, not a patch to the line below.

Until that is settled this class emits a ``FutureWarning`` on construction and is
withheld from ``pathway_subtyping.qc.__all__``. It remains importable so existing
code does not break. The unrelated
``KnowledgeGraph.get_pathway_crosstalk`` / ``get_shared_genes`` topology helpers
are **not** affected by this notice. Note also that "F9" is overloaded in this
project: this is F9 of the *molecular QC layer*, not F9 of the *v0.6 release*
(``qc.offtarget_sequence``, Evo 2 off-target scoring), which is unaffected.

Research use only. Not for clinical decision-making.
"""

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class InterferenceEdge:
    """Interference relationship between two pathways."""

    pathway_a: str
    pathway_b: str
    competition_score: float
    shared_genes: List[str] = field(default_factory=list)
    # ⚠️ EXPERIMENTAL: labelled "which pathway dominates shared nodes", but it is
    # currently derived from the pathways' EXCLUSIVE genes and is insensitive to
    # shared-node expression. See the module docstring. Do not interpret.
    dominant: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway_a": self.pathway_a,
            "pathway_b": self.pathway_b,
            "competition_score": round(self.competition_score, 4),
            "n_shared_genes": len(self.shared_genes),
            "dominant": self.dominant,
        }


@dataclass
class CrosstalkResult:
    """Result of crosstalk detection."""

    interference_edges: List[InterferenceEdge]
    n_significant: int = 0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_significant == 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_significant": self.n_significant,
            "summary": self.summary,
            "edges": [e.to_dict() for e in self.interference_edges],
        }


class CrosstalkDetector:
    """Detects pathway crosstalk via shared node competition.

    Usage::

        detector = CrosstalkDetector(kg)
        result = detector.detect(expression, pathways=["PW_A", "PW_B"])
    """

    def __init__(
        self,
        kg: Any,
        competition_threshold: float = 0.3,
        seed: Optional[int] = None,
    ):
        # FutureWarning, not DeprecationWarning: DeprecationWarning is hidden by
        # default outside __main__, and this needs to reach anyone who runs the
        # detector. The category is also honest -- the score's definition WILL
        # change once the intended formula is settled.
        warnings.warn(
            "CrosstalkDetector (molecular-QC F9) is EXPERIMENTAL and its competition_score does "
            "not currently measure competition at shared nodes: the score is "
            "computed from each pathway's EXCLUSIVE genes and is provably "
            "insensitive to shared-node expression. Do not interpret "
            "competition_score, dominant, or the PASS/FAIL summary. The scoring "
            "definition will change. See the module docstring in "
            "pathway_subtyping/qc/crosstalk.py for the demonstration and the "
            "candidate fixes.",
            FutureWarning,
            stacklevel=2,
        )
        self.kg = kg
        self.competition_threshold = competition_threshold

    def detect(
        self,
        expression: pd.DataFrame,
        pathways: List[str],
    ) -> CrosstalkResult:
        """Detect crosstalk between all pathway pairs.

        Args:
            expression: Expression matrix (cells x genes).
            pathways: List of pathway node IDs to check.

        Returns:
            CrosstalkResult with interference edges.
        """
        edges: List[InterferenceEdge] = []
        n_sig = 0

        for i in range(len(pathways)):
            for j in range(i + 1, len(pathways)):
                pw_a, pw_b = pathways[i], pathways[j]
                shared = self.kg.get_shared_genes(pw_a, pw_b)
                if not shared:
                    continue

                score, dominant = self._compute_competition(expression, pw_a, pw_b, shared)

                edge = InterferenceEdge(
                    pathway_a=pw_a,
                    pathway_b=pw_b,
                    competition_score=score,
                    shared_genes=shared,
                    dominant=dominant,
                )
                edges.append(edge)

                if abs(score) > self.competition_threshold:
                    n_sig += 1

        summary = (
            f"Crosstalk: {len(edges)} pathway pairs with shared nodes, "
            f"{n_sig} significant. "
            f"{'PASS' if n_sig == 0 else 'FAIL'}."
        )

        logger.info("[QC Crosstalk] %s", summary)

        return CrosstalkResult(
            interference_edges=edges,
            n_significant=n_sig,
            summary=summary,
        )

    def _compute_competition(
        self,
        expression: pd.DataFrame,
        pw_a: str,
        pw_b: str,
        shared: List[str],
    ) -> Tuple[float, str]:
        """Compute the (currently misnamed) competition score for a pathway pair.

        ⚠️ Despite the name, this does NOT measure competition at shared nodes.
        ``shared`` is subtracted from both gene sets below and then never enters
        the arithmetic; it acts only as a presence gate. The returned value is
        ``mean(A-exclusive) - mean(B-exclusive)``. See the module docstring for
        the demonstration and the candidate corrected formulas.
        """
        genes_a = set(self.kg.get_pathway_genes(pw_a)) - set(shared)
        genes_b = set(self.kg.get_pathway_genes(pw_b)) - set(shared)

        shared_present = [g for g in shared if g in expression.columns]
        a_present = [g for g in genes_a if g in expression.columns]
        b_present = [g for g in genes_b if g in expression.columns]

        if not shared_present or (not a_present and not b_present):
            return 0.0, ""

        a_mean = expression[a_present].values.mean() if a_present else 0.0
        b_mean = expression[b_present].values.mean() if b_present else 0.0

        # ⚠️ The comment this replaces claimed "positive = shared nodes favor A".
        # That is not what is computed: a_mean/b_mean are means over genes
        # EXCLUSIVE to each pathway, so no shared node influences the result.
        score = a_mean - b_mean
        dominant = pw_a if score > 0 else pw_b if score < 0 else ""
        return float(score), dominant
