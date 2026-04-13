"""
Neurosymbolic Combiner.

Merges neural (GNN) predictions with symbolic rule conclusions using
configurable combination methods. Pure Python implementation — PyTorch
integration deferred to [gnn] extra.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


class CombinationMethod(Enum):
    """Methods for combining neural and symbolic scores."""

    WEIGHTED_SUM = "weighted_sum"
    MAX = "max"
    PRODUCT = "product"
    RULE_GUIDED = "rule_guided"


@dataclass
class CombinerConfig:
    """Configuration for the neurosymbolic combiner."""

    method: CombinationMethod = CombinationMethod.WEIGHTED_SUM
    neural_weight: float = 0.6
    symbolic_weight: float = 0.4
    use_rule_confidence: bool = True
    normalize_outputs: bool = True

    def __post_init__(self) -> None:
        assert 0 <= self.neural_weight <= 1
        assert 0 <= self.symbolic_weight <= 1


@dataclass
class CombinedOutput:
    """Output from combining neural and symbolic predictions."""

    combined_scores: Dict[str, float]
    neural_contribution: Dict[str, float]
    symbolic_contribution: Dict[str, float]
    method: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "method": self.method,
            "n_genes": len(self.combined_scores),
            "top_genes": dict(
                sorted(self.combined_scores.items(), key=lambda x: -x[1])[:10]
            ),
            "metadata": self.metadata,
        }


class NeurosymbolicCombiner:
    """Combines GNN gene scores with symbolic rule-derived scores.

    Usage::

        combiner = NeurosymbolicCombiner()
        result = combiner.combine(
            neural_scores={"CHD8": 0.9, "SHANK3": 0.7},
            symbolic_scores={"CHD8": 0.92, "NRXN1": 0.85},
        )
    """

    def __init__(self, config: Optional[CombinerConfig] = None):
        self.config = config or CombinerConfig()

    def combine(
        self,
        neural_scores: Dict[str, float],
        symbolic_scores: Dict[str, float],
    ) -> CombinedOutput:
        """Combine neural and symbolic gene scores.

        Args:
            neural_scores: Gene -> score from GNN predictions.
            symbolic_scores: Gene -> score from rule engine.

        Returns:
            CombinedOutput with merged scores.
        """
        all_genes = set(neural_scores.keys()) | set(symbolic_scores.keys())
        combined: Dict[str, float] = {}
        neural_contrib: Dict[str, float] = {}
        symbolic_contrib: Dict[str, float] = {}

        for gene in all_genes:
            n_score = neural_scores.get(gene, 0.0)
            s_score = symbolic_scores.get(gene, 0.0)

            if self.config.method == CombinationMethod.WEIGHTED_SUM:
                c = self.config.neural_weight * n_score + self.config.symbolic_weight * s_score
            elif self.config.method == CombinationMethod.MAX:
                c = max(n_score, s_score)
            elif self.config.method == CombinationMethod.PRODUCT:
                c = n_score * s_score
            elif self.config.method == CombinationMethod.RULE_GUIDED:
                # If symbolic score exists, boost neural; otherwise use neural alone
                if s_score > 0:
                    c = n_score * (1.0 + s_score)
                else:
                    c = n_score * 0.5
            else:
                c = (n_score + s_score) / 2.0

            combined[gene] = c
            neural_contrib[gene] = n_score
            symbolic_contrib[gene] = s_score

        if self.config.normalize_outputs and combined:
            max_val = max(combined.values())
            if max_val > 0:
                combined = {g: v / max_val for g, v in combined.items()}

        return CombinedOutput(
            combined_scores=combined,
            neural_contribution=neural_contrib,
            symbolic_contribution=symbolic_contrib,
            method=self.config.method.value,
        )


def create_symbolic_score_vector(
    fired_rules: List[Any],
    gene_list: List[str],
    use_confidence: bool = True,
) -> np.ndarray:
    """Convert fired rules to a gene score vector.

    Args:
        fired_rules: List of FiredRule instances.
        gene_list: Ordered list of gene IDs.
        use_confidence: Weight by rule confidence.

    Returns:
        1D array of per-gene symbolic scores.
    """
    scores = np.zeros(len(gene_list))
    gene_idx = {g: i for i, g in enumerate(gene_list)}

    for fr in fired_rules:
        if fr.gene and fr.gene in gene_idx:
            idx = gene_idx[fr.gene]
            score = fr.confidence if use_confidence else 1.0
            scores[idx] = max(scores[idx], score)

    return scores
