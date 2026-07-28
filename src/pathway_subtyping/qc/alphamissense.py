"""
AlphaMissense-modulated cascade weighting (v0.6 F4).

PSF's cascade scoring (F1) currently treats every gene in a pathway as
either expressed or not. For carriers of rare missense variants, whether
that variant is pathogenic should modulate how much weight the gene
contributes to cascade completion. AlphaMissense (Cheng et al. 2023)
provides a proteome-wide pathogenicity prediction per missense variant
that we use to down-weight carrier genes proportionally to the
predicted pathogenicity.

Usage::

    from pathway_subtyping.qc.alphamissense import AlphaMissenseScorer

    scorer = AlphaMissenseScorer.from_table(pd.DataFrame({
        "variant_id": ["MYC:p.P72R", "TP53:p.R175H"],
        "am_score":  [0.12, 0.95],
    }))

    weights = scorer.weights_from_carriers(
        carriers=carriers_df,   # per (cell, gene, variant_id)
        cells=expr.index,
        genes=expr.columns,
    )
    cascade = CascadeAnalyzer(kg).analyze(expression, gene_weights=weights)

When ``carriers`` is empty or ``gene_weights=None`` is passed to
``CascadeAnalyzer.analyze``, downstream scoring is bit-identical to the
existing variant-naive behaviour.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


REQUIRED_SCORE_COLS = ("variant_id", "am_score")
REQUIRED_CARRIER_COLS = ("cell_id", "gene", "variant_id")


@dataclass
class AlphaMissenseScorer:
    """Look up AlphaMissense scores and produce per-cell per-gene weights.

    Attributes:
        table: Long-format DataFrame with columns ``variant_id`` and
            ``am_score``. Extra columns are preserved. A missing variant
            produces a :meth:`lookup` result of ``nan``.
    """

    table: pd.DataFrame

    # ------------------------------------------------------------ factories ---
    @classmethod
    def from_table(cls, table: pd.DataFrame) -> "AlphaMissenseScorer":
        return cls(table=table)

    @classmethod
    def empty(cls) -> "AlphaMissenseScorer":
        return cls(table=pd.DataFrame(columns=list(REQUIRED_SCORE_COLS)))

    # -------------------------------------------------------------- validate ---
    def __post_init__(self) -> None:
        for col in REQUIRED_SCORE_COLS:
            if col not in self.table.columns:
                raise ValueError(
                    f"AlphaMissenseScorer table missing required column '{col}'; "
                    f"got {list(self.table.columns)}"
                )
        # Normalise am_score into [0, 1] floats
        self.table = self.table.copy()
        self.table["am_score"] = pd.to_numeric(self.table["am_score"], errors="coerce").clip(
            0.0, 1.0
        )
        self._by_variant: Dict[str, float] = dict(
            zip(self.table["variant_id"], self.table["am_score"])
        )
        logger.info("[AlphaMissenseScorer] loaded %d variant scores", len(self._by_variant))

    # -------------------------------------------------------------- lookup ---
    def lookup(self, variant_id: str) -> float:
        """Return AM score in [0, 1] for a variant, or ``nan`` if unknown."""
        return float(self._by_variant.get(variant_id, np.nan))

    def lookup_many(self, variant_ids: Iterable[str]) -> np.ndarray:
        return np.asarray([self.lookup(v) for v in variant_ids], dtype=float)

    # --------------------------------------------------- carrier -> weights ---
    def weights_from_carriers(
        self,
        carriers: pd.DataFrame,
        cells: Iterable[Any],
        genes: Iterable[str],
        damage_floor: float = 0.0,
    ) -> pd.DataFrame:
        """Produce a (n_cells, n_genes) weight matrix for CascadeAnalyzer.

        Default weight is 1.0 (no carrier). For every (cell, gene) in the
        ``carriers`` table with a known AM score, the weight becomes
        ``max(1 - am_score, damage_floor)`` — so a highly pathogenic
        variant (AM close to 1) drives the weight near zero, and a benign
        variant (AM close to 0) leaves the weight near 1.

        If a single cell carries multiple pathogenic variants in the same
        gene, the most-damaging (lowest resulting weight) is retained.

        Args:
            carriers: DataFrame with columns ``cell_id``, ``gene``,
                ``variant_id``. Extra columns are ignored.
            cells: Row index of the resulting weight matrix. Typically
                ``expression.index``.
            genes: Column index of the weight matrix. Typically
                ``expression.columns``.
            damage_floor: Minimum possible weight per (cell, gene). The
                default of ``0.0`` means a fully-pathogenic variant zeroes
                out that gene's contribution; raise to, e.g., ``0.25``
                for a more conservative down-weighting.

        Returns:
            Weight DataFrame (float) indexed by ``cells`` with columns ``genes``.
        """
        for col in REQUIRED_CARRIER_COLS:
            if col not in carriers.columns:
                raise ValueError(
                    f"carriers DataFrame missing required column '{col}'; "
                    f"got {list(carriers.columns)}"
                )
        if not 0.0 <= damage_floor <= 1.0:
            raise ValueError("damage_floor must be in [0, 1]")

        cells = list(cells)
        genes = list(genes)
        weights = pd.DataFrame(1.0, index=cells, columns=genes)
        if carriers.empty:
            return weights

        cell_set = set(cells)
        gene_set = set(genes)
        for _, row in carriers.iterrows():
            cell_id = row["cell_id"]
            gene = row["gene"]
            am = self.lookup(row["variant_id"])
            if (cell_id not in cell_set) or (gene not in gene_set):
                continue
            if np.isnan(am):
                continue
            new_weight = max(1.0 - float(am), damage_floor)
            # keep the most damaging (lowest) weight across variants in this gene
            current = weights.at[cell_id, gene]
            weights.at[cell_id, gene] = min(current, new_weight)

        return weights

    # ---------------------------------------------------------- summary ---
    def summary(self) -> Dict[str, Any]:
        scores = self.table["am_score"].dropna()
        return {
            "n_variants": int(len(self._by_variant)),
            "mean_am_score": float(scores.mean()) if len(scores) else float("nan"),
            "n_likely_pathogenic": int((scores >= 0.564).sum()),
            "n_likely_benign": int((scores <= 0.34).sum()),
        }
