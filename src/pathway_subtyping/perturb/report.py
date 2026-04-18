"""
Perturbation impact reports.

Wraps a ``ScreenResult`` with human-oriented diagnostics:
    - top-K gene ranking (from the screen)
    - directional-expectation check against a user-supplied signature
      (e.g., "knocking out MYC should *decrease* proliferation pathway
      scores") — the roadmap acceptance test for master-regulator
      perturbations
    - optional conformal-interval integration for uncertainty-aware
      impact ranking (uses F1 ``ConformalPathwayPredictor``)

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional

import numpy as np
import pandas as pd

from .screen import ScreenResult

logger = logging.getLogger(__name__)


DirectionalSignature = Mapping[str, Mapping[str, int]]
"""Mapping gene -> {pathway -> expected sign of delta (+1, -1)}.

Used by :meth:`PerturbationReport.check_directional_signature` to test
whether a knockdown/overexpression screen produces the biologically
expected direction of change on master-regulator pathways.
"""


@dataclass
class PerturbationReport:
    """Diagnostic summary over a ``ScreenResult``.

    Usage::

        report = PerturbationReport.from_screen(screen_result)
        print(report.summary())
        hits = report.check_directional_signature({
            "MYC": {"HALLMARK_MYC_TARGETS_V1": -1},
        })
    """

    screen: ScreenResult

    # ------------------------------------------------------- factories ---
    @classmethod
    def from_screen(cls, screen: ScreenResult) -> "PerturbationReport":
        return cls(screen=screen)

    # ------------------------------------------------------- summaries ---
    @property
    def gene_panel(self) -> List[str]:
        return list(self.screen.gene_panel)

    def top_k(self, k: int = 5) -> pd.DataFrame:
        return self.screen.rank(k)

    def dominant_pathway_per_gene(self) -> pd.Series:
        """For each gene, return the pathway with the largest |delta|."""
        abs_delta = self.screen.delta_msv_by_gene.abs()
        return abs_delta.idxmax(axis=1).rename("dominant_pathway")

    def summary(self) -> str:
        top_gene = self.screen.l2_by_gene.idxmax() if len(self.screen.l2_by_gene) else None
        return (
            f"PerturbationReport(n_genes={len(self.screen.gene_panel)}, "
            f"mode={self.screen.mode.value}, "
            f"top_gene={top_gene!r}, "
            f"max_l2={float(self.screen.l2_by_gene.max()):.4f})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_genes": len(self.screen.gene_panel),
            "mode": self.screen.mode.value,
            "top5": self.top_k(5).to_dict()["l2_delta_msv"],
            "dominant_pathway_per_gene": self.dominant_pathway_per_gene().to_dict(),
        }

    # -------------------------------------- directional signature check ---
    def check_directional_signature(
        self,
        signature: DirectionalSignature,
        min_magnitude: float = 0.0,
    ) -> pd.DataFrame:
        """Check expected sign of delta-MSV per (gene, pathway) assertion.

        Returns a DataFrame indexed by (gene, pathway) with columns
        ``observed_sign``, ``expected_sign``, ``delta_value``, ``passed``.

        ``passed`` is True when ``observed_sign == expected_sign`` and
        ``|delta_value| >= min_magnitude``. A gene absent from the screen
        or a pathway absent from the MSV frame is reported with
        ``observed_sign=nan`` and ``passed=False``.
        """
        records: List[Dict[str, Any]] = []
        delta = self.screen.delta_msv_by_gene

        for gene, pathway_map in signature.items():
            for pathway, expected in pathway_map.items():
                if gene not in delta.index or pathway not in delta.columns:
                    records.append({
                        "gene": gene, "pathway": pathway,
                        "expected_sign": int(expected),
                        "observed_sign": float("nan"),
                        "delta_value": float("nan"),
                        "passed": False,
                        "reason": "missing gene or pathway",
                    })
                    continue
                value = float(delta.loc[gene, pathway])
                observed = int(np.sign(value))
                passes = (
                    observed == int(np.sign(expected))
                    and abs(value) >= float(min_magnitude)
                )
                records.append({
                    "gene": gene,
                    "pathway": pathway,
                    "expected_sign": int(np.sign(expected)),
                    "observed_sign": observed,
                    "delta_value": value,
                    "passed": bool(passes),
                    "reason": "" if passes else (
                        "sign mismatch" if observed != int(np.sign(expected))
                        else "below min_magnitude"
                    ),
                })

        df = pd.DataFrame.from_records(records)
        if not df.empty:
            df = df.set_index(["gene", "pathway"])
        return df

    # ------------------------------------------------------------ plotting ---
    def plot_top_k(self, k: int = 10, ax: Optional[Any] = None):
        """Horizontal bar chart of top-K gene impact."""
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "matplotlib required for PerturbationReport.plot_top_k(); "
                "install pathway-subtyping[viz]"
            ) from exc
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, max(3, 0.35 * k)))
        else:
            fig = ax.figure
        top = self.top_k(k).iloc[::-1]  # largest at top
        ax.barh(top.index, top["l2_delta_msv"])
        ax.set_xlabel("|delta-MSV| L2")
        ax.set_ylabel("gene")
        ax.set_title(f"Top {k} perturbation impact")
        ax.grid(True, alpha=0.3, axis="x")
        return fig
