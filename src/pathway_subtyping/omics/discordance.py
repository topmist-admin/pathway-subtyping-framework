"""
RNA vs protein discordance flagging.

Pathways where RNA and protein-level scores disagree carry information:
post-transcriptional regulation, protein-level feedback, secretion
losses, degradation rate differences. The roadmap notes this explicitly
— flag discordant pathways, don't smooth them away.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)


@dataclass
class DiscordanceReport:
    """Per-pathway RNA-vs-protein agreement diagnostic.

    Attributes:
        pathway_stats: DataFrame indexed by pathway with columns
            ``rho``, ``mean_absolute_diff``, ``discordant``.
        threshold_rho: The rho floor below which a pathway is flagged.
        threshold_abs_diff: The mean-abs-diff ceiling above which a
            pathway is flagged (used in combination with rho).
    """

    pathway_stats: pd.DataFrame
    threshold_rho: float
    threshold_abs_diff: float

    @property
    def discordant_pathways(self) -> List[str]:
        return list(self.pathway_stats.index[self.pathway_stats["discordant"]])

    @property
    def concordant_pathways(self) -> List[str]:
        return list(self.pathway_stats.index[~self.pathway_stats["discordant"]])

    def summary(self) -> str:
        n = len(self.pathway_stats)
        n_disc = int(self.pathway_stats["discordant"].sum())
        return (
            f"DiscordanceReport(n_pathways={n}, "
            f"discordant={n_disc} ({n_disc / max(n, 1):.0%}), "
            f"rho<{self.threshold_rho} or |diff|>{self.threshold_abs_diff})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_pathways": int(len(self.pathway_stats)),
            "n_discordant": int(self.pathway_stats["discordant"].sum()),
            "threshold_rho": float(self.threshold_rho),
            "threshold_abs_diff": float(self.threshold_abs_diff),
            "discordant_pathways": self.discordant_pathways,
        }


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

def flag_discordant_pathways(
    rna: pd.DataFrame,
    protein: pd.DataFrame,
    threshold_rho: float = 0.3,
    threshold_abs_diff: float = 1.0,
) -> DiscordanceReport:
    """Flag pathways whose RNA and protein scores disagree.

    Discordance criterion: Spearman rho across samples is below
    ``threshold_rho`` **and** mean absolute difference exceeds
    ``threshold_abs_diff``. Both conditions must hold to flag — a
    pathway with low rho but tiny absolute differences is noisy, not
    discordant, and vice versa.

    Args:
        rna: Samples × pathways RNA-derived score frame (typically
            z-normalised by the upstream scorer).
        protein: Samples × pathways protein-derived score frame.
            Index and columns must intersect non-trivially with ``rna``.
        threshold_rho: Spearman rho cutoff. Lower values → more
            pathways flagged.
        threshold_abs_diff: Mean absolute difference cutoff.

    Returns:
        DiscordanceReport with per-pathway stats and the discordant flag.
    """
    shared_samples = rna.index.intersection(protein.index)
    shared_pathways = rna.columns.intersection(protein.columns)
    if len(shared_samples) == 0 or len(shared_pathways) == 0:
        raise ValueError(
            "rna and protein must share at least one sample and pathway"
        )

    rna_sub = rna.loc[shared_samples, shared_pathways]
    prot_sub = protein.loc[shared_samples, shared_pathways]

    rows: List[Dict[str, Any]] = []
    for pw in shared_pathways:
        rna_v = rna_sub[pw].to_numpy()
        prot_v = prot_sub[pw].to_numpy()
        if np.std(rna_v) == 0 or np.std(prot_v) == 0:
            rho = 0.0
        else:
            rho_value, _ = spearmanr(rna_v, prot_v)
            rho = float(rho_value) if not np.isnan(rho_value) else 0.0
        mad = float(np.abs(rna_v - prot_v).mean())
        discordant = (rho < threshold_rho) and (mad > threshold_abs_diff)
        rows.append({
            "pathway": pw,
            "rho": rho,
            "mean_absolute_diff": mad,
            "discordant": bool(discordant),
        })

    stats = pd.DataFrame.from_records(rows).set_index("pathway")
    logger.info(
        "[Discordance] n_pathways=%d n_discordant=%d",
        len(stats), int(stats["discordant"].sum()),
    )
    return DiscordanceReport(
        pathway_stats=stats,
        threshold_rho=float(threshold_rho),
        threshold_abs_diff=float(threshold_abs_diff),
    )
