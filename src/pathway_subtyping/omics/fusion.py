"""
Multi-omics pathway-score fusion.

Takes per-modality pathway-score matrices (RNA, ATAC, proteomics), each
a sample-by-pathway DataFrame, and returns a fused per-pathway score.
Fusion weights can be supplied directly, or learned from a paired
reference cohort by minimising classification error against known cell-
type labels (the roadmap acceptance target is a >= 3% classification
accuracy improvement over RNA-only).

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


MODALITIES = ("rna", "atac", "protein")


# --------------------------------------------------------------------------- #
# Types
# --------------------------------------------------------------------------- #


@dataclass
class FusionWeights:
    """Per-modality scalar weights. Auto-normalises to sum to 1.0.

    Zero weights are permitted — useful when a modality is unavailable.
    """

    rna: float = 1.0
    atac: float = 0.0
    protein: float = 0.0

    def __post_init__(self) -> None:
        for name in MODALITIES:
            val = getattr(self, name)
            if val < 0:
                raise ValueError(f"fusion weight {name!r}={val} must be >= 0")
        s = self.rna + self.atac + self.protein
        if s <= 0:
            raise ValueError("at least one fusion weight must be positive")
        self.rna = self.rna / s
        self.atac = self.atac / s
        self.protein = self.protein / s

    def as_dict(self) -> Dict[str, float]:
        return {"rna": self.rna, "atac": self.atac, "protein": self.protein}


@dataclass
class FusionResult:
    """Fused pathway-score matrix plus provenance.

    Attributes:
        fused: Samples × pathways fused score matrix.
        weights: FusionWeights used to produce ``fused``.
        per_modality: Dict mapping modality name → the (aligned)
            per-modality score frame that was fed into the fusion.
        union_pathways: Pathways present in the fused frame (intersection
            of all supplied modalities).
    """

    fused: pd.DataFrame
    weights: FusionWeights
    per_modality: Dict[str, pd.DataFrame] = field(default_factory=dict)
    union_pathways: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "weights": self.weights.as_dict(),
            "n_samples": int(self.fused.shape[0]),
            "n_pathways": int(self.fused.shape[1]),
            "modalities_used": sorted(self.per_modality.keys()),
        }


# --------------------------------------------------------------------------- #
# Fusion class
# --------------------------------------------------------------------------- #


class MultiOmicsFusion:
    """Combine per-modality pathway-score matrices into a single MSV.

    Usage::

        fusion = MultiOmicsFusion()
        result = fusion.fuse(
            rna=rna_scores,
            atac=atac_scores,
            protein=protein_scores,
            weights=FusionWeights(rna=0.6, atac=0.2, protein=0.2),
        )
        # result.fused is (n_samples, n_pathways)

    Or learn weights from labeled reference data::

        weights = fusion.learn_weights(
            rna=rna_scores, atac=atac_scores, protein=protein_scores,
            labels=cell_type_labels,
        )
        result = fusion.fuse(rna=rna_scores, atac=atac_scores, protein=protein_scores,
                             weights=weights)
    """

    # ------------------------------------------------------------ fuse ---
    def fuse(
        self,
        rna: Optional[pd.DataFrame] = None,
        atac: Optional[pd.DataFrame] = None,
        protein: Optional[pd.DataFrame] = None,
        weights: Optional[FusionWeights] = None,
    ) -> FusionResult:
        modalities = {"rna": rna, "atac": atac, "protein": protein}
        supplied = {k: v for k, v in modalities.items() if v is not None}
        if not supplied:
            raise ValueError("at least one modality must be supplied")

        # Align on the intersection of sample indices and pathway columns
        sample_idx: Optional[pd.Index] = None
        pathway_idx: Optional[pd.Index] = None
        for df in supplied.values():
            sample_idx = df.index if sample_idx is None else sample_idx.intersection(df.index)
            pathway_idx = (
                df.columns if pathway_idx is None else pathway_idx.intersection(df.columns)
            )
        if sample_idx is None or pathway_idx is None:
            raise ValueError("no shared samples or pathways across modalities")
        if len(sample_idx) == 0:
            raise ValueError("no shared samples across supplied modalities")
        if len(pathway_idx) == 0:
            raise ValueError("no shared pathways across supplied modalities")

        aligned: Dict[str, pd.DataFrame] = {
            name: df.loc[sample_idx, pathway_idx] for name, df in supplied.items()
        }

        if weights is None:
            weights = FusionWeights(
                rna=1.0 if "rna" in supplied else 0.0,
                atac=1.0 if "atac" in supplied else 0.0,
                protein=1.0 if "protein" in supplied else 0.0,
            )

        fused = pd.DataFrame(
            0.0,
            index=sample_idx,
            columns=pathway_idx,
            dtype=float,
        )
        for name in MODALITIES:
            if name not in aligned:
                continue
            w = getattr(weights, name)
            if w == 0:
                continue
            fused = fused + w * aligned[name]

        logger.info(
            "[MultiOmicsFusion] fused n_samples=%d n_pathways=%d modalities=%s",
            len(sample_idx),
            len(pathway_idx),
            sorted(aligned.keys()),
        )

        return FusionResult(
            fused=fused,
            weights=weights,
            per_modality=aligned,
            union_pathways=list(pathway_idx),
        )

    # --------------------------------------------------- learn_weights ---
    def learn_weights(
        self,
        labels: pd.Series,
        rna: Optional[pd.DataFrame] = None,
        atac: Optional[pd.DataFrame] = None,
        protein: Optional[pd.DataFrame] = None,
        grid_step: float = 0.1,
    ) -> FusionWeights:
        """Grid-search per-modality weights that maximise 1-NN accuracy.

        Light-weight: iterates over weight triples on the simplex with
        step ``grid_step``, scoring each candidate with leave-one-out
        nearest-neighbour accuracy on the fused matrix. Returns the
        highest-scoring weight vector.
        """
        supplied = {
            k: v for k, v in [("rna", rna), ("atac", atac), ("protein", protein)] if v is not None
        }
        if not supplied:
            raise ValueError("at least one modality must be supplied")

        # Grid of weight triples summing to 1 on the simplex
        steps = [round(i * grid_step, 4) for i in range(int(1 / grid_step) + 1)]
        best_w: Optional[FusionWeights] = None
        best_acc = -1.0

        for w_rna in steps:
            for w_atac in steps:
                w_prot = 1.0 - w_rna - w_atac
                if w_prot < -1e-9 or w_prot > 1.0 + 1e-9:
                    continue
                w_prot = max(0.0, min(1.0, w_prot))
                # Skip triples touching a modality that wasn't supplied
                if "rna" not in supplied and w_rna > 0:
                    continue
                if "atac" not in supplied and w_atac > 0:
                    continue
                if "protein" not in supplied and w_prot > 0:
                    continue
                try:
                    wts = FusionWeights(rna=w_rna, atac=w_atac, protein=w_prot)
                except ValueError:
                    continue
                result = self.fuse(rna=rna, atac=atac, protein=protein, weights=wts)
                acc = _loo_nn_accuracy(result.fused, labels)
                if acc > best_acc:
                    best_acc = acc
                    best_w = wts

        if best_w is None:  # pragma: no cover - defensive
            raise RuntimeError("weight search produced no candidates")
        logger.info(
            "[MultiOmicsFusion] learned weights %s (LOO-NN acc=%.3f)",
            best_w.as_dict(),
            best_acc,
        )
        return best_w


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _loo_nn_accuracy(features: pd.DataFrame, labels: pd.Series) -> float:
    """Leave-one-out 1-NN accuracy (cosine similarity)."""
    X = features.loc[labels.index].to_numpy(dtype=float)
    y = labels.to_numpy()
    n = len(X)
    if n < 2:
        return 0.0
    # Normalise rows to unit length for cosine
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Xn = X / norms
    sim = Xn @ Xn.T
    # Mask self-similarity
    np.fill_diagonal(sim, -np.inf)
    predictions = y[sim.argmax(axis=1)]
    return float((predictions == y).mean())
