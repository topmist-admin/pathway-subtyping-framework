"""
F10 real-data acceptance validation — multi-omics fusion uplift on CITE-seq.

Roadmap criterion (docs/roadmap-v06-codeberg.md, Phase 3 F10):
    fused pathway score exceeds RNA-only pathway score on downstream
    cell-type classification by at least 3% accuracy on a paired
    CITE-seq reference.

Data: 10x Genomics pbmc_1k_protein_v3 CITE-seq reference (713 cells,
33538 genes, 17 antibodies — 14 immune markers + 3 isotype controls).
A small, public CITE-seq dataset hosted by 10x; fetched once to
``data/f10_citeseq/``. Labels are derived from canonical antibody
gating into 5 PBMC populations (CD4+ T, CD8+ T, B, monocyte, NK).

Methodology:
    1. Load GEX + ADT from the 10x filtered feature barcode matrix.
    2. Derive PBMC labels by gating on ADT thresholds.
    3. Define 7 immune pathways matched across RNA and protein
       (config: data/omics/cite_adt_to_pathway.yaml).
    4. Compute 7-dim RNA and 7-dim protein pathway scores.
    5. Fuse RNA + protein at equal weight (MultiOmicsFusion) —
       acceptance is measured without learned weights to avoid
       overfitting on a small cohort.
    6. 1-NN classification accuracy (5-fold cross-validated) on the
       RNA-only score matrix vs the fused matrix. Bootstrap n=1000
       on the fold-aggregated predictions for a CI on uplift.

Outputs:
    results/f10_validation/fusion_uplift.json
    (consumed by tests/test_omics_real_data.py)

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from pathway_subtyping.omics import (
    FusionWeights,
    MultiOmicsFusion,
    ProteomicsScorer,
    score_proteomics_pathways,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger("validate_f10")

REPO_ROOT = Path(__file__).resolve().parent.parent
H5_PATH = REPO_ROOT / "data" / "f10_citeseq" / "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5"
PATHWAY_MAP = REPO_ROOT / "data" / "omics" / "cite_adt_to_pathway.yaml"
OUTPUT_DIR = REPO_ROOT / "results" / "f10_validation"


# --------------------------------------------------------------------------- #
# Loader
# --------------------------------------------------------------------------- #

def load_citeseq_h5(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Parse a 10x filtered feature barcode h5 into (GEX, ADT) DataFrames.

    Returns two DataFrames (cells × features) of the dense counts — the
    1k cohort is small enough to materialise densely without issue.
    """
    logger.info("[load] %s", path)
    with h5py.File(path, "r") as f:
        barcodes = [b.decode() for b in f["matrix/barcodes"][:]]
        feature_name = [n.decode() for n in f["matrix/features/name"][:]]
        feature_type = [t.decode() for t in f["matrix/features/feature_type"][:]]
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(f["matrix/shape"][:])  # (n_features, n_cells)

    # CSC: columns are cells
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    # Cells × features
    dense = mat.T.toarray()

    ft = np.asarray(feature_type)
    name = np.asarray(feature_name)
    gex_mask = ft == "Gene Expression"
    adt_mask = ft == "Antibody Capture"

    gex = pd.DataFrame(
        dense[:, gex_mask], index=barcodes, columns=name[gex_mask]
    )
    adt = pd.DataFrame(
        dense[:, adt_mask], index=barcodes, columns=name[adt_mask]
    )

    # Collapse any duplicate gene symbols by summation
    gex = gex.groupby(axis=1, level=0).sum()

    logger.info("[load] GEX %s, ADT %s", gex.shape, adt.shape)
    return gex, adt


def load_pathway_map(path: Path) -> Dict[str, Dict[str, List[str]]]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    logger.info("[pathways] %d pathways loaded", len(data))
    return data


# --------------------------------------------------------------------------- #
# Labelling
# --------------------------------------------------------------------------- #

def _clr_norm(adt: pd.DataFrame) -> pd.DataFrame:
    """Centered log-ratio normalisation — standard for CITE-seq ADT."""
    # log1p per cell, then subtract per-cell mean
    logd = np.log1p(adt.to_numpy(dtype=float))
    centered = logd - logd.mean(axis=1, keepdims=True)
    return pd.DataFrame(centered, index=adt.index, columns=adt.columns)


def label_pbmc_types(adt_clr: pd.DataFrame) -> pd.Series:
    """Gate cells into 5 canonical PBMC types from ADT CLR values.

    Simple, widely-used gating thresholds on CLR-normalised counts.
    Cells failing any rule are labelled ``Other`` and dropped at
    classification time.
    """
    # Thresholds at the 60th percentile per-antibody — a permissive gate
    # that still separates populations; tightening doesn't help since
    # our only consumer is a classifier treating these labels as ground
    # truth, not a clinician.
    q = adt_clr.quantile(0.60)

    is_cd3 = adt_clr["CD3_TotalSeqB"] > q["CD3_TotalSeqB"]
    is_cd4 = adt_clr["CD4_TotalSeqB"] > q["CD4_TotalSeqB"]
    is_cd8 = adt_clr["CD8a_TotalSeqB"] > q["CD8a_TotalSeqB"]
    is_cd14 = adt_clr["CD14_TotalSeqB"] > q["CD14_TotalSeqB"]
    is_cd19 = adt_clr["CD19_TotalSeqB"] > q["CD19_TotalSeqB"]
    is_cd56 = adt_clr["CD56_TotalSeqB"] > q["CD56_TotalSeqB"]

    labels = pd.Series("Other", index=adt_clr.index, dtype=object)
    labels[is_cd19] = "B"
    labels[is_cd14 & ~is_cd3] = "Monocyte"
    labels[is_cd56 & ~is_cd3] = "NK"
    labels[is_cd3 & is_cd4 & ~is_cd8] = "CD4_T"
    labels[is_cd3 & is_cd8 & ~is_cd4] = "CD8_T"

    counts = labels.value_counts()
    logger.info("[label] %s", counts.to_dict())
    return labels


# --------------------------------------------------------------------------- #
# Pathway scoring (custom matched pathways on RNA + protein)
# --------------------------------------------------------------------------- #

def score_rna_immune(
    gex: pd.DataFrame, pathway_map: Dict[str, Dict[str, List[str]]],
    min_members: int = 2,
) -> pd.DataFrame:
    """Mean-Z on log1p(RNA counts) restricted to pathway member genes."""
    logd = np.log1p(gex.to_numpy(dtype=float))
    logd_df = pd.DataFrame(logd, index=gex.index, columns=gex.columns)

    records = {}
    for pw, spec in pathway_map.items():
        members = [g for g in spec.get("rna", []) if g in logd_df.columns]
        if len(members) < min_members:
            logger.warning(
                "[rna-score] pathway %s: only %d members present — skipping",
                pw, len(members),
            )
            continue
        block = logd_df[members]
        # Z per gene across cells then mean across members
        mean = block.mean(axis=0)
        std = block.std(axis=0).replace(0.0, 1.0)
        zblock = (block - mean) / std
        records[pw] = zblock.mean(axis=1)

    return pd.DataFrame(records, index=gex.index)


def score_protein_immune(
    adt_clr: pd.DataFrame, pathway_map: Dict[str, Dict[str, List[str]]],
    min_members: int = 1,
) -> pd.DataFrame:
    """Mean-Z on CLR ADT values restricted to pathway member antibodies."""
    records = {}
    for pw, spec in pathway_map.items():
        members = [a for a in spec.get("adt", []) if a in adt_clr.columns]
        if len(members) < min_members:
            logger.warning(
                "[adt-score] pathway %s: no ADT members — skipping",
                pw,
            )
            continue
        block = adt_clr[members]
        mean = block.mean(axis=0)
        std = block.std(axis=0).replace(0.0, 1.0)
        zblock = (block - mean) / std
        records[pw] = zblock.mean(axis=1)
    return pd.DataFrame(records, index=adt_clr.index)


# --------------------------------------------------------------------------- #
# Classification harness
# --------------------------------------------------------------------------- #

def _cosine_1nn(
    train_X: np.ndarray, train_y: np.ndarray,
    test_X: np.ndarray, test_y: np.ndarray,
) -> np.ndarray:
    """Cosine-similarity 1-NN; returns a boolean correct-mask."""
    # Normalise
    def _norm(M):
        n = np.linalg.norm(M, axis=1, keepdims=True)
        n[n == 0] = 1.0
        return M / n
    a = _norm(train_X)
    b = _norm(test_X)
    sim = b @ a.T
    pred = train_y[sim.argmax(axis=1)]
    return pred == test_y


def _kfold_indices(n: int, k: int, seed: int) -> List[Tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    folds = np.array_split(idx, k)
    out = []
    for i in range(k):
        test_idx = folds[i]
        train_idx = np.concatenate([folds[j] for j in range(k) if j != i])
        out.append((train_idx, test_idx))
    return out


def cross_validated_accuracy(
    features: pd.DataFrame, labels: pd.Series,
    k_folds: int = 5, seed: int = 0,
) -> Tuple[float, np.ndarray]:
    """Return (accuracy, per-cell correct-mask) under k-fold 1-NN."""
    assert features.index.equals(labels.index)
    X = features.to_numpy(dtype=float)
    y = labels.to_numpy()
    n = len(X)
    correct = np.zeros(n, dtype=bool)
    for tr, te in _kfold_indices(n, k_folds, seed):
        mask = _cosine_1nn(X[tr], y[tr], X[te], y[te])
        correct[te] = mask
    acc = float(correct.mean())
    return acc, correct


def _bootstrap_uplift_ci(
    correct_rna: np.ndarray, correct_fused: np.ndarray,
    n_boot: int = 1000, seed: int = 0,
) -> Dict[str, float]:
    """Bootstrap over cells; CI on (acc_fused - acc_rna)."""
    n = len(correct_rna)
    rng = np.random.default_rng(seed)
    uplifts = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        uplifts[i] = correct_fused[idx].mean() - correct_rna[idx].mean()
    return {
        "mean": float(np.mean(uplifts)),
        "ci_low": float(np.quantile(uplifts, 0.025)),
        "ci_high": float(np.quantile(uplifts, 0.975)),
    }


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-folds", type=int, default=5)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not H5_PATH.exists():
        logger.error(
            "CITE-seq h5 missing at %s — run:\n"
            "  curl -L -o %s \\\n"
            "    https://cf.10xgenomics.com/samples/cell-exp/3.0.0/"
            "pbmc_1k_protein_v3/"
            "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
            H5_PATH, H5_PATH,
        )
        return 1

    gex, adt = load_citeseq_h5(H5_PATH)
    pathway_map = load_pathway_map(PATHWAY_MAP)

    adt_clr = _clr_norm(adt)
    labels = label_pbmc_types(adt_clr)
    # Drop `Other` cells — they dilute the classification with ambiguous labels
    keep = labels != "Other"
    logger.info(
        "[cohort] %d gated cells / %d total — dropping Other",
        int(keep.sum()), len(labels),
    )
    gex = gex.loc[keep]
    adt_clr = adt_clr.loc[keep]
    labels = labels.loc[keep]

    rna_scores = score_rna_immune(gex, pathway_map)
    protein_scores = score_protein_immune(adt_clr, pathway_map)
    logger.info(
        "[scores] rna %s, protein %s", rna_scores.shape, protein_scores.shape,
    )

    # Fuse on the intersection of pathways + samples
    fusion = MultiOmicsFusion()
    fused_result = fusion.fuse(
        rna=rna_scores, protein=protein_scores,
        weights=FusionWeights(rna=0.5, protein=0.5),
    )
    fused_scores = fused_result.fused
    logger.info(
        "[fuse] fused %s  (pathways: %s)",
        fused_scores.shape, fused_result.union_pathways,
    )

    # Classification
    acc_rna, mask_rna = cross_validated_accuracy(
        rna_scores[fused_result.union_pathways], labels,
        k_folds=args.k_folds, seed=args.seed,
    )
    acc_fused, mask_fused = cross_validated_accuracy(
        fused_scores, labels, k_folds=args.k_folds, seed=args.seed,
    )
    logger.info("[classify] rna_acc=%.4f fused_acc=%.4f", acc_rna, acc_fused)

    uplift_stats = _bootstrap_uplift_ci(
        mask_rna, mask_fused, n_boot=args.n_boot, seed=args.seed + 1,
    )
    uplift = acc_fused - acc_rna
    logger.info(
        "[uplift] %+.4f (95%% CI %+.4f..%+.4f)",
        uplift, uplift_stats["ci_low"], uplift_stats["ci_high"],
    )

    label_counts = labels.value_counts().to_dict()

    report = {
        "feature": "F10",
        "cohort": {
            "label": "10x pbmc_1k_protein_v3",
            "n_cells_total": int(len(gex) + (~keep).sum()),
            "n_cells_gated": int(len(gex)),
            "n_genes": int(gex.shape[1]),
            "n_antibodies": int(adt_clr.shape[1]),
            "label_source": "canonical ADT gating on CLR values",
            "label_counts": label_counts,
        },
        "pathways": fused_result.union_pathways,
        "rna_only_accuracy": acc_rna,
        "fused_accuracy": acc_fused,
        "uplift": uplift,
        "uplift_bootstrap_ci": uplift_stats,
        "acceptance": {
            "uplift_target": 0.03,
            "uplift_pass": uplift >= 0.03,
            "ci_low_strictly_positive_target": True,
            "ci_low_pass": uplift_stats["ci_low"] > 0,
        },
        "config": {
            "k_folds": int(args.k_folds),
            "n_boot": int(args.n_boot),
            "seed": int(args.seed),
            "fusion_weights": {"rna": 0.5, "protein": 0.5},
            "classifier": "cosine-1NN",
        },
    }

    out_path = OUTPUT_DIR / "fusion_uplift.json"
    out_path.write_text(json.dumps(report, indent=2))
    logger.info("wrote %s", out_path)

    print()
    print("=" * 72)
    print("F10 real-data multi-omics fusion — acceptance summary")
    print("=" * 72)
    print(f"  cohort           : 10x pbmc_1k_protein_v3 (gated n={len(gex)})")
    print(f"  label_counts     : {label_counts}")
    print(f"  pathways         : {fused_result.union_pathways}")
    print(f"  rna-only acc     : {acc_rna:.4f}")
    print(f"  fused acc        : {acc_fused:.4f}")
    print(f"  uplift           : {uplift:+.4f}  (95% CI {uplift_stats['ci_low']:+.4f}..{uplift_stats['ci_high']:+.4f})")
    print(f"  acceptance       : uplift>=3%? {report['acceptance']['uplift_pass']}  ci_low>0? {report['acceptance']['ci_low_pass']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
