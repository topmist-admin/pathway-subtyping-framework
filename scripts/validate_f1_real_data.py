"""
F1 real-data acceptance validation — conformal coverage on TCGA-COAD + GSE28521.

Roadmap criterion (docs/roadmap-v06-codeberg.md, Phase 1 F1):
    Conformal intervals achieve nominal coverage (e.g., 90%) on held-out
    TCGA-COAD and autism benchmark datasets within +/- 2%.

Methodology:
    For each dataset, compute a pathway-score matrix (samples x pathways).
    Hold one pathway out as the target y; use the remaining pathways as the
    feature matrix X. Split into fit/cal/test, fit a linear regressor, wrap
    with ConformalPathwayPredictor, measure empirical coverage on the test
    split. Iterate across multiple target pathways and seeds; report the
    mean empirical deviation from nominal.

Outputs:
    results/f1_validation/TCGA-COAD_conformal_coverage.json
    results/f1_validation/GSE28521_conformal_coverage.json

Both JSONs are consumed by tests/test_uncertainty_real_data.py.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from pathway_subtyping.config import validate_gmt_file
from pathway_subtyping.expression import (
    ExpressionScoringMethod,
    score_pathways_from_expression,
)
from pathway_subtyping.uncertainty import ConformalPathwayPredictor

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger("validate_f1")

REPO_ROOT = Path(__file__).resolve().parent.parent
COAD_DIR = REPO_ROOT / "tcga_data" / "TCGA-COAD"
GSE28521_DIR = REPO_ROOT / "research-results" / "GSE28521"
HALLMARK_GMT = REPO_ROOT / "data" / "pathways" / "hallmark_200genes.gmt"
OUTPUT_DIR = REPO_ROOT / "results" / "f1_validation"


# --------------------------------------------------------------------------- #
# TCGA-COAD loader
# --------------------------------------------------------------------------- #

def load_tcga_coad_expression(coad_dir: Path) -> pd.DataFrame:
    """Load 57 TCGA-COAD STAR TSVs into a (samples x genes) log-TPM matrix."""
    tsv_files = sorted(coad_dir.glob("*.rna_seq.augmented_star_gene_counts.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"no COAD TSVs under {coad_dir}")

    logger.info("[COAD] loading %d TSVs", len(tsv_files))
    per_sample: List[pd.Series] = []
    sample_ids: List[str] = []
    for path in tsv_files:
        df = pd.read_csv(path, sep="\t", comment="#")
        df = df[~df["gene_id"].astype(str).str.startswith("N_")]
        df = df.dropna(subset=["gene_name"])
        # TPM column is numeric; drop any residual non-numeric rows
        df["tpm_unstranded"] = pd.to_numeric(
            df["tpm_unstranded"], errors="coerce"
        )
        df = df.dropna(subset=["tpm_unstranded"])
        # Collapse duplicate gene_name by max TPM
        s = df.groupby("gene_name")["tpm_unstranded"].max()
        per_sample.append(s)
        sample_ids.append(path.stem.split(".")[0])

    matrix = pd.concat(per_sample, axis=1, keys=sample_ids).T
    matrix.index.name = "sample"
    matrix = matrix.fillna(0.0)
    # log1p on TPM to stabilize variance
    matrix = np.log1p(matrix)
    logger.info("[COAD] expression matrix shape: %s", matrix.shape)
    return matrix


def compute_coad_pathway_scores(expr: pd.DataFrame) -> pd.DataFrame:
    """Hallmark pathway scores via mean-Z scoring (no extra deps)."""
    pathways = validate_gmt_file(str(HALLMARK_GMT))
    result = score_pathways_from_expression(
        gene_expression=expr,
        pathways=pathways,
        method=ExpressionScoringMethod.MEAN_Z,
        min_genes_per_pathway=20,
        show_progress=False,
    )
    scores = result.pathway_scores
    logger.info(
        "[COAD] pathway score matrix: %s (skipped: %d)",
        scores.shape, len(result.skipped_pathways),
    )
    return scores


# --------------------------------------------------------------------------- #
# GSE28521 loader (pre-computed)
# --------------------------------------------------------------------------- #

def load_gse28521_pathway_scores(gse_dir: Path) -> pd.DataFrame:
    """79 samples x 15 autism pathway scores, already Z-normalised."""
    path = gse_dir / "pathway_scores.csv"
    df = pd.read_csv(path, index_col=0)
    logger.info("[GSE28521] pre-computed pathway scores: %s", df.shape)
    return df


# --------------------------------------------------------------------------- #
# Conformal coverage harness
# --------------------------------------------------------------------------- #

def _linear_score_fn(coef: np.ndarray, intercept: float):
    def score_fn(X: np.ndarray) -> np.ndarray:
        return X @ coef + intercept
    return score_fn


def _fit_linear(X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, float]:
    """Closed-form OLS with a bias term. Ridge fallback if rank-deficient."""
    ones = np.ones((X.shape[0], 1))
    Xa = np.hstack([X, ones])
    try:
        beta, *_ = np.linalg.lstsq(Xa, y, rcond=None)
    except np.linalg.LinAlgError:
        lam = 1e-3 * np.eye(Xa.shape[1])
        beta = np.linalg.solve(Xa.T @ Xa + lam, Xa.T @ y)
    return beta[:-1], float(beta[-1])


def run_conformal_coverage(
    scores: pd.DataFrame,
    dataset_name: str,
    target_coverages: Tuple[float, ...] = (0.80, 0.90, 0.95),
    n_seeds: int = 10,
    split_frac: Tuple[float, float] = (0.35, 0.80),
) -> Dict:
    """Run the full conformal-coverage protocol on a pathway score matrix.

    For each target pathway and seed:
        - permute samples; split into fit / cal / test via split_frac
        - fit a linear map (other pathways -> target pathway) on fit split
        - calibrate ConformalPathwayPredictor on cal split
        - measure empirical coverage on test split

    Returns a dict keyed by target_coverage with mean/std deviation statistics.
    """
    scores = scores.dropna(axis=0, how="any")
    n_samples, n_pathways = scores.shape
    logger.info(
        "[%s] coverage protocol: n_samples=%d n_pathways=%d n_seeds=%d",
        dataset_name, n_samples, n_pathways, n_seeds,
    )

    pathway_names = list(scores.columns)
    score_values = scores.to_numpy()

    deviations: Dict[float, List[float]] = {t: [] for t in target_coverages}
    oracle_devs: Dict[float, List[float]] = {t: [] for t in target_coverages}
    per_pathway: Dict[str, Dict[float, List[float]]] = {
        p: {t: [] for t in target_coverages} for p in pathway_names
    }

    for pathway_idx, pathway_name in enumerate(pathway_names):
        y = score_values[:, pathway_idx]
        X = np.delete(score_values, pathway_idx, axis=1)

        for seed in range(n_seeds):
            rng = np.random.default_rng(1000 * pathway_idx + seed)
            perm = rng.permutation(n_samples)
            fit_end = int(split_frac[0] * n_samples)
            cal_end = int(split_frac[1] * n_samples)
            fit_idx = perm[:fit_end]
            cal_idx = perm[fit_end:cal_end]
            te_idx = perm[cal_end:]

            if len(cal_idx) < 10 or len(te_idx) < 5:
                continue  # skip degenerate splits

            coef, intercept = _fit_linear(X[fit_idx], y[fit_idx])
            score_fn = _linear_score_fn(coef, intercept)

            n_cal = len(cal_idx)
            for target in target_coverages:
                predictor = ConformalPathwayPredictor(
                    score_fn=score_fn, coverage=target
                )
                predictor.calibrate(X[cal_idx], y[cal_idx])
                empirical = predictor.coverage_on(X[te_idx], y[te_idx])
                dev = empirical - target
                # Oracle-adjusted: account for the finite-sample discretization
                # of split conformal. With n calibration points, the smallest
                # achievable upper coverage bound is ceil((n+1)*(1-alpha))/(n+1)
                # instead of (1-alpha). The oracle deviation reports how close
                # we are to that theoretical floor, removing the discretization
                # bias that dominates at small n_cal.
                alpha = 1.0 - target
                oracle = min(
                    1.0,
                    np.ceil((n_cal + 1) * (1.0 - alpha)) / (n_cal + 1),
                )
                oracle_dev = empirical - oracle
                deviations[target].append(dev)
                oracle_devs[target].append(oracle_dev)
                per_pathway[pathway_name][target].append(dev)

    report = {
        "dataset": dataset_name,
        "n_samples": int(n_samples),
        "n_pathways": int(n_pathways),
        "n_seeds_per_pathway": int(n_seeds),
        "target_coverages": list(target_coverages),
        "summary": {},
        "per_pathway_mean_dev": {},
    }

    for target, devs in deviations.items():
        arr = np.asarray(devs)
        oarr = np.asarray(oracle_devs[target])
        report["summary"][f"{target:.2f}"] = {
            "n_runs": int(len(arr)),
            "mean_deviation": float(arr.mean()) if arr.size else None,
            "std_deviation": float(arr.std(ddof=1)) if arr.size > 1 else None,
            "mean_abs_deviation": float(np.abs(arr).mean()) if arr.size else None,
            "mean_oracle_deviation": float(oarr.mean()) if oarr.size else None,
            "mean_abs_oracle_deviation": (
                float(np.abs(oarr).mean()) if oarr.size else None
            ),
            "within_2pct_fraction": (
                float((np.abs(arr) < 0.02).mean()) if arr.size else None
            ),
        }

    for pname, by_target in per_pathway.items():
        report["per_pathway_mean_dev"][pname] = {
            f"{t:.2f}": float(np.mean(v)) if v else None
            for t, v in by_target.items()
        }

    return report


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-coad", action="store_true",
        help="Skip TCGA-COAD (useful if STAR TSVs not available).",
    )
    parser.add_argument(
        "--skip-gse28521", action="store_true",
        help="Skip GSE28521 autism benchmark.",
    )
    parser.add_argument("--n-seeds", type=int, default=10)
    args = parser.parse_args(argv)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    reports: Dict[str, Dict] = {}

    if not args.skip_coad:
        logger.info("=== TCGA-COAD ===")
        expr = load_tcga_coad_expression(COAD_DIR)
        scores = compute_coad_pathway_scores(expr)
        reports["TCGA-COAD"] = run_conformal_coverage(
            scores, dataset_name="TCGA-COAD", n_seeds=args.n_seeds,
        )
        out_path = OUTPUT_DIR / "TCGA-COAD_conformal_coverage.json"
        out_path.write_text(json.dumps(reports["TCGA-COAD"], indent=2))
        logger.info("[COAD] wrote %s", out_path)

    if not args.skip_gse28521:
        logger.info("=== GSE28521 autism benchmark ===")
        scores = load_gse28521_pathway_scores(GSE28521_DIR)
        reports["GSE28521"] = run_conformal_coverage(
            scores, dataset_name="GSE28521", n_seeds=args.n_seeds,
        )
        out_path = OUTPUT_DIR / "GSE28521_conformal_coverage.json"
        out_path.write_text(json.dumps(reports["GSE28521"], indent=2))
        logger.info("[GSE28521] wrote %s", out_path)

    # Console summary
    print()
    print("=" * 72)
    print("F1 real-data conformal coverage — acceptance summary")
    print("=" * 72)
    for name, rep in reports.items():
        print(f"\n{name}  (n={rep['n_samples']}, pathways={rep['n_pathways']})")
        print(
            f"{'target':>8}  {'mean_dev':>10}  {'oracle_dev':>12}  "
            f"{'|mean_dev|':>10}  n_runs"
        )
        for target_key, stats in rep["summary"].items():
            if stats["mean_deviation"] is None:
                continue
            print(
                f"{target_key:>8}  "
                f"{stats['mean_deviation']:>+10.4f}  "
                f"{stats['mean_oracle_deviation']:>+12.4f}  "
                f"{stats['mean_abs_deviation']:>10.4f}  "
                f"{stats['n_runs']}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
