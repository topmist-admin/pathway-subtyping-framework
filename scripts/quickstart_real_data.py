"""
Quickstart — PSF on a real GEO cohort in 15 minutes.

Runs the full score → cluster → report flow on GSE28521 (79 post-mortem
frontal-cortex samples, autism vs control; Affymetrix U133). Written
for a beginner: no YAML config, no VCF, just pandas + the high-level PSF
APIs. Designed to be copied and adapted to another cohort.

What you get:
    outputs/quickstart_real_data/
        pathway_scores.csv   — 79 samples × 50 Hallmark pathways, Z-scored
        subtype_assignments.csv — GMM cluster ID per sample
        summary.md           — human-readable one-page summary

Run:
    python scripts/quickstart_real_data.py

Expected runtime: 1–3 minutes after the data is cached locally.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from pathway_subtyping.clustering import run_clustering
from pathway_subtyping.config import validate_gmt_file
from pathway_subtyping.expression import (
    ExpressionScoringMethod,
    score_pathways_from_expression,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger("quickstart_real_data")

REPO_ROOT = Path(__file__).resolve().parent.parent
EXPR_PATH = REPO_ROOT / "research-results" / "GSE28521" / "gene_expression_processed.csv"
HALLMARK_GMT = REPO_ROOT / "data" / "pathways" / "hallmark_200genes.gmt"
OUTPUT_DIR = REPO_ROOT / "outputs" / "quickstart_real_data"


def ensure_expression_matrix() -> pd.DataFrame:
    """Load the GSE28521 expression matrix; suggest the git-clone path if absent."""
    if not EXPR_PATH.exists():
        raise FileNotFoundError(
            f"Expected the GSE28521 expression matrix at {EXPR_PATH}. "
            f"This file ships with the cloned repo under research-results/. "
            f"If you installed from PyPI, clone the repo for sample data:\n"
            f"  git clone https://codeberg.org/pathways/pathway-subtyping-framework.git\n"
            f"  cd pathway-subtyping-framework\n"
            f"  python scripts/quickstart_real_data.py"
        )
    df = pd.read_csv(EXPR_PATH, index_col=0)
    logger.info("[load] GSE28521 expression: %s", df.shape)
    return df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--k", type=int, default=3,
        help="Number of subtypes to discover (default: 3).",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility.",
    )
    args = parser.parse_args(argv)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load real expression data (samples × genes, pre-normalised log2)
    expr = ensure_expression_matrix()

    # 2. Score Hallmark pathways via mean-Z (fastest option; see
    # docs/api/expression.md for SSGSEA / GSVA alternatives).
    pathways = validate_gmt_file(str(HALLMARK_GMT))
    result = score_pathways_from_expression(
        gene_expression=expr,
        pathways=pathways,
        method=ExpressionScoringMethod.MEAN_Z,
        min_genes_per_pathway=20,
        show_progress=False,
    )
    scores = result.pathway_scores
    logger.info("[score] pathway scores: %s (skipped: %d)",
                scores.shape, len(result.skipped_pathways))

    scores_path = OUTPUT_DIR / "pathway_scores.csv"
    scores.to_csv(scores_path)
    logger.info("[wrote] %s", scores_path)

    # 3. Discover subtypes via GMM. In production you'd sweep k via BIC
    # (see docs/api/clustering.md); here we pin k for reproducibility.
    clustering = run_clustering(
        data=scores.to_numpy(),
        n_clusters=args.k,
        seed=args.seed,
    )
    labels = pd.Series(
        clustering.labels, index=scores.index, name="subtype",
    )
    assignments_path = OUTPUT_DIR / "subtype_assignments.csv"
    labels.to_csv(assignments_path)
    logger.info("[wrote] %s", assignments_path)

    # 4. Summary report
    counts = labels.value_counts().sort_index()
    top_pathways_per_subtype = {}
    for subtype_id in sorted(labels.unique()):
        mask = labels == subtype_id
        pathway_means = scores.loc[mask].mean(axis=0).sort_values(ascending=False)
        top_pathways_per_subtype[int(subtype_id)] = pathway_means.head(3).index.tolist()

    summary_lines = [
        "# PSF Quickstart — Real Data (GSE28521)",
        "",
        f"Ran `scripts/quickstart_real_data.py` with k={args.k}, seed={args.seed}.",
        "",
        "## Input",
        "",
        f"- Cohort: GSE28521 (post-mortem frontal cortex, autism + control)",
        f"- Samples: {scores.shape[0]}",
        f"- Pathways: {scores.shape[1]} (MSigDB Hallmark)",
        f"- Scoring: mean-Z on log2 expression",
        "",
        "## Subtype sizes",
        "",
        "| Subtype | n samples |",
        "|---------|-----------|",
    ]
    for sid, n in counts.items():
        summary_lines.append(f"| {sid} | {n} |")
    summary_lines += [
        "",
        "## Top 3 pathways by mean score, per subtype",
        "",
    ]
    for sid, top in sorted(top_pathways_per_subtype.items()):
        summary_lines.append(f"- **Subtype {sid}**: " + ", ".join(top))
    summary_lines += [
        "",
        "## Files",
        "",
        f"- `{scores_path.relative_to(REPO_ROOT)}` — sample × pathway score matrix",
        f"- `{assignments_path.relative_to(REPO_ROOT)}` — per-sample subtype IDs",
        "",
        "## Where to go next",
        "",
        "- **[Notebook 10](../../examples/notebooks/10_geo_autism_bulk.ipynb)** —"
        " the same cohort walked through cell-by-cell with full validation.",
        "- **[docs/api/index.md](../../docs/api/index.md)** — v0.5 and v0.6 module"
        " catalogue (pathway scoring, clustering, validation gates, uncertainty,"
        " harmonize, perturb, multi-omics, etc).",
        "- **[docs/guides/validation-gates.md](../../docs/guides/validation-gates.md)** —"
        " run the 5 validation gates on these subtypes before trusting them.",
        "- **Adapt this script** — the two inputs you'd change for your own data"
        " are `EXPR_PATH` (samples-by-genes CSV) and `HALLMARK_GMT` (any GMT).",
        "",
        "Research use only. Not for clinical decision-making.",
        "",
    ]
    summary_path = OUTPUT_DIR / "summary.md"
    summary_path.write_text("\n".join(summary_lines))
    logger.info("[wrote] %s", summary_path)

    print()
    print("=" * 72)
    print("PSF Quickstart (Real Data) — GSE28521 autism cortex")
    print("=" * 72)
    print(f"  samples     : {scores.shape[0]}")
    print(f"  pathways    : {scores.shape[1]}")
    print(f"  subtypes (k): {args.k}")
    print(f"  sizes       : {counts.to_dict()}")
    print(f"  outputs     : {OUTPUT_DIR.relative_to(REPO_ROOT)}/")
    print()
    print("Open outputs/quickstart_real_data/summary.md for top pathways per subtype.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
