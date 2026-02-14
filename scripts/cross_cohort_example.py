#!/usr/bin/env python3
"""
Cross-cohort validation example.

Generates synthetic cohort pairs, runs cross-cohort comparison,
and saves a markdown report with interpretation guidance.

Usage:
    python scripts/cross_cohort_example.py [--output outputs/cross_cohort/]
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pathway_subtyping import (
    compare_cohorts,
    generate_cross_cohort_report,
    generate_synthetic_cohort_pair,
)


def main():
    parser = argparse.ArgumentParser(description="Cross-cohort validation example")
    parser.add_argument(
        "--output",
        default="outputs/cross_cohort/",
        help="Output directory for report (default: outputs/cross_cohort/)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Cross-Cohort Validation Example")
    print("=" * 40)

    # Step 1: Generate synthetic cohort pair
    print("\n1. Generating synthetic cohort pair...")
    cohort_a, cohort_b = generate_synthetic_cohort_pair(
        n_samples_a=100,
        n_samples_b=80,
        n_subtypes=2,
        n_pathways=10,
        effect_size=1.5,
        seed=args.seed,
    )
    print(f"   Cohort A: {cohort_a.n_samples} samples, {cohort_a.n_clusters} subtypes")
    print(f"   Cohort B: {cohort_b.n_samples} samples, {cohort_b.n_clusters} subtypes")

    # Step 2: Compare cohorts
    print("\n2. Comparing cohorts...")
    result = compare_cohorts(cohort_a, cohort_b, seed=args.seed)

    # Step 3: Display results
    print("\n3. Results:")
    print(result.format_report())

    # Step 4: Generate markdown report
    report_path = output_dir / "cross_cohort_report.md"
    generate_cross_cohort_report([result], str(report_path))
    print(f"\n4. Report saved to: {report_path}")

    # Step 5: Interpretation guidance
    print("\n" + "=" * 40)
    print("Interpretation Guide")
    print("=" * 40)
    print("""
Transfer ARI (Adjusted Rand Index):
  > 0.5  : Good replication — subtypes are consistent across cohorts
  0.3-0.5: Moderate — partially shared subtype structure
  < 0.3  : Weak — subtypes may not replicate

Projection ARI:
  > 0.5  : Cluster geometry is preserved in shared PCA space
  < 0.3  : Cluster structure differs between cohorts

Pathway Correlation:
  > 0.7  : Similar pathway importance across cohorts
  0.3-0.7: Partially shared pathway effects
  < 0.3  : Different pathways drive variation in each cohort

For publication, report all three metrics and note that cross-cohort
validation is essential for distinguishing real subtypes from
dataset-specific artifacts.
""")

    # Citations
    print("Citations:")
    for citation in result.get_citations():
        print(f"  - {citation}")


if __name__ == "__main__":
    main()
