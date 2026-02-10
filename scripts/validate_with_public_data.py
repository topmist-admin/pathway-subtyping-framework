#!/usr/bin/env python3
"""
Validate pathway-subtyping framework against public genomic datasets.

Downloads ClinVar and Reactome data, validates pathway definitions,
generates disease-realistic synthetic data, and produces a validation report.

Usage:
    python scripts/validate_with_public_data.py
    python scripts/validate_with_public_data.py --disease autism --output outputs/validation/
    python scripts/validate_with_public_data.py --offline  # Skip downloads, use cache only

Research use only. Not for clinical decision-making.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

# Add src to path if running from repo root
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pathway_subtyping.validation_datasets import run_full_validation


def main():
    parser = argparse.ArgumentParser(
        description="Validate pathway-subtyping framework against public datasets"
    )
    parser.add_argument(
        "--disease",
        default="autism",
        help="Disease name for pathway GMT lookup (default: autism)",
    )
    parser.add_argument(
        "--gmt-path",
        default=None,
        help="Explicit path to GMT file (overrides --disease lookup)",
    )
    parser.add_argument(
        "--output",
        default="outputs/validation",
        help="Output directory for reports (default: outputs/validation/)",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Cache directory for downloads (default: data/validation_cache/)",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Number of synthetic samples (default: 100)",
    )
    parser.add_argument(
        "--n-subtypes",
        type=int,
        default=3,
        help="Number of planted subtypes (default: 3)",
    )
    parser.add_argument(
        "--effect-size",
        type=float,
        default=1.0,
        help="Effect size for synthetic data (default: 1.0)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Skip downloads, use cached data only",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download even if cached",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    cache_dir = Path(args.cache_dir) if args.cache_dir else None

    # Run validation
    print(f"Running validation for disease: {args.disease}")
    print(f"  Samples: {args.n_samples}, Subtypes: {args.n_subtypes}, "
          f"Effect size: {args.effect_size}")
    print(f"  Offline mode: {args.offline}")
    print()

    report = run_full_validation(
        gmt_path=args.gmt_path,
        disease_name=args.disease,
        n_samples=args.n_samples,
        n_subtypes=args.n_subtypes,
        effect_size=args.effect_size,
        cache_dir=cache_dir,
        seed=args.seed,
        skip_download=args.offline,
    )

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write JSON report
    json_path = output_dir / "validation_report.json"
    with open(json_path, "w") as f:
        json.dump(report.to_dict(), f, indent=2)
    print(f"JSON report: {json_path}")

    # Write Markdown report
    md_path = output_dir / "validation_report.md"
    with open(md_path, "w") as f:
        f.write(report.format_report())
    print(f"Markdown report: {md_path}")

    # Print summary
    print()
    print(report.format_report())

    # Exit code
    if report.overall_pass:
        print("\nValidation PASSED")
        sys.exit(0)
    else:
        print("\nValidation NEEDS REVIEW")
        for w in report.warnings:
            print(f"  WARNING: {w}")
        sys.exit(1)


if __name__ == "__main__":
    main()
