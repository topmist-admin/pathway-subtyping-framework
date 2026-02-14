#!/usr/bin/env python3
"""
Generate calibration lookup tables for threshold_calibration.py.

This script runs simulations to produce the pre-computed lookup tables
embedded in the threshold_calibration module. It is long-running
(~10-30 minutes depending on hardware) and only needs to be run when
regenerating the tables.

Usage:
    python scripts/generate_calibration_table.py [--n_simulations 500] [--seed 42]
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pathway_subtyping.threshold_calibration import generate_calibration_table


def main():
    parser = argparse.ArgumentParser(
        description="Generate calibration lookup tables"
    )
    parser.add_argument(
        "--n_simulations", type=int, default=500,
        help="Number of simulations per grid point (default: 500)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)",
    )
    args = parser.parse_args()

    print(f"Generating calibration tables with n_simulations={args.n_simulations}, seed={args.seed}")
    print("This may take a while...")

    tables = generate_calibration_table(
        n_simulations=args.n_simulations,
        seed=args.seed,
    )

    # Print as Python dict literals for copy-paste
    print("\n# Null ARI table (95th percentile)")
    print("_NULL_ARI_TABLE = {")
    for (n, k), v in sorted(tables["null_ari"].items()):
        print(f"    ({n}, {k}): {v:.4f},")
    print("}")

    print("\n# Stability table (5th percentile)")
    print("_STABILITY_TABLE = {")
    for (n, k), v in sorted(tables["stability"].items()):
        print(f"    ({n}, {k}): {v:.4f},")
    print("}")


if __name__ == "__main__":
    main()
