"""
Command-line interface for pathway subtyping framework.
"""

import argparse
import sys
from .config import load_config
from .pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Pathway Subtyping Framework - Disease-agnostic molecular subtyping"
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        help="Path to YAML configuration file"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output"
    )

    args = parser.parse_args()

    try:
        config = load_config(args.config)
        if args.verbose:
            config["pipeline"]["verbose"] = True

        results = run_pipeline(config)

        # Print summary
        if results.get("validation_gates", {}).get("all_passed"):
            print("✓ All validation gates passed")
        else:
            print("✗ Some validation gates failed - review results carefully")

    except NotImplementedError as e:
        print(f"Error: {e}")
        print("\nThis is a template repository. To use it:")
        print("1. Copy core modules from autism-pathway-framework")
        print("2. Adapt for disease-agnostic use")
        print("3. Run with your disease-specific config")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
