#!/usr/bin/env python3
"""
Performance benchmark for the Pathway Subtyping Framework.

Generates synthetic data at scale, runs clustering and validation,
and reports wall-clock time and peak memory per step.

Usage:
    python scripts/benchmark_performance.py
    python scripts/benchmark_performance.py --n-samples 10000
    python scripts/benchmark_performance.py --n-samples 5000 --skip-validation

Targets (Issue #8):
    - 10K samples in < 30 minutes total
    - Peak memory < 8 GB
"""

import argparse
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pathway_subtyping.simulation import SimulationConfig, generate_synthetic_data
from pathway_subtyping.validation import ValidationGates


def format_time(seconds: float) -> str:
    """Format seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    secs = seconds % 60
    return f"{minutes}m {secs:.1f}s"


def format_memory(bytes_val: float) -> str:
    """Format bytes into human-readable string."""
    mb = bytes_val / (1024 * 1024)
    if mb > 1024:
        return f"{mb / 1024:.2f} GB"
    return f"{mb:.1f} MB"


def run_benchmark(
    n_samples: int = 10000,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    seed: int = 42,
    skip_validation: bool = False,
    n_permutations: int = 50,
    n_bootstrap: int = 25,
):
    """Run the full benchmark."""
    print("=" * 60)
    print("Pathway Subtyping Framework — Performance Benchmark")
    print("=" * 60)
    print(f"  Samples:      {n_samples:,}")
    print(f"  Pathways:     {n_pathways}")
    print(f"  Subtypes:     {n_subtypes}")
    print(f"  Seed:         {seed}")
    print(f"  Validation:   {'skip' if skip_validation else f'{n_permutations} perms / {n_bootstrap} bootstrap'}")
    print()

    tracemalloc.start()
    results = []

    # Step 1: Generate synthetic data
    print("[1/3] Generating synthetic data...")
    t0 = time.time()
    config = SimulationConfig(
        n_samples=n_samples,
        n_pathways=n_pathways,
        n_subtypes=n_subtypes,
        effect_size=1.5,
        seed=seed,
    )
    data = generate_synthetic_data(config)
    t_synth = time.time() - t0
    mem_synth = tracemalloc.get_traced_memory()[1]
    results.append(("Synthetic data generation", t_synth, mem_synth))
    print(f"      Done in {format_time(t_synth)} | Peak: {format_memory(mem_synth)}")
    print(f"      Shape: {data.pathway_scores.shape}")

    # Step 2: GMM clustering
    print("[2/3] Running GMM clustering...")
    t0 = time.time()
    from sklearn.mixture import GaussianMixture

    gmm = GaussianMixture(
        n_components=n_subtypes,
        random_state=seed,
        n_init=3,
    )
    labels = gmm.fit_predict(data.pathway_scores.values)
    t_gmm = time.time() - t0
    mem_gmm = tracemalloc.get_traced_memory()[1]
    results.append(("GMM clustering", t_gmm, mem_gmm))
    print(f"      Done in {format_time(t_gmm)} | Peak: {format_memory(mem_gmm)}")

    unique, counts = np.unique(labels, return_counts=True)
    for u, c in zip(unique, counts):
        print(f"      Cluster {u}: {c:,} samples")

    # Step 3: Validation gates
    if skip_validation:
        print("[3/3] Validation gates — SKIPPED")
        t_val = 0.0
    else:
        print(f"[3/3] Running validation gates ({n_permutations} perms, {n_bootstrap} bootstrap)...")
        t0 = time.time()
        vg = ValidationGates(
            seed=seed,
            n_permutations=n_permutations,
            n_bootstrap=n_bootstrap,
            show_progress=True,
        )
        val_result = vg.run_all(
            pathway_scores=data.pathway_scores,
            cluster_labels=labels,
            pathways=data.pathways,
            gene_burdens=data.gene_burdens,
            n_clusters=n_subtypes,
            gmm_seed=seed,
        )
        t_val = time.time() - t0
        mem_val = tracemalloc.get_traced_memory()[1]
        results.append(("Validation gates", t_val, mem_val))
        print(f"      Done in {format_time(t_val)} | Peak: {format_memory(mem_val)}")
        print(f"      All passed: {val_result.all_passed}")
        for r in val_result.results:
            status = "PASS" if r.passed else "FAIL"
            print(f"        [{status}] {r.name}: {r.metric_value:.4f} ({r.comparison} {r.threshold})")

    # Summary
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    total_time = sum(r[1] for r in results)

    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"{'Step':<30} {'Time':>10} {'Peak Memory':>15}")
    print("-" * 60)
    for name, t, m in results:
        print(f"{name:<30} {format_time(t):>10} {format_memory(m):>15}")
    print("-" * 60)
    print(f"{'TOTAL':<30} {format_time(total_time):>10} {format_memory(peak_mem):>15}")
    print()

    # Pass/fail against targets
    time_target = 30 * 60  # 30 minutes
    mem_target = 8 * 1024 * 1024 * 1024  # 8 GB

    time_ok = total_time < time_target
    mem_ok = peak_mem < mem_target

    print("Targets:")
    status = "PASS" if time_ok else "FAIL"
    print(f"  [{status}] Total time < 30 min: {format_time(total_time)}")
    status = "PASS" if mem_ok else "FAIL"
    print(f"  [{status}] Peak memory < 8 GB:  {format_memory(peak_mem)}")
    print()

    if time_ok and mem_ok:
        print("All targets met.")
    else:
        print("Some targets not met. Consider:")
        if not time_ok:
            print("  - Reducing n_permutations / n_bootstrap")
            print("  - Using fewer pathways")
        if not mem_ok:
            print("  - Using chunked processing (use_chunked_processing=True)")
            print("  - Reducing sample count")

    return time_ok and mem_ok


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark the Pathway Subtyping Framework"
    )
    parser.add_argument(
        "--n-samples", type=int, default=10000,
        help="Number of synthetic samples (default: 10000)",
    )
    parser.add_argument(
        "--n-pathways", type=int, default=15,
        help="Number of pathways (default: 15)",
    )
    parser.add_argument(
        "--n-subtypes", type=int, default=3,
        help="Number of planted subtypes (default: 3)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--skip-validation", action="store_true",
        help="Skip validation gates (benchmark data generation + clustering only)",
    )
    parser.add_argument(
        "--n-permutations", type=int, default=50,
        help="Number of permutations for validation (default: 50)",
    )
    parser.add_argument(
        "--n-bootstrap", type=int, default=25,
        help="Number of bootstrap iterations (default: 25)",
    )
    args = parser.parse_args()

    success = run_benchmark(
        n_samples=args.n_samples,
        n_pathways=args.n_pathways,
        n_subtypes=args.n_subtypes,
        seed=args.seed,
        skip_validation=args.skip_validation,
        n_permutations=args.n_permutations,
        n_bootstrap=args.n_bootstrap,
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
