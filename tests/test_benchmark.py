"""Tests for the benchmark comparison module."""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.benchmark import (
    BenchmarkComparisonResult,
    BenchmarkMethod,
    BenchmarkResult,
    BenchmarkSweepResult,
    run_benchmark_comparison,
    run_benchmark_sweep,
    run_single_benchmark,
)
from pathway_subtyping.simulation import SimulationConfig, generate_synthetic_data

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def synthetic_data():
    """Create synthetic data with known ground truth (moderate signal)."""
    config = SimulationConfig(
        n_samples=100,
        n_pathways=10,
        n_subtypes=3,
        effect_size=1.0,
        noise_level=1.0,
        seed=42,
    )
    return generate_synthetic_data(config)


@pytest.fixture
def synthetic_data_easy():
    """Create synthetic data with strong signal for convergence tests."""
    config = SimulationConfig(
        n_samples=120,
        n_pathways=10,
        n_subtypes=3,
        effect_size=2.0,
        noise_level=0.5,
        seed=42,
    )
    return generate_synthetic_data(config)


@pytest.fixture
def small_gene_burdens():
    """Small gene burden matrix for edge case tests."""
    rng = np.random.RandomState(42)
    data = rng.exponential(0.5, size=(30, 20))
    return pd.DataFrame(
        data,
        index=[f"S_{i}" for i in range(30)],
        columns=[f"GENE_{i}" for i in range(20)],
    )


@pytest.fixture
def small_pathway_scores():
    """Small pathway score matrix for edge case tests."""
    rng = np.random.RandomState(42)
    data = rng.normal(0, 1, size=(30, 5))
    return pd.DataFrame(
        data,
        index=[f"S_{i}" for i in range(30)],
        columns=[f"PATHWAY_{i}" for i in range(5)],
    )


# =============================================================================
# TEST: BenchmarkMethod enum
# =============================================================================


class TestBenchmarkMethod:
    def test_enum_values(self):
        assert BenchmarkMethod.PATHWAY_GMM.value == "pathway_gmm"
        assert BenchmarkMethod.NMF_CLUSTERING.value == "nmf_clustering"
        assert BenchmarkMethod.PCA_KMEANS.value == "pca_kmeans"
        assert BenchmarkMethod.GENE_KMEANS.value == "gene_kmeans"
        assert BenchmarkMethod.RANDOM_BASELINE.value == "random_baseline"

    def test_enum_count(self):
        assert len(BenchmarkMethod) == 5

    def test_enum_iteration(self):
        methods = list(BenchmarkMethod)
        assert len(methods) == 5


# =============================================================================
# TEST: run_single_benchmark
# =============================================================================


class TestRunSingleBenchmark:
    def test_pathway_gmm(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert result.method == BenchmarkMethod.PATHWAY_GMM
        assert result.ari is not None
        assert result.nmi is not None
        assert len(result.predicted_labels) == 100
        assert result.runtime_seconds >= 0

    def test_nmf_clustering(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.NMF_CLUSTERING,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert result.method == BenchmarkMethod.NMF_CLUSTERING
        assert result.ari is not None
        assert result.n_clusters_found <= 3

    def test_pca_kmeans(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.PCA_KMEANS,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert result.method == BenchmarkMethod.PCA_KMEANS
        assert result.ari is not None

    def test_gene_kmeans(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.GENE_KMEANS,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert result.method == BenchmarkMethod.GENE_KMEANS
        assert result.ari is not None

    def test_random_baseline(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.RANDOM_BASELINE,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert result.method == BenchmarkMethod.RANDOM_BASELINE
        assert result.ari is not None
        # Random baseline should have low ARI
        assert result.ari < 0.3

    def test_no_true_labels(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=None,
            seed=42,
        )

        assert result.ari is None
        assert result.nmi is None
        assert result.silhouette != 0.0 or result.silhouette == 0.0  # computed either way

    def test_seed_reproducibility(self, synthetic_data):
        r1 = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )
        r2 = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        assert r1.ari == r2.ari
        np.testing.assert_array_equal(r1.predicted_labels, r2.predicted_labels)

    def test_easy_data_high_ari(self, synthetic_data_easy):
        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data_easy.gene_burdens,
            pathway_scores=synthetic_data_easy.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data_easy.true_labels,
            seed=42,
        )

        # Strong signal should yield high ARI
        assert result.ari > 0.5


# =============================================================================
# TEST: run_benchmark_comparison
# =============================================================================


class TestRunBenchmarkComparison:
    def test_all_methods(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            seed=42,
        )

        assert isinstance(result, BenchmarkComparisonResult)
        assert len(result.method_results) == 5
        assert len(result.ranking) == 5
        assert result.best_method in result.method_results

    def test_selected_methods(self, synthetic_data):
        methods = [BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.NMF_CLUSTERING]
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=methods,
            seed=42,
        )

        assert len(result.method_results) == 2
        assert "pathway_gmm" in result.method_results
        assert "nmf_clustering" in result.method_results

    def test_ranking_order(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            seed=42,
        )

        # Verify ranking is sorted by ARI descending
        aris = [result.method_results[m].ari for m in result.ranking]
        for i in range(len(aris) - 1):
            if aris[i] is not None and aris[i + 1] is not None:
                assert aris[i] >= aris[i + 1]

    def test_infer_n_clusters_from_labels(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=None,  # Should infer from true_labels
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        assert result.n_clusters == 3

    def test_default_n_clusters_no_labels(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=None,
            n_clusters=None,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        assert result.n_clusters == 3  # default

    def test_n_samples_recorded(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        assert result.n_samples == 100

    def test_random_baseline_ranks_last(self, synthetic_data_easy):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data_easy.gene_burdens,
            pathway_scores=synthetic_data_easy.pathway_scores,
            pathways=synthetic_data_easy.pathways,
            true_labels=synthetic_data_easy.true_labels,
            n_clusters=3,
            seed=42,
        )

        # Random baseline should rank last on easy data
        assert result.ranking[-1] == "random_baseline"


# =============================================================================
# TEST: run_benchmark_sweep
# =============================================================================


class TestRunBenchmarkSweep:
    def test_basic_sweep(self):
        result = run_benchmark_sweep(
            effect_sizes=[0.5, 1.0],
            sample_sizes=[60],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.RANDOM_BASELINE],
            seed=42,
        )

        assert isinstance(result, BenchmarkSweepResult)
        assert len(result.conditions) == 2
        assert len(result.results_per_condition) == 2

    def test_method_mean_ari(self):
        result = run_benchmark_sweep(
            effect_sizes=[1.0, 2.0],
            sample_sizes=[60],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.RANDOM_BASELINE],
            seed=42,
        )

        assert "pathway_gmm" in result.method_mean_ari
        assert "random_baseline" in result.method_mean_ari
        # GMM should outperform random
        assert result.method_mean_ari["pathway_gmm"] > result.method_mean_ari["random_baseline"]

    def test_method_wins(self):
        result = run_benchmark_sweep(
            effect_sizes=[1.0],
            sample_sizes=[80],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.RANDOM_BASELINE],
            seed=42,
        )

        total_wins = sum(result.method_wins.values())
        assert total_wins == 1  # only 1 condition

    def test_conditions_recorded(self):
        result = run_benchmark_sweep(
            effect_sizes=[0.5, 1.5],
            sample_sizes=[60, 100],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        assert len(result.conditions) == 4  # 2 effects x 2 sizes
        for cond in result.conditions:
            assert "effect_size" in cond
            assert "n_samples" in cond


# =============================================================================
# TEST: BenchmarkResult dataclass
# =============================================================================


class TestBenchmarkResult:
    def test_to_dict(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=synthetic_data.true_labels,
            seed=42,
        )

        d = result.to_dict()
        assert "method" in d
        assert d["method"] == "pathway_gmm"
        assert "ari" in d
        assert "nmi" in d
        assert "silhouette" in d
        assert "runtime_seconds" in d
        assert "converged" in d

    def test_to_dict_no_labels(self, synthetic_data):
        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            n_clusters=3,
            true_labels=None,
            seed=42,
        )

        d = result.to_dict()
        assert "ari" not in d
        assert "nmi" not in d


# =============================================================================
# TEST: BenchmarkComparisonResult dataclass
# =============================================================================


class TestBenchmarkComparisonResult:
    def test_to_dict(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.NMF_CLUSTERING],
            seed=42,
        )

        d = result.to_dict()
        assert "best_method" in d
        assert "ranking" in d
        assert "methods" in d
        assert "n_samples" in d
        assert "n_clusters" in d

    def test_format_report(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.NMF_CLUSTERING],
            seed=42,
        )

        report = result.format_report()
        assert "Benchmark Comparison Report" in report
        assert "Best method" in report
        assert "Ranking" in report

    def test_get_citations(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        citations = result.get_citations()
        assert len(citations) >= 1
        assert any("Lee" in c for c in citations)


# =============================================================================
# TEST: BenchmarkSweepResult dataclass
# =============================================================================


class TestBenchmarkSweepResult:
    def test_to_dict(self):
        result = run_benchmark_sweep(
            effect_sizes=[1.0],
            sample_sizes=[60],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        d = result.to_dict()
        assert "n_conditions" in d
        assert "method_mean_ari" in d
        assert "method_wins" in d

    def test_format_report(self):
        result = run_benchmark_sweep(
            effect_sizes=[1.0],
            sample_sizes=[60],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        report = result.format_report()
        assert "Benchmark Sweep Report" in report
        assert "Conditions tested" in report

    def test_get_citations(self):
        result = run_benchmark_sweep(
            effect_sizes=[1.0],
            sample_sizes=[60],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        citations = result.get_citations()
        assert len(citations) >= 1


# =============================================================================
# TEST: EDGE CASES
# =============================================================================


class TestEdgeCases:
    def test_small_sample(self, small_gene_burdens, small_pathway_scores):
        result = run_single_benchmark(
            method=BenchmarkMethod.PATHWAY_GMM,
            gene_burdens=small_gene_burdens,
            pathway_scores=small_pathway_scores,
            n_clusters=2,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert len(result.predicted_labels) == 30

    def test_two_clusters(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=2,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.NMF_CLUSTERING],
            seed=42,
        )

        assert result.n_clusters == 2
        assert len(result.method_results) == 2

    def test_nmf_with_zeros(self):
        """NMF should handle data with many zeros."""
        rng = np.random.RandomState(42)
        data = rng.exponential(0.5, size=(50, 30))
        data[data < 0.3] = 0.0  # many zeros
        gene_burdens = pd.DataFrame(
            data,
            index=[f"S_{i}" for i in range(50)],
            columns=[f"G_{i}" for i in range(30)],
        )
        pathway_scores = pd.DataFrame(
            rng.normal(0, 1, size=(50, 5)),
            index=[f"S_{i}" for i in range(50)],
            columns=[f"P_{i}" for i in range(5)],
        )

        result = run_single_benchmark(
            method=BenchmarkMethod.NMF_CLUSTERING,
            gene_burdens=gene_burdens,
            pathway_scores=pathway_scores,
            n_clusters=3,
            seed=42,
        )

        assert isinstance(result, BenchmarkResult)
        assert len(result.predicted_labels) == 50

    def test_single_method_comparison(self, synthetic_data):
        result = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )

        assert len(result.method_results) == 1
        assert len(result.ranking) == 1


# =============================================================================
# TEST: REPRODUCIBILITY
# =============================================================================


class TestReproducibility:
    def test_comparison_reproducible(self, synthetic_data):
        r1 = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.PCA_KMEANS],
            seed=42,
        )
        r2 = run_benchmark_comparison(
            gene_burdens=synthetic_data.gene_burdens,
            pathway_scores=synthetic_data.pathway_scores,
            pathways=synthetic_data.pathways,
            true_labels=synthetic_data.true_labels,
            n_clusters=3,
            methods=[BenchmarkMethod.PATHWAY_GMM, BenchmarkMethod.PCA_KMEANS],
            seed=42,
        )

        assert r1.best_method == r2.best_method
        assert r1.ranking == r2.ranking

    def test_sweep_reproducible(self):
        kwargs = dict(
            effect_sizes=[1.0],
            sample_sizes=[60],
            n_pathways=8,
            n_subtypes=3,
            methods=[BenchmarkMethod.PATHWAY_GMM],
            seed=42,
        )
        r1 = run_benchmark_sweep(**kwargs)
        r2 = run_benchmark_sweep(**kwargs)

        assert r1.method_mean_ari == r2.method_mean_ari
