"""Tests for multi-omic pathway score fusion module."""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.multi_omic import (
    FusionStrategy,
    MissingStrategy,
    ModalityInput,
    ModalityType,
    MultiOmicFusionResult,
    MultiOmicQualityReport,
    SampleOverlapStats,
    compute_sample_overlap,
    correlation_analysis,
    fuse_modalities,
    prepare_modality,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def vcf_pathway_scores():
    """Simulated VCF-based pathway scores (60 samples, 10 pathways)."""
    rng = np.random.RandomState(42)
    n_samples, n_pathways = 60, 10
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    pathway_names = [f"PATHWAY_{i}" for i in range(1, n_pathways + 1)]
    data = rng.normal(0, 1, (n_samples, n_pathways))
    return pd.DataFrame(data, index=sample_ids, columns=pathway_names)


@pytest.fixture
def expr_pathway_scores():
    """Simulated expression-based pathway scores (50 samples, 10 pathways).

    First 40 overlap with vcf_pathway_scores, last 10 are expression-only.
    """
    rng = np.random.RandomState(43)
    n_pathways = 10
    # 40 shared + 10 expression-only
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, 41)] + [f"EXPR_{i:03d}" for i in range(1, 11)]
    pathway_names = [f"PATHWAY_{i}" for i in range(1, n_pathways + 1)]
    data = rng.normal(0, 1, (50, n_pathways))
    return pd.DataFrame(data, index=sample_ids, columns=pathway_names)


@pytest.fixture
def partial_pathway_scores():
    """Expression scores with different pathways (5 shared + 5 unique)."""
    rng = np.random.RandomState(44)
    n_samples = 40
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    pathway_names = [f"PATHWAY_{i}" for i in range(1, 6)] + [
        f"EXPR_PATHWAY_{i}" for i in range(1, 6)
    ]
    data = rng.normal(0, 1, (n_samples, len(pathway_names)))
    return pd.DataFrame(data, index=sample_ids, columns=pathway_names)


@pytest.fixture
def vcf_gene_burdens():
    """Simulated gene burden data for VCF modality."""
    rng = np.random.RandomState(42)
    n_samples, n_genes = 60, 50
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    gene_names = [f"GENE{i}" for i in range(1, n_genes + 1)]
    data = rng.normal(0, 1, (n_samples, n_genes))
    return pd.DataFrame(data, index=sample_ids, columns=gene_names)


@pytest.fixture
def vcf_modality(vcf_pathway_scores, vcf_gene_burdens):
    """ModalityInput for VCF data."""
    return prepare_modality(
        ModalityType.VCF,
        vcf_pathway_scores,
        gene_data=vcf_gene_burdens,
        label="WES",
    )


@pytest.fixture
def expr_modality(expr_pathway_scores):
    """ModalityInput for expression data."""
    return prepare_modality(
        ModalityType.EXPRESSION,
        expr_pathway_scores,
        label="RNA-seq",
    )


@pytest.fixture
def sc_pathway_scores():
    """Simulated single-cell pathway scores (8 cell types, 10 pathways)."""
    rng = np.random.RandomState(45)
    cell_types = [f"CellType_{i}" for i in range(1, 9)]
    pathway_names = [f"PATHWAY_{i}" for i in range(1, 11)]
    data = rng.normal(0, 1, (8, 10))
    return pd.DataFrame(data, index=cell_types, columns=pathway_names)


# ---------------------------------------------------------------------------
# Enum tests
# ---------------------------------------------------------------------------


class TestModalityType:
    def test_values(self):
        assert ModalityType.VCF.value == "vcf"
        assert ModalityType.EXPRESSION.value == "expression"
        assert ModalityType.SINGLE_CELL.value == "single_cell"

    def test_from_string(self):
        assert ModalityType("vcf") == ModalityType.VCF
        assert ModalityType("expression") == ModalityType.EXPRESSION


class TestFusionStrategy:
    def test_values(self):
        assert FusionStrategy.CONCATENATE.value == "concatenate"
        assert FusionStrategy.WEIGHTED_AVERAGE.value == "weighted_average"
        assert FusionStrategy.INTERSECTION_ONLY.value == "intersection_only"


class TestMissingStrategy:
    def test_values(self):
        assert MissingStrategy.IMPUTE_ZERO.value == "impute_zero"
        assert MissingStrategy.IMPUTE_MEAN.value == "impute_mean"
        assert MissingStrategy.DROP.value == "drop"


# ---------------------------------------------------------------------------
# prepare_modality tests
# ---------------------------------------------------------------------------


class TestPrepareModality:
    def test_basic(self, vcf_pathway_scores):
        m = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="WES")
        assert m.modality_type == ModalityType.VCF
        assert m.label == "WES"
        assert m.n_samples == 60
        assert m.n_pathways == 10

    def test_with_gene_data(self, vcf_pathway_scores, vcf_gene_burdens):
        m = prepare_modality(ModalityType.VCF, vcf_pathway_scores, gene_data=vcf_gene_burdens)
        assert m.gene_data is not None
        assert m.gene_data.shape == (60, 50)

    def test_properties(self, vcf_pathway_scores):
        m = prepare_modality(ModalityType.VCF, vcf_pathway_scores)
        assert len(m.sample_ids) == 60
        assert len(m.pathway_names) == 10
        assert m.sample_ids[0] == "SAMPLE_001"
        assert m.pathway_names[0] == "PATHWAY_1"

    def test_to_dict(self, vcf_pathway_scores):
        m = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="WES")
        d = m.to_dict()
        assert d["modality_type"] == "vcf"
        assert d["label"] == "WES"
        assert d["n_samples"] == 60
        assert d["n_pathways"] == 10

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty or None"):
            prepare_modality(ModalityType.VCF, pd.DataFrame())

    def test_none_raises(self):
        with pytest.raises(ValueError, match="empty or None"):
            prepare_modality(ModalityType.VCF, None)

    def test_single_sample_raises(self):
        scores = pd.DataFrame({"P1": [1.0]}, index=["S1"])
        with pytest.raises(ValueError, match="at least 2"):
            prepare_modality(ModalityType.VCF, scores)

    def test_default_label(self, vcf_pathway_scores):
        m = prepare_modality(ModalityType.VCF, vcf_pathway_scores)
        assert m.label == "vcf"


# ---------------------------------------------------------------------------
# compute_sample_overlap tests
# ---------------------------------------------------------------------------


class TestSampleOverlap:
    def test_full_overlap(self, vcf_pathway_scores):
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, vcf_pathway_scores.copy(), label="B")
        overlap = compute_sample_overlap([m1, m2])
        assert overlap.n_shared == 60
        assert overlap.n_union == 60
        assert overlap.overlap_fraction == 1.0

    def test_partial_overlap(self, vcf_modality, expr_modality):
        overlap = compute_sample_overlap([vcf_modality, expr_modality])
        assert overlap.n_shared == 40  # SAMPLE_001..040
        assert overlap.n_union == 70  # 60 VCF + 10 EXPR-only
        assert 0.5 < overlap.overlap_fraction < 0.6

    def test_no_overlap(self, vcf_pathway_scores):
        rng = np.random.RandomState(99)
        disjoint = pd.DataFrame(
            rng.normal(0, 1, (10, 10)),
            index=[f"OTHER_{i}" for i in range(10)],
            columns=vcf_pathway_scores.columns,
        )
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, disjoint, label="B")
        overlap = compute_sample_overlap([m1, m2])
        assert overlap.n_shared == 0
        assert overlap.overlap_fraction == 0.0

    def test_three_modalities(self, vcf_modality, expr_modality, sc_pathway_scores):
        sc_mod = prepare_modality(ModalityType.SINGLE_CELL, sc_pathway_scores, label="SC")
        overlap = compute_sample_overlap([vcf_modality, expr_modality, sc_mod])
        assert overlap.n_modalities == 3
        assert overlap.n_shared == 0  # SC cell types don't match sample IDs
        assert len(overlap.pairwise_overlap) == 3

    def test_to_dict(self, vcf_modality, expr_modality):
        overlap = compute_sample_overlap([vcf_modality, expr_modality])
        d = overlap.to_dict()
        assert "n_modalities" in d
        assert "n_shared_samples" in d
        assert "overlap_fraction" in d
        assert "pairwise_overlap" in d


# ---------------------------------------------------------------------------
# fuse_modalities — CONCATENATE tests
# ---------------------------------------------------------------------------


class TestFuseModalitiesConcatenate:
    def test_basic(self, vcf_modality, expr_modality):
        result = fuse_modalities([vcf_modality, expr_modality], strategy=FusionStrategy.CONCATENATE)
        assert isinstance(result, MultiOmicFusionResult)
        assert result.strategy == FusionStrategy.CONCATENATE
        assert result.n_modalities == 2

    def test_prefixed_columns(self, vcf_modality, expr_modality):
        result = fuse_modalities([vcf_modality, expr_modality], strategy=FusionStrategy.CONCATENATE)
        cols = list(result.fused_pathway_scores.columns)
        assert any(c.startswith("WES:") for c in cols)
        assert any(c.startswith("RNA-seq:") for c in cols)
        # 10 pathways per modality = 20 total
        assert len(cols) == 20

    def test_partial_samples_impute_zero(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.IMPUTE_ZERO,
            renormalize=False,
        )
        # Union: 60 + 10 expr-only = 70
        assert len(result.fused_pathway_scores) == 70
        # Check that VCF-only samples have 0 for RNA-seq columns (before renorm)
        vcf_only = [s for s in result.fused_pathway_scores.index if s.startswith("SAMPLE_04")]
        rnaseq_cols = [c for c in result.fused_pathway_scores.columns if c.startswith("RNA-seq:")]
        for s in vcf_only:
            if s not in expr_modality.sample_ids:
                for c in rnaseq_cols:
                    assert result.fused_pathway_scores.loc[s, c] == 0.0

    def test_partial_samples_impute_mean(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.IMPUTE_MEAN,
            renormalize=False,
        )
        assert len(result.fused_pathway_scores) == 70
        # No NaN values
        assert not result.fused_pathway_scores.isna().any().any()

    def test_drop_missing(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.DROP,
        )
        # Intersection: 40 shared samples
        assert len(result.fused_pathway_scores) == 40

    def test_preserves_all_pathways(self, vcf_pathway_scores, partial_pathway_scores):
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="VCF")
        m2 = prepare_modality(ModalityType.EXPRESSION, partial_pathway_scores, label="EXPR")
        result = fuse_modalities([m1, m2], strategy=FusionStrategy.CONCATENATE)
        # 10 VCF + 10 EXPR pathways = 20 (with prefixes)
        assert len(result.fused_pathway_scores.columns) == 20

    def test_three_modalities(self, vcf_modality, expr_modality, sc_pathway_scores):
        sc_mod = prepare_modality(ModalityType.SINGLE_CELL, sc_pathway_scores, label="SC")
        result = fuse_modalities(
            [vcf_modality, expr_modality, sc_mod],
            strategy=FusionStrategy.CONCATENATE,
        )
        assert result.n_modalities == 3
        # 10 + 10 + 10 = 30 columns
        assert len(result.fused_pathway_scores.columns) == 30

    def test_with_renormalize(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.DROP,
            renormalize=True,
        )
        # After Z-normalization, mean should be ~0
        means = result.fused_pathway_scores.mean()
        assert all(abs(m) < 0.1 for m in means)

    def test_without_renormalize(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.DROP,
            renormalize=False,
        )
        assert isinstance(result.fused_pathway_scores, pd.DataFrame)
        assert len(result.fused_pathway_scores) == 40


# ---------------------------------------------------------------------------
# fuse_modalities — WEIGHTED_AVERAGE tests
# ---------------------------------------------------------------------------


class TestFuseModalitiesWeightedAverage:
    def test_equal_weights(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.WEIGHTED_AVERAGE,
            missing_strategy=MissingStrategy.DROP,
        )
        assert result.strategy == FusionStrategy.WEIGHTED_AVERAGE
        # All 10 pathways shared
        assert len(result.fused_pathway_scores.columns) == 10
        # Only shared samples (intersection)
        assert len(result.fused_pathway_scores) == 40
        # Equal weights
        assert abs(result.modality_weights["WES"] - 0.5) < 0.01
        assert abs(result.modality_weights["RNA-seq"] - 0.5) < 0.01

    def test_custom_weights(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.WEIGHTED_AVERAGE,
            weights={"WES": 0.7, "RNA-seq": 0.3},
            missing_strategy=MissingStrategy.DROP,
        )
        assert abs(result.modality_weights["WES"] - 0.7) < 0.01
        assert abs(result.modality_weights["RNA-seq"] - 0.3) < 0.01

    def test_shared_pathways_only(self, vcf_pathway_scores, partial_pathway_scores):
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="VCF")
        m2 = prepare_modality(ModalityType.EXPRESSION, partial_pathway_scores, label="EXPR")
        result = fuse_modalities(
            [m1, m2],
            strategy=FusionStrategy.WEIGHTED_AVERAGE,
            missing_strategy=MissingStrategy.DROP,
        )
        # Only 5 shared pathways (PATHWAY_1..5)
        assert len(result.fused_pathway_scores.columns) == 5
        assert all(c.startswith("PATHWAY_") for c in result.fused_pathway_scores.columns)

    def test_renormalize(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.WEIGHTED_AVERAGE,
            missing_strategy=MissingStrategy.DROP,
            renormalize=True,
        )
        means = result.fused_pathway_scores.mean()
        assert all(abs(m) < 0.1 for m in means)

    def test_no_shared_pathways_raises(self, vcf_pathway_scores):
        rng = np.random.RandomState(99)
        disjoint = pd.DataFrame(
            rng.normal(0, 1, (60, 5)),
            index=vcf_pathway_scores.index,
            columns=[f"OTHER_PATHWAY_{i}" for i in range(5)],
        )
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, disjoint, label="B")
        with pytest.raises(ValueError, match="No shared pathways"):
            fuse_modalities([m1, m2], strategy=FusionStrategy.WEIGHTED_AVERAGE)

    def test_invalid_weights_raises(self, vcf_modality, expr_modality):
        with pytest.raises(ValueError, match="sum to 1.0"):
            fuse_modalities(
                [vcf_modality, expr_modality],
                strategy=FusionStrategy.WEIGHTED_AVERAGE,
                weights={"WES": 0.5, "RNA-seq": 0.8},
            )

    def test_missing_weight_raises(self, vcf_modality, expr_modality):
        with pytest.raises(ValueError, match="Weight missing"):
            fuse_modalities(
                [vcf_modality, expr_modality],
                strategy=FusionStrategy.WEIGHTED_AVERAGE,
                weights={"WES": 1.0},
            )

    def test_partial_sample_overlap(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.WEIGHTED_AVERAGE,
            missing_strategy=MissingStrategy.IMPUTE_ZERO,
        )
        # Union: 70 samples
        assert len(result.fused_pathway_scores) == 70


# ---------------------------------------------------------------------------
# fuse_modalities — INTERSECTION_ONLY tests
# ---------------------------------------------------------------------------


class TestFuseModalitiesIntersectionOnly:
    def test_basic(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.INTERSECTION_ONLY,
        )
        assert result.strategy == FusionStrategy.INTERSECTION_ONLY
        assert len(result.fused_pathway_scores) == 40  # shared samples
        assert len(result.fused_pathway_scores.columns) == 10  # shared pathways

    def test_partial_pathways(self, vcf_pathway_scores, partial_pathway_scores):
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="VCF")
        m2 = prepare_modality(ModalityType.EXPRESSION, partial_pathway_scores, label="EXPR")
        result = fuse_modalities([m1, m2], strategy=FusionStrategy.INTERSECTION_ONLY)
        # 5 shared pathways, 40 shared samples
        assert len(result.fused_pathway_scores.columns) == 5
        assert len(result.fused_pathway_scores) == 40

    def test_no_shared_samples_raises(self, vcf_pathway_scores):
        rng = np.random.RandomState(99)
        disjoint = pd.DataFrame(
            rng.normal(0, 1, (10, 10)),
            index=[f"OTHER_{i}" for i in range(10)],
            columns=vcf_pathway_scores.columns,
        )
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, disjoint, label="B")
        with pytest.raises(ValueError, match="shared sample"):
            fuse_modalities([m1, m2], strategy=FusionStrategy.INTERSECTION_ONLY)

    def test_no_shared_pathways_raises(self, vcf_pathway_scores):
        rng = np.random.RandomState(99)
        disjoint = pd.DataFrame(
            rng.normal(0, 1, (60, 5)),
            index=vcf_pathway_scores.index,
            columns=[f"OTHER_{i}" for i in range(5)],
        )
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, disjoint, label="B")
        with pytest.raises(ValueError, match="shared pathway"):
            fuse_modalities([m1, m2], strategy=FusionStrategy.INTERSECTION_ONLY)


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


class TestFuseModalitiesEdgeCases:
    def test_single_modality_raises(self, vcf_modality):
        with pytest.raises(ValueError, match="at least 2"):
            fuse_modalities([vcf_modality])

    def test_gene_data_fusion(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.DROP,
        )
        # VCF has gene_data, expression doesn't
        assert result.fused_gene_data is not None
        assert len(result.fused_gene_data) == 40  # shared samples
        assert any(c.startswith("WES:") for c in result.fused_gene_data.columns)

    def test_seed_reproducibility(self, vcf_modality, expr_modality):
        r1 = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            seed=42,
        )
        r2 = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            seed=42,
        )
        pd.testing.assert_frame_equal(r1.fused_pathway_scores, r2.fused_pathway_scores)

    def test_no_nan_in_output(self, vcf_modality, expr_modality):
        for strategy in FusionStrategy:
            if strategy == FusionStrategy.INTERSECTION_ONLY:
                result = fuse_modalities([vcf_modality, expr_modality], strategy=strategy)
            else:
                result = fuse_modalities(
                    [vcf_modality, expr_modality],
                    strategy=strategy,
                    missing_strategy=MissingStrategy.IMPUTE_ZERO,
                )
            assert (
                not result.fused_pathway_scores.isna().any().any()
            ), f"NaN found in {strategy.value} output"


# ---------------------------------------------------------------------------
# MultiOmicFusionResult tests
# ---------------------------------------------------------------------------


class TestMultiOmicFusionResult:
    def test_to_dict(self, vcf_modality, expr_modality):
        result = fuse_modalities([vcf_modality, expr_modality])
        d = result.to_dict()
        assert d["strategy"] == "concatenate"
        assert d["n_modalities"] == 2
        assert "fused_shape" in d
        assert d["fused_shape"]["n_samples"] > 0

    def test_format_report(self, vcf_modality, expr_modality):
        result = fuse_modalities([vcf_modality, expr_modality])
        report = result.format_report()
        assert "Multi-Omic Fusion Report" in report
        assert "WES" in report
        assert "RNA-seq" in report
        assert "Sample Overlap" in report

    def test_get_citations_concatenate(self, vcf_modality, expr_modality):
        result = fuse_modalities([vcf_modality, expr_modality], strategy=FusionStrategy.CONCATENATE)
        citations = result.get_citations()
        assert len(citations) >= 1
        assert "Ritchie" in citations[0]

    def test_get_citations_weighted(self, vcf_modality, expr_modality):
        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.WEIGHTED_AVERAGE,
            missing_strategy=MissingStrategy.DROP,
        )
        citations = result.get_citations()
        assert len(citations) >= 2
        assert any("Speicher" in c for c in citations)


# ---------------------------------------------------------------------------
# correlation_analysis tests
# ---------------------------------------------------------------------------


class TestCorrelationAnalysis:
    def test_pearson(self, vcf_modality, expr_modality):
        df = correlation_analysis([vcf_modality, expr_modality], method="pearson")
        assert isinstance(df, pd.DataFrame)
        assert "correlation" in df.columns
        assert "p_value" in df.columns
        assert "pathway" in df.columns
        assert len(df) == 10  # 10 shared pathways

    def test_spearman(self, vcf_modality, expr_modality):
        df = correlation_analysis([vcf_modality, expr_modality], method="spearman")
        assert len(df) == 10
        assert all(df["p_value"] >= 0)

    def test_no_shared_pathways(self, vcf_pathway_scores):
        rng = np.random.RandomState(99)
        disjoint = pd.DataFrame(
            rng.normal(0, 1, (60, 5)),
            index=vcf_pathway_scores.index,
            columns=[f"OTHER_{i}" for i in range(5)],
        )
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, disjoint, label="B")
        df = correlation_analysis([m1, m2])
        assert len(df) == 0

    def test_three_modalities(self, vcf_modality, expr_modality, sc_pathway_scores):
        sc_mod = prepare_modality(ModalityType.SINGLE_CELL, sc_pathway_scores, label="SC")
        df = correlation_analysis([vcf_modality, expr_modality, sc_mod])
        # Only WES vs RNA-seq have shared samples AND shared pathways
        modality_pairs = set()
        for _, row in df.iterrows():
            modality_pairs.add((row["modality_a"], row["modality_b"]))
        assert ("WES", "RNA-seq") in modality_pairs

    def test_invalid_method_raises(self, vcf_modality, expr_modality):
        with pytest.raises(ValueError, match="pearson.*spearman"):
            correlation_analysis([vcf_modality, expr_modality], method="kendall")


# ---------------------------------------------------------------------------
# MultiOmicQualityReport tests
# ---------------------------------------------------------------------------


class TestMultiOmicQualityReport:
    def test_to_dict(self, vcf_modality, expr_modality):
        result = fuse_modalities([vcf_modality, expr_modality])
        d = result.quality_report.to_dict()
        assert "per_modality_reports" in d
        assert "sample_overlap" in d
        assert "pathway_overlap_count" in d
        assert "is_usable" in d

    def test_warnings_low_overlap(self, vcf_pathway_scores):
        rng = np.random.RandomState(99)
        # Only 5 shared out of 65 total
        few_shared = pd.DataFrame(
            rng.normal(0, 1, (10, 10)),
            index=[f"SAMPLE_{i:03d}" for i in range(1, 6)] + [f"OTHER_{i}" for i in range(5)],
            columns=vcf_pathway_scores.columns,
        )
        m1 = prepare_modality(ModalityType.VCF, vcf_pathway_scores, label="A")
        m2 = prepare_modality(ModalityType.EXPRESSION, few_shared, label="B")
        result = fuse_modalities([m1, m2], strategy=FusionStrategy.CONCATENATE)
        assert any("overlap" in w.lower() for w in result.quality_report.warnings)


# ---------------------------------------------------------------------------
# Downstream compatibility tests
# ---------------------------------------------------------------------------


class TestDownstreamCompatibility:
    def test_clustering_compatible(self, vcf_modality, expr_modality):
        """Fused output works with run_clustering."""
        from pathway_subtyping.clustering import ClusteringAlgorithm, run_clustering

        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.DROP,
        )
        data = result.fused_pathway_scores.values
        cluster_result = run_clustering(
            data, n_clusters=3, algorithm=ClusteringAlgorithm.GMM, seed=42
        )
        assert len(cluster_result.labels) == len(result.fused_pathway_scores)
        assert cluster_result.n_clusters == 3

    def test_characterization_compatible(self, vcf_modality, expr_modality):
        """Fused output works with characterize_subtypes."""
        from pathway_subtyping.characterization import characterize_subtypes
        from pathway_subtyping.clustering import ClusteringAlgorithm, run_clustering

        result = fuse_modalities(
            [vcf_modality, expr_modality],
            strategy=FusionStrategy.CONCATENATE,
            missing_strategy=MissingStrategy.DROP,
        )
        data = result.fused_pathway_scores.values
        cluster_result = run_clustering(
            data, n_clusters=3, algorithm=ClusteringAlgorithm.GMM, seed=42
        )
        char_result = characterize_subtypes(
            pathway_scores=result.fused_pathway_scores,
            cluster_labels=cluster_result.labels,
            seed=42,
        )
        assert char_result.n_subtypes == 3
        assert char_result.n_samples == len(result.fused_pathway_scores)
