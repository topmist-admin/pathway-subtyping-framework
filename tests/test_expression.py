"""
Tests for the expression pathway scoring module.
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.expression import (
    ExpressionDataQualityReport,
    ExpressionInputType,
    ExpressionScoringMethod,
    ExpressionScoringResult,
    _detect_orientation,
    _looks_like_gene_symbols,
    _preprocess_expression,
    _score_gsva,
    _score_mean_z,
    _score_ssgsea,
    load_expression_matrix,
    score_pathways_from_expression,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def expression_matrix():
    """Generate a sample gene expression matrix (samples x genes), log2-scale."""
    rng = np.random.RandomState(42)
    n_samples = 60
    n_genes = 100
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    gene_names = [f"GENE{i}" for i in range(1, n_genes + 1)]
    expression = rng.normal(6.0, 2.0, (n_samples, n_genes))
    expression = np.clip(expression, 0, None)
    return pd.DataFrame(expression, index=sample_ids, columns=gene_names)


@pytest.fixture
def expression_pathways():
    """Pathway definitions mapping to GENE1..GENE100."""
    return {
        "PATHWAY_A": [f"GENE{i}" for i in range(1, 21)],
        "PATHWAY_B": [f"GENE{i}" for i in range(21, 41)],
        "PATHWAY_C": [f"GENE{i}" for i in range(41, 61)],
        "PATHWAY_D": [f"GENE{i}" for i in range(61, 81)],
        "PATHWAY_E": [f"GENE{i}" for i in range(81, 101)],
    }


@pytest.fixture
def enriched_expression(expression_pathways):
    """Expression matrix with clear pathway enrichment patterns."""
    rng = np.random.RandomState(42)
    n_per = 20
    n_samples = n_per * 3
    all_genes = []
    for genes in expression_pathways.values():
        all_genes.extend(genes)
    all_genes = sorted(set(all_genes))

    sample_ids = [f"S_{i:03d}" for i in range(n_samples)]
    data = rng.normal(6.0, 1.0, (n_samples, len(all_genes)))
    df = pd.DataFrame(data, index=sample_ids, columns=all_genes)

    # Cluster 0: upregulate PATHWAY_A genes
    for g in expression_pathways["PATHWAY_A"]:
        df.loc[sample_ids[:n_per], g] += 4.0

    # Cluster 1: upregulate PATHWAY_B genes
    for g in expression_pathways["PATHWAY_B"]:
        df.loc[sample_ids[n_per : 2 * n_per], g] += 4.0

    # Cluster 2: upregulate PATHWAY_C genes
    for g in expression_pathways["PATHWAY_C"]:
        df.loc[sample_ids[2 * n_per :], g] += 4.0

    return df


@pytest.fixture
def expression_csv_path(expression_matrix, tmp_path):
    """Write expression matrix to CSV and return path."""
    path = tmp_path / "expression.csv"
    expression_matrix.to_csv(path)
    return str(path)


@pytest.fixture
def expression_tsv_path(expression_matrix, tmp_path):
    """Write expression matrix to TSV and return path."""
    path = tmp_path / "expression.tsv"
    expression_matrix.to_csv(path, sep="\t")
    return str(path)


@pytest.fixture
def transposed_csv_path(expression_matrix, tmp_path):
    """Write genes-as-rows expression matrix to CSV."""
    path = tmp_path / "expression_transposed.csv"
    expression_matrix.T.to_csv(path)
    return str(path)


@pytest.fixture
def counts_csv_path(tmp_path):
    """Write a raw counts matrix to CSV."""
    rng = np.random.RandomState(42)
    n_samples = 20
    n_genes = 50
    sample_ids = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]
    gene_names = [f"GENE{i}" for i in range(1, n_genes + 1)]
    counts = rng.poisson(lam=100, size=(n_samples, n_genes))
    df = pd.DataFrame(counts, index=sample_ids, columns=gene_names)
    path = tmp_path / "counts.csv"
    df.to_csv(path)
    return str(path)


# ---------------------------------------------------------------------------
# Enum Tests
# ---------------------------------------------------------------------------


class TestExpressionScoringMethod:
    def test_enum_values(self):
        assert ExpressionScoringMethod.MEAN_Z.value == "mean_z"
        assert ExpressionScoringMethod.SSGSEA.value == "ssgsea"
        assert ExpressionScoringMethod.GSVA.value == "gsva"

    def test_enum_from_string(self):
        assert ExpressionScoringMethod("mean_z") == ExpressionScoringMethod.MEAN_Z
        assert ExpressionScoringMethod("ssgsea") == ExpressionScoringMethod.SSGSEA
        assert ExpressionScoringMethod("gsva") == ExpressionScoringMethod.GSVA

    def test_invalid_value(self):
        with pytest.raises(ValueError):
            ExpressionScoringMethod("invalid")


class TestExpressionInputType:
    def test_enum_values(self):
        assert ExpressionInputType.COUNTS.value == "counts"
        assert ExpressionInputType.TPM.value == "tpm"
        assert ExpressionInputType.FPKM.value == "fpkm"
        assert ExpressionInputType.LOG2.value == "log2"

    def test_enum_from_string(self):
        assert ExpressionInputType("counts") == ExpressionInputType.COUNTS
        assert ExpressionInputType("log2") == ExpressionInputType.LOG2


# ---------------------------------------------------------------------------
# Orientation Detection Tests
# ---------------------------------------------------------------------------


class TestOrientationDetection:
    def test_genes_as_columns_detected(self, expression_matrix):
        result = _detect_orientation(expression_matrix)
        assert result == "genes_as_columns"

    def test_genes_as_rows_detected(self, expression_matrix):
        transposed = expression_matrix.T
        result = _detect_orientation(transposed)
        assert result == "genes_as_rows"

    def test_explicit_gene_column(self, expression_matrix):
        # Add a gene_name column
        df = expression_matrix.T.reset_index()
        df = df.rename(columns={"index": "gene_name"})
        result = _detect_orientation(df, gene_column="gene_name")
        assert result == "genes_as_rows"

    def test_missing_gene_column_fallback(self, expression_matrix):
        # Specified column doesn't exist — falls back to auto-detection
        result = _detect_orientation(expression_matrix, gene_column="nonexistent")
        assert result in ("genes_as_columns", "genes_as_rows")

    def test_looks_like_gene_symbols(self):
        assert _looks_like_gene_symbols(["SHANK3", "CHD8", "NRXN1"])
        assert _looks_like_gene_symbols(["TP53", "BRCA1", "EGFR"])
        assert not _looks_like_gene_symbols(["sample_001", "sample_002"])
        assert not _looks_like_gene_symbols(["1", "2", "3"])
        assert not _looks_like_gene_symbols([])


# ---------------------------------------------------------------------------
# Preprocessing Tests
# ---------------------------------------------------------------------------


class TestPreprocessExpression:
    def test_counts_log_transform(self):
        df = pd.DataFrame({"G1": [0, 100, 1000], "G2": [50, 200, 0]})
        result = _preprocess_expression(df, ExpressionInputType.COUNTS)
        # log2(100 + 1) ≈ 6.66
        assert result.loc[1, "G1"] == pytest.approx(np.log2(101), rel=1e-6)
        # log2(0 + 1) = 0
        assert result.loc[0, "G1"] == pytest.approx(0.0)

    def test_log2_no_transform(self):
        df = pd.DataFrame({"G1": [5.0, 6.0, 7.0]})
        result = _preprocess_expression(df, ExpressionInputType.LOG2)
        pd.testing.assert_frame_equal(result, df)

    def test_tpm_high_values_transformed(self):
        # Values > 20 trigger log transform
        df = pd.DataFrame({"G1": [100.0, 200.0, 300.0]})
        result = _preprocess_expression(df, ExpressionInputType.TPM)
        assert result.loc[0, "G1"] == pytest.approx(np.log2(101), rel=1e-6)

    def test_tpm_low_values_not_transformed(self):
        # Values <= 20 are assumed already log-scaled
        df = pd.DataFrame({"G1": [5.0, 8.0, 12.0]})
        result = _preprocess_expression(df, ExpressionInputType.TPM)
        pd.testing.assert_frame_equal(result, df)


# ---------------------------------------------------------------------------
# Load Expression Matrix Tests
# ---------------------------------------------------------------------------


class TestLoadExpressionMatrix:
    def test_load_genes_as_columns(self, expression_csv_path):
        df, report = load_expression_matrix(
            expression_csv_path, input_type=ExpressionInputType.LOG2
        )
        assert df.shape[0] == 60
        assert df.shape[1] > 0
        assert report.is_usable
        assert not report.was_transposed

    def test_load_genes_as_rows(self, transposed_csv_path):
        df, report = load_expression_matrix(
            transposed_csv_path, input_type=ExpressionInputType.LOG2
        )
        assert df.shape[0] == 60  # samples
        assert report.was_transposed

    def test_load_tsv(self, expression_tsv_path):
        df, report = load_expression_matrix(
            expression_tsv_path, input_type=ExpressionInputType.LOG2
        )
        assert df.shape[0] == 60
        assert report.is_usable

    def test_load_counts(self, counts_csv_path):
        df, report = load_expression_matrix(
            counts_csv_path, input_type=ExpressionInputType.COUNTS
        )
        # Values should be log-transformed
        assert df.values.max() < 20  # log2(counts) should be much smaller than raw

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_expression_matrix("/nonexistent/path.csv")

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty.csv"
        path.write_text("")
        with pytest.raises(Exception):
            load_expression_matrix(str(path))

    def test_quality_report_populated(self, expression_csv_path):
        _, report = load_expression_matrix(
            expression_csv_path, input_type=ExpressionInputType.LOG2
        )
        assert report.n_samples == 60
        assert report.n_genes > 0
        assert report.n_genes_before_filter == 100
        assert isinstance(report.orientation_detected, str)
        assert report.input_type == "log2"

    def test_quality_report_to_dict(self, expression_csv_path):
        _, report = load_expression_matrix(
            expression_csv_path, input_type=ExpressionInputType.LOG2
        )
        d = report.to_dict()
        assert "n_samples" in d
        assert "n_genes" in d
        assert "is_usable" in d
        assert isinstance(d["warnings"], list)

    def test_zero_gene_removal(self, tmp_path):
        """Genes with all-zero expression should be removed."""
        df = pd.DataFrame(
            {
                "G1": [1.0, 2.0, 3.0],
                "G2": [0.0, 0.0, 0.0],
                "G3": [4.0, 5.0, 6.0],
            },
            index=["S1", "S2", "S3"],
        )
        path = tmp_path / "test.csv"
        df.to_csv(path)
        result, report = load_expression_matrix(
            str(path),
            input_type=ExpressionInputType.LOG2,
            min_genes_per_sample=1,
        )
        assert "G2" not in result.columns
        assert report.n_zero_genes == 1


# ---------------------------------------------------------------------------
# Mean-Z Scoring Tests
# ---------------------------------------------------------------------------


class TestMeanZScoring:
    def test_basic_scoring(self, expression_matrix, expression_pathways):
        scores, skipped = _score_mean_z(expression_matrix, expression_pathways)
        assert scores.shape[0] == 60
        assert scores.shape[1] == 5
        assert len(skipped) == 0

    def test_known_values(self):
        """Manual verification of mean-Z calculation."""
        df = pd.DataFrame(
            {
                "G1": [10.0, 0.0, 0.0],
                "G2": [10.0, 0.0, 0.0],
                "G3": [0.0, 10.0, 0.0],
            },
            index=["S1", "S2", "S3"],
        )
        pathways = {"P1": ["G1", "G2"], "P2": ["G3"]}
        scores, skipped = _score_mean_z(df, pathways, min_genes=1)

        # S1 should have high P1 score (both genes high)
        assert scores.loc["S1", "P1"] > scores.loc["S2", "P1"]

    def test_skip_small_pathways(self, expression_matrix):
        pathways = {"SMALL": ["GENE1"]}  # Only 1 gene
        scores, skipped = _score_mean_z(expression_matrix, pathways, min_genes=2)
        assert "SMALL" in skipped
        assert scores.shape[1] == 0

    def test_no_overlap(self, expression_matrix):
        pathways = {"MISSING": ["NONEXISTENT1", "NONEXISTENT2"]}
        scores, skipped = _score_mean_z(expression_matrix, pathways)
        assert "MISSING" in skipped

    def test_zero_variance_gene_handling(self):
        """All-same expression should not crash."""
        df = pd.DataFrame(
            {"G1": [5.0, 5.0, 5.0], "G2": [5.0, 5.0, 5.0]},
            index=["S1", "S2", "S3"],
        )
        pathways = {"P1": ["G1", "G2"]}
        scores, skipped = _score_mean_z(df, pathways)
        assert scores.shape == (3, 1)


# ---------------------------------------------------------------------------
# ssGSEA Scoring Tests
# ---------------------------------------------------------------------------


class TestSsGSEAScoring:
    def test_basic_scoring(self, expression_matrix, expression_pathways):
        scores, skipped = _score_ssgsea(expression_matrix, expression_pathways)
        assert scores.shape[0] == 60
        assert scores.shape[1] == 5
        assert len(skipped) == 0

    def test_known_enrichment(self, enriched_expression, expression_pathways):
        """Pathway genes highly expressed should yield high scores."""
        scores, _ = _score_ssgsea(enriched_expression, expression_pathways)

        # Cluster 0 samples (first 20) should have highest PATHWAY_A score
        c0_mean = scores.iloc[:20]["PATHWAY_A"].mean()
        c1_mean = scores.iloc[20:40]["PATHWAY_A"].mean()
        c2_mean = scores.iloc[40:60]["PATHWAY_A"].mean()
        assert c0_mean > c1_mean
        assert c0_mean > c2_mean

    def test_alpha_parameter(self, expression_matrix, expression_pathways):
        scores_025, _ = _score_ssgsea(
            expression_matrix, expression_pathways, alpha=0.25
        )
        scores_100, _ = _score_ssgsea(
            expression_matrix, expression_pathways, alpha=1.0
        )
        # Different alpha should produce different scores
        assert not np.allclose(scores_025.values, scores_100.values)

    def test_skip_small_pathways(self, expression_matrix):
        pathways = {"TINY": ["GENE1"]}
        _, skipped = _score_ssgsea(expression_matrix, pathways, min_genes=2)
        assert "TINY" in skipped

    def test_all_genes_in_pathway(self, expression_matrix):
        """Edge case: pathway contains all genes in expression matrix."""
        all_genes = list(expression_matrix.columns)
        pathways = {"ALL": all_genes}
        scores, skipped = _score_ssgsea(expression_matrix, pathways)
        # Should handle gracefully (all zeros since no non-pathway genes)
        assert scores.shape[1] == 1
        assert np.allclose(scores["ALL"].values, 0.0)


# ---------------------------------------------------------------------------
# GSVA Scoring Tests
# ---------------------------------------------------------------------------


class TestGSVAScoring:
    def test_basic_scoring(self, expression_matrix, expression_pathways):
        scores, skipped = _score_gsva(expression_matrix, expression_pathways)
        assert scores.shape[0] == 60
        assert scores.shape[1] == 5
        assert len(skipped) == 0

    def test_known_enrichment(self, enriched_expression, expression_pathways):
        """Pathway genes highly expressed should yield high GSVA scores."""
        scores, _ = _score_gsva(enriched_expression, expression_pathways)

        c0_mean = scores.iloc[:20]["PATHWAY_A"].mean()
        c1_mean = scores.iloc[20:40]["PATHWAY_A"].mean()
        assert c0_mean > c1_mean

    def test_skip_small_pathways(self, expression_matrix):
        pathways = {"TINY": ["GENE1"]}
        _, skipped = _score_gsva(expression_matrix, pathways, min_genes=2)
        assert "TINY" in skipped

    def test_output_shape(self, expression_matrix, expression_pathways):
        scores, _ = _score_gsva(expression_matrix, expression_pathways)
        assert scores.index.tolist() == expression_matrix.index.tolist()


# ---------------------------------------------------------------------------
# Main Scoring Function Tests
# ---------------------------------------------------------------------------


class TestScorePathwaysFromExpression:
    def test_mean_z_method(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.MEAN_Z,
        )
        assert isinstance(result, ExpressionScoringResult)
        assert result.n_pathways_scored == 5
        assert result.n_pathways_skipped == 0
        assert result.method == ExpressionScoringMethod.MEAN_Z

    def test_ssgsea_method(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.SSGSEA,
        )
        assert result.n_pathways_scored == 5
        assert result.method == ExpressionScoringMethod.SSGSEA

    def test_gsva_method(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.GSVA,
        )
        assert result.n_pathways_scored == 5
        assert result.method == ExpressionScoringMethod.GSVA

    def test_output_is_z_normalized(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.MEAN_Z,
        )
        # Mean should be approximately 0
        means = result.pathway_scores.mean()
        for m in means:
            assert abs(m) < 0.01

        # Std should be approximately 1
        stds = result.pathway_scores.std()
        for s in stds:
            assert abs(s - 1.0) < 0.01

    def test_min_genes_filter(self, expression_matrix):
        pathways = {
            "GOOD": [f"GENE{i}" for i in range(1, 11)],
            "BAD": ["GENE1"],
        }
        result = score_pathways_from_expression(
            expression_matrix,
            pathways,
            min_genes_per_pathway=2,
        )
        assert result.n_pathways_scored == 1
        assert "BAD" in result.skipped_pathways

    def test_no_overlapping_genes(self, expression_matrix):
        pathways = {"MISSING": ["NONEXIST1", "NONEXIST2", "NONEXIST3"]}
        result = score_pathways_from_expression(expression_matrix, pathways)
        assert result.n_pathways_scored == 0
        assert result.n_pathways_skipped == 1

    def test_empty_pathways(self, expression_matrix):
        result = score_pathways_from_expression(expression_matrix, {})
        assert result.n_pathways_scored == 0

    def test_reproducibility_with_seed(self, expression_matrix, expression_pathways):
        r1 = score_pathways_from_expression(
            expression_matrix, expression_pathways, seed=42
        )
        r2 = score_pathways_from_expression(
            expression_matrix, expression_pathways, seed=42
        )
        pd.testing.assert_frame_equal(r1.pathway_scores, r2.pathway_scores)

    def test_quality_report_in_result(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix, expression_pathways
        )
        assert result.quality_report.n_samples == 60
        assert result.quality_report.n_pathways_total == 5
        assert result.quality_report.n_pathways_covered == 5
        assert result.quality_report.mean_pathway_gene_coverage > 0.0


# ---------------------------------------------------------------------------
# Result Serialization Tests
# ---------------------------------------------------------------------------


class TestExpressionScoringResult:
    def test_to_dict(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix, expression_pathways
        )
        d = result.to_dict()
        assert d["method"] == "ssgsea"
        assert d["n_samples"] == 60
        assert d["n_pathways_scored"] == 5
        assert "quality_report" in d
        assert isinstance(d["pathway_names"], list)

    def test_format_report(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix, expression_pathways
        )
        report = result.format_report()
        assert "Expression Pathway Scoring Report" in report
        assert "ssgsea" in report
        assert "5" in report  # pathways scored

    def test_get_citations_ssgsea(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.SSGSEA,
        )
        citations = result.get_citations()
        assert len(citations) >= 1
        assert "Barbie" in citations[0]

    def test_get_citations_gsva(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.GSVA,
        )
        citations = result.get_citations()
        assert "Hanzelmann" in citations[0]

    def test_get_citations_mean_z(self, expression_matrix, expression_pathways):
        result = score_pathways_from_expression(
            expression_matrix,
            expression_pathways,
            method=ExpressionScoringMethod.MEAN_Z,
        )
        citations = result.get_citations()
        assert "Lee" in citations[0]


class TestExpressionDataQualityReport:
    def test_to_dict(self):
        report = ExpressionDataQualityReport(
            n_samples=100,
            n_genes=5000,
            is_usable=True,
        )
        d = report.to_dict()
        assert d["n_samples"] == 100
        assert d["n_genes"] == 5000
        assert d["is_usable"] is True

    def test_defaults(self):
        report = ExpressionDataQualityReport()
        assert report.n_samples == 0
        assert report.is_usable is True
        assert report.warnings == []
