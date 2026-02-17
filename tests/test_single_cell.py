"""
Tests for the single-cell pathway scoring module.
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

anndata = pytest.importorskip("anndata")

from pathway_subtyping.single_cell import (  # noqa: E402
    SingleCellInputType,
    SingleCellQualityReport,
    SingleCellScoringMethod,
    SingleCellScoringResult,
    _compute_pseudobulk,
    _normalize_counts,
    _score_cells_mean_z,
    load_single_cell_data,
    score_single_cell_pathways,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_adata():
    """Small AnnData: 90 cells, 50 genes, 3 cell types (30 each)."""
    rng = np.random.RandomState(42)
    n_cells, n_genes = 90, 50
    X = rng.normal(6.0, 2.0, (n_cells, n_genes)).clip(0)

    cell_types = ["Type_A"] * 30 + ["Type_B"] * 30 + ["Type_C"] * 30
    obs = pd.DataFrame(
        {"cell_type": cell_types},
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)])

    return anndata.AnnData(X=X.astype(np.float64), obs=obs, var=var)


@pytest.fixture
def sparse_adata(simple_adata):
    """Same as simple_adata but with sparse X (90% zeros)."""
    rng = np.random.RandomState(42)
    X = simple_adata.X.copy()
    mask = rng.random(X.shape) < 0.9
    X[mask] = 0.0
    adata = simple_adata.copy()
    adata.X = scipy.sparse.csr_matrix(X)
    return adata


@pytest.fixture
def single_cell_pathways():
    """5 pathways mapping to GENE0..GENE49."""
    return {
        "PATHWAY_0": [f"GENE{i}" for i in range(0, 10)],
        "PATHWAY_1": [f"GENE{i}" for i in range(10, 20)],
        "PATHWAY_2": [f"GENE{i}" for i in range(20, 30)],
        "PATHWAY_3": [f"GENE{i}" for i in range(30, 40)],
        "PATHWAY_4": [f"GENE{i}" for i in range(40, 50)],
    }


@pytest.fixture
def enriched_adata(single_cell_pathways):
    """AnnData with clear cell-type/pathway enrichment patterns."""
    rng = np.random.RandomState(42)
    n_per = 30
    n_cells = n_per * 3
    n_genes = 50

    X = rng.normal(6.0, 1.0, (n_cells, n_genes))

    # Type_A: upregulate PATHWAY_0 genes (0-9)
    X[:n_per, 0:10] += 4.0
    # Type_B: upregulate PATHWAY_1 genes (10-19)
    X[n_per : 2 * n_per, 10:20] += 4.0
    # Type_C: upregulate PATHWAY_2 genes (20-29)
    X[2 * n_per :, 20:30] += 4.0

    X = X.clip(0)
    cell_types = ["Type_A"] * n_per + ["Type_B"] * n_per + ["Type_C"] * n_per
    obs = pd.DataFrame(
        {"cell_type": cell_types},
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)])

    return anndata.AnnData(X=X.astype(np.float64), obs=obs, var=var)


@pytest.fixture
def h5ad_path(simple_adata, tmp_path):
    """Write simple_adata to h5ad and return path."""
    path = tmp_path / "test_data.h5ad"
    simple_adata.write_h5ad(path)
    return str(path)


@pytest.fixture
def csv_path(simple_adata, tmp_path):
    """Write simple_adata expression to CSV (cells x genes) with cell_type in obs."""
    path = tmp_path / "test_cells.csv"
    df = pd.DataFrame(
        simple_adata.X,
        index=simple_adata.obs_names,
        columns=simple_adata.var_names,
    )
    df.to_csv(path)
    return str(path)


# ---------------------------------------------------------------------------
# Test: Enums
# ---------------------------------------------------------------------------


class TestSingleCellScoringMethod:
    def test_enum_values(self):
        assert SingleCellScoringMethod.MEAN_Z.value == "mean_z"
        assert SingleCellScoringMethod.PSEUDOBULK_MEAN_Z.value == "pseudobulk_mean_z"
        assert SingleCellScoringMethod.PSEUDOBULK_SSGSEA.value == "pseudobulk_ssgsea"
        assert SingleCellScoringMethod.PSEUDOBULK_GSVA.value == "pseudobulk_gsva"

    def test_from_string(self):
        assert SingleCellScoringMethod("mean_z") == SingleCellScoringMethod.MEAN_Z
        assert (
            SingleCellScoringMethod("pseudobulk_ssgsea")
            == SingleCellScoringMethod.PSEUDOBULK_SSGSEA
        )

    def test_invalid_value(self):
        with pytest.raises(ValueError):
            SingleCellScoringMethod("invalid_method")


class TestSingleCellInputType:
    def test_enum_values(self):
        assert SingleCellInputType.RAW_COUNTS.value == "raw_counts"
        assert SingleCellInputType.LOG_NORMALIZED.value == "log_normalized"
        assert SingleCellInputType.H5AD.value == "h5ad"


# ---------------------------------------------------------------------------
# Test: Quality Report
# ---------------------------------------------------------------------------


class TestSingleCellQualityReport:
    def test_to_dict(self):
        report = SingleCellQualityReport(
            n_cells=1000,
            n_genes=500,
            n_cell_types=5,
            sparsity=0.9,
        )
        d = report.to_dict()
        assert d["n_cells"] == 1000
        assert d["n_genes"] == 500
        assert d["n_cell_types"] == 5
        assert d["sparsity"] == 0.9
        assert isinstance(d["warnings"], list)
        assert isinstance(d["cell_type_counts"], dict)

    def test_defaults(self):
        report = SingleCellQualityReport()
        assert report.n_cells == 0
        assert report.is_usable
        assert report.warnings == []
        assert report.excluded_cell_types == []


# ---------------------------------------------------------------------------
# Test: _normalize_counts
# ---------------------------------------------------------------------------


class TestNormalizeCounts:
    def test_dense_normalization(self):
        X = np.array([[100, 200, 300], [50, 50, 0]], dtype=np.float64)
        result = _normalize_counts(X, scale_factor=1e4)
        assert result.shape == (2, 3)
        # Row 1: totals=600, first cell = log1p(100/600 * 10000) = log1p(1666.67)
        expected_00 = np.log1p(100 / 600 * 1e4)
        assert abs(result[0, 0] - expected_00) < 1e-6

    def test_sparse_normalization(self):
        X_dense = np.array([[100, 200, 300], [50, 50, 0]], dtype=np.float64)
        X_sparse = scipy.sparse.csr_matrix(X_dense)
        result_dense = _normalize_counts(X_dense, scale_factor=1e4)
        result_sparse = _normalize_counts(X_sparse, scale_factor=1e4)
        np.testing.assert_allclose(result_dense, result_sparse, atol=1e-10)

    def test_zero_total_handling(self):
        X = np.array([[0, 0, 0], [1, 2, 3]], dtype=np.float64)
        result = _normalize_counts(X)
        # Row with all zeros should produce zeros (no NaN/inf)
        assert np.all(np.isfinite(result))
        assert result[0, 0] == 0.0

    def test_output_shape(self):
        X = np.random.rand(100, 50)
        result = _normalize_counts(X)
        assert result.shape == (100, 50)


# ---------------------------------------------------------------------------
# Test: _compute_pseudobulk
# ---------------------------------------------------------------------------


class TestComputePseudobulk:
    def test_basic_shape(self, simple_adata):
        pb, counts, excluded = _compute_pseudobulk(simple_adata, "cell_type", min_cells_per_type=5)
        assert pb.shape == (3, 50)  # 3 cell types, 50 genes
        assert len(excluded) == 0
        assert set(pb.index) == {"Type_A", "Type_B", "Type_C"}

    def test_min_cells_filter(self, simple_adata):
        # Set threshold higher than any cell type has
        pb, counts, excluded = _compute_pseudobulk(simple_adata, "cell_type", min_cells_per_type=50)
        assert pb.shape[0] == 0
        assert len(excluded) == 3

    def test_sparse_input(self, sparse_adata):
        pb, counts, excluded = _compute_pseudobulk(sparse_adata, "cell_type", min_cells_per_type=5)
        assert pb.shape == (3, 50)
        assert len(excluded) == 0

    def test_cell_counts(self, simple_adata):
        _, counts, _ = _compute_pseudobulk(simple_adata, "cell_type", min_cells_per_type=5)
        assert counts["Type_A"] == 30
        assert counts["Type_B"] == 30
        assert counts["Type_C"] == 30

    def test_single_cell_type(self):
        """One cell type should work (returned as single-row DataFrame)."""
        rng = np.random.RandomState(42)
        X = rng.normal(5.0, 1.0, (20, 10))
        obs = pd.DataFrame(
            {"cell_type": ["Only_Type"] * 20},
            index=[f"c{i}" for i in range(20)],
        )
        var = pd.DataFrame(index=[f"G{i}" for i in range(10)])
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        pb, counts, excluded = _compute_pseudobulk(adata, "cell_type", min_cells_per_type=5)
        assert pb.shape == (1, 10)
        assert "Only_Type" in pb.index


# ---------------------------------------------------------------------------
# Test: _score_cells_mean_z
# ---------------------------------------------------------------------------


class TestScoreCellsMeanZ:
    def test_basic_scoring(self, simple_adata, single_cell_pathways):
        scores, skipped = _score_cells_mean_z(
            simple_adata.X,
            list(simple_adata.var_names),
            single_cell_pathways,
            show_progress=False,
        )
        assert scores.shape == (90, 5)
        assert len(skipped) == 0

    def test_known_enrichment(self, enriched_adata, single_cell_pathways):
        scores, _ = _score_cells_mean_z(
            enriched_adata.X,
            list(enriched_adata.var_names),
            single_cell_pathways,
            show_progress=False,
        )
        # Type_A (cells 0-29) should have highest PATHWAY_0 scores
        type_a_mean = scores.iloc[:30]["PATHWAY_0"].mean()
        type_b_mean = scores.iloc[30:60]["PATHWAY_0"].mean()
        type_c_mean = scores.iloc[60:90]["PATHWAY_0"].mean()
        assert type_a_mean > type_b_mean
        assert type_a_mean > type_c_mean

    def test_chunk_size_invariance(self, simple_adata, single_cell_pathways):
        scores_small, _ = _score_cells_mean_z(
            simple_adata.X,
            list(simple_adata.var_names),
            single_cell_pathways,
            chunk_size=10,
            show_progress=False,
        )
        scores_large, _ = _score_cells_mean_z(
            simple_adata.X,
            list(simple_adata.var_names),
            single_cell_pathways,
            chunk_size=1000,
            show_progress=False,
        )
        pd.testing.assert_frame_equal(scores_small, scores_large)

    def test_skip_small_pathways(self, simple_adata):
        pathways = {"TINY": ["GENE0"], "OK": [f"GENE{i}" for i in range(10)]}
        scores, skipped = _score_cells_mean_z(
            simple_adata.X,
            list(simple_adata.var_names),
            pathways,
            min_genes=2,
            show_progress=False,
        )
        assert "TINY" in skipped
        assert "OK" in scores.columns

    def test_sparse_input(self, sparse_adata, single_cell_pathways):
        scores, _ = _score_cells_mean_z(
            sparse_adata.X,
            list(sparse_adata.var_names),
            single_cell_pathways,
            show_progress=False,
        )
        assert scores.shape == (90, 5)
        assert np.all(np.isfinite(scores.values))

    def test_no_gene_overlap(self, simple_adata):
        pathways = {"NO_MATCH": ["NONEXISTENT1", "NONEXISTENT2"]}
        scores, skipped = _score_cells_mean_z(
            simple_adata.X,
            list(simple_adata.var_names),
            pathways,
            show_progress=False,
        )
        assert "NO_MATCH" in skipped
        assert scores.empty or len(scores.columns) == 0


# ---------------------------------------------------------------------------
# Test: load_single_cell_data
# ---------------------------------------------------------------------------


class TestLoadSingleCellData:
    def test_load_h5ad(self, h5ad_path):
        adata, report = load_single_cell_data(
            h5ad_path,
            cell_type_column="cell_type",
            input_type=SingleCellInputType.H5AD,
        )
        assert adata.n_obs == 90
        assert adata.n_vars <= 50  # may filter some genes
        assert report.n_cells == 90
        assert report.n_cell_types == 3

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_single_cell_data("/nonexistent/path.h5ad")

    def test_missing_cell_type_column(self, h5ad_path):
        with pytest.raises(ValueError, match="Cell type column"):
            load_single_cell_data(
                h5ad_path,
                cell_type_column="nonexistent_column",
            )

    def test_quality_report_populated(self, h5ad_path):
        _, report = load_single_cell_data(h5ad_path, cell_type_column="cell_type")
        assert report.n_cells > 0
        assert report.n_genes > 0
        assert report.n_cell_types == 3
        assert report.median_genes_per_cell > 0
        assert report.median_umi_per_cell > 0
        assert 0.0 <= report.sparsity <= 1.0
        assert len(report.cell_type_counts) == 3

    def test_quality_report_to_dict(self, h5ad_path):
        _, report = load_single_cell_data(h5ad_path, cell_type_column="cell_type")
        d = report.to_dict()
        assert isinstance(d, dict)
        assert "n_cells" in d
        assert "sparsity" in d
        assert "cell_type_counts" in d


# ---------------------------------------------------------------------------
# Test: score_single_cell_pathways
# ---------------------------------------------------------------------------


class TestScoreSingleCellPathways:
    def test_pseudobulk_mean_z(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        assert result.pathway_scores.shape == (3, 5)
        assert set(result.pathway_scores.index) == {"Type_A", "Type_B", "Type_C"}
        assert result.per_cell_scores is None

    def test_pseudobulk_ssgsea(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
            show_progress=False,
        )
        assert result.pathway_scores.shape == (3, 5)
        assert result.n_pathways_scored == 5

    def test_pseudobulk_gsva(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_GSVA,
            show_progress=False,
        )
        assert result.pathway_scores.shape == (3, 5)

    def test_per_cell_mean_z(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.MEAN_Z,
            show_progress=False,
        )
        assert result.pathway_scores.shape == (3, 5)  # aggregated per cell type
        assert result.per_cell_scores is not None
        assert result.per_cell_scores.shape == (90, 5)

    def test_output_shape(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
            show_progress=False,
        )
        # Index = cell types, columns = pathways
        assert list(sorted(result.pathway_scores.index)) == [
            "Type_A",
            "Type_B",
            "Type_C",
        ]
        assert list(result.pathway_scores.columns) == [
            "PATHWAY_0",
            "PATHWAY_1",
            "PATHWAY_2",
            "PATHWAY_3",
            "PATHWAY_4",
        ]

    def test_z_normalized_output(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        scores = result.pathway_scores
        # With 3 cell types, Z-normalization should give mean ~0
        for col in scores.columns:
            col_mean = scores[col].mean()
            assert abs(col_mean) < 0.1, f"Mean of {col} = {col_mean}"

    def test_known_enrichment_recovery(self, enriched_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            enriched_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        scores = result.pathway_scores
        # Type_A should have highest PATHWAY_0
        assert scores.loc["Type_A", "PATHWAY_0"] > scores.loc["Type_B", "PATHWAY_0"]
        assert scores.loc["Type_A", "PATHWAY_0"] > scores.loc["Type_C", "PATHWAY_0"]
        # Type_B should have highest PATHWAY_1
        assert scores.loc["Type_B", "PATHWAY_1"] > scores.loc["Type_A", "PATHWAY_1"]
        assert scores.loc["Type_B", "PATHWAY_1"] > scores.loc["Type_C", "PATHWAY_1"]
        # Type_C should have highest PATHWAY_2
        assert scores.loc["Type_C", "PATHWAY_2"] > scores.loc["Type_A", "PATHWAY_2"]
        assert scores.loc["Type_C", "PATHWAY_2"] > scores.loc["Type_B", "PATHWAY_2"]

    def test_sparse_input(self, sparse_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            sparse_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        assert result.pathway_scores.shape == (3, 5)
        assert np.all(np.isfinite(result.pathway_scores.values))

    def test_compute_per_cell_flag(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
            compute_per_cell=True,
            show_progress=False,
        )
        assert result.per_cell_scores is not None
        assert result.per_cell_scores.shape == (90, 5)
        # Primary output is still cell-type level
        assert result.pathway_scores.shape == (3, 5)

    def test_min_cells_per_type(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            min_cells_per_type=50,  # Higher than any type (30 each)
            show_progress=False,
        )
        # All types excluded
        assert result.pathway_scores.empty or len(result.pathway_scores) == 0
        assert len(result.quality_report.excluded_cell_types) == 3

    def test_seed_reproducibility(self, simple_adata, single_cell_pathways):
        result1 = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
            seed=42,
            show_progress=False,
        )
        result2 = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA,
            seed=42,
            show_progress=False,
        )
        pd.testing.assert_frame_equal(result1.pathway_scores, result2.pathway_scores)

    def test_no_gene_overlap(self, simple_adata):
        pathways = {"PW1": ["NONEXIST1", "NONEXIST2", "NONEXIST3"]}
        result = score_single_cell_pathways(
            simple_adata,
            pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        assert result.n_pathways_scored == 0
        assert result.n_pathways_skipped == 1

    def test_missing_cell_type_column(self, simple_adata, single_cell_pathways):
        with pytest.raises(ValueError, match="Cell type column"):
            score_single_cell_pathways(
                simple_adata,
                single_cell_pathways,
                cell_type_column="nonexistent",
            )

    def test_quality_report(self, simple_adata, single_cell_pathways):
        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        qr = result.quality_report
        assert qr.n_cells == 90
        assert qr.n_cell_types == 3
        assert qr.n_pathways_covered == 5
        assert qr.n_pathways_total == 5
        assert qr.mean_pathway_gene_coverage == 1.0


# ---------------------------------------------------------------------------
# Test: SingleCellScoringResult
# ---------------------------------------------------------------------------


class TestSingleCellScoringResult:
    def _make_result(self, method=SingleCellScoringMethod.PSEUDOBULK_SSGSEA):
        scores = pd.DataFrame(
            {"PW_A": [0.5, -0.5], "PW_B": [-0.3, 0.3]},
            index=["TypeX", "TypeY"],
        )
        report = SingleCellQualityReport(
            n_cells=100,
            n_genes=50,
            n_cell_types=2,
            sparsity=0.85,
        )
        return SingleCellScoringResult(
            pathway_scores=scores,
            per_cell_scores=None,
            pseudobulk_expression=pd.DataFrame(),
            method=method,
            quality_report=report,
            n_pathways_scored=2,
            n_pathways_skipped=0,
            skipped_pathways=[],
            n_cells_per_type={"TypeX": 50, "TypeY": 50},
        )

    def test_to_dict(self):
        result = self._make_result()
        d = result.to_dict()
        assert d["method"] == "pseudobulk_ssgsea"
        assert d["n_cell_types"] == 2
        assert d["n_pathways_scored"] == 2
        assert d["cell_types"] == ["TypeX", "TypeY"]
        assert d["has_per_cell_scores"] is False
        assert "quality_report" in d

    def test_format_report(self):
        result = self._make_result()
        report_str = result.format_report()
        assert "Single-Cell Pathway Scoring Report" in report_str
        assert "pseudobulk_ssgsea" in report_str
        assert "TypeX" in report_str
        assert "TypeY" in report_str

    def test_format_report_with_excluded(self):
        result = self._make_result()
        result.quality_report.excluded_cell_types = ["Rare_Type"]
        report_str = result.format_report()
        assert "Excluded" in report_str
        assert "Rare_Type" in report_str

    def test_format_report_with_skipped_pathways(self):
        result = self._make_result()
        result.skipped_pathways = ["SKIP_1", "SKIP_2"]
        report_str = result.format_report()
        assert "Skipped Pathways" in report_str
        assert "SKIP_1" in report_str

    def test_get_citations_ssgsea(self):
        result = self._make_result(SingleCellScoringMethod.PSEUDOBULK_SSGSEA)
        citations = result.get_citations()
        assert len(citations) == 2
        assert "Satija" in citations[0]
        assert "Barbie" in citations[1]

    def test_get_citations_gsva(self):
        result = self._make_result(SingleCellScoringMethod.PSEUDOBULK_GSVA)
        citations = result.get_citations()
        assert len(citations) == 2
        assert "Satija" in citations[0]
        assert "Hanzelmann" in citations[1]

    def test_get_citations_mean_z(self):
        result = self._make_result(SingleCellScoringMethod.MEAN_Z)
        citations = result.get_citations()
        assert len(citations) == 2
        assert "Lee" in citations[1]


# ---------------------------------------------------------------------------
# Test: Integration
# ---------------------------------------------------------------------------


class TestIntegration:
    def test_full_workflow(self, enriched_adata, single_cell_pathways):
        """End-to-end: score -> cluster -> verify compatibility."""
        from pathway_subtyping.clustering import ClusteringAlgorithm, run_clustering

        result = score_single_cell_pathways(
            enriched_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )

        # Verify output is compatible with run_clustering
        clustering = run_clustering(
            result.pathway_scores.values,
            n_clusters=3,
            algorithm=ClusteringAlgorithm.GMM,
            seed=42,
        )
        assert len(clustering.labels) == len(result.pathway_scores)
        assert clustering.n_clusters == 3

    def test_to_dict_serializable(self, simple_adata, single_cell_pathways):
        """Verify to_dict produces JSON-safe output."""
        import json

        result = score_single_cell_pathways(
            simple_adata,
            single_cell_pathways,
            method=SingleCellScoringMethod.PSEUDOBULK_MEAN_Z,
            show_progress=False,
        )
        d = result.to_dict()
        # Should not raise
        json_str = json.dumps(d)
        assert len(json_str) > 0
