"""Tests for the bulk deconvolution module."""

import json

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.deconvolution import (
    DeconvolutionMethod,
    DeconvolutionQualityReport,
    DeconvolutionResult,
    _nnls_deconvolve,
    build_reference_profile,
    combine_features,
    deconvolve_bulk,
    generate_synthetic_bulk,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def reference_profile():
    """5 cell types x 100 genes reference profile (deterministic)."""
    rng = np.random.RandomState(42)
    cell_types = ["Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "Endothelial"]
    genes = [f"GENE_{i}" for i in range(100)]
    # Each cell type has a distinct expression pattern
    data = np.abs(rng.normal(5, 2, (5, 100)))
    for i in range(5):
        # Boost 20 marker genes per cell type
        data[i, i * 20 : (i + 1) * 20] += 10
    return pd.DataFrame(data, index=cell_types, columns=genes)


@pytest.fixture
def synthetic_data(reference_profile):
    """Synthetic bulk + true proportions + labels from known reference."""
    return generate_synthetic_bulk(
        reference_profile,
        n_samples=60,
        n_subtypes=3,
        noise_level=0.1,
        seed=42,
    )


@pytest.fixture
def simple_bulk(reference_profile):
    """Small 10-sample bulk expression for edge case tests."""
    rng = np.random.RandomState(99)
    genes = list(reference_profile.columns)
    data = np.abs(rng.normal(5, 2, (10, len(genes))))
    sample_ids = [f"S_{i}" for i in range(10)]
    return pd.DataFrame(data, index=sample_ids, columns=genes)


@pytest.fixture
def single_cell_reference():
    """Raw single-cell expression (100 cells, 100 genes, 5 cell types)."""
    rng = np.random.RandomState(42)
    n_cells = 100
    n_genes = 100
    labels = np.repeat(["Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "Endothelial"], 20)
    data = np.abs(rng.normal(5, 2, (n_cells, n_genes)))
    # Add marker genes per cell type
    for i in range(5):
        start = i * 20
        end = (i + 1) * 20
        data[start:end, i * 20 : (i + 1) * 20] += 10
    genes = [f"GENE_{i}" for i in range(n_genes)]
    cells = [f"CELL_{i}" for i in range(n_cells)]
    return pd.DataFrame(data, index=cells, columns=genes), labels


# ---------------------------------------------------------------------------
# Enum tests
# ---------------------------------------------------------------------------


class TestDeconvolutionMethod:
    def test_nnls_value(self):
        assert DeconvolutionMethod.NNLS.value == "nnls"

    def test_from_string(self):
        assert DeconvolutionMethod("nnls") == DeconvolutionMethod.NNLS

    def test_invalid_value(self):
        with pytest.raises(ValueError):
            DeconvolutionMethod("invalid_method")


# ---------------------------------------------------------------------------
# Quality report tests
# ---------------------------------------------------------------------------


class TestDeconvolutionQualityReport:
    def test_to_dict(self):
        report = DeconvolutionQualityReport(
            n_samples=50,
            n_genes_bulk=200,
            n_genes_reference=100,
            n_genes_shared=80,
            n_cell_types=5,
            gene_coverage=0.8,
            proportion_stats={"Neuron": {"mean": 0.3, "min": 0.1, "max": 0.5, "std": 0.1}},
            warnings=["test warning"],
            is_usable=True,
        )
        d = report.to_dict()
        assert isinstance(d, dict)
        assert d["n_samples"] == 50
        assert d["gene_coverage"] == 0.8
        assert "Neuron" in d["proportion_stats"]
        assert d["warnings"] == ["test warning"]

    def test_default_values(self):
        report = DeconvolutionQualityReport()
        d = report.to_dict()
        assert d["n_samples"] == 0
        assert d["is_usable"]

    def test_json_serializable(self):
        report = DeconvolutionQualityReport(
            n_samples=10,
            proportion_stats={"A": {"mean": 0.5, "min": 0.1, "max": 0.9, "std": 0.2}},
        )
        json_str = json.dumps(report.to_dict())
        assert len(json_str) > 0


# ---------------------------------------------------------------------------
# Result tests
# ---------------------------------------------------------------------------


class TestDeconvolutionResult:
    def test_to_dict(self, synthetic_data, reference_profile):
        bulk_expr, true_proportions, labels = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        d = result.to_dict()
        assert isinstance(d, dict)
        assert d["method"] == "nnls"
        assert d["n_samples"] == 60
        assert d["n_cell_types"] == 5

    def test_format_report(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        report = result.format_report()
        assert "Bulk Deconvolution Report" in report
        assert "NNLS" in report
        assert "Neuron" in report
        assert "Gene Overlap" in report

    def test_get_citations(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        citations = result.get_citations()
        assert len(citations) >= 2
        assert any("Lawson" in c for c in citations)
        assert any("Newman" in c for c in citations)

    def test_json_serializable(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        json_str = json.dumps(result.to_dict())
        assert len(json_str) > 0


# ---------------------------------------------------------------------------
# build_reference_profile tests
# ---------------------------------------------------------------------------


class TestBuildReferenceProfile:
    def test_shape(self, single_cell_reference):
        expr, labels = single_cell_reference
        profile = build_reference_profile(expr, labels)
        assert profile.shape == (5, 100)

    def test_cell_types_present(self, single_cell_reference):
        expr, labels = single_cell_reference
        profile = build_reference_profile(expr, labels)
        assert "Neuron" in profile.index
        assert "Astrocyte" in profile.index

    def test_min_cells_filter(self, single_cell_reference):
        expr, labels = single_cell_reference
        # With min_cells=25, all types (20 cells each) get excluded → should raise
        with pytest.raises(ValueError, match="No cell types passed"):
            build_reference_profile(expr, labels, min_cells_per_type=25)

    def test_min_cells_excludes_all_raises(self, single_cell_reference):
        expr, labels = single_cell_reference
        with pytest.raises(ValueError, match="No cell types passed"):
            build_reference_profile(expr, labels, min_cells_per_type=100)

    def test_empty_reference_raises(self):
        with pytest.raises(ValueError, match="empty or None"):
            build_reference_profile(pd.DataFrame(), np.array([]))

    def test_mismatched_lengths_raises(self, single_cell_reference):
        expr, labels = single_cell_reference
        with pytest.raises(ValueError, match="does not match"):
            build_reference_profile(expr, labels[:10])

    def test_values_are_means(self, single_cell_reference):
        expr, labels = single_cell_reference
        profile = build_reference_profile(expr, labels)
        # Neuron mean should match manual computation
        neuron_mask = labels == "Neuron"
        expected_mean = expr.values[neuron_mask].mean(axis=0)
        np.testing.assert_allclose(profile.loc["Neuron"].values, expected_mean, rtol=1e-10)


# ---------------------------------------------------------------------------
# deconvolve_bulk tests
# ---------------------------------------------------------------------------


class TestDeconvolveBulk:
    def test_proportions_shape(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        assert result.cell_type_proportions.shape == (60, 5)

    def test_proportions_sum_to_one(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        row_sums = result.cell_type_proportions.sum(axis=1)
        np.testing.assert_allclose(row_sums.values, 1.0, atol=1e-10)

    def test_proportions_non_negative(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        assert (result.cell_type_proportions.values >= 0).all()

    def test_synthetic_recovery(self, synthetic_data, reference_profile):
        """Estimated proportions should correlate with true proportions."""
        bulk_expr, true_proportions, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)

        # Check per-cell-type correlation
        for ct in reference_profile.index:
            true_col = true_proportions[ct].values
            est_col = result.cell_type_proportions[ct].values
            corr = np.corrcoef(true_col, est_col)[0, 1]
            # With low noise, correlation should be positive
            assert corr > 0.3, f"Low correlation for {ct}: {corr:.3f}"

    def test_quality_report_populated(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        qr = result.quality_report
        assert qr.n_samples == 60
        assert qr.n_genes_shared > 0
        assert qr.n_cell_types == 5
        assert qr.gene_coverage > 0

    def test_min_genes_guard(self, reference_profile):
        """Raises when too few shared genes."""
        rng = np.random.RandomState(42)
        # Bulk with completely different gene names
        bulk = pd.DataFrame(
            rng.normal(0, 1, (10, 20)),
            columns=[f"OTHER_GENE_{i}" for i in range(20)],
        )
        with pytest.raises(ValueError, match="shared genes"):
            deconvolve_bulk(bulk, reference_profile, min_genes=50, show_progress=False)

    def test_empty_bulk_raises(self, reference_profile):
        with pytest.raises(ValueError, match="empty or None"):
            deconvolve_bulk(pd.DataFrame(), reference_profile, show_progress=False)

    def test_empty_reference_raises(self):
        rng = np.random.RandomState(42)
        bulk = pd.DataFrame(rng.normal(0, 1, (10, 20)))
        with pytest.raises(ValueError, match="empty or None"):
            deconvolve_bulk(bulk, pd.DataFrame(), show_progress=False)

    def test_method_nnls(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(
            bulk_expr, reference_profile, method=DeconvolutionMethod.NNLS, show_progress=False
        )
        assert result.method == DeconvolutionMethod.NNLS

    def test_reference_cell_types(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        assert result.reference_cell_types == list(reference_profile.index)

    def test_sample_ids_preserved(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        assert list(result.cell_type_proportions.index) == list(bulk_expr.index)


# ---------------------------------------------------------------------------
# NNLS helper tests
# ---------------------------------------------------------------------------


class TestNNLSDeconvolve:
    def test_non_negative(self):
        rng = np.random.RandomState(42)
        ref = np.abs(rng.normal(5, 2, (3, 50)))
        bulk = ref[0] * 0.5 + ref[1] * 0.3 + ref[2] * 0.2
        result = _nnls_deconvolve(bulk, ref)
        assert (result >= 0).all()

    def test_sums_to_one(self):
        rng = np.random.RandomState(42)
        ref = np.abs(rng.normal(5, 2, (3, 50)))
        bulk = ref[0] * 0.5 + ref[1] * 0.3 + ref[2] * 0.2
        result = _nnls_deconvolve(bulk, ref)
        assert abs(result.sum() - 1.0) < 1e-10

    def test_known_solution_recovery(self):
        """With noiseless data, should recover true proportions."""
        rng = np.random.RandomState(42)
        ref = np.abs(rng.normal(10, 3, (3, 100)))
        # Make reference profiles more distinct
        ref[0, :33] += 20
        ref[1, 33:66] += 20
        ref[2, 66:] += 20

        true_props = np.array([0.5, 0.3, 0.2])
        bulk = true_props @ ref

        estimated = _nnls_deconvolve(bulk, ref)
        np.testing.assert_allclose(estimated, true_props, atol=0.05)

    def test_zero_bulk_returns_uniform(self):
        ref = np.ones((3, 50))
        bulk = np.zeros(50)
        result = _nnls_deconvolve(bulk, ref)
        # Should return uniform when NNLS finds no signal
        expected = np.ones(3) / 3
        np.testing.assert_allclose(result, expected, atol=1e-10)


# ---------------------------------------------------------------------------
# combine_features tests
# ---------------------------------------------------------------------------


class TestCombineFeatures:
    def test_shape(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        pathway_scores = pd.DataFrame(
            np.random.RandomState(42).normal(0, 1, (60, 10)),
            index=list(bulk_expr.index),
            columns=[f"PATHWAY_{i}" for i in range(10)],
        )
        combined = combine_features(pathway_scores, result.cell_type_proportions)
        assert combined.shape == (60, 15)  # 10 pathways + 5 cell types

    def test_column_names_prefixed(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        pathway_scores = pd.DataFrame(
            np.random.RandomState(42).normal(0, 1, (60, 10)),
            index=list(bulk_expr.index),
            columns=[f"PATHWAY_{i}" for i in range(10)],
        )
        combined = combine_features(pathway_scores, result.cell_type_proportions)
        ct_cols = [c for c in combined.columns if c.startswith("CELLTYPE:")]
        assert len(ct_cols) == 5
        assert "CELLTYPE:Neuron_proportion" in combined.columns

    def test_weight_scaling(self):
        rng = np.random.RandomState(42)
        n = 20
        ps = pd.DataFrame(rng.normal(0, 1, (n, 3)), columns=["P1", "P2", "P3"])
        ct = pd.DataFrame(rng.normal(0.5, 0.1, (n, 2)), columns=["A", "B"])

        # With weight=0, cell types should be zeroed out
        combined_0 = combine_features(ps, ct, proportion_weight=0.0)
        ct_cols = [c for c in combined_0.columns if c.startswith("CELLTYPE:")]
        assert np.allclose(combined_0[ct_cols].values, 0.0)

        # With weight=1, pathways should be zeroed out
        combined_1 = combine_features(ps, ct, proportion_weight=1.0)
        pw_cols = [c for c in combined_1.columns if not c.startswith("CELLTYPE:")]
        assert np.allclose(combined_1[pw_cols].values, 0.0)

    def test_partial_sample_overlap(self):
        rng = np.random.RandomState(42)
        ps = pd.DataFrame(
            rng.normal(0, 1, (20, 3)),
            index=[f"S_{i}" for i in range(20)],
            columns=["P1", "P2", "P3"],
        )
        ct = pd.DataFrame(
            rng.normal(0.5, 0.1, (15, 2)),
            index=[f"S_{i}" for i in range(5, 20)],  # overlap: S_5 to S_19
            columns=["A", "B"],
        )
        combined = combine_features(ps, ct)
        assert len(combined) == 15

    def test_no_shared_samples_raises(self):
        ps = pd.DataFrame(
            np.ones((5, 2)),
            index=["A1", "A2", "A3", "A4", "A5"],
            columns=["P1", "P2"],
        )
        ct = pd.DataFrame(
            np.ones((5, 2)),
            index=["B1", "B2", "B3", "B4", "B5"],
            columns=["X", "Y"],
        )
        with pytest.raises(ValueError, match="No shared samples"):
            combine_features(ps, ct)

    def test_invalid_weight_raises(self):
        ps = pd.DataFrame(np.ones((5, 2)), columns=["P1", "P2"])
        ct = pd.DataFrame(np.ones((5, 2)), columns=["A", "B"])
        with pytest.raises(ValueError, match="proportion_weight"):
            combine_features(ps, ct, proportion_weight=1.5)

    def test_no_nan_in_output(self, synthetic_data, reference_profile):
        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)
        pathway_scores = pd.DataFrame(
            np.random.RandomState(42).normal(0, 1, (60, 10)),
            index=list(bulk_expr.index),
            columns=[f"PATHWAY_{i}" for i in range(10)],
        )
        combined = combine_features(pathway_scores, result.cell_type_proportions)
        assert not combined.isna().any().any()


# ---------------------------------------------------------------------------
# generate_synthetic_bulk tests
# ---------------------------------------------------------------------------


class TestGenerateSyntheticBulk:
    def test_shapes(self, reference_profile):
        bulk, proportions, labels = generate_synthetic_bulk(
            reference_profile, n_samples=90, n_subtypes=3, seed=42
        )
        assert bulk.shape == (90, 100)
        assert proportions.shape == (90, 5)
        assert labels.shape == (90,)

    def test_labels_planted(self, reference_profile):
        _, _, labels = generate_synthetic_bulk(
            reference_profile, n_samples=90, n_subtypes=3, seed=42
        )
        unique = np.unique(labels)
        assert len(unique) == 3
        for label in unique:
            assert np.sum(labels == label) == 30

    def test_distinct_compositions(self, reference_profile):
        """Each subtype should have a distinct dominant cell type."""
        _, proportions, labels = generate_synthetic_bulk(
            reference_profile, n_samples=90, n_subtypes=3, seed=42
        )
        cell_types = list(proportions.columns)

        for subtype in range(3):
            mask = labels == subtype
            mean_props = proportions.loc[mask].mean()
            dominant_idx = mean_props.values.argmax()
            # The dominant cell type should correspond to the subtype index
            assert dominant_idx == subtype, (
                f"Subtype {subtype} dominant type is {cell_types[dominant_idx]}, "
                f"expected {cell_types[subtype]}"
            )

    def test_proportions_sum_to_one(self, reference_profile):
        _, proportions, _ = generate_synthetic_bulk(reference_profile, n_samples=60, seed=42)
        row_sums = proportions.sum(axis=1)
        np.testing.assert_allclose(row_sums.values, 1.0, atol=1e-10)

    def test_proportions_non_negative(self, reference_profile):
        _, proportions, _ = generate_synthetic_bulk(reference_profile, n_samples=60, seed=42)
        assert (proportions.values >= 0).all()

    def test_bulk_non_negative(self, reference_profile):
        bulk, _, _ = generate_synthetic_bulk(
            reference_profile, n_samples=60, noise_level=0.1, seed=42
        )
        assert (bulk.values >= 0).all()

    def test_reproducible_with_seed(self, reference_profile):
        b1, p1, l1 = generate_synthetic_bulk(reference_profile, seed=42)
        b2, p2, l2 = generate_synthetic_bulk(reference_profile, seed=42)
        pd.testing.assert_frame_equal(b1, b2)
        pd.testing.assert_frame_equal(p1, p2)
        np.testing.assert_array_equal(l1, l2)

    def test_different_seeds_vary(self, reference_profile):
        b1, _, _ = generate_synthetic_bulk(reference_profile, seed=42)
        b2, _, _ = generate_synthetic_bulk(reference_profile, seed=99)
        assert not np.allclose(b1.values, b2.values)

    def test_too_many_subtypes_raises(self, reference_profile):
        with pytest.raises(ValueError, match="n_subtypes.*cannot exceed"):
            generate_synthetic_bulk(reference_profile, n_subtypes=10, seed=42)


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_single_cell_type(self):
        """Deconvolution with one cell type → all proportions = 1.0."""
        rng = np.random.RandomState(42)
        ref = pd.DataFrame(
            np.abs(rng.normal(5, 2, (1, 60))),
            index=["OnlyType"],
            columns=[f"GENE_{i}" for i in range(60)],
        )
        bulk = pd.DataFrame(
            np.abs(rng.normal(5, 2, (10, 60))),
            index=[f"S_{i}" for i in range(10)],
            columns=[f"GENE_{i}" for i in range(60)],
        )
        result = deconvolve_bulk(bulk, ref, show_progress=False)
        # All proportions should be 1.0 for the single type
        np.testing.assert_allclose(result.cell_type_proportions["OnlyType"].values, 1.0, atol=1e-10)

    def test_partial_gene_overlap(self, reference_profile):
        """Works when bulk has extra genes not in reference."""
        rng = np.random.RandomState(42)
        # Bulk has 100 shared + 50 extra genes
        shared_genes = list(reference_profile.columns)
        extra_genes = [f"EXTRA_{i}" for i in range(50)]
        all_genes = shared_genes + extra_genes
        bulk = pd.DataFrame(
            np.abs(rng.normal(5, 2, (20, 150))),
            index=[f"S_{i}" for i in range(20)],
            columns=all_genes,
        )
        result = deconvolve_bulk(bulk, reference_profile, show_progress=False)
        assert result.cell_type_proportions.shape == (20, 5)

    def test_low_gene_coverage_warning(self):
        """Warning when gene overlap is low."""
        rng = np.random.RandomState(42)
        # Reference with 200 genes, bulk shares only 55
        ref_genes = [f"GENE_{i}" for i in range(200)]
        ref = pd.DataFrame(
            np.abs(rng.normal(5, 2, (3, 200))),
            index=["A", "B", "C"],
            columns=ref_genes,
        )
        bulk_genes = [f"GENE_{i}" for i in range(55)] + [f"OTHER_{i}" for i in range(45)]
        bulk = pd.DataFrame(
            np.abs(rng.normal(5, 2, (10, 100))),
            index=[f"S_{i}" for i in range(10)],
            columns=bulk_genes,
        )
        result = deconvolve_bulk(bulk, ref, min_genes=50, show_progress=False)
        assert any("Low gene overlap" in w for w in result.quality_report.warnings)

    def test_combine_zero_variance_column(self):
        """combine_features handles constant columns gracefully."""
        rng = np.random.RandomState(42)
        n = 20
        ps = pd.DataFrame(rng.normal(0, 1, (n, 3)), columns=["P1", "P2", "P3"])
        # One cell type has constant value
        ct = pd.DataFrame(
            {
                "TypeA": np.ones(n) * 0.5,
                "TypeB": rng.normal(0.3, 0.1, n),
            }
        )
        combined = combine_features(ps, ct)
        assert not combined.isna().any().any()


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------


class TestIntegration:
    def test_deconvolve_combine_cluster(self, synthetic_data, reference_profile):
        """Full workflow: deconvolve → combine → cluster."""
        from pathway_subtyping.clustering import run_clustering

        bulk_expr, true_proportions, labels = synthetic_data
        n_clusters = 3

        # Deconvolve
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)

        # Create pathway scores (using some of the expression data as proxy)
        rng = np.random.RandomState(42)
        pathway_scores = pd.DataFrame(
            rng.normal(0, 1, (60, 10)),
            index=list(bulk_expr.index),
            columns=[f"PATHWAY_{i}" for i in range(10)],
        )
        # Add subtype signal to pathway scores
        for s in range(n_clusters):
            mask = labels == s
            pathway_scores.iloc[mask, s * 3 : (s + 1) * 3] += 3

        # Combine
        combined = combine_features(pathway_scores, result.cell_type_proportions)
        assert combined.shape[1] == 15  # 10 pathways + 5 cell types

        # Cluster
        clustering = run_clustering(combined.values, n_clusters=n_clusters, seed=42)
        assert len(clustering.labels) == 60
        assert len(set(clustering.labels)) == n_clusters

    def test_multi_omic_integration(self, synthetic_data, reference_profile):
        """Deconvolution output is compatible with prepare_modality."""
        from pathway_subtyping.multi_omic import ModalityType, prepare_modality

        bulk_expr, _, _ = synthetic_data
        result = deconvolve_bulk(bulk_expr, reference_profile, show_progress=False)

        modality = prepare_modality(
            ModalityType.DECONVOLUTION,
            result.cell_type_proportions,
            label="deconv",
        )
        assert modality.modality_type == ModalityType.DECONVOLUTION
        assert modality.n_samples == 60
        assert modality.label == "deconv"
