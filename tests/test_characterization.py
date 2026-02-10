"""
Tests for the Subtype Characterization Module.
"""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.characterization import (
    CharacterizationResult,
    GeneContribution,
    PathwayEnrichment,
    SubtypeProfile,
    _cohens_d,
    characterize_subtypes,
    export_characterization,
    gene_contribution_scores,
    generate_gene_heatmap,
    generate_subtype_heatmap,
    pathway_enrichment_analysis,
)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def synthetic_data():
    """Generate well-separated synthetic data for testing."""
    rng = np.random.RandomState(42)

    n_per_cluster = 30
    n_pathways = 5

    # Cluster 0: high in pathway 0 and 1
    c0 = rng.normal(loc=[2.0, 1.5, 0.0, 0.0, 0.0], scale=0.5, size=(n_per_cluster, n_pathways))
    # Cluster 1: high in pathway 2 and 3
    c1 = rng.normal(loc=[0.0, 0.0, 2.0, 1.5, 0.0], scale=0.5, size=(n_per_cluster, n_pathways))
    # Cluster 2: high in pathway 4
    c2 = rng.normal(loc=[0.0, 0.0, 0.0, 0.0, 2.0], scale=0.5, size=(n_per_cluster, n_pathways))

    data = np.vstack([c0, c1, c2])
    labels = np.array([0] * n_per_cluster + [1] * n_per_cluster + [2] * n_per_cluster)

    pathway_names = [f"Pathway_{i}" for i in range(n_pathways)]
    sample_ids = [f"S{i:03d}" for i in range(len(labels))]

    pathway_scores = pd.DataFrame(data, columns=pathway_names, index=sample_ids)

    return pathway_scores, labels


@pytest.fixture
def synthetic_data_with_genes(synthetic_data):
    """Add gene burden data to synthetic data."""
    pathway_scores, labels = synthetic_data
    rng = np.random.RandomState(42)

    n_samples = len(labels)
    genes = [f"GENE{i}" for i in range(20)]

    # Gene burdens correlated with pathway scores
    gene_data = rng.exponential(scale=0.5, size=(n_samples, len(genes)))
    # Make first 4 genes high in cluster 0
    gene_data[labels == 0, :4] += 2.0
    # Make genes 4-8 high in cluster 1
    gene_data[labels == 1, 4:8] += 2.0
    # Make genes 8-12 high in cluster 2
    gene_data[labels == 2, 8:12] += 2.0

    gene_burdens = pd.DataFrame(gene_data, columns=genes, index=pathway_scores.index)

    pathways = {
        "Pathway_0": ["GENE0", "GENE1", "GENE2", "GENE3"],
        "Pathway_1": ["GENE4", "GENE5", "GENE6", "GENE7"],
        "Pathway_2": ["GENE8", "GENE9", "GENE10", "GENE11"],
        "Pathway_3": ["GENE12", "GENE13", "GENE14", "GENE15"],
        "Pathway_4": ["GENE16", "GENE17", "GENE18", "GENE19"],
    }

    return pathway_scores, labels, gene_burdens, pathways


@pytest.fixture
def confidence_scores(synthetic_data):
    """Generate confidence scores."""
    _, labels = synthetic_data
    rng = np.random.RandomState(42)
    return rng.uniform(0.7, 1.0, size=len(labels))


# =============================================================================
# TEST: Cohen's d
# =============================================================================


class TestCohensD:
    def test_identical_groups(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0, 3.0])
        assert _cohens_d(a, b) == pytest.approx(0.0)

    def test_separated_groups(self):
        a = np.array([10.0, 11.0, 12.0])
        b = np.array([0.0, 1.0, 2.0])
        d = _cohens_d(a, b)
        assert d > 5.0  # Very large effect

    def test_positive_when_group_higher(self):
        a = np.array([5.0, 6.0, 7.0])
        b = np.array([1.0, 2.0, 3.0])
        assert _cohens_d(a, b) > 0

    def test_negative_when_group_lower(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([5.0, 6.0, 7.0])
        assert _cohens_d(a, b) < 0

    def test_small_group(self):
        a = np.array([1.0])
        b = np.array([2.0, 3.0, 4.0])
        assert _cohens_d(a, b) == 0.0  # n1 < 2

    def test_zero_variance(self):
        a = np.array([5.0, 5.0, 5.0])
        b = np.array([5.0, 5.0, 5.0])
        assert _cohens_d(a, b) == 0.0


# =============================================================================
# TEST: Pathway Enrichment
# =============================================================================


class TestPathwayEnrichment:
    def test_basic_enrichment(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = pathway_enrichment_analysis(pathway_scores, labels)

        assert len(result) == 3  # 3 clusters
        for cid in [0, 1, 2]:
            assert cid in result
            assert len(result[cid]) == 5  # 5 pathways

    def test_significant_pathways_found(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = pathway_enrichment_analysis(pathway_scores, labels)

        # Cluster 0 should have Pathway_0 as top enriched
        c0_top = result[0][0]  # Sorted by abs(effect_size)
        assert c0_top.pathway in ["Pathway_0", "Pathway_1"]
        assert c0_top.effect_size > 0

    def test_fdr_correction_applied(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = pathway_enrichment_analysis(pathway_scores, labels)

        # All enrichments should have q_value >= p_value
        for cid, enrichments in result.items():
            for e in enrichments:
                assert e.q_value >= e.p_value or np.isclose(
                    e.q_value, e.p_value
                )

    def test_custom_fdr_alpha(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        strict = pathway_enrichment_analysis(
            pathway_scores, labels, fdr_alpha=0.001
        )
        lenient = pathway_enrichment_analysis(
            pathway_scores, labels, fdr_alpha=0.5
        )

        # Lenient should have at least as many significant results
        n_sig_strict = sum(
            sum(1 for e in enrs if e.significant) for enrs in strict.values()
        )
        n_sig_lenient = sum(
            sum(1 for e in enrs if e.significant) for enrs in lenient.values()
        )
        assert n_sig_lenient >= n_sig_strict

    def test_sorted_by_effect_size(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = pathway_enrichment_analysis(pathway_scores, labels)

        for cid, enrichments in result.items():
            effect_sizes = [abs(e.effect_size) for e in enrichments]
            assert effect_sizes == sorted(effect_sizes, reverse=True)

    def test_to_dict(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = pathway_enrichment_analysis(pathway_scores, labels)

        d = result[0][0].to_dict()
        assert "pathway" in d
        assert "mean_score" in d
        assert "effect_size" in d
        assert "p_value" in d
        assert "q_value" in d
        assert "significant" in d

    def test_single_cluster(self, synthetic_data):
        pathway_scores, _ = synthetic_data
        single_labels = np.zeros(len(pathway_scores), dtype=int)
        result = pathway_enrichment_analysis(pathway_scores, single_labels)

        assert len(result) == 1
        assert 0 in result
        for e in result[0]:
            assert not e.significant
            assert e.p_value == 1.0


# =============================================================================
# TEST: Gene Contribution Scores
# =============================================================================


class TestGeneContributions:
    def test_basic_contributions(self, synthetic_data_with_genes):
        _, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = gene_contribution_scores(gene_burdens, labels, pathways)

        assert len(result) == 3
        for cid in [0, 1, 2]:
            assert cid in result
            assert len(result[cid]) <= 20  # top_n default

    def test_cluster0_top_genes(self, synthetic_data_with_genes):
        _, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = gene_contribution_scores(gene_burdens, labels, pathways)

        # Cluster 0 should have GENE0-3 as top contributors
        c0_genes = {g.gene for g in result[0][:4]}
        expected = {"GENE0", "GENE1", "GENE2", "GENE3"}
        assert len(c0_genes & expected) >= 2  # At least 2 of expected in top 4

    def test_custom_top_n(self, synthetic_data_with_genes):
        _, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = gene_contribution_scores(
            gene_burdens, labels, pathways, top_n=5
        )

        for cid in result:
            assert len(result[cid]) <= 5

    def test_genes_have_pathway_labels(self, synthetic_data_with_genes):
        _, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = gene_contribution_scores(gene_burdens, labels, pathways)

        for cid, contribs in result.items():
            for g in contribs:
                assert g.pathway != ""
                assert g.pathway.startswith("Pathway_")

    def test_to_dict(self, synthetic_data_with_genes):
        _, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = gene_contribution_scores(gene_burdens, labels, pathways)

        d = result[0][0].to_dict()
        assert "gene" in d
        assert "pathway" in d
        assert "effect_size" in d
        assert "fold_change" in d

    def test_no_matching_genes(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        gene_burdens = pd.DataFrame(
            np.random.rand(len(labels), 3),
            columns=["UNRELATED1", "UNRELATED2", "UNRELATED3"],
            index=pathway_scores.index,
        )
        pathways = {"Pathway_0": ["MISSING_GENE1", "MISSING_GENE2"]}

        result = gene_contribution_scores(gene_burdens, labels, pathways)
        for cid in result:
            assert result[cid] == []

    def test_sorted_by_effect_size(self, synthetic_data_with_genes):
        _, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = gene_contribution_scores(gene_burdens, labels, pathways)

        for cid, contribs in result.items():
            effect_sizes = [abs(g.effect_size) for g in contribs]
            assert effect_sizes == sorted(effect_sizes, reverse=True)


# =============================================================================
# TEST: characterize_subtypes (main entry)
# =============================================================================


class TestCharacterizeSubtypes:
    def test_basic_characterization(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        assert isinstance(result, CharacterizationResult)
        assert result.n_subtypes == 3
        assert result.n_samples == 90
        assert result.n_pathways == 5

    def test_with_gene_data(self, synthetic_data_with_genes):
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )

        assert result.n_genes == 20
        for p in result.subtype_profiles:
            assert len(p.top_genes) > 0

    def test_with_confidence_scores(self, synthetic_data, confidence_scores):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(
            pathway_scores, labels,
            confidence_scores=confidence_scores,
        )

        for p in result.subtype_profiles:
            assert p.mean_confidence > 0.0

    def test_with_cluster_names(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        names = {0: "Synaptic", 1: "Chromatin", 2: "Immune"}
        result = characterize_subtypes(
            pathway_scores, labels, cluster_names=names
        )

        labels_set = {p.subtype_label for p in result.subtype_profiles}
        assert labels_set == {"Synaptic", "Chromatin", "Immune"}

    def test_default_cluster_names(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        for p in result.subtype_profiles:
            assert p.subtype_label.startswith("Subtype_")

    def test_subtype_fractions_sum_to_one(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        total = sum(p.fraction for p in result.subtype_profiles)
        assert total == pytest.approx(1.0)

    def test_sample_counts_correct(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        total = sum(p.n_samples for p in result.subtype_profiles)
        assert total == len(labels)

    def test_mismatched_lengths_raises(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        with pytest.raises(ValueError, match="samples"):
            characterize_subtypes(pathway_scores, labels[:10])

    def test_gene_burdens_without_pathways_skips(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        gene_burdens = pd.DataFrame(
            np.random.rand(len(labels), 5),
            columns=[f"G{i}" for i in range(5)],
            index=pathway_scores.index,
        )
        result = characterize_subtypes(
            pathway_scores, labels, gene_burdens=gene_burdens
        )
        # Should complete without error, but no gene contributions
        for p in result.subtype_profiles:
            assert p.top_genes == []

    def test_custom_fdr_alpha(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(
            pathway_scores, labels, fdr_alpha=0.001
        )
        assert result.fdr_alpha == 0.001

    def test_custom_top_n_genes(self, synthetic_data_with_genes):
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
            top_n_genes=5,
        )
        for p in result.subtype_profiles:
            assert len(p.top_genes) <= 5


# =============================================================================
# TEST: CharacterizationResult methods
# =============================================================================


class TestCharacterizationResult:
    def test_to_dict(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        d = result.to_dict()

        assert "summary" in d
        assert "subtype_profiles" in d
        assert d["summary"]["n_subtypes"] == 3
        assert len(d["subtype_profiles"]) == 3

    def test_format_report(self, synthetic_data_with_genes):
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )
        report = result.format_report()

        assert "## Subtype Characterization" in report
        assert "Subtypes discovered" in report
        assert "Subtype_0" in report or "Subtype_1" in report

    def test_format_report_with_enriched_pathways(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        report = result.format_report()

        assert "Pathway" in report
        assert "Effect Size" in report

    def test_get_citations(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        citations = result.get_citations()

        assert len(citations) == 3
        assert any("Kruskal" in c for c in citations)
        assert any("Cohen" in c for c in citations)
        assert any("Benjamini" in c for c in citations)


# =============================================================================
# TEST: SubtypeProfile
# =============================================================================


class TestSubtypeProfile:
    def test_to_dict(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        profile = result.subtype_profiles[0]
        d = profile.to_dict()

        assert "subtype_id" in d
        assert "subtype_label" in d
        assert "n_samples" in d
        assert "fraction" in d
        assert "enriched_pathways" in d
        assert "pathway_score_means" in d

    def test_pathway_score_means_populated(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        for profile in result.subtype_profiles:
            assert len(profile.pathway_score_means) == 5
            for pw, val in profile.pathway_score_means.items():
                assert isinstance(val, float)


# =============================================================================
# TEST: Heatmaps
# =============================================================================


class TestHeatmaps:
    def test_subtype_heatmap_returns_figure(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        fig = generate_subtype_heatmap(result)

        # Should return figure (matplotlib available in test env)
        assert fig is not None

        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_subtype_heatmap_saves_file(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "heatmap.png")
            fig = generate_subtype_heatmap(result, output_path=path)
            assert os.path.exists(path)
            assert os.path.getsize(path) > 0

    def test_gene_heatmap_returns_figure(self, synthetic_data_with_genes):
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )
        fig = generate_gene_heatmap(result)
        assert fig is not None

        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_gene_heatmap_saves_file(self, synthetic_data_with_genes):
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "gene_heatmap.png")
            fig = generate_gene_heatmap(result, output_path=path)
            assert os.path.exists(path)
            assert os.path.getsize(path) > 0

    def test_gene_heatmap_no_genes(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        fig = generate_gene_heatmap(result)
        assert fig is None  # No gene data, should return None

    def test_heatmap_custom_figsize(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)
        fig = generate_subtype_heatmap(result, figsize=(8, 4))
        assert fig is not None

        import matplotlib.pyplot as plt
        plt.close(fig)


# =============================================================================
# TEST: Export
# =============================================================================


class TestExport:
    def test_csv_export(self, synthetic_data_with_genes):
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            files = export_characterization(result, tmpdir)

            assert len(files) == 4
            for f in files:
                assert os.path.exists(f)
                assert f.endswith(".csv")

            # Check summary CSV
            summary = pd.read_csv(
                os.path.join(tmpdir, "subtype_summary.csv")
            )
            assert len(summary) == 3
            assert "subtype_id" in summary.columns
            assert "n_samples" in summary.columns

            # Check enrichment CSV
            enrichment = pd.read_csv(
                os.path.join(tmpdir, "pathway_enrichment.csv")
            )
            assert "effect_size" in enrichment.columns
            assert "q_value" in enrichment.columns

            # Check gene contributions CSV
            genes = pd.read_csv(
                os.path.join(tmpdir, "gene_contributions.csv")
            )
            assert "gene" in genes.columns
            assert "pathway" in genes.columns

            # Check matrix CSV
            matrix = pd.read_csv(
                os.path.join(tmpdir, "pathway_scores_matrix.csv"),
                index_col=0,
            )
            assert matrix.shape[0] == 3  # 3 subtypes

    def test_csv_export_no_genes(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        with tempfile.TemporaryDirectory() as tmpdir:
            files = export_characterization(result, tmpdir)
            assert len(files) == 4

            # Gene contributions CSV should exist (may be empty or header-only)
            gene_path = os.path.join(tmpdir, "gene_contributions.csv")
            assert os.path.exists(gene_path)
            assert os.path.getsize(gene_path) < 100  # Very small file

    def test_excel_export(self, synthetic_data_with_genes):
        """Test Excel export (requires openpyxl)."""
        pytest.importorskip("openpyxl")

        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            files = export_characterization(
                result, tmpdir, formats=["excel"]
            )
            assert len(files) == 1
            assert files[0].endswith(".xlsx")
            assert os.path.exists(files[0])

    def test_both_formats(self, synthetic_data):
        """Test exporting both CSV and Excel."""
        pytest.importorskip("openpyxl")

        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        with tempfile.TemporaryDirectory() as tmpdir:
            files = export_characterization(
                result, tmpdir, formats=["csv", "excel"]
            )
            csv_files = [f for f in files if f.endswith(".csv")]
            xlsx_files = [f for f in files if f.endswith(".xlsx")]
            assert len(csv_files) == 4
            assert len(xlsx_files) == 1

    def test_creates_output_dir(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        with tempfile.TemporaryDirectory() as tmpdir:
            new_dir = os.path.join(tmpdir, "nested", "output")
            files = export_characterization(result, new_dir)
            assert os.path.isdir(new_dir)
            assert len(files) == 4

    def test_default_format_is_csv(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        result = characterize_subtypes(pathway_scores, labels)

        with tempfile.TemporaryDirectory() as tmpdir:
            files = export_characterization(result, tmpdir)
            assert all(f.endswith(".csv") for f in files)


# =============================================================================
# TEST: Edge Cases
# =============================================================================


class TestEdgeCases:
    def test_two_clusters(self):
        rng = np.random.RandomState(42)
        data = np.vstack([
            rng.normal(2.0, 0.5, (20, 3)),
            rng.normal(-2.0, 0.5, (20, 3)),
        ])
        labels = np.array([0] * 20 + [1] * 20)
        pathway_scores = pd.DataFrame(
            data, columns=["P0", "P1", "P2"]
        )

        result = characterize_subtypes(pathway_scores, labels)
        assert result.n_subtypes == 2

    def test_single_cluster(self):
        rng = np.random.RandomState(42)
        data = rng.normal(0, 1, (30, 3))
        labels = np.zeros(30, dtype=int)
        pathway_scores = pd.DataFrame(
            data, columns=["P0", "P1", "P2"]
        )

        result = characterize_subtypes(pathway_scores, labels)
        assert result.n_subtypes == 1
        assert result.subtype_profiles[0].fraction == pytest.approx(1.0)

    def test_unequal_cluster_sizes(self):
        rng = np.random.RandomState(42)
        data = np.vstack([
            rng.normal(2.0, 0.5, (50, 3)),
            rng.normal(-2.0, 0.5, (10, 3)),
        ])
        labels = np.array([0] * 50 + [1] * 10)
        pathway_scores = pd.DataFrame(
            data, columns=["P0", "P1", "P2"]
        )

        result = characterize_subtypes(pathway_scores, labels)
        fracs = {p.subtype_id: p.fraction for p in result.subtype_profiles}
        assert fracs[0] == pytest.approx(50 / 60)
        assert fracs[1] == pytest.approx(10 / 60)

    def test_many_clusters(self):
        rng = np.random.RandomState(42)
        n_clusters = 8
        data = np.vstack([
            rng.normal(i * 3, 0.5, (10, 4)) for i in range(n_clusters)
        ])
        labels = np.repeat(range(n_clusters), 10)
        pathway_scores = pd.DataFrame(
            data, columns=[f"P{i}" for i in range(4)]
        )

        result = characterize_subtypes(pathway_scores, labels)
        assert result.n_subtypes == 8

    def test_single_pathway(self):
        rng = np.random.RandomState(42)
        data = np.vstack([
            rng.normal(2, 0.5, (20, 1)),
            rng.normal(-2, 0.5, (20, 1)),
        ])
        labels = np.array([0] * 20 + [1] * 20)
        pathway_scores = pd.DataFrame(data, columns=["P0"])

        result = characterize_subtypes(pathway_scores, labels)
        assert result.n_pathways == 1

    def test_reproducibility_with_seed(self, synthetic_data):
        pathway_scores, labels = synthetic_data
        r1 = characterize_subtypes(pathway_scores, labels, seed=42)
        r2 = characterize_subtypes(pathway_scores, labels, seed=42)

        # Results should be identical
        d1 = r1.to_dict()
        d2 = r2.to_dict()
        assert d1 == d2

    def test_no_samples(self):
        """Empty labels returns empty result without crashing."""
        pathway_scores = pd.DataFrame()
        labels = np.array([], dtype=int)
        result = characterize_subtypes(pathway_scores, labels)
        assert result.n_subtypes == 0
        assert result.n_samples == 0
        assert result.subtype_profiles == []

    def test_empty_pathway_scores(self):
        """Empty pathway scores returns empty result."""
        pathway_scores = pd.DataFrame(index=[f"S{i}" for i in range(10)])
        labels = np.zeros(10, dtype=int)
        result = characterize_subtypes(pathway_scores, labels)
        assert result.n_subtypes == 0
        assert result.n_pathways == 0

    def test_none_pathway_scores(self):
        """None pathway scores returns empty result."""
        labels = np.array([0, 0, 1, 1])
        result = characterize_subtypes(None, labels)
        assert result.n_subtypes == 0

    def test_to_dict_is_json_serializable(self, synthetic_data_with_genes):
        """Ensure to_dict output is JSON serializable (no numpy types)."""
        import json
        pathway_scores, labels, gene_burdens, pathways = synthetic_data_with_genes
        result = characterize_subtypes(
            pathway_scores, labels,
            gene_burdens=gene_burdens,
            pathways=pathways,
        )
        # This will raise TypeError if numpy types leak through
        json.dumps(result.to_dict())
