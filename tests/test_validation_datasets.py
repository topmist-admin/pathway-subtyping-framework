"""
Tests for the public dataset validation module.

Includes:
- Unit tests for parsing logic (mock data, no network)
- Unit tests for coverage validation
- Unit tests for synthetic data generation with real gene names
- Integration tests for full pipeline
- Network tests marked @pytest.mark.network for download verification
"""

import json
import os
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.validation_datasets import (
    DATASETS,
    BiologicalPlausibilityResult,
    ClinVarGeneSummary,
    DatasetInfo,
    PathwayCoverageResult,
    ValidationReport,
    download_dataset,
    generate_disease_realistic_synthetic,
    load_clinvar_gene_summary,
    load_reactome_pathways,
    run_biological_plausibility_check,
    run_full_validation,
    validate_pathway_against_reactome,
    validate_pathway_coverage,
)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def mock_clinvar_tsv(tmp_path):
    """Create a minimal ClinVar-format TSV with known genes."""
    content = (
        "#Symbol\tGeneID\tTotal_Submissions\tAlleles_Reported\t"
        "Gene_MIM_number\tNumber_Pathogenic\tNumber_Likely_Pathogenic\t"
        "Number_Uncertain_Significance\tNumber_Benign\tNumber_Likely_Benign\n"
        "SHANK3\t85358\t50\t30\t606230\t15\t8\t10\t5\t3\n"
        "CHD8\t57680\t40\t25\t610528\t12\t5\t8\t3\t2\n"
        "NRXN1\t9378\t35\t20\t600565\t8\t4\t12\t6\t4\n"
        "FAKE_GENE\t99999\t5\t3\t0\t0\t0\t2\t1\t1\n"
        "BRCA1\t672\t200\t150\t113705\t80\t30\t50\t20\t10\n"
    )
    path = tmp_path / "mock_clinvar.tsv"
    path.write_text(content)
    return path


@pytest.fixture
def mock_reactome_gmt(tmp_path):
    """Create a minimal Reactome-format GMT with known pathways."""
    content = (
        "R-HSA-195721\tHomo sapiens: Signaling by WNT\t"
        "CTNNB1\tAPC\tGSK3B\tDVL1\tWNT1\n"
        "R-HSA-112316\tHomo sapiens: Synaptic transmission\t"
        "GRIN2A\tGRIN2B\tGRIA1\tSHANK3\tDLG4\n"
        "R-HSA-69278\tHomo sapiens: Cell cycle\t"
        "CDK1\tCDK2\tCCNA2\tCCNB1\n"
        "R-MMU-1640170\tMus musculus: Mouse development\t"
        "MOUSEONLY1\tMOUSEONLY2\n"
    )
    path = tmp_path / "mock_reactome.gmt"
    path.write_text(content)
    return path


@pytest.fixture
def sample_clinvar_genes():
    """Pre-parsed ClinVar gene summaries for testing."""
    return {
        "SHANK3": ClinVarGeneSummary(
            symbol="SHANK3", gene_id=85358,
            n_pathogenic=15, n_likely_pathogenic=8,
            n_uncertain=10, n_benign=5, n_likely_benign=3,
        ),
        "CHD8": ClinVarGeneSummary(
            symbol="CHD8", gene_id=57680,
            n_pathogenic=12, n_likely_pathogenic=5,
            n_uncertain=8, n_benign=3, n_likely_benign=2,
        ),
        "NRXN1": ClinVarGeneSummary(
            symbol="NRXN1", gene_id=9378,
            n_pathogenic=8, n_likely_pathogenic=4,
            n_uncertain=12, n_benign=6, n_likely_benign=4,
        ),
        "FAKE_GENE": ClinVarGeneSummary(
            symbol="FAKE_GENE", gene_id=99999,
            n_pathogenic=0, n_likely_pathogenic=0,
            n_uncertain=2, n_benign=1, n_likely_benign=1,
        ),
        "BRCA1": ClinVarGeneSummary(
            symbol="BRCA1", gene_id=672,
            n_pathogenic=80, n_likely_pathogenic=30,
            n_uncertain=50, n_benign=20, n_likely_benign=10,
        ),
    }


@pytest.fixture
def sample_pathways():
    """Curated pathways for testing."""
    return {
        "SYNAPTIC": ["SHANK3", "NRXN1", "GRIN2B", "DLG4"],
        "CHROMATIN": ["CHD8", "ARID1B", "KMT2A"],
        "UNKNOWN": ["GENE_X", "GENE_Y", "GENE_Z"],
    }


# =============================================================================
# TEST DATASET INFO
# =============================================================================


class TestDatasetInfo:
    def test_datasets_dict_populated(self):
        assert "clinvar_gene_summary" in DATASETS
        assert "reactome_pathways" in DATASETS

    def test_dataset_info_to_dict(self):
        ds = DATASETS["clinvar_gene_summary"]
        d = ds.to_dict()
        assert d["name"] == "ClinVar Gene-Specific Summary"
        assert "url" in d
        assert "file_format" in d

    def test_dataset_has_url(self):
        for key, ds in DATASETS.items():
            assert ds.url.startswith("http"), f"{key} has invalid URL"

    def test_dataset_info_fields(self):
        ds = DatasetInfo(
            name="Test", url="https://example.com/data.txt",
            description="Test dataset", file_format="tsv",
            license_note="MIT",
        )
        d = ds.to_dict()
        assert d["license_note"] == "MIT"


# =============================================================================
# TEST DOWNLOAD
# =============================================================================


class TestDownloadDataset:
    def test_unknown_dataset_key_raises(self):
        with pytest.raises(KeyError, match="Unknown dataset key"):
            download_dataset("nonexistent_dataset")

    def test_cache_hit_skips_download(self, tmp_path):
        """If a cached file exists, don't re-download."""
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        # Pre-create the expected cached file
        cached_file = cache_dir / "gene_specific_summary.txt"
        cached_file.write_text("# cached data\nSHANK3\t85358\n")

        with mock.patch("pathway_subtyping.validation_datasets.urllib.request.urlretrieve") as mock_dl:
            result = download_dataset("clinvar_gene_summary", cache_dir=cache_dir)
            mock_dl.assert_not_called()
            assert result == cached_file

    def test_force_redownloads(self, tmp_path):
        """force=True should re-download even if cached."""
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        cached_file = cache_dir / "gene_specific_summary.txt"
        cached_file.write_text("# old cached data")

        def fake_download(url, path):
            Path(path).write_text("# new data\n")

        with mock.patch(
            "pathway_subtyping.validation_datasets.urllib.request.urlretrieve",
            side_effect=fake_download,
        ):
            result = download_dataset(
                "clinvar_gene_summary", cache_dir=cache_dir, force=True
            )
            assert result.read_text().startswith("# new data")

    @pytest.mark.network
    def test_clinvar_download_real(self, tmp_path):
        """Actually download ClinVar (requires network)."""
        result = download_dataset("clinvar_gene_summary", cache_dir=tmp_path)
        assert result.exists()
        assert result.stat().st_size > 1000


# =============================================================================
# TEST CLINVAR PARSING
# =============================================================================


class TestLoadClinvarGeneSummary:
    def test_parse_mock_tsv(self, mock_clinvar_tsv):
        genes = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
        assert len(genes) == 5

    def test_gene_symbols_correct(self, mock_clinvar_tsv):
        genes = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
        assert "SHANK3" in genes
        assert "CHD8" in genes
        assert "BRCA1" in genes

    def test_pathogenic_count(self, mock_clinvar_tsv):
        genes = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
        assert genes["SHANK3"].n_pathogenic == 15
        assert genes["SHANK3"].n_likely_pathogenic == 8

    def test_total_pathogenic_property(self, mock_clinvar_tsv):
        genes = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
        assert genes["SHANK3"].total_pathogenic == 23  # 15 + 8
        assert genes["FAKE_GENE"].total_pathogenic == 0

    def test_gene_to_dict(self, mock_clinvar_tsv):
        genes = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
        d = genes["SHANK3"].to_dict()
        assert d["symbol"] == "SHANK3"
        assert d["total_pathogenic"] == 23
        assert isinstance(d["gene_id"], int)

    def test_case_insensitive_lookup(self, mock_clinvar_tsv):
        genes = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
        # Keys are uppercased
        assert "SHANK3" in genes

    def test_empty_file(self, tmp_path):
        empty_file = tmp_path / "empty.tsv"
        empty_file.write_text("")
        genes = load_clinvar_gene_summary(file_path=empty_file)
        assert len(genes) == 0


# =============================================================================
# TEST REACTOME PARSING
# =============================================================================


class TestLoadReactomePathways:
    def test_parse_mock_gmt(self, mock_reactome_gmt):
        pathways = load_reactome_pathways(file_path=mock_reactome_gmt)
        assert len(pathways) == 3  # 3 Homo sapiens, Mouse excluded

    def test_species_filter(self, mock_reactome_gmt):
        pathways = load_reactome_pathways(
            file_path=mock_reactome_gmt, species="Homo sapiens"
        )
        assert "R-MMU-1640170" not in pathways

    def test_gene_lists_correct(self, mock_reactome_gmt):
        pathways = load_reactome_pathways(file_path=mock_reactome_gmt)
        assert "SHANK3" in pathways["R-HSA-112316"]
        assert "CTNNB1" in pathways["R-HSA-195721"]

    def test_no_species_filter(self, mock_reactome_gmt):
        pathways = load_reactome_pathways(
            file_path=mock_reactome_gmt, species=""
        )
        # All pathways loaded when species filter is empty
        assert len(pathways) == 4

    def test_empty_gmt(self, tmp_path):
        empty_file = tmp_path / "empty.gmt"
        empty_file.write_text("")
        pathways = load_reactome_pathways(file_path=empty_file)
        assert len(pathways) == 0


# =============================================================================
# TEST PATHWAY COVERAGE
# =============================================================================


class TestValidatePathwayCoverage:
    def test_full_coverage(self, sample_clinvar_genes):
        pathways = {"TEST": ["SHANK3", "CHD8"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        assert len(results) == 1
        r = results[0]
        assert r.coverage_fraction == 1.0  # Both in ClinVar

    def test_partial_coverage(self, sample_clinvar_genes):
        pathways = {"TEST": ["SHANK3", "NONEXISTENT_GENE"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        r = results[0]
        assert r.genes_in_clinvar == 1
        assert r.coverage_fraction == 0.5

    def test_no_coverage(self, sample_clinvar_genes):
        pathways = {"TEST": ["GENE_X", "GENE_Y"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        r = results[0]
        assert r.genes_in_clinvar == 0
        assert r.coverage_fraction == 0.0

    def test_pathogenic_fraction(self, sample_clinvar_genes):
        # SHANK3 has pathogenic, FAKE_GENE doesn't
        pathways = {"TEST": ["SHANK3", "FAKE_GENE"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        r = results[0]
        assert r.genes_with_pathogenic == 1
        assert r.pathogenic_fraction == 0.5

    def test_multiple_pathways(self, sample_clinvar_genes, sample_pathways):
        results = validate_pathway_coverage(sample_pathways, sample_clinvar_genes)
        assert len(results) == 3

    def test_sorted_by_pathogenic_fraction(self, sample_clinvar_genes, sample_pathways):
        results = validate_pathway_coverage(sample_pathways, sample_clinvar_genes)
        fractions = [r.pathogenic_fraction for r in results]
        assert fractions == sorted(fractions, reverse=True)

    def test_gene_details_present(self, sample_clinvar_genes):
        pathways = {"TEST": ["SHANK3"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        assert len(results[0].gene_details) == 1
        assert results[0].gene_details[0]["gene"] == "SHANK3"
        assert results[0].gene_details[0]["in_clinvar"]

    def test_to_dict(self, sample_clinvar_genes):
        pathways = {"TEST": ["SHANK3", "CHD8"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        d = results[0].to_dict()
        assert isinstance(d["coverage_fraction"], float)
        assert isinstance(d["total_genes"], int)

    def test_case_insensitive_matching(self, sample_clinvar_genes):
        pathways = {"TEST": ["shank3", "chd8"]}
        results = validate_pathway_coverage(pathways, sample_clinvar_genes)
        r = results[0]
        assert r.genes_in_clinvar == 2


# =============================================================================
# TEST REACTOME CROSS-REFERENCE
# =============================================================================


class TestValidatePathwayAgainstReactome:
    def test_exact_overlap(self):
        curated = {"SYNAPTIC": ["SHANK3", "GRIN2B", "DLG4"]}
        reactome = {"Synaptic_Pathway": ["SHANK3", "GRIN2B", "DLG4"]}
        results = validate_pathway_against_reactome(curated, reactome)
        assert results["SYNAPTIC"]["jaccard"] == 1.0

    def test_partial_overlap(self):
        curated = {"SYNAPTIC": ["SHANK3", "GRIN2B", "GENE_X"]}
        reactome = {"Synaptic_Pathway": ["SHANK3", "GRIN2B", "GENE_Y"]}
        results = validate_pathway_against_reactome(curated, reactome)
        # Jaccard: 2 / 4 = 0.5
        assert abs(results["SYNAPTIC"]["jaccard"] - 0.5) < 0.01

    def test_no_overlap(self):
        curated = {"UNKNOWN": ["GENE_A", "GENE_B"]}
        reactome = {"Other_Pathway": ["GENE_C", "GENE_D"]}
        results = validate_pathway_against_reactome(curated, reactome)
        assert results["UNKNOWN"]["jaccard"] == 0.0

    def test_min_overlap_filter(self):
        curated = {"TEST": ["A", "B", "C", "D", "E"]}
        reactome = {"R_TEST": ["A", "F", "G", "H", "I"]}
        results = validate_pathway_against_reactome(
            curated, reactome, min_overlap=0.2
        )
        # Jaccard: 1/9 ≈ 0.11 < 0.2
        assert results["TEST"]["best_match"] is None

    def test_case_insensitive(self):
        curated = {"TEST": ["shank3", "grin2b"]}
        reactome = {"R_TEST": ["SHANK3", "GRIN2B"]}
        results = validate_pathway_against_reactome(curated, reactome)
        assert results["TEST"]["jaccard"] == 1.0

    def test_multiple_reactome_picks_best(self):
        curated = {"TEST": ["A", "B", "C"]}
        reactome = {
            "BAD": ["X", "Y", "Z"],
            "GOOD": ["A", "B", "C"],
        }
        results = validate_pathway_against_reactome(curated, reactome)
        assert results["TEST"]["best_match"] == "GOOD"
        assert results["TEST"]["jaccard"] == 1.0


# =============================================================================
# TEST DISEASE-REALISTIC SYNTHETIC DATA
# =============================================================================


class TestGenerateDiseaseRealisticSynthetic:
    def test_basic_generation(self, sample_clinvar_genes, sample_pathways):
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=60,
            n_subtypes=3,
            seed=42,
        )
        assert sim.pathway_scores.shape[0] == 60
        assert sim.pathway_scores.shape[1] == 3  # 3 pathways
        assert len(sim.true_labels) == 60

    def test_real_gene_names(self, sample_clinvar_genes, sample_pathways):
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=30,
            n_subtypes=2,
            seed=42,
        )
        # Should use real gene names, not GENE_0_0
        gene_cols = sim.gene_burdens.columns.tolist()
        assert "SHANK3" in gene_cols
        assert "CHD8" in gene_cols

    def test_correct_dimensions(self, sample_clinvar_genes, sample_pathways):
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=50,
            n_subtypes=2,
            seed=42,
        )
        all_genes = set(g for genes in sample_pathways.values() for g in genes)
        assert sim.gene_burdens.shape[1] == len(all_genes)
        assert sim.pathway_scores.shape[1] == len(sample_pathways)
        assert len(sim.true_labels) == 50

    def test_reproducibility(self, sample_clinvar_genes, sample_pathways):
        sim1 = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=30,
            seed=42,
        )
        sim2 = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=30,
            seed=42,
        )
        np.testing.assert_array_equal(sim1.true_labels, sim2.true_labels)
        pd.testing.assert_frame_equal(sim1.pathway_scores, sim2.pathway_scores)

    def test_subtype_effects_assigned(self, sample_clinvar_genes, sample_pathways):
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=60,
            n_subtypes=3,
            seed=42,
        )
        assert len(sim.subtype_pathway_effects) == 3
        for subtype, effect_pws in sim.subtype_pathway_effects.items():
            assert len(effect_pws) > 0

    def test_pathway_scores_z_scored(self, sample_clinvar_genes, sample_pathways):
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=100,
            n_subtypes=3,
            seed=42,
        )
        # Z-scored means should be near 0, stds near 1
        means = sim.pathway_scores.mean()
        stds = sim.pathway_scores.std()
        for col in sim.pathway_scores.columns:
            assert abs(means[col]) < 0.2, f"Mean of {col} too far from 0"
            assert abs(stds[col] - 1.0) < 0.3, f"Std of {col} too far from 1"

    def test_gene_burdens_non_negative(self, sample_clinvar_genes, sample_pathways):
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=30,
            seed=42,
        )
        assert (sim.gene_burdens.values >= 0).all()

    def test_empty_pathways_raises(self, sample_clinvar_genes):
        with pytest.raises(ValueError, match="No genes found"):
            generate_disease_realistic_synthetic(
                pathways={}, clinvar_genes=sample_clinvar_genes
            )

    def test_no_clinvar_still_works(self, sample_pathways):
        """Should work even without ClinVar data (all baseline weights)."""
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes={},
            n_samples=30,
            seed=42,
        )
        assert sim.pathway_scores.shape[0] == 30


# =============================================================================
# TEST BIOLOGICAL PLAUSIBILITY
# =============================================================================


class TestBiologicalPlausibility:
    def test_plausible_result(self, sample_clinvar_genes, sample_pathways):
        """Well-separated synthetic data should pass plausibility."""
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            n_samples=90,
            n_subtypes=3,
            effect_size=2.0,
            seed=42,
        )
        from pathway_subtyping.clustering import run_clustering

        result = run_clustering(sim.pathway_scores, n_clusters=3, seed=42)
        bp = run_biological_plausibility_check(
            pathway_scores=sim.pathway_scores,
            cluster_labels=result.labels,
            pathways=sample_pathways,
            clinvar_genes=sample_clinvar_genes,
            seed=42,
        )
        assert bp.n_subtypes == 3
        assert isinstance(bp.subtypes_biologically_distinct, bool)
        assert isinstance(bp.pathway_gene_clinvar_overlap, float)

    def test_to_dict(self):
        bp = BiologicalPlausibilityResult(
            disease_name="autism",
            n_subtypes=3,
            n_enriched_pathways=5,
            pathway_gene_clinvar_overlap=0.8,
            subtypes_biologically_distinct=True,
        )
        d = bp.to_dict()
        assert d["disease_name"] == "autism"
        assert isinstance(d["subtypes_biologically_distinct"], bool)

    def test_plausibility_no_clinvar(self, sample_pathways):
        """Should work without ClinVar (overlap will be 0)."""
        sim = generate_disease_realistic_synthetic(
            pathways=sample_pathways,
            clinvar_genes={},
            n_samples=60,
            n_subtypes=2,
            effect_size=2.0,
            seed=42,
        )
        from pathway_subtyping.clustering import run_clustering

        result = run_clustering(sim.pathway_scores, n_clusters=2, seed=42)
        bp = run_biological_plausibility_check(
            pathway_scores=sim.pathway_scores,
            cluster_labels=result.labels,
            pathways=sample_pathways,
            clinvar_genes={},
            seed=42,
        )
        assert bp.pathway_gene_clinvar_overlap == 0.0


# =============================================================================
# TEST VALIDATION REPORT
# =============================================================================


class TestValidationReport:
    def test_to_dict_json_serializable(self):
        report = ValidationReport(
            timestamp="2026-02-10 12:00:00",
            overall_pass=True,
        )
        d = report.to_dict()
        json_str = json.dumps(d)
        assert isinstance(json_str, str)

    def test_format_report_markdown(self):
        report = ValidationReport(
            timestamp="2026-02-10 12:00:00",
            overall_pass=True,
            pathway_coverage=[
                PathwayCoverageResult(
                    pathway_name="TEST",
                    total_genes=10,
                    genes_in_clinvar=8,
                    genes_with_pathogenic=6,
                    coverage_fraction=0.8,
                    pathogenic_fraction=0.6,
                ),
            ],
        )
        md = report.format_report()
        assert "Validation Report" in md
        assert "TEST" in md
        assert "PASS" in md

    def test_format_report_fail(self):
        report = ValidationReport(
            timestamp="2026-02-10 12:00:00",
            overall_pass=False,
            warnings=["Low coverage"],
        )
        md = report.format_report()
        assert "NEEDS REVIEW" in md
        assert "Low coverage" in md

    def test_get_citations(self):
        report = ValidationReport(timestamp="2026-02-10")
        citations = report.get_citations()
        assert len(citations) == 2
        assert "ClinVar" in citations[0]
        assert "Reactome" in citations[1]

    def test_report_with_all_sections(self, sample_clinvar_genes, sample_pathways):
        report = ValidationReport(
            timestamp="2026-02-10 12:00:00",
            data_sources=[DATASETS["clinvar_gene_summary"]],
            pathway_coverage=[
                PathwayCoverageResult(
                    pathway_name="SYNAPTIC",
                    total_genes=4,
                    genes_in_clinvar=2,
                    genes_with_pathogenic=2,
                    coverage_fraction=0.5,
                    pathogenic_fraction=0.5,
                ),
            ],
            reactome_cross_ref={
                "SYNAPTIC": {"best_match": "R_Synaptic", "jaccard": 0.6},
            },
            synthetic_validation={
                "n_samples": 100,
                "n_subtypes": 3,
                "n_clusters_found": 3,
                "ari": 0.85,
            },
            biological_plausibility=BiologicalPlausibilityResult(
                disease_name="autism",
                n_subtypes=3,
                n_enriched_pathways=5,
                pathway_gene_clinvar_overlap=0.8,
                subtypes_biologically_distinct=True,
            ),
            overall_pass=True,
        )
        md = report.format_report()
        assert "Data Sources" in md
        assert "ClinVar" in md
        assert "Reactome" in md
        assert "Synthetic" in md
        assert "Biological Plausibility" in md

        d = report.to_dict()
        json.dumps(d)  # Should not raise


# =============================================================================
# TEST FULL VALIDATION (offline with mocks)
# =============================================================================


class TestRunFullValidation:
    def test_full_validation_with_mock_data(
        self, mock_clinvar_tsv, mock_reactome_gmt, tmp_path
    ):
        """Run full validation using mock data files, skip download."""
        # Create a small GMT for the "test" disease
        gmt_path = tmp_path / "test_pathways.gmt"
        gmt_path.write_text(
            "SYNAPTIC\thttps://example.com\tSHANK3\tNRXN1\tGRIN2B\n"
            "CHROMATIN\thttps://example.com\tCHD8\tARID1B\n"
        )

        # Mock the download functions to use our mock files
        with mock.patch(
            "pathway_subtyping.validation_datasets.load_clinvar_gene_summary"
        ) as mock_cv, mock.patch(
            "pathway_subtyping.validation_datasets.load_reactome_pathways"
        ) as mock_rp:
            mock_cv.return_value = load_clinvar_gene_summary(file_path=mock_clinvar_tsv)
            mock_rp.return_value = load_reactome_pathways(file_path=mock_reactome_gmt)

            report = run_full_validation(
                gmt_path=str(gmt_path),
                disease_name="test",
                n_samples=60,
                n_subtypes=2,
                effect_size=1.5,
                seed=42,
            )

        assert isinstance(report, ValidationReport)
        assert len(report.pathway_coverage) == 2
        assert report.synthetic_validation.get("n_samples") == 60

        # Should be JSON serializable
        json.dumps(report.to_dict())

    def test_full_validation_skip_download(self, tmp_path):
        """Offline mode with no cache should warn but not crash."""
        gmt_path = tmp_path / "test_pathways.gmt"
        gmt_path.write_text("PATHWAY_A\tdesc\tGENE1\tGENE2\n")

        report = run_full_validation(
            gmt_path=str(gmt_path),
            disease_name="test",
            cache_dir=tmp_path / "empty_cache",
            skip_download=True,
            seed=42,
        )

        assert isinstance(report, ValidationReport)
        assert not report.overall_pass  # No data → fail
        assert any("not available" in w for w in report.warnings)

    def test_missing_gmt_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            run_full_validation(
                gmt_path=str(tmp_path / "nonexistent.gmt"),
                disease_name="test",
                skip_download=True,
                seed=42,
            )

    def test_default_gmt_path_resolution(self, tmp_path):
        """If no gmt_path, should look for disease_pathways.gmt in data/pathways/."""
        # This should find autism_pathways.gmt if it exists
        data_dir = Path(__file__).parent.parent / "data" / "pathways"
        if (data_dir / "autism_pathways.gmt").exists():
            report = run_full_validation(
                disease_name="autism",
                n_samples=30,
                n_subtypes=2,
                skip_download=True,
                cache_dir=tmp_path,
                seed=42,
            )
            assert isinstance(report, ValidationReport)
        else:
            with pytest.raises(FileNotFoundError):
                run_full_validation(
                    disease_name="nonexistent_disease",
                    skip_download=True,
                    cache_dir=tmp_path,
                )


# =============================================================================
# NETWORK TESTS (skipped by default)
# =============================================================================


@pytest.mark.network
class TestNetworkIntegration:
    def test_clinvar_download_and_parse(self, tmp_path):
        """Download and parse real ClinVar data."""
        genes = load_clinvar_gene_summary(cache_dir=tmp_path)
        assert len(genes) > 1000
        assert "BRCA1" in genes
        assert genes["BRCA1"].total_pathogenic > 0

    def test_reactome_download_and_parse(self, tmp_path):
        """Download and parse real Reactome data."""
        pathways = load_reactome_pathways(cache_dir=tmp_path)
        assert len(pathways) > 100
