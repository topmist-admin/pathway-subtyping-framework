"""
Tests for the cross_cohort module.

Tests cover:
- CohortResult dataclass
- CrossCohortResult dataclass
- load_cohort_result function
- compare_cohorts function
- Transfer learning validation
- Projection validation
- Pathway correlation
- Shared subtype detection
- Batch comparison
- Report generation
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.cross_cohort import (
    CohortResult,
    CrossCohortResult,
    _find_shared_subtypes,
    _pathway_correlation,
    _projection_validation,
    _transfer_learning_validation,
    batch_compare_cohorts,
    compare_cohorts,
    generate_cross_cohort_report,
    load_cohort_result,
)


class TestCohortResult:
    """Tests for CohortResult dataclass."""

    def test_create_cohort_result(self):
        """Test creating a CohortResult."""
        scores = pd.DataFrame(
            {"pathway1": [1.0, 2.0, 3.0], "pathway2": [0.5, 1.5, 2.5]},
            index=["S1", "S2", "S3"],
        )
        result = CohortResult(
            name="test_cohort",
            pathway_scores=scores,
            cluster_labels=np.array([0, 0, 1]),
            cluster_names={0: "synaptic", 1: "chromatin"},
            n_samples=3,
            n_clusters=2,
        )

        assert result.name == "test_cohort"
        assert result.n_samples == 3
        assert result.n_clusters == 2
        assert len(result.cluster_labels) == 3


class TestCrossCohortResult:
    """Tests for CrossCohortResult dataclass."""

    def test_create_cross_cohort_result(self):
        """Test creating a CrossCohortResult."""
        result = CrossCohortResult(
            cohort_a="cohort1",
            cohort_b="cohort2",
            transfer_ari=0.75,
            projection_ari=0.68,
            pathway_correlation=0.85,
            shared_subtypes=["synaptic", "chromatin"],
            details={"common_pathways": 10},
        )

        assert result.cohort_a == "cohort1"
        assert result.transfer_ari == 0.75
        assert len(result.shared_subtypes) == 2


class TestLoadCohortResult:
    """Tests for load_cohort_result function."""

    @pytest.fixture
    def mock_pipeline_output(self, tmp_path):
        """Create a mock pipeline output directory."""
        # Create pathway scores
        scores = pd.DataFrame(
            {"pathway1": [1.0, 2.0, 3.0], "pathway2": [0.5, 1.5, 2.5]},
            index=["S1", "S2", "S3"],
        )
        scores.to_csv(tmp_path / "pathway_scores.csv")

        # Create cluster assignments
        assignments = pd.DataFrame(
            {
                "sample_id": ["S1", "S2", "S3"],
                "cluster_id": [0, 0, 1],
                "cluster_label": ["synaptic", "synaptic", "chromatin"],
            }
        )
        assignments.to_csv(tmp_path / "subtype_assignments.csv", index=False)

        # Create report
        report = {"pipeline_name": "test_run", "summary": {"n_clusters": 2}}
        with open(tmp_path / "report.json", "w") as f:
            json.dump(report, f)

        return tmp_path

    def test_load_cohort_result(self, mock_pipeline_output):
        """Test loading a cohort result from directory."""
        result = load_cohort_result(str(mock_pipeline_output))

        assert result.name == "test_run"
        assert result.n_samples == 3
        assert result.n_clusters == 2
        assert len(result.pathway_scores.columns) == 2

    def test_load_cohort_result_missing_scores(self, tmp_path):
        """Test error when pathway scores are missing."""
        with pytest.raises(FileNotFoundError):
            load_cohort_result(str(tmp_path))

    def test_load_cohort_result_missing_assignments(self, tmp_path):
        """Test error when assignments are missing."""
        # Create only pathway scores
        scores = pd.DataFrame({"pathway1": [1.0]}, index=["S1"])
        scores.to_csv(tmp_path / "pathway_scores.csv")

        with pytest.raises(FileNotFoundError):
            load_cohort_result(str(tmp_path))


class TestTransferLearningValidation:
    """Tests for transfer learning validation."""

    def test_transfer_perfect_match(self):
        """Test transfer learning with identical cohorts."""
        np.random.seed(42)

        # Create two cohorts with identical structure
        scores_a = pd.DataFrame(
            np.random.randn(30, 5),
            columns=[f"p{i}" for i in range(5)],
        )
        # Create clear cluster structure
        scores_a.iloc[:15, :] += 3  # Group 1
        labels_a = np.array([0] * 15 + [1] * 15)

        # Identical cohort B
        scores_b = scores_a.copy()
        labels_b = labels_a.copy()

        ari = _transfer_learning_validation(
            scores_a, labels_a, scores_b, labels_b, n_clusters=2, seed=42
        )

        # Should have high ARI for identical cohorts
        assert ari > 0.8

    def test_transfer_random_labels(self):
        """Test transfer learning with random labels (should have low ARI)."""
        np.random.seed(42)

        scores_a = pd.DataFrame(
            np.random.randn(30, 5),
            columns=[f"p{i}" for i in range(5)],
        )
        labels_a = np.array([0] * 15 + [1] * 15)

        scores_b = pd.DataFrame(
            np.random.randn(30, 5),
            columns=[f"p{i}" for i in range(5)],
        )
        # Random labels for B
        labels_b = np.random.randint(0, 2, 30)

        ari = _transfer_learning_validation(
            scores_a, labels_a, scores_b, labels_b, n_clusters=2, seed=42
        )

        # Random labels should have ARI close to 0
        assert ari < 0.5


class TestProjectionValidation:
    """Tests for projection validation."""

    def test_projection_similar_structure(self):
        """Test projection with similar cluster structure."""
        np.random.seed(42)

        # Create cohort A with clear clusters
        scores_a = pd.DataFrame(
            np.vstack(
                [
                    np.random.randn(15, 5) + [3, 0, 0, 0, 0],
                    np.random.randn(15, 5) + [0, 3, 0, 0, 0],
                ]
            ),
            columns=[f"p{i}" for i in range(5)],
        )
        labels_a = np.array([0] * 15 + [1] * 15)

        # Create cohort B with similar structure
        scores_b = pd.DataFrame(
            np.vstack(
                [
                    np.random.randn(10, 5) + [3, 0, 0, 0, 0],
                    np.random.randn(10, 5) + [0, 3, 0, 0, 0],
                ]
            ),
            columns=[f"p{i}" for i in range(5)],
        )
        labels_b = np.array([0] * 10 + [1] * 10)

        ari = _projection_validation(scores_a, labels_a, scores_b, labels_b, seed=42)

        # Should recover similar structure
        assert ari > 0.5


class TestPathwayCorrelation:
    """Tests for pathway correlation."""

    def test_pathway_correlation_identical(self):
        """Test correlation with identical variance patterns."""
        scores_a = pd.DataFrame(
            {"p1": [1, 2, 3, 4], "p2": [0.1, 0.2, 0.3, 0.4]},
        )
        scores_b = pd.DataFrame(
            {"p1": [2, 4, 6, 8], "p2": [0.2, 0.4, 0.6, 0.8]},
        )

        corr = _pathway_correlation(scores_a, scores_b)
        assert corr == 1.0

    def test_pathway_correlation_different(self):
        """Test correlation with different variance patterns."""
        np.random.seed(42)
        scores_a = pd.DataFrame(
            {"p1": np.random.randn(20), "p2": np.random.randn(20) * 0.1},
        )
        scores_b = pd.DataFrame(
            {"p1": np.random.randn(20) * 0.1, "p2": np.random.randn(20)},
        )

        corr = _pathway_correlation(scores_a, scores_b)
        # Inverted variance pattern should have negative correlation
        assert corr < 0


class TestFindSharedSubtypes:
    """Tests for shared subtype detection."""

    def test_find_shared_subtypes_overlap(self):
        """Test finding shared subtypes with overlap."""
        names_a = {0: "synaptic", 1: "chromatin", 2: "immune"}
        names_b = {0: "synaptic", 1: "mtor", 2: "chromatin"}

        shared = _find_shared_subtypes(names_a, names_b)
        assert set(shared) == {"synaptic", "chromatin"}

    def test_find_shared_subtypes_no_overlap(self):
        """Test finding shared subtypes with no overlap."""
        names_a = {0: "synaptic", 1: "chromatin"}
        names_b = {0: "mtor", 1: "immune"}

        shared = _find_shared_subtypes(names_a, names_b)
        assert shared == []


class TestCompareCohorts:
    """Tests for compare_cohorts function."""

    @pytest.fixture
    def two_cohorts(self):
        """Create two test cohorts with similar structure."""
        np.random.seed(42)

        # Cohort A
        scores_a = pd.DataFrame(
            np.vstack(
                [
                    np.random.randn(20, 5) + [3, 0, 0, 0, 0],
                    np.random.randn(20, 5) + [0, 3, 0, 0, 0],
                ]
            ),
            columns=["pathway1", "pathway2", "pathway3", "pathway4", "pathway5"],
            index=[f"A{i}" for i in range(40)],
        )
        cohort_a = CohortResult(
            name="cohort_a",
            pathway_scores=scores_a,
            cluster_labels=np.array([0] * 20 + [1] * 20),
            cluster_names={0: "synaptic", 1: "chromatin"},
            n_samples=40,
            n_clusters=2,
        )

        # Cohort B with similar structure
        scores_b = pd.DataFrame(
            np.vstack(
                [
                    np.random.randn(15, 5) + [3, 0, 0, 0, 0],
                    np.random.randn(15, 5) + [0, 3, 0, 0, 0],
                ]
            ),
            columns=["pathway1", "pathway2", "pathway3", "pathway4", "pathway5"],
            index=[f"B{i}" for i in range(30)],
        )
        cohort_b = CohortResult(
            name="cohort_b",
            pathway_scores=scores_b,
            cluster_labels=np.array([0] * 15 + [1] * 15),
            cluster_names={0: "synaptic", 1: "immune"},
            n_samples=30,
            n_clusters=2,
        )

        return cohort_a, cohort_b

    def test_compare_cohorts(self, two_cohorts):
        """Test comparing two cohorts."""
        cohort_a, cohort_b = two_cohorts
        result = compare_cohorts(cohort_a, cohort_b, seed=42)

        assert result.cohort_a == "cohort_a"
        assert result.cohort_b == "cohort_b"
        assert -1 <= result.transfer_ari <= 1
        assert -1 <= result.projection_ari <= 1
        assert -1 <= result.pathway_correlation <= 1
        assert "synaptic" in result.shared_subtypes

    def test_compare_cohorts_insufficient_pathways(self):
        """Test error when cohorts have no common pathways."""
        scores_a = pd.DataFrame({"p1": [1, 2], "p2": [3, 4]})
        scores_b = pd.DataFrame({"p3": [1, 2], "p4": [3, 4]})

        cohort_a = CohortResult(
            name="a",
            pathway_scores=scores_a,
            cluster_labels=np.array([0, 1]),
            cluster_names={0: "x"},
            n_samples=2,
            n_clusters=1,
        )
        cohort_b = CohortResult(
            name="b",
            pathway_scores=scores_b,
            cluster_labels=np.array([0, 1]),
            cluster_names={0: "y"},
            n_samples=2,
            n_clusters=1,
        )

        with pytest.raises(ValueError, match="Insufficient common pathways"):
            compare_cohorts(cohort_a, cohort_b)


class TestBatchCompareCohorts:
    """Tests for batch_compare_cohorts function."""

    @pytest.fixture
    def three_pipeline_outputs(self, tmp_path):
        """Create three mock pipeline outputs."""
        dirs = []
        for i in range(3):
            cohort_dir = tmp_path / f"cohort_{i}"
            cohort_dir.mkdir()

            # Pathway scores with same column names
            np.random.seed(i)
            scores = pd.DataFrame(
                np.random.randn(20, 3),
                columns=["pathway1", "pathway2", "pathway3"],
                index=[f"S{j}" for j in range(20)],
            )
            scores.to_csv(cohort_dir / "pathway_scores.csv")

            # Assignments
            assignments = pd.DataFrame(
                {
                    "sample_id": [f"S{j}" for j in range(20)],
                    "cluster_id": [j % 2 for j in range(20)],
                    "cluster_label": ["type_a" if j % 2 == 0 else "type_b" for j in range(20)],
                }
            )
            assignments.to_csv(cohort_dir / "subtype_assignments.csv", index=False)

            # Report
            report = {"pipeline_name": f"cohort_{i}"}
            with open(cohort_dir / "report.json", "w") as f:
                json.dump(report, f)

            dirs.append(str(cohort_dir))

        return dirs

    def test_batch_compare_cohorts(self, three_pipeline_outputs):
        """Test batch comparison of multiple cohorts."""
        results = batch_compare_cohorts(three_pipeline_outputs, seed=42)

        # 3 cohorts = 3 pairwise comparisons
        assert len(results) == 3

    def test_batch_compare_insufficient_cohorts(self, tmp_path):
        """Test error with fewer than 2 cohorts."""
        with pytest.raises(ValueError, match="at least 2 cohorts"):
            batch_compare_cohorts([str(tmp_path)])


class TestGenerateCrossCohortReport:
    """Tests for report generation."""

    def test_generate_report(self, tmp_path):
        """Test generating a cross-cohort report."""
        results = [
            CrossCohortResult(
                cohort_a="cohort1",
                cohort_b="cohort2",
                transfer_ari=0.75,
                projection_ari=0.68,
                pathway_correlation=0.85,
                shared_subtypes=["synaptic"],
                details={
                    "common_pathways": 10,
                    "cohort_a_samples": 50,
                    "cohort_b_samples": 40,
                    "cohort_a_clusters": 3,
                    "cohort_b_clusters": 3,
                },
            )
        ]

        report_path = tmp_path / "cross_cohort_report.md"
        generate_cross_cohort_report(results, str(report_path))

        assert report_path.exists()
        content = report_path.read_text()
        assert "Cross-Cohort Validation Report" in content
        assert "cohort1" in content
        assert "0.75" in content
        assert "synaptic" in content
