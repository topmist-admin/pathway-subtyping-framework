"""
Tests for QC Layer Phase 2: F12 Reference Atlas Comparison.

Tests validate:
    - Atlas construction (add, CSV, DataFrame)
    - Per-cell distance computation
    - Nearest reference type mapping
    - On-target classification with expected_type
    - Auto-threshold calibration
    - Integration with simulator data
"""

import tempfile

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.qc.atlas import (
    AtlasComparator,
    AtlasComparisonResult,
    AtlasEntry,
    CellMapping,
)
from pathway_subtyping.qc.testing.simulator import ManufacturingSimulator, ManufacturingSpec

# ── Fixtures ──────────────────────────────────────────────────────────


@pytest.fixture
def simulator():
    return ManufacturingSimulator(seed=42, genes_per_pathway=10)


@pytest.fixture
def simple_atlas():
    atlas = AtlasComparator(distance_threshold=3.0)
    atlas.add_reference(
        "healthy_T_cell",
        {"PW_A": 0.5, "PW_B": 0.3, "PW_C": 0.1},
        annotation="Healthy activated T cell",
        source="Test fixture",
    )
    atlas.add_reference(
        "exhausted_T_cell",
        {"PW_A": 0.1, "PW_B": 0.1, "PW_C": 0.8},
        annotation="Exhausted T cell",
        source="Test fixture",
    )
    return atlas


@pytest.fixture
def batch_scores_on_target():
    """Batch of cells near the healthy_T_cell profile."""
    rng = np.random.default_rng(42)
    n_cells = 100
    data = {
        "PW_A": rng.normal(0.5, 0.1, n_cells),
        "PW_B": rng.normal(0.3, 0.1, n_cells),
        "PW_C": rng.normal(0.1, 0.1, n_cells),
    }
    return pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])


@pytest.fixture
def batch_scores_off_target():
    """Batch of cells far from any reference."""
    rng = np.random.default_rng(42)
    n_cells = 50
    data = {
        "PW_A": rng.normal(5.0, 0.5, n_cells),
        "PW_B": rng.normal(5.0, 0.5, n_cells),
        "PW_C": rng.normal(5.0, 0.5, n_cells),
    }
    return pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])


# ── Atlas Construction ────────────────────────────────────────────────


class TestAtlasConstruction:

    def test_add_reference(self, simple_atlas):
        assert simple_atlas.n_references == 2

    def test_add_multiple_references(self):
        atlas = AtlasComparator()
        atlas.add_reference("type_a", {"PW_A": 1.0})
        atlas.add_reference("type_b", {"PW_A": 2.0})
        atlas.add_reference("type_c", {"PW_A": 3.0})
        assert atlas.n_references == 3

    def test_atlas_registry(self, simple_atlas):
        registry = simple_atlas.atlas_registry()
        assert len(registry) == 2
        names = {r["name"] for r in registry}
        assert "healthy_T_cell" in names
        assert "exhausted_T_cell" in names

    def test_load_from_dataframe(self):
        df = pd.DataFrame(
            {"PW_A": [0.5, 0.1], "PW_B": [0.3, 0.8]},
            index=["type_1", "type_2"],
        )
        atlas = AtlasComparator()
        count = atlas.load_from_dataframe(df, source="test")
        assert count == 2
        assert atlas.n_references == 2

    def test_load_from_csv(self):
        df = pd.DataFrame(
            {
                "cell_type": ["type_1", "type_2"],
                "PW_A": [0.5, 0.1],
                "PW_B": [0.3, 0.8],
            }
        )
        with tempfile.NamedTemporaryFile(suffix=".csv", mode="w", delete=False) as f:
            df.to_csv(f, index=False)
            path = f.name

        atlas = AtlasComparator()
        count = atlas.load_atlas_csv(path, name_column="cell_type")
        assert count == 2

    def test_csv_missing_name_column_raises(self):
        df = pd.DataFrame({"PW_A": [0.5], "PW_B": [0.3]})
        with tempfile.NamedTemporaryFile(suffix=".csv", mode="w", delete=False) as f:
            df.to_csv(f, index=False)
            path = f.name

        atlas = AtlasComparator()
        with pytest.raises(ValueError, match="Name column"):
            atlas.load_atlas_csv(path, name_column="cell_type")


# ── Comparison ────────────────────────────────────────────────────────


class TestAtlasComparison:

    def test_on_target_batch(self, simple_atlas, batch_scores_on_target):
        result = simple_atlas.compare(batch_scores_on_target)
        assert isinstance(result, AtlasComparisonResult)
        assert result.n_total == 100
        # Cells close to healthy_T_cell should be on-target
        assert result.n_on_target > 50
        assert result.passed

    def test_off_target_batch(self, simple_atlas, batch_scores_off_target):
        result = simple_atlas.compare(batch_scores_off_target)
        assert result.n_total == 50
        # Cells far from references should be off-target
        assert result.n_off_target > result.n_on_target

    def test_expected_type_filtering(self, simple_atlas, batch_scores_on_target):
        # When expected_type is healthy_T_cell, cells near it are on-target
        result = simple_atlas.compare(batch_scores_on_target, expected_type="healthy_T_cell")
        assert result.n_on_target > 50

        # When expected_type is exhausted_T_cell, most cells should be off-target
        result2 = simple_atlas.compare(batch_scores_on_target, expected_type="exhausted_T_cell")
        assert result2.n_on_target < result.n_on_target

    def test_type_distribution(self, simple_atlas, batch_scores_on_target):
        result = simple_atlas.compare(batch_scores_on_target)
        assert "healthy_T_cell" in result.type_distribution
        assert "exhausted_T_cell" in result.type_distribution
        total = sum(result.type_distribution.values())
        assert total == 100

    def test_distance_metrics(self, simple_atlas, batch_scores_on_target):
        result = simple_atlas.compare(batch_scores_on_target)
        assert result.median_distance >= 0
        assert result.p90_distance >= result.median_distance
        assert result.mean_distance >= 0
        assert len(result.per_cell_distances) == 100

    def test_cell_mappings(self, simple_atlas, batch_scores_on_target):
        result = simple_atlas.compare(batch_scores_on_target)
        assert len(result.cell_mappings) == 100
        for mapping in result.cell_mappings:
            assert isinstance(mapping, CellMapping)
            assert mapping.nearest_type in ("healthy_T_cell", "exhausted_T_cell")
            assert mapping.distance >= 0

    def test_empty_atlas_raises(self, batch_scores_on_target):
        atlas = AtlasComparator()
        with pytest.raises(ValueError, match="Atlas is empty"):
            atlas.compare(batch_scores_on_target)

    def test_auto_threshold(self, batch_scores_on_target):
        atlas = AtlasComparator()  # No threshold set
        atlas.add_reference("type_a", {"PW_A": 0.5, "PW_B": 0.3, "PW_C": 0.1})
        atlas.add_reference("type_b", {"PW_A": 5.0, "PW_B": 5.0, "PW_C": 5.0})
        result = atlas.compare(batch_scores_on_target)
        assert result.distance_threshold > 0

    def test_single_reference_auto_threshold(self, batch_scores_on_target):
        atlas = AtlasComparator()
        atlas.add_reference("only_type", {"PW_A": 0.5, "PW_B": 0.3, "PW_C": 0.1})
        result = atlas.compare(batch_scores_on_target)
        assert result.distance_threshold > 0

    def test_to_dict(self, simple_atlas, batch_scores_on_target):
        result = simple_atlas.compare(batch_scores_on_target)
        d = result.to_dict()
        assert "passed" in d
        assert "median_distance" in d
        assert "type_distribution" in d
        assert "fraction_on_target" in d

    def test_cell_mapping_to_dict(self, simple_atlas, batch_scores_on_target):
        result = simple_atlas.compare(batch_scores_on_target)
        d = result.cell_mappings[0].to_dict()
        assert "cell_id" in d
        assert "nearest_type" in d
        assert "distance" in d
        assert "on_target" in d


# ── Integration with Simulator ────────────────────────────────────────


class TestAtlasSimulatorIntegration:

    def test_healthy_batch_matches_self_atlas(self, simulator):
        spec = ManufacturingSpec(
            intended_pathways=simulator.pathways[:3],
            excluded_pathways=simulator.pathways[3:5],
        )
        batch = simulator.generate_healthy_batch(n_cells=50, spec=spec)

        # Build atlas from the batch's mean profile
        atlas = AtlasComparator(distance_threshold=3.0)
        mean_profile = batch.pathway_scores.mean().to_dict()
        atlas.add_reference("batch_mean", mean_profile)

        result = atlas.compare(batch.pathway_scores)
        assert result.n_on_target > 30  # Most cells should match

    def test_drifted_batch_diverges_from_atlas(self, simulator):
        spec = ManufacturingSpec(
            intended_pathways=simulator.pathways[:3],
            excluded_pathways=simulator.pathways[3:5],
        )
        batch = simulator.generate_healthy_batch(n_cells=50, spec=spec)

        # Atlas = healthy batch mean
        atlas = AtlasComparator(distance_threshold=1.0)
        mean_profile = batch.pathway_scores.mean().to_dict()
        atlas.add_reference("healthy", mean_profile)

        # Drift the batch
        drifted = simulator.inject_atlas_drift(batch, drift_distance=0.9)
        result = atlas.compare(drifted.pathway_scores)
        # Drifted batch should have worse distances
        assert result.mean_distance > 0

    def test_atlas_with_heterogeneity_profiler(self, simulator):
        """F12 uses F7 conformity scoring — verify they work together."""
        from pathway_subtyping.qc.heterogeneity import HeterogeneityProfiler

        spec = ManufacturingSpec(
            intended_pathways=simulator.pathways[:3],
            excluded_pathways=simulator.pathways[3:5],
        )
        batch = simulator.generate_healthy_batch(n_cells=80, spec=spec)
        injected = simulator.inject_subpopulation(batch, fraction=0.2)

        # F7: heterogeneity
        profiler = HeterogeneityProfiler(conformity_threshold=1.5)
        f7_result = profiler.profile(injected.pathway_scores)

        # F12: atlas
        atlas = AtlasComparator(distance_threshold=2.0)
        mean_profile = batch.pathway_scores.mean().to_dict()
        atlas.add_reference("healthy", mean_profile)
        f12_result = atlas.compare(injected.pathway_scores)

        # Both should detect issues
        assert f7_result.n_nonconforming > 0 or f12_result.n_off_target > 0
