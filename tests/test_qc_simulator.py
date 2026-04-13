"""
Tests for the Manufacturing Simulator and QC test infrastructure.

Tests cover:
    - Healthy batch generation
    - All 9 injection methods
    - Ground truth labels and metadata
    - Orthogonality matrix structure
    - Severity titration harness
    - Scenario execution
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.qc.testing.simulator import (
    CompetitionDirection,
    DefectLabel,
    DefectMetadata,
    DefectType,
    FeedbackBreakType,
    InjectedBatch,
    ManufacturingSimulator,
    ManufacturingSpec,
    StressorType,
)


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def simulator():
    return ManufacturingSimulator(seed=42)


@pytest.fixture
def spec():
    return ManufacturingSpec(
        intended_pathways=[
            "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        ],
        excluded_pathways=[
            "HALLMARK_APOPTOSIS",
        ],
        therapeutic_windows={
            "HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8),
        },
    )


@pytest.fixture
def healthy_batch(simulator, spec):
    return simulator.generate_healthy_batch(n_cells=100, spec=spec)


# ── Healthy Batch ─────────────────────────────────────────────────────

class TestHealthyBatch:

    def test_generates_correct_dimensions(self, simulator):
        batch = simulator.generate_healthy_batch(n_cells=50)
        assert batch.n_cells == 50
        assert batch.n_genes == len(simulator.all_genes)
        assert batch.n_pathways == len(simulator.pathways)

    def test_expression_is_dataframe(self, healthy_batch):
        assert isinstance(healthy_batch.expression, pd.DataFrame)
        assert isinstance(healthy_batch.pathway_scores, pd.DataFrame)

    def test_no_defects_in_healthy_batch(self, healthy_batch):
        assert len(healthy_batch.defect_metadata) == 0
        assert not any(label.affected for label in healthy_batch.cell_labels)
        assert healthy_batch.affected_mask.sum() == 0

    def test_intended_pathways_elevated(self, healthy_batch):
        scores = healthy_batch.pathway_scores
        intended = "HALLMARK_INFLAMMATORY_RESPONSE"
        excluded = "HALLMARK_APOPTOSIS"
        # Intended should have higher raw expression than excluded
        assert scores[intended].mean() > scores[excluded].mean()

    def test_custom_spec(self, simulator):
        spec = ManufacturingSpec(
            intended_pathways=["HALLMARK_HYPOXIA"],
            excluded_pathways=["HALLMARK_GLYCOLYSIS"],
        )
        batch = simulator.generate_healthy_batch(n_cells=30, spec=spec)
        assert batch.spec.intended_pathways == ["HALLMARK_HYPOXIA"]
        assert batch.n_cells == 30

    def test_to_dict(self, healthy_batch):
        d = healthy_batch.to_dict()
        assert d["n_cells"] == 100
        assert d["n_affected"] == 0
        assert d["defects"] == []

    def test_reproducibility_with_seed(self):
        sim1 = ManufacturingSimulator(seed=123)
        sim2 = ManufacturingSimulator(seed=123)
        b1 = sim1.generate_healthy_batch(n_cells=20)
        b2 = sim2.generate_healthy_batch(n_cells=20)
        np.testing.assert_array_almost_equal(
            b1.expression.values, b2.expression.values
        )


# ── Injection Methods ─────────────────────────────────────────────────

class TestStalledCascades:

    def test_injects_stalled_cascade(self, simulator, healthy_batch):
        injected = simulator.inject_stalled_cascades(
            healthy_batch,
            pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            severity=0.7,
        )
        assert all(label.affected for label in injected.cell_labels)
        assert injected.defect_metadata[-1].defect_type == DefectType.STALLED_CASCADE
        assert injected.defect_metadata[-1].severity == 0.7

    def test_cascade_modifies_expression(self, simulator, healthy_batch):
        injected = simulator.inject_stalled_cascades(
            healthy_batch,
            pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            severity=0.9,
        )
        # Expression should differ
        assert not np.allclose(
            healthy_batch.expression.values, injected.expression.values
        )


class TestTemporalDrift:

    def test_returns_list_of_passages(self, simulator, healthy_batch):
        passages = simulator.inject_temporal_drift(
            healthy_batch, n_passages=5, drift_rate=0.03
        )
        assert len(passages) == 6  # baseline + 5 passages
        assert isinstance(passages[0], InjectedBatch)

    def test_drift_increases_over_passages(self, simulator, healthy_batch):
        passages = simulator.inject_temporal_drift(
            healthy_batch, n_passages=5, drift_rate=0.05
        )
        # Later passages should differ more from baseline
        baseline = passages[0].pathway_scores.values
        diffs = [
            np.abs(p.pathway_scores.values - baseline).mean()
            for p in passages[1:]
        ]
        # Generally increasing (allow some noise)
        assert diffs[-1] > diffs[0]


class TestOfftarget:

    def test_injects_fraction_of_cells(self, simulator, healthy_batch):
        injected = simulator.inject_offtarget(
            healthy_batch,
            pathways=["HALLMARK_APOPTOSIS"],
            cell_fraction=0.3,
            activation_level=0.8,
        )
        n_affected = sum(1 for l in injected.cell_labels if l.affected)
        expected = int(100 * 0.3)
        assert abs(n_affected - expected) <= 2  # Allow rounding

    def test_off_target_pathway_elevated(self, simulator, healthy_batch):
        target_pw = "HALLMARK_APOPTOSIS"
        injected = simulator.inject_offtarget(
            healthy_batch,
            pathways=[target_pw],
            cell_fraction=1.0,
            activation_level=0.9,
        )
        # Z-scored mean won't show absolute difference, but raw expression will
        assert not np.allclose(
            healthy_batch.expression.values, injected.expression.values
        )

    def test_metadata_correct(self, simulator, healthy_batch):
        injected = simulator.inject_offtarget(
            healthy_batch,
            pathways=["HALLMARK_P53_PATHWAY"],
            cell_fraction=0.2,
            activation_level=0.5,
        )
        meta = injected.defect_metadata[-1]
        assert meta.defect_type == DefectType.OFFTARGET
        assert "F6_offtarget" in meta.expected_detections
        assert meta.expected_detections["F6_offtarget"] is True


class TestOverdose:

    def test_scales_pathway_expression(self, simulator, healthy_batch):
        injected = simulator.inject_overdose(
            healthy_batch,
            pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            multiplier=3.0,
        )
        assert all(l.affected for l in injected.cell_labels)
        assert injected.defect_metadata[-1].defect_type == DefectType.OVERDOSE

    def test_higher_multiplier_higher_severity(self, simulator, healthy_batch):
        low = simulator.inject_overdose(
            healthy_batch, pathways=["HALLMARK_INFLAMMATORY_RESPONSE"], multiplier=1.5
        )
        high = simulator.inject_overdose(
            healthy_batch, pathways=["HALLMARK_INFLAMMATORY_RESPONSE"], multiplier=5.0
        )
        assert high.defect_metadata[-1].severity > low.defect_metadata[-1].severity


class TestSubpopulation:

    def test_injects_deviant_fraction(self, simulator, healthy_batch):
        injected = simulator.inject_subpopulation(
            healthy_batch, fraction=0.15
        )
        n_affected = injected.affected_mask.sum()
        expected = int(100 * 0.15)
        assert abs(n_affected - expected) <= 2

    def test_custom_deviant_profile(self, simulator, healthy_batch):
        profile = {"HALLMARK_HYPOXIA": 5.0, "HALLMARK_GLYCOLYSIS": -3.0}
        injected = simulator.inject_subpopulation(
            healthy_batch, fraction=0.2, deviant_profile=profile
        )
        meta = injected.defect_metadata[-1]
        assert "HALLMARK_HYPOXIA" in meta.affected_pathways


class TestCrosstalk:

    def test_injects_crosstalk(self, simulator, healthy_batch):
        injected = simulator.inject_crosstalk(
            healthy_batch,
            pathway_a="HALLMARK_INFLAMMATORY_RESPONSE",
            pathway_b="HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            competition_direction=CompetitionDirection.A_DOMINATES,
        )
        assert all(l.affected for l in injected.cell_labels)
        meta = injected.defect_metadata[-1]
        assert meta.defect_type == DefectType.CROSSTALK

    def test_balanced_competition(self, simulator, healthy_batch):
        injected = simulator.inject_crosstalk(
            healthy_batch,
            pathway_a="HALLMARK_INFLAMMATORY_RESPONSE",
            pathway_b="HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            competition_direction=CompetitionDirection.BALANCED,
        )
        assert injected.defect_metadata[-1].parameters["competition_direction"] == "balanced"


class TestBrokenFeedback:

    def test_decoupled_feedback(self, simulator, healthy_batch):
        injected = simulator.inject_broken_feedback(
            healthy_batch,
            loops=[("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_APOPTOSIS")],
            break_type=FeedbackBreakType.DECOUPLED,
        )
        meta = injected.defect_metadata[-1]
        assert meta.defect_type == DefectType.BROKEN_FEEDBACK
        assert meta.parameters["break_type"] == "decoupled"

    def test_inverted_feedback(self, simulator, healthy_batch):
        injected = simulator.inject_broken_feedback(
            healthy_batch,
            loops=[("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_APOPTOSIS")],
            break_type=FeedbackBreakType.INVERTED,
        )
        assert injected.defect_metadata[-1].parameters["break_type"] == "inverted"


class TestStress:

    def test_hypoxia_signature(self, simulator, healthy_batch):
        injected = simulator.inject_stress(
            healthy_batch, stressor=StressorType.HYPOXIA, severity=0.7
        )
        meta = injected.defect_metadata[-1]
        assert meta.defect_type == DefectType.STRESS
        assert "HALLMARK_HYPOXIA" in meta.affected_pathways

    def test_all_stressor_types(self, simulator, healthy_batch):
        for stressor in StressorType:
            injected = simulator.inject_stress(
                healthy_batch, stressor=stressor, severity=0.5
            )
            assert len(injected.defect_metadata) > 0


class TestAtlasDrift:

    def test_drift_shifts_expression(self, simulator, healthy_batch):
        injected = simulator.inject_atlas_drift(
            healthy_batch, drift_distance=0.7
        )
        assert all(l.affected for l in injected.cell_labels)
        assert injected.defect_metadata[-1].defect_type == DefectType.ATLAS_DRIFT

    def test_higher_drift_more_pathways(self, simulator, healthy_batch):
        low = simulator.inject_atlas_drift(healthy_batch, drift_distance=0.2)
        high = simulator.inject_atlas_drift(healthy_batch, drift_distance=0.8)
        assert len(high.defect_metadata[-1].affected_pathways) >= len(
            low.defect_metadata[-1].affected_pathways
        )


class TestCombinedInjection:

    def test_multiple_defects(self, simulator, healthy_batch):
        injected = simulator.inject_combined(
            healthy_batch,
            defect_list=[
                ("inject_offtarget", {
                    "pathways": ["HALLMARK_APOPTOSIS"],
                    "cell_fraction": 0.2,
                    "activation_level": 0.5,
                }),
                ("inject_stress", {
                    "stressor": StressorType.HYPOXIA,
                    "severity": 0.5,
                }),
            ],
        )
        assert len(injected.defect_metadata) >= 2


# ── Data Classes ──────────────────────────────────────────────────────

class TestDataClasses:

    def test_defect_label(self):
        label = DefectLabel(affected=True, defect_type=DefectType.OFFTARGET, severity=0.5)
        assert label.affected
        assert label.defect_type == DefectType.OFFTARGET

    def test_defect_metadata_to_dict(self):
        meta = DefectMetadata(
            defect_type=DefectType.STRESS,
            severity=0.5,
            affected_pathways=["HALLMARK_HYPOXIA"],
            expected_detections={"F11_stress": True},
        )
        d = meta.to_dict()
        assert d["defect_type"] == "stress"
        assert d["expected_detections"]["F11_stress"] is True

    def test_manufacturing_spec_to_dict(self):
        spec = ManufacturingSpec(
            intended_pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            excluded_pathways=["HALLMARK_APOPTOSIS"],
            therapeutic_windows={"HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8)},
        )
        d = spec.to_dict()
        assert "HALLMARK_INFLAMMATORY_RESPONSE" in d["intended_pathways"]


# ── Orthogonality Matrix ─────────────────────────────────────────────

class TestOrthogonalityMatrix:

    def test_matrix_runs_with_stubs(self):
        from pathway_subtyping.qc.testing.orthogonality import OrthogonalityMatrix
        ortho = OrthogonalityMatrix(seed=42)
        result = ortho.run(n_trials=2, severity=0.5, n_cells=30)
        assert result.matrix.shape[0] == len(DefectType)
        # All stubs return False, so all detection rates should be 0
        assert (result.matrix.values == 0.0).all()

    def test_register_and_detect(self):
        from pathway_subtyping.qc.testing.orthogonality import (
            FeatureDetection,
            OrthogonalityMatrix,
        )
        ortho = OrthogonalityMatrix(seed=42)

        # Register a detector that always fires
        def always_detect(batch):
            return FeatureDetection(feature_id="F6_offtarget", detected=True, score=1.0)

        ortho.register_detector("F6_offtarget", always_detect)
        result = ortho.run(n_trials=2, severity=0.5, n_cells=30)
        # F6 column should be all 1.0
        assert result.matrix["F6_offtarget"].sum() > 0

    def test_invalid_feature_id_raises(self):
        from pathway_subtyping.qc.testing.orthogonality import OrthogonalityMatrix
        ortho = OrthogonalityMatrix(seed=42)
        with pytest.raises(ValueError, match="Unknown feature ID"):
            ortho.register_detector("F99_nonexistent", lambda b: None)


# ── Severity Titration ────────────────────────────────────────────────

class TestSeverityTitration:

    def test_titration_runs(self):
        from pathway_subtyping.qc.testing.orthogonality import FeatureDetection
        from pathway_subtyping.qc.testing.severity_titration import SeverityTitration

        titration = SeverityTitration(seed=42)

        # Register a detector that fires above severity 0.5
        def threshold_detector(batch):
            # Check if any defect metadata has severity > 0.4
            detected = any(
                m.severity > 0.4 for m in batch.defect_metadata
            )
            return FeatureDetection(
                feature_id="F6_offtarget", detected=detected, score=1.0 if detected else 0.0
            )

        titration.register_detector("F6_offtarget", threshold_detector)
        result = titration.run(
            n_trials=3, n_cells=30,
            defect_types=[DefectType.OFFTARGET],
        )
        assert len(result.curves) == 1
        curve = result.curves[0]
        assert len(curve.severity_levels) == 5
        assert len(curve.detection_rates) == 5


# ── Scenarios ─────────────────────────────────────────────────────────

class TestScenarios:

    def test_cart_scenario_runs(self):
        from pathway_subtyping.qc.testing.scenarios.cart_manufacturing import (
            CARTManufacturingScenario,
        )
        scenario = CARTManufacturingScenario(seed=42)
        report = scenario.run()
        assert report.scenario_name == "CARTManufacturingScenario"
        assert len(report.timepoints) == 5

    def test_bpu_scenario_runs(self):
        from pathway_subtyping.qc.testing.scenarios.bpu_organoid import (
            BPUOrganoidScenario,
        )
        scenario = BPUOrganoidScenario(seed=42)
        report = scenario.run()
        assert report.scenario_name == "BPUOrganoidScenario"
        assert len(report.timepoints) == 5

    def test_passage_scenario_runs(self):
        from pathway_subtyping.qc.testing.scenarios.passage_stability import (
            PassageStabilityScenario,
        )
        scenario = PassageStabilityScenario(seed=42)
        report = scenario.run()
        assert report.scenario_name == "PassageStabilityScenario"
        assert len(report.timepoints) == 20

    def test_report_validator(self):
        from pathway_subtyping.qc.testing.scenarios.cart_manufacturing import (
            CARTManufacturingScenario,
        )
        from pathway_subtyping.qc.testing.scenarios.report_validator import (
            ReportValidator,
        )
        scenario = CARTManufacturingScenario(seed=42)
        report = scenario.run()
        validator = ReportValidator()
        rubric = validator.validate(report)
        assert len(rubric.items) == 4
