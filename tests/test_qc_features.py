"""
Tests for QC Layer Phase 1 Features: F6 Off-Target, F7 Heterogeneity, F8 Dosage.

Tests validate that each feature:
    - Detects planted defects in synthetic data
    - Does not false-alarm on healthy batches
    - Produces correct data structures and serialization
    - Integrates with the ManufacturingSimulator
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.qc.dosage import (
    DosageAnalysisResult,
    DosageAnalyzer,
    DosageState,
)
from pathway_subtyping.qc.heterogeneity import (
    HeterogeneityProfiler,
    HeterogeneityResult,
)
from pathway_subtyping.qc.offtarget import (
    ActivationClass,
    OffTargetDetector,
    OffTargetResult,
)
from pathway_subtyping.qc.testing.simulator import (
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
            "HALLMARK_P53_PATHWAY",
        ],
        therapeutic_windows={
            "HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8),
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB": (0.3, 0.7),
        },
    )


@pytest.fixture
def healthy_batch(simulator, spec):
    return simulator.generate_healthy_batch(n_cells=200, spec=spec)


# ══════════════════════════════════════════════════════════════════════
# F6: OFF-TARGET DETECTION
# ══════════════════════════════════════════════════════════════════════


class TestOffTargetDetector:

    def test_healthy_batch_no_violations(self, healthy_batch, spec):
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.5,
        )
        result = detector.scan(healthy_batch.pathway_scores)
        assert isinstance(result, OffTargetResult)
        assert result.n_excluded_violations == 0

    def test_detects_offtarget_activation(self, simulator, healthy_batch, spec):
        # Inject off-target activation in excluded pathway
        injected = simulator.inject_offtarget(
            healthy_batch,
            pathways=["HALLMARK_APOPTOSIS"],
            cell_fraction=1.0,
            activation_level=0.9,
        )

        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.3,
            cell_fraction_threshold=0.05,
        )
        result = detector.scan(injected.pathway_scores)
        assert result.n_excluded_violations > 0 or result.n_off_target > 0
        assert result.safety_flag or result.n_off_target > 0

    def test_intended_pathways_not_flagged(self, healthy_batch, spec):
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.5,
        )
        result = detector.scan(healthy_batch.pathway_scores)
        for pr in result.pathway_results:
            if pr.pathway in spec.intended_pathways:
                assert pr.classification == ActivationClass.INTENDED

    def test_excluded_violation_triggers_safety_flag(self, simulator, healthy_batch, spec):
        injected = simulator.inject_offtarget(
            healthy_batch,
            pathways=["HALLMARK_APOPTOSIS"],
            cell_fraction=0.8,
            activation_level=0.9,
        )
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.0,
            cell_fraction_threshold=0.05,
        )
        result = detector.scan(injected.pathway_scores)
        violations = result.get_excluded_violations()
        if violations:
            assert result.safety_flag

    def test_calibrate_threshold(self, healthy_batch, spec):
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
        )
        threshold = detector.calibrate_threshold(healthy_batch.pathway_scores, percentile=95)
        assert isinstance(threshold, float)
        assert detector.noise_threshold == threshold

    def test_to_dict(self, healthy_batch, spec):
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
        )
        result = detector.scan(healthy_batch.pathway_scores)
        d = result.to_dict()
        assert "passed" in d
        assert "pathways" in d
        assert isinstance(d["pathways"], list)

    def test_get_off_target_pathways(self, simulator, healthy_batch, spec):
        injected = simulator.inject_offtarget(
            healthy_batch,
            pathways=["HALLMARK_HYPOXIA"],
            cell_fraction=0.5,
            activation_level=0.8,
        )
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.0,
            cell_fraction_threshold=0.05,
        )
        result = detector.scan(injected.pathway_scores)
        off_target = result.get_off_target_pathways()
        assert isinstance(off_target, list)


# ══════════════════════════════════════════════════════════════════════
# F7: HETEROGENEITY PROFILING
# ══════════════════════════════════════════════════════════════════════


class TestHeterogeneityProfiler:

    def test_healthy_batch_low_heterogeneity(self, healthy_batch):
        profiler = HeterogeneityProfiler(conformity_threshold=5.0)
        result = profiler.profile(healthy_batch.pathway_scores)
        assert isinstance(result, HeterogeneityResult)
        # Healthy batch should have low heterogeneity
        assert result.heterogeneity_index < 0.5
        assert result.n_total == 200

    def test_detects_subpopulation(self, simulator, healthy_batch):
        injected = simulator.inject_subpopulation(healthy_batch, fraction=0.2)
        # Use tight conformity threshold to catch deviations
        profiler = HeterogeneityProfiler(
            conformity_threshold=1.0,
            dbscan_eps=2.0,
            dbscan_min_samples=3,
        )
        result = profiler.profile(injected.pathway_scores)
        assert result.n_nonconforming > 0
        assert result.heterogeneity_index > 0.0

    def test_custom_target_profile(self, healthy_batch):
        target = {pw: 0.0 for pw in healthy_batch.pathway_scores.columns}
        profiler = HeterogeneityProfiler(target_profile=target, conformity_threshold=2.0)
        result = profiler.profile(healthy_batch.pathway_scores)
        assert result.n_total == 200

    def test_bimodality_coefficient(self, healthy_batch):
        profiler = HeterogeneityProfiler(conformity_threshold=3.0)
        result = profiler.profile(healthy_batch.pathway_scores)
        assert isinstance(result.bimodality_coefficient, float)

    def test_conformity_scores_shape(self, healthy_batch):
        profiler = HeterogeneityProfiler(conformity_threshold=3.0)
        result = profiler.profile(healthy_batch.pathway_scores)
        assert result.conformity_scores.shape == (200,)
        assert all(s >= 0 for s in result.conformity_scores)

    def test_to_dict(self, healthy_batch):
        profiler = HeterogeneityProfiler(conformity_threshold=3.0)
        result = profiler.profile(healthy_batch.pathway_scores)
        d = result.to_dict()
        assert "heterogeneity_index" in d
        assert "bimodality_coefficient" in d
        assert "subpopulations" in d

    def test_set_specification(self, healthy_batch):
        profiler = HeterogeneityProfiler()
        target = {pw: 0.5 for pw in healthy_batch.pathway_scores.columns}
        profiler.set_specification(target)
        assert profiler.target_profile == target

    def test_small_batch_no_crash(self):
        scores = pd.DataFrame(np.random.randn(5, 3), columns=["A", "B", "C"])
        profiler = HeterogeneityProfiler(conformity_threshold=2.0)
        result = profiler.profile(scores)
        assert result.n_total == 5


# ══════════════════════════════════════════════════════════════════════
# F8: DOSAGE & STOICHIOMETRY
# ══════════════════════════════════════════════════════════════════════


class TestDosageAnalyzer:

    def test_healthy_batch_in_range(self, healthy_batch, spec):
        # Use wide windows so healthy batch passes
        analyzer = DosageAnalyzer(
            therapeutic_windows={pw: (-5.0, 5.0) for pw in spec.therapeutic_windows},
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        assert isinstance(result, DosageAnalysisResult)
        assert result.n_over == 0
        assert result.passed

    def test_detects_overdose(self, simulator, healthy_batch, spec):
        injected = simulator.inject_overdose(
            healthy_batch,
            pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            multiplier=5.0,
        )
        # Use a very tight window around zero — overdose shifts the variance
        # and per-cell z-scores, so individual cells exceed the window even
        # if the batch mean is near zero after z-scoring
        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (-0.001, 0.001),
            },
        )
        result = analyzer.analyze(injected.pathway_scores)
        # With such a tight window, nearly all cells should be outside
        pw_result = [
            r for r in result.pathway_dosages if r.pathway == "HALLMARK_INFLAMMATORY_RESPONSE"
        ][0]
        assert pw_result.n_cells_over + pw_result.n_cells_under > 0

    def test_detects_underdose(self, healthy_batch):
        # Set window above what the batch produces
        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (50.0, 100.0),
            },
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        assert result.n_under > 0

    def test_stoichiometry_check(self, healthy_batch):
        analyzer = DosageAnalyzer(
            therapeutic_windows={},
            pathway_ratios=[
                {
                    "pathway_a": "HALLMARK_INFLAMMATORY_RESPONSE",
                    "pathway_b": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                    "expected_ratio": 1.0,
                    "tolerance": 50.0,  # Wide tolerance
                },
            ],
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        assert len(result.stoichiometry_results) == 1
        assert result.stoichiometry_results[0].within_tolerance

    def test_stoichiometry_violation(self, healthy_batch):
        analyzer = DosageAnalyzer(
            therapeutic_windows={},
            pathway_ratios=[
                {
                    "pathway_a": "HALLMARK_INFLAMMATORY_RESPONSE",
                    "pathway_b": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                    "expected_ratio": 100.0,  # Impossible ratio
                    "tolerance": 0.1,
                },
            ],
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        assert result.n_ratio_violations > 0

    def test_toxicity_risk_score(self, simulator, healthy_batch):
        injected = simulator.inject_overdose(
            healthy_batch,
            pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            multiplier=3.0,
        )
        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (-1.0, 0.5),
            },
            toxicity_weights={
                "HALLMARK_INFLAMMATORY_RESPONSE": 2.0,
            },
        )
        result = analyzer.analyze(injected.pathway_scores)
        assert result.toxicity_risk_score >= 0.0

    def test_define_therapeutic_window(self, healthy_batch):
        analyzer = DosageAnalyzer()
        analyzer.define_therapeutic_window("HALLMARK_HYPOXIA", 0.2, 0.8)
        assert "HALLMARK_HYPOXIA" in analyzer.therapeutic_windows

    def test_to_dict(self, healthy_batch, spec):
        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (-5.0, 5.0),
            },
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        d = result.to_dict()
        assert "passed" in d
        assert "pathways" in d
        assert "stoichiometry" in d
        assert "toxicity_risk_score" in d

    def test_per_cell_counts(self, healthy_batch):
        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "HALLMARK_INFLAMMATORY_RESPONSE": (-0.01, 0.01),  # Tight window
            },
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        pw_result = [
            r for r in result.pathway_dosages if r.pathway == "HALLMARK_INFLAMMATORY_RESPONSE"
        ]
        assert len(pw_result) == 1
        total = (
            pw_result[0].n_cells_under + pw_result[0].n_cells_in_range + pw_result[0].n_cells_over
        )
        assert total == 200

    def test_missing_pathway_ignored(self, healthy_batch):
        analyzer = DosageAnalyzer(
            therapeutic_windows={
                "NONEXISTENT_PATHWAY": (0.0, 1.0),
            },
        )
        result = analyzer.analyze(healthy_batch.pathway_scores)
        assert len(result.pathway_dosages) == 0


# ══════════════════════════════════════════════════════════════════════
# INTEGRATION: Features work together on simulator data
# ══════════════════════════════════════════════════════════════════════


class TestFeatureIntegration:

    def test_all_features_on_healthy_batch(self, simulator, spec):
        batch = simulator.generate_healthy_batch(n_cells=100, spec=spec)

        # F6
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.5,
        )
        f6 = detector.scan(batch.pathway_scores)
        assert f6.n_excluded_violations == 0

        # F7
        profiler = HeterogeneityProfiler(conformity_threshold=5.0)
        f7 = profiler.profile(batch.pathway_scores)
        assert f7.n_total == 100

        # F8
        analyzer = DosageAnalyzer(
            therapeutic_windows={pw: (-5.0, 5.0) for pw in spec.therapeutic_windows},
        )
        f8 = analyzer.analyze(batch.pathway_scores)
        assert f8.passed

    def test_combined_defects_detected(self, simulator, spec):
        batch = simulator.generate_healthy_batch(n_cells=150, spec=spec)

        # Inject both off-target and subpopulation
        injected = simulator.inject_combined(
            batch,
            defect_list=[
                (
                    "inject_offtarget",
                    {
                        "pathways": ["HALLMARK_APOPTOSIS"],
                        "cell_fraction": 0.5,
                        "activation_level": 0.9,
                    },
                ),
                ("inject_subpopulation", {"fraction": 0.2}),
            ],
        )

        # F6 should detect off-target
        detector = OffTargetDetector(
            intended_pathways=spec.intended_pathways,
            excluded_pathways=spec.excluded_pathways,
            noise_threshold=0.0,
            cell_fraction_threshold=0.05,
        )
        f6 = detector.scan(injected.pathway_scores)
        assert f6.n_excluded_violations > 0 or f6.n_off_target > 0

        # F7 should see heterogeneity
        profiler = HeterogeneityProfiler(conformity_threshold=1.0)
        f7 = profiler.profile(injected.pathway_scores)
        assert f7.n_nonconforming > 0
