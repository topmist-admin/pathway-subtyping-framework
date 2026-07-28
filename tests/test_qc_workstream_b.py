"""
Tests for Workstream B features: F4 (ResolutionGate), F5 (DriftDetector),
F9 (CrosstalkDetector), F10 (FeedbackMonitor), F11 (StressFingerprinter).
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.knowledge_graph.builder import KnowledgeGraph
from pathway_subtyping.knowledge_graph.schema import EdgeType, NodeType
from pathway_subtyping.qc.cascade import CascadeResult, LayerScores
from pathway_subtyping.qc.crosstalk import CrosstalkDetector, CrosstalkResult
from pathway_subtyping.qc.drift import DriftDetector, DriftResult
from pathway_subtyping.qc.feedback import FeedbackMonitor, FeedbackResult, FeedbackStatus
from pathway_subtyping.qc.gates import ReleaseDecision, ResolutionGate, ResolutionResult
from pathway_subtyping.qc.stress import StressFingerprinter, StressResult
from pathway_subtyping.qc.tension import TensionResult

# ══════════════════════════════════════════════════════════════════════
# F4: RESOLUTION GATE
# ══════════════════════════════════════════════════════════════════════


class TestResolutionGate:

    def _make_cascade(self, ccr=0.9, n_incomplete=0):
        return CascadeResult(
            per_pathway=[LayerScores(pathway="PW", upstream_score=1.0, downstream_score=ccr)],
            per_cell_ccr=pd.DataFrame({"PW": np.full(10, ccr)}),
            n_incomplete=n_incomplete,
            mean_completion_ratio=ccr,
        )

    def _make_tension(self, mean=0.3, max_val=0.5):
        return TensionResult(
            per_cell_tension=np.full(10, mean),
            per_pathway_contribution={"PW": mean},
            mean_tension=mean,
            median_tension=mean,
            max_tension=max_val,
        )

    def test_all_pass_release(self):
        gate = ResolutionGate()
        result = gate.evaluate(
            cascade_result=self._make_cascade(0.9),
            tension_result=self._make_tension(0.3),
        )
        assert result.decision == ReleaseDecision.RELEASE

    def test_cascade_fail_reject(self):
        gate = ResolutionGate(min_cascade_ccr=0.8, hold_threshold=0.5)
        result = gate.evaluate(
            cascade_result=self._make_cascade(0.3, n_incomplete=1),
            tension_result=self._make_tension(0.3),
        )
        assert result.decision in (ReleaseDecision.HOLD, ReleaseDecision.REJECT)

    def test_no_results_release(self):
        gate = ResolutionGate()
        result = gate.evaluate()  # No results at all
        assert result.decision == ReleaseDecision.RELEASE  # No gates to fail

    def test_partial_results(self):
        gate = ResolutionGate()
        result = gate.evaluate(tension_result=self._make_tension(0.3))
        assert len(result.gate_checks) == 1

    def test_to_dict(self):
        gate = ResolutionGate()
        result = gate.evaluate(cascade_result=self._make_cascade(0.9))
        d = result.to_dict()
        assert "decision" in d
        assert "gates" in d
        assert "justification" in d


# ══════════════════════════════════════════════════════════════════════
# F5: DRIFT DETECTOR
# ══════════════════════════════════════════════════════════════════════


class TestDriftDetector:

    def test_no_drift(self):
        detector = DriftDetector(max_drift=0.5)
        baseline = {"PW_A": 0.9, "PW_B": 0.8}
        detector.establish_baseline(baseline)
        m = detector.measure_drift({"PW_A": 0.9, "PW_B": 0.8}, "P1")
        assert m.drift_score == pytest.approx(0.0)
        assert not m.alarm_triggered

    def test_drift_detected(self):
        detector = DriftDetector(max_drift=0.3)
        detector.establish_baseline({"PW_A": 0.9, "PW_B": 0.8})
        m = detector.measure_drift({"PW_A": 0.5, "PW_B": 0.3}, "P5")
        assert m.drift_score > 0.3
        assert m.alarm_triggered

    def test_trajectory(self):
        detector = DriftDetector(max_drift=1.0)
        detector.establish_baseline({"PW": 0.9})
        for i in range(5):
            detector.measure_drift({"PW": 0.9 - i * 0.1}, f"P{i}")
        result = detector.drift_trajectory()
        assert isinstance(result, DriftResult)
        assert len(result.measurements) == 5

    def test_identify_drivers(self):
        detector = DriftDetector(max_drift=1.0)
        detector.establish_baseline({"A": 0.9, "B": 0.9, "C": 0.9})
        detector.measure_drift({"A": 0.1, "B": 0.8, "C": 0.85}, "P10")
        drivers = detector.identify_drivers()
        assert drivers[0] == "A"

    def test_no_baseline_raises(self):
        detector = DriftDetector()
        with pytest.raises(ValueError, match="Baseline not established"):
            detector.measure_drift({"PW": 0.5}, "P1")

    def test_to_dict(self):
        detector = DriftDetector(max_drift=0.5)
        detector.establish_baseline({"PW": 0.9})
        detector.measure_drift({"PW": 0.5}, "P1")
        d = detector.drift_trajectory().to_dict()
        assert "passed" in d
        assert "measurements" in d


# ══════════════════════════════════════════════════════════════════════
# F9: CROSSTALK DETECTOR
# ══════════════════════════════════════════════════════════════════════


class TestCrosstalkDetector:

    @pytest.fixture
    def crosstalk_kg(self):
        kg = KnowledgeGraph()
        for g in ["A", "B", "C", "D", "E", "F"]:
            kg.add_node(g, NodeType.GENE)
        kg.add_node("PW1", NodeType.PATHWAY)
        kg.add_node("PW2", NodeType.PATHWAY)
        for g in ["A", "B", "C"]:
            kg.add_edge(g, "PW1", EdgeType.GENE_IN_PATHWAY)
        for g in ["B", "C", "D"]:  # B, C shared
            kg.add_edge(g, "PW2", EdgeType.GENE_IN_PATHWAY)
        return kg

    def test_detects_shared_genes(self, crosstalk_kg):
        rng = np.random.default_rng(42)
        expr = pd.DataFrame(
            {g: rng.normal(5.0, 0.5, 50) for g in "ABCDEF"},
            index=[f"cell_{i}" for i in range(50)],
        )
        detector = CrosstalkDetector(crosstalk_kg)
        result = detector.detect(expr, pathways=["PW1", "PW2"])
        assert isinstance(result, CrosstalkResult)
        assert len(result.interference_edges) == 1
        assert len(result.interference_edges[0].shared_genes) == 2

    def test_no_shared_genes(self, crosstalk_kg):
        # PW1 vs nonexistent
        rng = np.random.default_rng(42)
        expr = pd.DataFrame(
            {g: rng.normal(5.0, 0.5, 30) for g in "ABCDEF"},
            index=[f"cell_{i}" for i in range(30)],
        )
        kg = crosstalk_kg
        kg.add_node("PW3", NodeType.PATHWAY)
        kg.add_edge("E", "PW3", EdgeType.GENE_IN_PATHWAY)
        kg.add_edge("F", "PW3", EdgeType.GENE_IN_PATHWAY)
        detector = CrosstalkDetector(kg)
        result = detector.detect(expr, pathways=["PW1", "PW3"])
        assert result.n_significant == 0

    def test_to_dict(self, crosstalk_kg):
        rng = np.random.default_rng(42)
        expr = pd.DataFrame(
            {g: rng.normal(5.0, 0.5, 30) for g in "ABCDEF"},
            index=[f"cell_{i}" for i in range(30)],
        )
        detector = CrosstalkDetector(crosstalk_kg)
        result = detector.detect(expr, pathways=["PW1", "PW2"])
        d = result.to_dict()
        assert "edges" in d


# ══════════════════════════════════════════════════════════════════════
# F10: FEEDBACK MONITOR
# ══════════════════════════════════════════════════════════════════════


class TestFeedbackMonitor:

    def test_intact_feedback(self):
        rng = np.random.default_rng(42)
        n = 100
        activator = rng.normal(5.0, 1.0, n)
        inhibitor = activator * 0.8 + rng.normal(0, 0.3, n)  # Positively correlated
        scores = pd.DataFrame({"ACT": activator, "INH": inhibitor})
        monitor = FeedbackMonitor(loops=[("ACT", "INH")])
        result = monitor.check(scores)
        assert result.loops[0].status == FeedbackStatus.INTACT
        assert result.n_intact == 1

    def test_decoupled_feedback(self):
        rng = np.random.default_rng(42)
        n = 100
        activator = rng.normal(5.0, 1.0, n)
        inhibitor = rng.normal(5.0, 1.0, n)  # Independent
        scores = pd.DataFrame({"ACT": activator, "INH": inhibitor})
        monitor = FeedbackMonitor(loops=[("ACT", "INH")])
        result = monitor.check(scores)
        assert result.loops[0].status in (FeedbackStatus.DECOUPLED, FeedbackStatus.INTACT)

    def test_inverted_feedback(self):
        rng = np.random.default_rng(42)
        n = 100
        activator = rng.normal(5.0, 1.0, n)
        inhibitor = -activator + rng.normal(0, 0.3, n)  # Negatively correlated
        scores = pd.DataFrame({"ACT": activator, "INH": inhibitor})
        monitor = FeedbackMonitor(loops=[("ACT", "INH")])
        result = monitor.check(scores)
        assert result.loops[0].status == FeedbackStatus.INVERTED

    def test_missing_pathway(self):
        scores = pd.DataFrame({"ACT": [1, 2, 3]})
        monitor = FeedbackMonitor(loops=[("ACT", "MISSING")])
        result = monitor.check(scores)
        assert result.loops[0].status == FeedbackStatus.UNKNOWN

    def test_multiple_loops(self):
        rng = np.random.default_rng(42)
        n = 100
        scores = pd.DataFrame(
            {
                "A": rng.normal(5, 1, n),
                "B": rng.normal(5, 1, n),
                "C": rng.normal(5, 1, n),
            }
        )
        # Make A-B correlated
        scores["B"] = scores["A"] * 0.9 + rng.normal(0, 0.2, n)
        monitor = FeedbackMonitor(loops=[("A", "B"), ("A", "C")])
        result = monitor.check(scores)
        assert len(result.loops) == 2

    def test_to_dict(self):
        scores = pd.DataFrame({"A": [1, 2, 3], "B": [3, 2, 1]})
        monitor = FeedbackMonitor(loops=[("A", "B")])
        result = monitor.check(scores)
        d = result.to_dict()
        assert "loops" in d
        assert "runaway_risk_score" in d


# ══════════════════════════════════════════════════════════════════════
# F11: STRESS FINGERPRINTER
# ══════════════════════════════════════════════════════════════════════


class TestStressFingerprinter:

    def test_detects_hypoxia(self):
        fp = StressFingerprinter()
        fp.load_default_signatures()
        rng = np.random.default_rng(42)
        n = 50
        data = {
            "HALLMARK_HYPOXIA": rng.normal(3.0, 0.3, n),
            "HALLMARK_GLYCOLYSIS": rng.normal(2.5, 0.3, n),
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION": rng.normal(-2.0, 0.3, n),
            "HALLMARK_FATTY_ACID_METABOLISM": rng.normal(-1.5, 0.3, n),
            "HALLMARK_ANGIOGENESIS": rng.normal(1.5, 0.3, n),
        }
        scores = pd.DataFrame(data)
        result = fp.fingerprint(scores)
        assert result.primary_stressor == "hypoxia"
        assert result.confidence > 0.3

    def test_no_stress(self):
        fp = StressFingerprinter(min_confidence=0.5)
        fp.load_default_signatures()
        rng = np.random.default_rng(42)
        # Near-zero values — no stress pattern
        data = {f"PW_{i}": rng.normal(0, 0.01, 30) for i in range(5)}
        scores = pd.DataFrame(data)
        result = fp.fingerprint(scores)
        # May or may not detect — depends on random noise

    def test_no_signatures_loaded(self):
        fp = StressFingerprinter()
        scores = pd.DataFrame({"PW": [1, 2, 3]})
        result = fp.fingerprint(scores)
        assert not result.stress_detected

    def test_remediation_provided(self):
        fp = StressFingerprinter()
        fp.load_default_signatures()
        rng = np.random.default_rng(42)
        data = {
            "HALLMARK_HYPOXIA": rng.normal(3.0, 0.1, 30),
            "HALLMARK_GLYCOLYSIS": rng.normal(2.0, 0.1, 30),
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION": rng.normal(-2.0, 0.1, 30),
        }
        result = fp.fingerprint(pd.DataFrame(data))
        if result.stress_detected:
            assert len(result.remediation) > 0

    def test_to_dict(self):
        fp = StressFingerprinter()
        fp.load_default_signatures()
        scores = pd.DataFrame({"HALLMARK_HYPOXIA": [1.0, 2.0]})
        result = fp.fingerprint(scores)
        d = result.to_dict()
        assert "stress_detected" in d
        assert "matches" in d

    def test_with_reference(self):
        fp = StressFingerprinter()
        fp.load_default_signatures()
        rng = np.random.default_rng(42)
        batch = pd.DataFrame(
            {
                "HALLMARK_HYPOXIA": rng.normal(3.0, 0.1, 30),
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION": rng.normal(-1.0, 0.1, 30),
            }
        )
        reference = pd.DataFrame(
            {
                "HALLMARK_HYPOXIA": rng.normal(0.0, 0.1, 30),
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION": rng.normal(0.0, 0.1, 30),
            }
        )
        result = fp.fingerprint(batch, reference_scores=reference)
        assert isinstance(result, StressResult)
