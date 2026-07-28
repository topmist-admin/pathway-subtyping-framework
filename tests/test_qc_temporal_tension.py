"""
Tests for F2 (TemporalTracker) and F3 (TensionScorer).
"""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.qc.temporal import (
    PathwayTrajectory,
    TemporalResult,
    TemporalTracker,
    TrajectoryType,
)
from pathway_subtyping.qc.tension import TensionResult, TensionScorer

# ══════════════════════════════════════════════════════════════════════
# F2: TEMPORAL TRACKER
# ══════════════════════════════════════════════════════════════════════


class TestTrajectoryClassification:

    def test_resolving(self):
        t = PathwayTrajectory(
            pathway="PW",
            timepoint_ids=["t0", "t1", "t2", "t3"],
            scores=[0.1, 0.3, 0.6, 0.9],
            trajectory_type=TrajectoryType.RESOLVING,
        )
        assert t.trend > 0
        assert t.n_timepoints == 4

    def test_to_dict(self):
        t = PathwayTrajectory(
            pathway="PW",
            timepoint_ids=["t0", "t1"],
            scores=[0.5, 0.8],
            trajectory_type=TrajectoryType.RESOLVING,
        )
        d = t.to_dict()
        assert d["pathway"] == "PW"
        assert "trajectory_type" in d
        assert "trend" in d


class TestTemporalTracker:

    @pytest.fixture
    def rng(self):
        return np.random.default_rng(42)

    def _make_scores(self, rng, pathways, n_cells, base_values):
        data = {}
        for pw, base in zip(pathways, base_values):
            data[pw] = rng.normal(base, 0.1, n_cells)
        return pd.DataFrame(data, index=[f"cell_{i}" for i in range(n_cells)])

    def test_register_timepoints(self, rng):
        tracker = TemporalTracker()
        for i in range(3):
            df = self._make_scores(rng, ["A", "B"], 50, [0.5, 0.3])
            tracker.register_timepoint(df, f"t{i}")
        assert tracker.n_timepoints == 3

    def test_resolving_trajectory(self, rng):
        tracker = TemporalTracker()
        # Increasing scores across timepoints
        for i, base in enumerate([0.2, 0.4, 0.6, 0.8]):
            df = self._make_scores(rng, ["PW"], 50, [base])
            tracker.register_timepoint(df, f"t{i}")
        result = tracker.compute_all_trajectories(pathways=["PW"])
        assert result.trajectories[0].trajectory_type == TrajectoryType.RESOLVING

    def test_stalled_trajectory(self, rng):
        tracker = TemporalTracker(
            expected_profile={"PW": 0.8},
            stall_threshold=0.05,
            min_timepoints_for_stall=2,
        )
        # Plateau at 0.3, below target of 0.8
        for i, base in enumerate([0.2, 0.3, 0.3, 0.3]):
            df = self._make_scores(rng, ["PW"], 50, [base])
            tracker.register_timepoint(df, f"t{i}")
        result = tracker.compute_all_trajectories(pathways=["PW"])
        assert result.trajectories[0].trajectory_type == TrajectoryType.STALLED
        assert result.n_stalled == 1

    def test_reversing_trajectory(self, rng):
        tracker = TemporalTracker(stall_threshold=0.05)
        # Rise then fall
        for i, base in enumerate([0.2, 0.5, 0.8, 0.4]):
            df = self._make_scores(rng, ["PW"], 50, [base])
            tracker.register_timepoint(df, f"t{i}")
        result = tracker.compute_all_trajectories(pathways=["PW"])
        assert result.trajectories[0].trajectory_type == TrajectoryType.REVERSING
        assert result.n_reversing == 1

    def test_oscillating_trajectory(self, rng):
        tracker = TemporalTracker(stall_threshold=0.05)
        # Zigzag
        for i, base in enumerate([0.3, 0.7, 0.3, 0.7, 0.3]):
            df = self._make_scores(rng, ["PW"], 50, [base])
            tracker.register_timepoint(df, f"t{i}")
        result = tracker.compute_all_trajectories(pathways=["PW"])
        assert result.trajectories[0].trajectory_type == TrajectoryType.OSCILLATING

    def test_insufficient_data(self, rng):
        tracker = TemporalTracker()
        df = self._make_scores(rng, ["PW"], 50, [0.5])
        tracker.register_timepoint(df, "t0")
        result = tracker.compute_all_trajectories(pathways=["PW"])
        assert result.trajectories[0].trajectory_type == TrajectoryType.INSUFFICIENT_DATA

    def test_multiple_pathways(self, rng):
        tracker = TemporalTracker(expected_profile={"A": 0.9, "B": 0.9})
        for i, (a, b) in enumerate([(0.2, 0.2), (0.5, 0.3), (0.8, 0.3), (0.9, 0.3)]):
            df = self._make_scores(rng, ["A", "B"], 50, [a, b])
            tracker.register_timepoint(df, f"t{i}")
        result = tracker.compute_all_trajectories()
        types = {t.pathway: t.trajectory_type for t in result.trajectories}
        assert types["A"] == TrajectoryType.RESOLVING

    def test_get_problem_pathways(self, rng):
        tracker = TemporalTracker(expected_profile={"PW": 0.8}, stall_threshold=0.05)
        for base in [0.3, 0.3, 0.3, 0.3]:
            df = self._make_scores(rng, ["PW"], 50, [base])
            tracker.register_timepoint(df, f"t{len(tracker._timepoints)}")
        result = tracker.compute_all_trajectories(pathways=["PW"])
        problems = result.get_problem_pathways()
        assert isinstance(problems, list)

    def test_to_dict(self, rng):
        tracker = TemporalTracker()
        for base in [0.3, 0.5]:
            df = self._make_scores(rng, ["PW"], 30, [base])
            tracker.register_timepoint(df, f"t{len(tracker._timepoints)}")
        result = tracker.compute_all_trajectories()
        d = result.to_dict()
        assert "passed" in d
        assert "trajectories" in d

    def test_no_timepoints(self):
        tracker = TemporalTracker()
        result = tracker.compute_all_trajectories()
        assert len(result.trajectories) == 0


# ══════════════════════════════════════════════════════════════════════
# F3: TENSION SCORER
# ══════════════════════════════════════════════════════════════════════


class TestTensionScorer:

    @pytest.fixture
    def complete_ccr(self):
        """All cascades complete (CCR = 1.0)."""
        return pd.DataFrame(
            {"PW_A": np.ones(50), "PW_B": np.ones(50)},
            index=[f"cell_{i}" for i in range(50)],
        )

    @pytest.fixture
    def stalled_ccr(self):
        """Some cascades incomplete."""
        rng = np.random.default_rng(42)
        return pd.DataFrame(
            {
                "PW_A": np.clip(rng.normal(0.3, 0.1, 50), 0, 1),  # Mostly stalled
                "PW_B": np.clip(rng.normal(0.9, 0.05, 50), 0, 1),  # Mostly complete
            },
            index=[f"cell_{i}" for i in range(50)],
        )

    def test_zero_tension_when_complete(self, complete_ccr):
        scorer = TensionScorer()
        tension = scorer.score_sample(complete_ccr)
        assert np.allclose(tension, 0.0, atol=1e-10)

    def test_high_tension_when_stalled(self, stalled_ccr):
        scorer = TensionScorer()
        tension = scorer.score_sample(stalled_ccr)
        assert np.mean(tension) > 0.3

    def test_pathway_weights(self, stalled_ccr):
        scorer_equal = TensionScorer(pathway_weights={"PW_A": 1.0, "PW_B": 1.0})
        scorer_weighted = TensionScorer(pathway_weights={"PW_A": 5.0, "PW_B": 1.0})
        t_equal = scorer_equal.score_sample(stalled_ccr)
        t_weighted = scorer_weighted.score_sample(stalled_ccr)
        # PW_A is stalled — weighting it 5x should increase tension
        assert np.mean(t_weighted) > np.mean(t_equal)

    def test_score_batch(self):
        from pathway_subtyping.qc.cascade import CascadeResult, LayerScores

        ccr = pd.DataFrame(
            {"PW": np.full(30, 0.4)},
            index=[f"cell_{i}" for i in range(30)],
        )
        cascade_result = CascadeResult(
            per_pathway=[LayerScores(pathway="PW", upstream_score=0.8, downstream_score=0.3)],
            per_cell_ccr=ccr,
            n_incomplete=1,
            mean_completion_ratio=0.4,
        )
        scorer = TensionScorer(tension_threshold=0.5)
        result = scorer.score_batch(cascade_result)
        assert isinstance(result, TensionResult)
        assert result.mean_tension > 0
        assert "PW" in result.per_pathway_contribution

    def test_tension_gate_pass(self):
        scorer = TensionScorer(tension_threshold=1.0)
        scores = np.full(50, 0.3)
        assert scorer.tension_gate(scores) is True

    def test_tension_gate_fail(self):
        scorer = TensionScorer(tension_threshold=0.1)
        scores = np.full(50, 0.5)
        assert scorer.tension_gate(scores) is False

    def test_empty_cascade_result(self):
        from pathway_subtyping.qc.cascade import CascadeResult

        empty = CascadeResult(
            per_pathway=[], per_cell_ccr=None, n_incomplete=0, mean_completion_ratio=1.0
        )
        scorer = TensionScorer()
        result = scorer.score_batch(empty)
        assert len(result.per_cell_tension) == 0

    def test_to_dict(self, stalled_ccr):
        from pathway_subtyping.qc.cascade import CascadeResult, LayerScores

        cascade_result = CascadeResult(
            per_pathway=[LayerScores(pathway="PW_A"), LayerScores(pathway="PW_B")],
            per_cell_ccr=stalled_ccr,
        )
        scorer = TensionScorer()
        result = scorer.score_batch(cascade_result)
        d = result.to_dict()
        assert "mean_tension" in d
        assert "per_pathway_contribution" in d
