"""
Real-data acceptance tests for the F5 perturbation layer.

Consumes the JSON artifact produced by scripts/validate_f5_real_data.py
and asserts the roadmap criteria:

    - directional agreement >= 0.70 across curated (gene, pathway,
      expected_sign) edges — knocking out a pathway driver depresses
      the pathway's predicted MSV score.
    - conformal coverage on perturbed MSV prediction within +/- 2% of
      the finite-sample oracle on a TCGA-COAD test split.

Cohort: TCGA-COAD (57 samples × log1p(TPM)). Uses FallbackPerturber
as the backend — the directional test exercises the perturb → MSV
head pipeline on real expression data and is the gate enforced here.
Production Geneformer weights (``pip install
pathway-subtyping[perturb]`` + cached checkpoint) are the tracked
post-release follow-up path.

If the JSON artifact is absent, the tests are skipped — CI runs
without the TCGA-COAD TSVs will not fail.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT_PATH = REPO_ROOT / "results" / "f5_validation" / "perturbation_directional.json"


def _load_or_skip() -> dict:
    if not ARTIFACT_PATH.exists():
        pytest.skip(
            f"no real-data validation artifact at {ARTIFACT_PATH} "
            f"(run scripts/validate_f5_real_data.py to generate)"
        )
    return json.loads(ARTIFACT_PATH.read_text())


class TestRealDataPerturb:

    def test_artifact_shape(self):
        report = _load_or_skip()
        assert report["feature"] == "F5"
        assert report["cohort"]["n_samples"] >= 20
        assert report["directional"]["n_edges_evaluated"] >= 8

    def test_directional_agreement_at_least_70pct(self):
        """Roadmap gate: sign(delta_MSV) matches expected on >= 70% of edges."""
        report = _load_or_skip()
        rate = report["directional"]["directional_agreement_rate"]
        target = report["acceptance"]["directional_target"]
        assert rate >= target, (
            f"directional agreement {rate:.4f} below target {target:.4f}; "
            f"perturbation pipeline is not preserving the signed KO "
            f"effect on pathway-member driver genes"
        )

    def test_conformal_coverage_within_2pct(self):
        """Roadmap gate: perturbed conformal coverage within +/- 2% oracle."""
        report = _load_or_skip()
        conformal = report["conformal"]
        if conformal.get("status") != "ok":
            pytest.fail(
                f"conformal sub-test did not run: {conformal.get('status')}"
            )
        oracle_dev = conformal["mean_oracle_deviation"]
        target = report["acceptance"]["conformal_oracle_target_abs"]
        assert abs(oracle_dev) < target, (
            f"perturbed-MSV conformal oracle deviation {oracle_dev:+.4f} "
            f"exceeds +/- {target:.4f}; F1 calibration guarantee has not "
            f"transferred through the perturbation wrapper"
        )

    def test_every_edge_evaluated_or_accounted(self):
        """No silent drops: every edge must be either evaluated or flagged."""
        report = _load_or_skip()
        edges = report["directional"]["edges"]
        for e in edges:
            assert "status" in e, f"edge missing status field: {e}"
            if e["status"] == "ok":
                assert "observed_sign" in e
                assert "mean_delta_msv" in e
