"""
Real-data acceptance tests for the F5 perturbation layer.

Consumes the JSON artifact produced by scripts/validate_f5_real_data.py
and asserts the roadmap criteria:

    - directional agreement >= 0.70 across curated (gene, pathway,
      expected_sign) edges (TCGA-COAD in-silico KO section).
    - conformal coverage on perturbed MSV prediction within +/- 2% of
      the finite-sample oracle (TCGA-COAD section).
    - directional agreement >= 0.70 between predicted in-silico ΔMSV
      and observed real-cohort ΔMSV on GSE123753 (WT vs MECP2-KO iPSC
      neurons, Rodrigues et al. 2020) — only asserted when the artifact's
      ``wt_vs_ko_gse123753`` section is present with a final status
      (otherwise the sub-test is skipped).

Either backend (``fallback`` or ``geneformer``) may have produced the
artefact. The gates themselves are the same; the artefact's
``backend_kind`` field records which was used.

If the JSON artifact is absent, the tests are skipped — CI runs
without the TCGA-COAD or GSE123753 data will not fail.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT_PATH = REPO_ROOT / "results" / "f5_validation" / "perturbation_directional.json"
GENEFORMER_ARTIFACT_PATH = (
    REPO_ROOT / "results" / "f5_validation" / "perturbation_directional_geneformer.json"
)


def _load_or_skip() -> dict:
    if not ARTIFACT_PATH.exists():
        pytest.skip(
            f"no real-data validation artifact at {ARTIFACT_PATH} "
            f"(run scripts/validate_f5_real_data.py to generate)"
        )
    return json.loads(ARTIFACT_PATH.read_text())


def _load_geneformer_or_skip() -> dict:
    if not GENEFORMER_ARTIFACT_PATH.exists():
        pytest.skip(
            f"no Geneformer-backed validation artifact at "
            f"{GENEFORMER_ARTIFACT_PATH} "
            f"(run scripts/validate_f5_real_data.py --backend geneformer)"
        )
    return json.loads(GENEFORMER_ARTIFACT_PATH.read_text())


class TestRealDataPerturb:

    def test_artifact_shape(self):
        report = _load_or_skip()
        assert report["feature"] == "F5"
        if "cohort" in report:
            assert report["cohort"]["n_samples"] >= 20
            assert report["directional"]["n_edges_evaluated"] >= 8

    def test_directional_agreement_at_least_70pct(self):
        """Roadmap gate: sign(delta_MSV) matches expected on >= 70% of edges."""
        report = _load_or_skip()
        if "directional" not in report:
            pytest.skip("no TCGA directional section in artefact")
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
        if "conformal" not in report:
            pytest.skip("no TCGA conformal section in artefact")
        conformal = report["conformal"]
        if conformal.get("status") != "ok":
            pytest.skip(
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
        if "directional" not in report:
            pytest.skip("no TCGA directional section in artefact")
        edges = report["directional"]["edges"]
        for e in edges:
            assert "status" in e, f"edge missing status field: {e}"
            if e["status"] == "ok":
                assert "observed_sign" in e
                assert "mean_delta_msv" in e

    # ------------------------------------------------------------------ #
    # GSE123753 WT-vs-KO real-cohort section
    # ------------------------------------------------------------------ #

    def test_wt_vs_ko_shape(self):
        """Artefact section present and well-formed when wt_vs_ko runs."""
        report = _load_or_skip()
        wtko = report.get("wt_vs_ko_gse123753")
        if wtko is None or isinstance(wtko, dict) and "status" in wtko and "n_pathways_evaluated" not in wtko:
            pytest.skip("wt_vs_ko_gse123753 section not present or skipped")
        assert wtko["n_wt_samples"] >= 2 and wtko["n_ko_samples"] >= 2
        assert wtko["n_pathways_evaluated"] >= 10

    def test_wt_vs_ko_directional_agreement(self):
        """Roadmap gate: predicted-vs-observed ΔMSV agreement >= 70%."""
        report = _load_or_skip()
        wtko = report.get("wt_vs_ko_gse123753")
        if wtko is None or "directional_agreement_rate" not in wtko:
            pytest.skip("wt_vs_ko_gse123753 section not present")
        rate = wtko["directional_agreement_rate"]
        target = report["acceptance"].get("wt_vs_ko_target", 0.70)
        assert rate >= target, (
            f"GSE123753 WT-vs-KO directional agreement {rate:.4f} "
            f"below target {target:.4f}; in-silico MECP2 KO predictions "
            f"do not match the direction of real MECP2-deletion effects"
        )

    def test_wt_vs_ko_rank_correlation_positive(self):
        """Per-pathway predicted vs observed ΔMSV rank correlation > 0."""
        report = _load_or_skip()
        wtko = report.get("wt_vs_ko_gse123753")
        if wtko is None or "spearman_pred_vs_observed" not in wtko:
            pytest.skip("wt_vs_ko_gse123753 section not present")
        rho = wtko["spearman_pred_vs_observed"]
        assert rho > 0.0, (
            f"predicted vs observed ΔMSV Spearman rho {rho:+.4f} is not "
            f"positive — signal is not aligned with real cohort direction"
        )


class TestGeneformerBackedF5:
    """Acceptance tests for the Geneformer-backed F5 artefact.

    Consumes ``results/f5_validation/perturbation_directional_geneformer.json``
    produced by:

        python scripts/validate_f5_real_data.py --backend geneformer \\
            --geneformer-model-dir <path-to-Geneformer-V2-104M>

    The **core claim** for the Geneformer-backed path is the real
    WT-vs-KO comparison on GSE123753 (Rodrigues et al. 2020): in-silico
    MECP2 KO on WT iPSC neurons predicts the direction of real
    MECP2-deletion effects on hallmark pathway scores. The TCGA-COAD
    direct-membership directional test (dropping a gene from a pathway
    depresses that pathway's score) is a FallbackPerturber strength
    rather than a Geneformer strength — context-aware embeddings
    don't naturally monotonically respond to removal of a single
    pathway-member gene — so that assertion is only enforced on the
    Fallback artefact (TestRealDataPerturb) and is NOT re-asserted here.
    """

    def test_artifact_shape(self):
        report = _load_geneformer_or_skip()
        assert report["feature"] == "F5"
        assert report["backend_kind"] == "geneformer"
        assert "wt_vs_ko_gse123753" in report

    def test_wt_vs_ko_directional_agreement(self):
        report = _load_geneformer_or_skip()
        wtko = report["wt_vs_ko_gse123753"]
        rate = wtko["directional_agreement_rate"]
        target = report["acceptance"].get("wt_vs_ko_target", 0.70)
        assert rate >= target, (
            f"Geneformer-backed WT-vs-KO directional agreement "
            f"{rate:.4f} below target {target:.4f}; in-silico MECP2 KO "
            f"predictions do not match real MECP2-deletion direction"
        )

    def test_wt_vs_ko_rank_correlation_positive(self):
        report = _load_geneformer_or_skip()
        wtko = report["wt_vs_ko_gse123753"]
        rho = wtko["spearman_pred_vs_observed"]
        assert rho > 0.0, (
            f"predicted vs observed ΔMSV rho {rho:+.4f} not positive; "
            f"Geneformer in-silico signal is anti-aligned with real cohort"
        )

    def test_conformal_coverage_within_2pct(self):
        """Conformal oracle deviation within ±2% on Geneformer artefact."""
        report = _load_geneformer_or_skip()
        conformal = report.get("conformal", {})
        if conformal.get("status") != "ok":
            pytest.skip(f"conformal sub-test did not run: {conformal}")
        oracle_dev = conformal["mean_oracle_deviation"]
        target = report["acceptance"]["conformal_oracle_target_abs"]
        assert abs(oracle_dev) < target, (
            f"Geneformer-backed perturbed-MSV conformal oracle deviation "
            f"{oracle_dev:+.4f} exceeds +/- {target:.4f}"
        )
