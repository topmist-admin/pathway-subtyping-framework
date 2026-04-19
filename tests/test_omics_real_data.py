"""
Real-data acceptance tests for the F10 multi-omics fusion layer.

Consumes the JSON artifact produced by scripts/validate_f10_real_data.py
and asserts the roadmap criterion: fused (RNA + protein) pathway score
matrix classifies PBMC cell types at least 3% more accurately than the
RNA-only score matrix on a paired CITE-seq reference, with the 95%
bootstrap CI lower bound strictly greater than zero.

Cohort: 10x Genomics pbmc_1k_protein_v3 (713 cells, 17-antibody panel).
Labels are derived from canonical ADT gating; the uplift reflects the
downstream signal the protein modality adds over RNA alone.

If the JSON artifact is absent, the tests are skipped — CI runs without
the CITE-seq data will not fail.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT_PATH = REPO_ROOT / "results" / "f10_validation" / "fusion_uplift.json"


def _load_or_skip() -> dict:
    if not ARTIFACT_PATH.exists():
        pytest.skip(
            f"no real-data validation artifact at {ARTIFACT_PATH} "
            f"(run scripts/validate_f10_real_data.py to generate)"
        )
    return json.loads(ARTIFACT_PATH.read_text())


class TestRealDataFusion:

    def test_artifact_shape(self):
        report = _load_or_skip()
        assert report["feature"] == "F10"
        assert "cohort" in report
        assert report["cohort"]["n_cells_gated"] >= 300
        assert report["cohort"]["n_antibodies"] >= 10
        assert len(report["pathways"]) >= 5

    def test_both_classifiers_above_chance(self):
        """Sanity check: both RNA-only and fused beat uniform chance."""
        report = _load_or_skip()
        n_classes = len(report["cohort"]["label_counts"])
        chance = 1.0 / n_classes
        assert report["rna_only_accuracy"] > chance, (
            f"RNA-only accuracy {report['rna_only_accuracy']:.4f} at or "
            f"below chance ({chance:.4f}) — feature pipeline is broken"
        )
        assert report["fused_accuracy"] > chance, (
            f"fused accuracy {report['fused_accuracy']:.4f} at or below "
            f"chance ({chance:.4f})"
        )

    def test_fusion_uplift_at_least_3pct(self):
        """Roadmap gate: fused - rna_only >= 0.03 on cell-type accuracy."""
        report = _load_or_skip()
        uplift = report["uplift"]
        target = report["acceptance"]["uplift_target"]
        assert uplift >= target, (
            f"fusion uplift {uplift:+.4f} below {target:+.4f}; protein "
            f"modality is not adding downstream classification signal"
        )

    def test_uplift_ci_lower_bound_positive(self):
        """Roadmap gate: bootstrap 95% CI lower bound > 0 on uplift."""
        report = _load_or_skip()
        ci_low = report["uplift_bootstrap_ci"]["ci_low"]
        assert ci_low > 0, (
            f"uplift 95% CI lower bound {ci_low:+.4f} does not exclude "
            f"zero — uplift is not statistically significant"
        )
