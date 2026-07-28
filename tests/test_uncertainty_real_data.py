"""
Real-data acceptance tests for the F1 uncertainty layer.

Consumes JSON artifacts produced by scripts/validate_f1_real_data.py and
asserts the roadmap criterion: conformal coverage within +/- 2% of the
(finite-sample oracle) target on TCGA-COAD and GSE28521.

If the JSON artifacts are absent, the tests are skipped — the validation
script is an offline step that depends on non-packaged datasets, so CI
runs without the data will not fail.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT_DIR = REPO_ROOT / "results" / "f1_validation"

DATASETS = [
    ("TCGA-COAD", ARTIFACT_DIR / "TCGA-COAD_conformal_coverage.json"),
    ("GSE28521", ARTIFACT_DIR / "GSE28521_conformal_coverage.json"),
]

TARGET_COVERAGES = ("0.80", "0.90", "0.95")


def _load_or_skip(path: Path) -> dict:
    if not path.exists():
        pytest.skip(
            f"no real-data validation artifact at {path} "
            f"(run scripts/validate_f1_real_data.py to generate)"
        )
    return json.loads(path.read_text())


@pytest.mark.parametrize("dataset,path", DATASETS)
class TestRealDataConformal:
    """Real-data acceptance — one class per dataset artifact."""

    def test_artifact_shape(self, dataset, path):
        report = _load_or_skip(path)
        assert report["dataset"] == dataset
        assert report["n_samples"] > 0
        assert report["n_pathways"] > 0
        assert set(TARGET_COVERAGES).issubset(report["summary"].keys())

    @pytest.mark.parametrize("target", TARGET_COVERAGES)
    def test_oracle_deviation_within_2pct(self, dataset, path, target):
        """Roadmap acceptance: coverage within +/- 2% of finite-sample oracle.

        The oracle accounts for split-conformal's discretization at small
        n_cal — the theoretical realized coverage is ceil((n+1)(1-alpha))/(n+1)
        rather than (1-alpha). That offset is a property of the algorithm,
        not a deficiency of the data or the implementation.
        """
        report = _load_or_skip(path)
        stats = report["summary"][target]
        assert (
            stats["mean_oracle_deviation"] is not None
        ), f"{dataset} target={target}: no runs recorded"
        assert abs(stats["mean_oracle_deviation"]) < 0.02, (
            f"{dataset} target={target}: oracle deviation "
            f"{stats['mean_oracle_deviation']:+.4f} exceeds +/- 2%"
        )

    @pytest.mark.parametrize("target", TARGET_COVERAGES)
    def test_enough_runs_for_statistical_confidence(self, dataset, path, target):
        report = _load_or_skip(path)
        stats = report["summary"][target]
        assert stats["n_runs"] >= 150, (
            f"{dataset} target={target}: only {stats['n_runs']} runs "
            f"— rerun scripts/validate_f1_real_data.py with more seeds"
        )
