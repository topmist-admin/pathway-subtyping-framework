"""
Tests for Gate 7 (somatic mode) — the cancer counterpart to feature-level
genetic anchoring.

Where the feature-level gate asks "are a subtype's defining genes enriched for
germline risk?", the somatic gate asks "do the tumors in a subtype carry a
somatic stratum (BRAF/KRAS/MSI, CNA, signature) more than the others?". It is the
confound gate's statistic (chi-square + Bergsma Cramér's V, BH-adjusted) with
inverted polarity: a strong association is the *desired* positive evidence.

These tests reconstruct a colorectal-style scenario synthetically and
deterministically — a subtype aligned with a BRAF-V600E stratum and a KRAS
stratum, versus a random stratum that should not anchor.
"""

import numpy as np
import pytest

from pathway_subtyping.genetics import somatic_alignment
from pathway_subtyping.validation import ValidationGates


class TestSomaticAlignmentCore:
    def _scenario(self, seed=20260716):
        rng = np.random.default_rng(seed)
        n = 200
        # partition: subtype 0 vs 1
        clusters = np.array([0] * 100 + [1] * 100)
        # BRAF carrier concentrated in subtype 0 (strong, noisy alignment)
        braf = np.where(clusters == 0, "mut", "wt")
        flip = rng.random(n) < 0.12
        braf = np.where(flip, np.where(braf == "mut", "wt", "mut"), braf)
        # a random stratum independent of the partition
        random_stratum = rng.integers(0, 2, n).astype(str)
        return clusters, braf, random_stratum

    def test_driver_anchors_random_does_not(self):
        clusters, braf, rnd = self._scenario()
        res = somatic_alignment(
            clusters,
            {"BRAF_V600E": braf, "random": rnd},
            cramers_v_min=0.30,
            alpha=0.05,
        )
        assert res["passed"] is True
        assert "BRAF_V600E" in res["anchored_strata"]
        assert "random" not in res["anchored_strata"]
        assert res["best_stratum"] == "BRAF_V600E"
        assert res["per_stratum"]["BRAF_V600E"]["cramers_v"] >= 0.5
        assert res["per_stratum"]["BRAF_V600E"]["significant"] is True
        assert res["per_stratum"]["random"]["anchored"] is False

    def test_missing_values_dropped_per_stratum(self):
        clusters = np.array([0] * 50 + [1] * 50)
        # MSI known for only half the tumors; aligned where known
        msi = np.array(["high" if c == 0 else "low" for c in clusters], dtype=object)
        msi[25:75] = None  # unknown for a middle block
        res = somatic_alignment(clusters, {"MSI": msi})
        assert res["per_stratum"]["MSI"]["n"] == 50  # only the non-missing tumors
        assert res["per_stratum"]["MSI"]["anchored"] is True

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            somatic_alignment(np.array([0, 1, 0]), {"BRAF": np.array(["mut", "wt"])})


class TestSomaticAnchoringGate:
    def _scenario(self, seed=20260716):
        rng = np.random.default_rng(seed)
        n = 200
        clusters = np.array([0] * 100 + [1] * 100)
        braf = np.where(clusters == 0, "mut", "wt")
        flip = rng.random(n) < 0.12
        braf = np.where(flip, np.where(braf == "mut", "wt", "mut"), braf)
        rnd = rng.integers(0, 2, n).astype(str)
        return clusters, braf, rnd

    def test_gate_passes_on_driver_alignment(self):
        clusters, braf, rnd = self._scenario()
        gates = ValidationGates(show_progress=False)
        result = gates.somatic_anchoring_gate(
            cluster_labels=clusters,
            somatic_strata={"BRAF_V600E": braf, "random": rnd},
        )
        assert result.passed is True
        assert result.name == "Gate 7 (somatic): Genetic Anchoring"
        assert result.metric_name == "max_somatic_cramers_v"
        assert result.details["anchored_strata"] == ["BRAF_V600E"]
        assert result.details["gate_polarity"] == "positive_evidence"
        assert result.metric_value >= 0.5

    def test_gate_fails_when_no_stratum_aligns(self):
        rng = np.random.default_rng(1)
        n = 200
        clusters = rng.integers(0, 2, n)
        # strata independent of the partition
        strata = {
            "BRAF": rng.integers(0, 2, n).astype(str),
            "KRAS": rng.integers(0, 2, n).astype(str),
        }
        result = ValidationGates(show_progress=False).somatic_anchoring_gate(
            cluster_labels=clusters, somatic_strata=strata
        )
        assert result.passed is False
        assert result.details["anchored_strata"] == []
        assert "low-power" in result.details["interpretation"]

    def test_run_all_includes_somatic_gate_when_supplied(self):
        import pandas as pd

        rng = np.random.default_rng(20260716)
        n = 60
        pathway_scores = pd.DataFrame(rng.normal(size=(n, 8)), columns=[f"PW{i}" for i in range(8)])
        gene_burdens = pd.DataFrame(rng.normal(size=(n, 12)), columns=[f"G{i}" for i in range(12)])
        pathways = {f"PW{i}": [f"G{i}", f"G{(i + 1) % 12}"] for i in range(8)}
        clusters = np.array([0] * 30 + [1] * 30)
        braf = np.where(clusters == 0, "mut", "wt")

        gates = ValidationGates(
            seed=20260716, n_permutations=10, n_bootstrap=10, show_progress=False
        )
        res = gates.run_all(
            pathway_scores=pathway_scores,
            cluster_labels=clusters,
            pathways=pathways,
            gene_burdens=gene_burdens,
            n_clusters=2,
            gmm_seed=20260716,
            somatic_anchoring={"somatic_strata": {"BRAF_V600E": braf}},
        )
        names = [r.name for r in res.results]
        assert "Gate 7 (somatic): Genetic Anchoring" in names
        gate = next(r for r in res.results if r.name == "Gate 7 (somatic): Genetic Anchoring")
        assert gate.passed is True
        assert gate.details["anchored_strata"] == ["BRAF_V600E"]


class TestSparseTableGuard:
    """Chi-square is invalid AND anticonservative on sparse tables.

    This is the positive-evidence gate for cancer, so a rare driver in a small
    cohort could otherwise fabricate a "somatic anchor". Concrete case: min
    expected 3.0, chi-square p=0.0079 vs Fisher exact p=0.0202 — a 2.6x
    overstatement of significance.
    """

    def test_sparse_2x2_falls_back_to_fisher_exact(self):
        labels = np.array([0] * 20 + [1] * 20)
        stratum = np.array(["MUT"] * 6 + ["WT"] * 14 + ["WT"] * 20)
        e = somatic_alignment(labels, {"DRIVER": stratum}, cramers_v_min=0.30)["per_stratum"][
            "DRIVER"
        ]
        assert e["sparse_table"] is True
        assert e["min_expected_count"] == pytest.approx(3.0, abs=1e-6)
        assert e["test"] == "fisher_exact"
        # the exact p, not the anticonservative chi-square value of 0.0079
        assert e["p_value"] == pytest.approx(0.0202, abs=1e-3)

    def test_well_powered_table_still_uses_chi2(self):
        labels = np.array([0] * 40 + [1] * 40)
        stratum = np.array(["MUT"] * 30 + ["WT"] * 10 + ["MUT"] * 10 + ["WT"] * 30)
        e = somatic_alignment(labels, {"DRIVER": stratum})["per_stratum"]["DRIVER"]
        assert e["sparse_table"] is False
        assert e["test"] == "chi2"
        assert e["min_expected_count"] >= 5.0

    def test_sparse_non_2x2_cannot_anchor(self):
        """No exact fallback exists for r x c, so it must not anchor."""
        labels = np.array([0] * 15 + [1] * 15 + [2] * 15)
        stratum = np.array(
            ["A"] * 13
            + ["B"] * 1
            + ["C"] * 1
            + ["A"] * 1
            + ["B"] * 13
            + ["C"] * 1
            + ["A"] * 1
            + ["B"] * 1
            + ["C"] * 13
        )
        e = somatic_alignment(labels, {"DRIVER": stratum}, cramers_v_min=0.30)["per_stratum"][
            "DRIVER"
        ]
        if e["sparse_table"] and e["test"] == "chi2":
            assert e["anchored"] is False, "invalid chi-square must not anchor"
            assert "not eligible to anchor" in e.get("note", "")
