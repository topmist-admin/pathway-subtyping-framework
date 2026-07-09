"""
Tests for the Confound Association Gate.

This is the gate whose absence let a brain-region artifact pass the full
validation battery on GSE80655: a k=3 partition passed every stability control
(bootstrap ARI ~0.92) yet was really a cortex-vs-accumbens region classifier
(Cramer's V ~0.67, chi2 ~125, p ~4e-26) and was independent of diagnosis
(p ~0.41). The gate must flag region as a confound while clearing diagnosis.
"""

import numpy as np
import pytest

from pathway_subtyping.validation import ValidationGates, cramers_v


class TestCramersV:
    def test_perfect_association(self):
        # Block-diagonal table -> V = 1.0 (no bias correction to keep it exact).
        table = np.array([[50, 0], [0, 50]])
        assert cramers_v(table, bias_correction=False) == pytest.approx(1.0)

    def test_no_association(self):
        # Independent table -> V ~ 0.
        table = np.array([[25, 25], [25, 25]])
        assert cramers_v(table, bias_correction=False) == pytest.approx(0.0, abs=1e-9)

    def test_degenerate_table_returns_zero(self):
        assert cramers_v(np.array([[10]])) == 0.0
        assert cramers_v(np.array([[10, 5]])) == 0.0  # single row


class TestConfoundGateGSE80655Like:
    """Synthetic reconstruction of the GSE80655 region-artifact scenario."""

    def _make_scenario(self, seed=20260708):
        rng = np.random.default_rng(seed)
        # 141 samples across 3 regions; partition tracks region almost exactly.
        regions = np.array(["AnCg"] * 47 + ["DLPFC"] * 55 + ["nAcc"] * 39)
        clusters = np.where(regions == "nAcc", 0, np.where(regions == "DLPFC", 1, 2))
        # inject a tiny bit of region/cluster noise so V < 1 (like the real 0.67)
        flip = rng.random(len(clusters)) < 0.12
        clusters = np.where(flip, rng.integers(0, 3, len(clusters)), clusters)
        # diagnosis independent of cluster (SCZ vs Control, ~balanced per region)
        diagnosis = rng.permutation(
            np.array(["SCZ"] * 70 + ["Control"] * 71)
        )
        return clusters, regions, diagnosis

    def test_region_flagged_diagnosis_cleared(self):
        clusters, regions, diagnosis = self._make_scenario()
        gates = ValidationGates(show_progress=False)
        result = gates.confound_association_gate(
            cluster_labels=clusters,
            confounds={"region": regions, "diagnosis": diagnosis},
        )
        pc = result.details["per_confound"]

        # Region: strong, significant association -> gate FAILS on it.
        assert pc["region"]["cramers_v"] >= 0.5
        assert pc["region"]["significant"] is True
        assert "region" in result.details["failing_confounds"]
        assert result.passed is False

        # Diagnosis: treated as biology-of-interest, never a failing confound.
        assert pc["diagnosis"]["is_biology_of_interest"] is True
        assert "diagnosis" not in result.details["failing_confounds"]
        # and it should be weakly associated (independent by construction)
        assert pc["diagnosis"]["cramers_v"] < 0.3

    def test_clean_partition_passes(self):
        # A partition independent of all nuisance confounds passes.
        rng = np.random.default_rng(1)
        n = 150
        clusters = rng.integers(0, 3, n)
        batch = rng.integers(0, 2, n).astype(str)
        region = rng.integers(0, 3, n).astype(str)
        gates = ValidationGates(show_progress=False)
        result = gates.confound_association_gate(
            cluster_labels=clusters,
            confounds={"batch": batch, "region": region},
        )
        assert result.passed is True
        assert result.details["failing_confounds"] == []

    def test_run_all_includes_gate_when_confounds_present(self):
        import pandas as pd

        rng = np.random.default_rng(20260708)
        n = 60
        pathway_scores = pd.DataFrame(
            rng.normal(size=(n, 8)), columns=[f"PW{i}" for i in range(8)]
        )
        gene_burdens = pd.DataFrame(
            rng.normal(size=(n, 12)), columns=[f"G{i}" for i in range(12)]
        )
        pathways = {f"PW{i}": [f"G{i}", f"G{(i+1) % 12}"] for i in range(8)}
        clusters = rng.integers(0, 2, n)
        region = np.where(clusters == 0, "A", "B")  # perfectly confounded

        gates = ValidationGates(
            seed=20260708, n_permutations=10, n_bootstrap=10, show_progress=False
        )
        res = gates.run_all(
            pathway_scores=pathway_scores,
            cluster_labels=clusters,
            pathways=pathways,
            gene_burdens=gene_burdens,
            n_clusters=2,
            gmm_seed=20260708,
            confounds={"region": region},
        )
        names = [r.name for r in res.results]
        assert "Confound Association Gate" in names
        confound = next(r for r in res.results if r.name == "Confound Association Gate")
        assert confound.passed is False  # perfect region confound
