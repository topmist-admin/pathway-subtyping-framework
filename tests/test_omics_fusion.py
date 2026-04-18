"""
Tests for v0.6 F10 multi-omics pathway-score fusion.

Covers:
    - ATACScorer: peak -> gene collapse, pathway scoring, aggregation modes
    - ProteomicsScorer: optional gene -> protein mapping, min members
    - MultiOmicsFusion: weight normalisation, alignment on sample and
      pathway intersections, learn_weights grid search
    - DiscordanceReport: rho + absolute-difference flagging logic
    - Roadmap acceptance: fused score improves 1-NN cell-type
      classification accuracy by >= 3% over RNA-only on a paired
      CITE-seq-style synthetic reference
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.omics import (
    ATACScorer,
    DiscordanceReport,
    FusionResult,
    FusionWeights,
    MultiOmicsFusion,
    ProteomicsScorer,
    flag_discordant_pathways,
    score_atac_pathways,
    score_proteomics_pathways,
)


# --------------------------------------------------------------------------- #
# ATAC scorer
# --------------------------------------------------------------------------- #

@pytest.fixture
def atac_cohort():
    rng = np.random.default_rng(0)
    n_samples = 40
    peaks = [f"peak_{i}" for i in range(30)]
    accessibility = pd.DataFrame(
        rng.standard_normal((n_samples, len(peaks))),
        columns=peaks,
        index=[f"sample_{i}" for i in range(n_samples)],
    )
    # Every 3 peaks link to a single gene; 10 genes total
    peak_to_gene = {p: f"GENE_{i // 3}" for i, p in enumerate(peaks)}
    pathways = {
        "PATH_A": ["GENE_0", "GENE_1", "GENE_2"],
        "PATH_B": ["GENE_3", "GENE_4", "GENE_5"],
    }
    return accessibility, peak_to_gene, pathways


class TestATACScorer:

    def test_score_shape(self, atac_cohort):
        acc, p2g, pathways = atac_cohort
        scored = ATACScorer(peak_to_gene=p2g).score(acc, pathways)
        assert scored.shape == (len(acc), 2)
        assert list(scored.columns) == ["PATH_A", "PATH_B"]

    def test_aggregation_modes(self, atac_cohort):
        acc, p2g, pathways = atac_cohort
        for agg in ("mean", "sum", "max"):
            scored = ATACScorer(peak_to_gene=p2g, aggregation=agg).score(acc, pathways)
            assert scored.shape == (len(acc), 2)

    def test_rejects_bad_aggregation(self):
        with pytest.raises(ValueError):
            ATACScorer(peak_to_gene={}, aggregation="median")

    def test_zero_normalised_output(self, atac_cohort):
        acc, p2g, pathways = atac_cohort
        scored = ATACScorer(peak_to_gene=p2g).score(acc, pathways)
        np.testing.assert_allclose(scored.mean(axis=0).values, 0.0, atol=1e-10)

    def test_missing_peaks_ignored(self):
        acc = pd.DataFrame(
            np.random.default_rng(0).standard_normal((10, 2)),
            columns=["real_peak_1", "real_peak_2"],
        )
        p2g = {"missing_peak": "GENE_X", "real_peak_1": "GENE_Y"}
        pathways = {"PATH": ["GENE_Y"]}
        scored = ATACScorer(peak_to_gene=p2g).score(acc, pathways, min_genes_per_pathway=1)
        assert scored.shape == (10, 1)

    def test_convenience_helper(self, atac_cohort):
        acc, p2g, pathways = atac_cohort
        a = ATACScorer(peak_to_gene=p2g).score(acc, pathways)
        b = score_atac_pathways(acc, p2g, pathways)
        pd.testing.assert_frame_equal(a, b)


# --------------------------------------------------------------------------- #
# Proteomics scorer
# --------------------------------------------------------------------------- #

class TestProteomicsScorer:

    def test_direct_protein_ids(self):
        rng = np.random.default_rng(0)
        abundance = pd.DataFrame(
            rng.standard_normal((20, 6)),
            columns=["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E", "GENE_F"],
        )
        pathways = {"PATH_1": ["GENE_A", "GENE_B", "GENE_C"]}
        scored = ProteomicsScorer().score(abundance, pathways)
        assert scored.shape == (20, 1)

    def test_gene_to_protein_translation(self):
        rng = np.random.default_rng(0)
        abundance = pd.DataFrame(
            rng.standard_normal((20, 3)),
            columns=["P01", "P02", "P03"],
        )
        pathways = {"PATH_1": ["GENE_A", "GENE_B"]}
        mapping = {"GENE_A": "P01", "GENE_B": "P02", "GENE_C": "P99"}
        scored = ProteomicsScorer(gene_to_protein=mapping).score(
            abundance, pathways
        )
        assert scored.shape == (20, 1)

    def test_min_members_filter(self):
        abundance = pd.DataFrame(
            np.zeros((5, 3)), columns=["A", "B", "C"],
        )
        pathways = {"PATH": ["A"]}  # only 1 member — below default 2
        scored = ProteomicsScorer().score(abundance, pathways)
        assert scored.empty

    def test_convenience_helper(self):
        rng = np.random.default_rng(0)
        abundance = pd.DataFrame(
            rng.standard_normal((10, 4)),
            columns=["A", "B", "C", "D"],
        )
        pathways = {"P1": ["A", "B"], "P2": ["C", "D"]}
        a = ProteomicsScorer().score(abundance, pathways)
        b = score_proteomics_pathways(abundance, pathways)
        pd.testing.assert_frame_equal(a, b)


# --------------------------------------------------------------------------- #
# FusionWeights
# --------------------------------------------------------------------------- #

class TestFusionWeights:

    def test_auto_normalisation(self):
        w = FusionWeights(rna=2.0, atac=1.0, protein=1.0)
        assert w.rna == pytest.approx(0.5)
        assert w.atac == pytest.approx(0.25)
        assert w.protein == pytest.approx(0.25)
        total = w.rna + w.atac + w.protein
        assert total == pytest.approx(1.0)

    def test_negative_weight_raises(self):
        with pytest.raises(ValueError, match="must be >= 0"):
            FusionWeights(rna=-1.0, atac=1.0)

    def test_all_zero_raises(self):
        with pytest.raises(ValueError, match="positive"):
            FusionWeights(rna=0.0, atac=0.0, protein=0.0)


# --------------------------------------------------------------------------- #
# MultiOmicsFusion
# --------------------------------------------------------------------------- #

def _synthetic_citeseq(
    n_per_cluster: int = 50, n_pathways: int = 8, seed: int = 0
):
    """Synthetic paired RNA + protein cohort with 3 cell types.

    True pathway profile per cluster; RNA and protein observations are
    independently noisy copies. Each modality alone has non-trivial
    classification error; fusing both reduces the error.
    """
    rng = np.random.default_rng(seed)
    cluster_profiles = rng.standard_normal((3, n_pathways)) * 0.9
    cluster_ids = np.repeat(np.arange(3), n_per_cluster)
    n = 3 * n_per_cluster

    true = cluster_profiles[cluster_ids]
    noise_scale = 1.6  # high enough that RNA-only misses some cells
    rna = pd.DataFrame(
        true + rng.normal(0, noise_scale, size=true.shape),
        columns=[f"PATH_{i}" for i in range(n_pathways)],
        index=[f"cell_{i}" for i in range(n)],
    )
    protein = pd.DataFrame(
        true + rng.normal(0, noise_scale, size=true.shape),
        columns=[f"PATH_{i}" for i in range(n_pathways)],
        index=[f"cell_{i}" for i in range(n)],
    )
    labels = pd.Series(cluster_ids, index=rna.index, name="cluster")
    return rna, protein, labels


class TestMultiOmicsFusion:

    def test_rejects_no_modalities(self):
        with pytest.raises(ValueError, match="at least one"):
            MultiOmicsFusion().fuse()

    def test_rejects_disjoint_samples(self):
        df_a = pd.DataFrame(np.zeros((3, 2)), index=["a", "b", "c"],
                            columns=["P1", "P2"])
        df_b = pd.DataFrame(np.zeros((3, 2)), index=["x", "y", "z"],
                            columns=["P1", "P2"])
        with pytest.raises(ValueError, match="shared samples"):
            MultiOmicsFusion().fuse(rna=df_a, protein=df_b)

    def test_default_weights_per_supplied_modality(self):
        df = pd.DataFrame(np.arange(6).reshape(3, 2), columns=["P1", "P2"])
        result = MultiOmicsFusion().fuse(rna=df, protein=df)
        assert result.weights.rna == pytest.approx(0.5)
        assert result.weights.protein == pytest.approx(0.5)
        assert result.weights.atac == 0.0

    def test_fused_shape_equals_intersection(self):
        rna, prot, _ = _synthetic_citeseq()
        result = MultiOmicsFusion().fuse(rna=rna, protein=prot)
        assert result.fused.shape == rna.shape
        assert isinstance(result, FusionResult)


# --------------------------------------------------------------------------- #
# Roadmap acceptance: fused >= RNA-only + 3%
# --------------------------------------------------------------------------- #

class TestRoadmapAccuracyUplift:

    def _accuracy(self, features: pd.DataFrame, labels: pd.Series) -> float:
        X = features.loc[labels.index].to_numpy(dtype=float)
        y = labels.to_numpy()
        norms = np.linalg.norm(X, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Xn = X / norms
        sim = Xn @ Xn.T
        np.fill_diagonal(sim, -np.inf)
        predictions = y[sim.argmax(axis=1)]
        return float((predictions == y).mean())

    def test_fusion_beats_rna_only_by_3pct(self):
        rna, prot, labels = _synthetic_citeseq()
        rna_acc = self._accuracy(rna, labels)
        fused = MultiOmicsFusion().fuse(
            rna=rna, protein=prot,
            weights=FusionWeights(rna=1.0, protein=1.0),
        )
        fused_acc = self._accuracy(fused.fused, labels)
        assert fused_acc >= rna_acc + 0.03, (
            f"fused 1-NN accuracy {fused_acc:.3f} did not beat RNA-only "
            f"{rna_acc:.3f} by >= 3% (delta = {fused_acc - rna_acc:+.3f})"
        )

    def test_learn_weights_finds_nontrivial_mix(self):
        rna, prot, labels = _synthetic_citeseq()
        fusion = MultiOmicsFusion()
        weights = fusion.learn_weights(
            labels=labels, rna=rna, protein=prot, grid_step=0.2,
        )
        # A non-degenerate mix should result (not all RNA nor all protein)
        assert 0.0 < weights.rna < 1.0
        assert 0.0 < weights.protein < 1.0


# --------------------------------------------------------------------------- #
# Discordance
# --------------------------------------------------------------------------- #

class TestDiscordance:

    def test_concordant_pathway_not_flagged(self):
        rng = np.random.default_rng(0)
        n = 50
        true = rng.standard_normal(n)
        rna = pd.DataFrame({"PATH_A": true + rng.normal(0, 0.1, n)})
        prot = pd.DataFrame({"PATH_A": true + rng.normal(0, 0.1, n)})
        report = flag_discordant_pathways(rna, prot)
        assert "PATH_A" not in report.discordant_pathways

    def test_discordant_pathway_flagged(self):
        rng = np.random.default_rng(0)
        n = 50
        rna = pd.DataFrame({"PATH_X": rng.standard_normal(n)})
        # Protein: uncorrelated AND large absolute difference (shift + noise)
        prot = pd.DataFrame({
            "PATH_X": 3.0 + rng.standard_normal(n) * 1.5,
        }, index=rna.index)
        report = flag_discordant_pathways(
            rna, prot, threshold_rho=0.3, threshold_abs_diff=1.0,
        )
        assert "PATH_X" in report.discordant_pathways

    def test_summary_renders(self):
        rna = pd.DataFrame({"P": np.random.default_rng(0).standard_normal(10)})
        prot = pd.DataFrame({"P": np.random.default_rng(1).standard_normal(10)})
        report = flag_discordant_pathways(rna, prot)
        assert "DiscordanceReport" in report.summary()

    def test_no_overlap_raises(self):
        rna = pd.DataFrame({"PATH_A": [1.0, 2.0]}, index=["a", "b"])
        prot = pd.DataFrame({"PATH_B": [1.0, 2.0]}, index=["x", "y"])
        with pytest.raises(ValueError, match="share"):
            flag_discordant_pathways(rna, prot)

    def test_to_dict_structure(self):
        rna = pd.DataFrame({"P": np.random.default_rng(0).standard_normal(20)})
        prot = pd.DataFrame({"P": np.random.default_rng(0).standard_normal(20)})
        d = flag_discordant_pathways(rna, prot).to_dict()
        assert {"n_pathways", "n_discordant", "discordant_pathways"}.issubset(d)
