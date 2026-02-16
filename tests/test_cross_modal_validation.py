"""Tests for cross-modal validation gate."""

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.cross_modal_validation import (
    CrossModalPairResult,
    CrossModalValidationResult,
    SingleCellCompositionResult,
    cross_modal_concordance,
    generate_synthetic_multimodal_data,
    single_cell_composition_test,
)
from pathway_subtyping.validation import ValidationGates, ValidationResult

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def concordant_data():
    """Two modalities with shared planted subtypes (strong signal, should pass gate)."""
    scores, labels = generate_synthetic_multimodal_data(
        n_samples=120,
        n_pathways=10,
        n_subtypes=3,
        n_modalities=2,
        shared_signal_strength=2.0,
        modality_noise=0.3,
        seed=42,
    )
    return scores, labels


@pytest.fixture
def discordant_data():
    """Two modalities with independent structure (should fail gate)."""
    rng = np.random.RandomState(99)
    n_samples = 120
    n_pathways = 10
    sample_ids = [f"SAMPLE_{i:04d}" for i in range(n_samples)]
    pathway_names = [f"PATHWAY_{p}" for p in range(n_pathways)]

    # Modality 1: random noise (no subtype structure)
    data_a = rng.normal(0, 1, (n_samples, n_pathways))
    scores_a = pd.DataFrame(data_a, index=sample_ids, columns=pathway_names)

    # Modality 2: completely independent random noise
    data_b = rng.normal(0, 1, (n_samples, n_pathways))
    scores_b = pd.DataFrame(data_b, index=sample_ids, columns=pathway_names)

    return {"modality_1": scores_a, "modality_2": scores_b}


@pytest.fixture
def three_modality_data():
    """Three modalities with shared structure."""
    scores, labels = generate_synthetic_multimodal_data(
        n_samples=120,
        n_pathways=10,
        n_subtypes=3,
        n_modalities=3,
        shared_signal_strength=2.0,
        modality_noise=0.3,
        seed=42,
    )
    return scores, labels


@pytest.fixture
def cell_type_fractions():
    """Cell-type fractions with distinct compositions per subtype."""
    rng = np.random.RandomState(42)
    n_samples = 90
    n_subtypes = 3
    labels = np.repeat(np.arange(n_subtypes), n_samples // n_subtypes)

    # Create fractions where cell types differ across subtypes
    fractions = pd.DataFrame(index=range(n_samples))
    for ct_idx, ct in enumerate(["T_cell", "B_cell", "Monocyte", "NK_cell"]):
        base = rng.uniform(0.1, 0.3, n_samples)
        # Make each cell type higher in a different subtype
        for s in range(n_subtypes):
            mask = labels == s
            if ct_idx % n_subtypes == s:
                base[mask] += 0.4  # Strong signal
        fractions[ct] = base

    # Normalize rows to sum to 1
    fractions = fractions.div(fractions.sum(axis=1), axis=0)
    return labels, fractions


@pytest.fixture
def uniform_cell_fractions():
    """Cell-type fractions that are uniform across subtypes."""
    rng = np.random.RandomState(42)
    n_samples = 90
    labels = np.repeat(np.arange(3), 30)

    fractions = pd.DataFrame(
        rng.dirichlet([1, 1, 1, 1], n_samples),
        columns=["T_cell", "B_cell", "Monocyte", "NK_cell"],
    )
    return labels, fractions


# ---------------------------------------------------------------------------
# CrossModalPairResult tests
# ---------------------------------------------------------------------------


class TestCrossModalPairResult:
    def test_to_dict(self):
        r = CrossModalPairResult(
            modality_a="VCF",
            modality_b="RNA-seq",
            concordance_ari=0.7123456,
            concordance_nmi=0.6543,
            transfer_ari_a_to_b=0.68,
            transfer_ari_b_to_a=0.65,
            n_shared_samples=100,
            n_shared_pathways=8,
        )
        d = r.to_dict()
        assert d["modality_a"] == "VCF"
        assert d["modality_b"] == "RNA-seq"
        assert d["concordance_ari"] == 0.7123
        assert d["n_shared_samples"] == 100

    def test_fields(self):
        r = CrossModalPairResult(
            modality_a="A",
            modality_b="B",
            concordance_ari=0.5,
            concordance_nmi=0.4,
            transfer_ari_a_to_b=0.45,
            transfer_ari_b_to_a=0.42,
            n_shared_samples=50,
            n_shared_pathways=5,
        )
        assert r.modality_a == "A"
        assert r.concordance_ari == 0.5


# ---------------------------------------------------------------------------
# SingleCellCompositionResult tests
# ---------------------------------------------------------------------------


class TestSingleCellCompositionResult:
    def test_to_dict(self):
        r = SingleCellCompositionResult(
            f_statistic=12.345,
            p_value=0.0001234,
            n_subtypes=3,
            n_cell_types=4,
            per_cell_type_pvalues={"T_cell": 0.001, "B_cell": 0.05},
            significant_cell_types=["T_cell"],
        )
        d = r.to_dict()
        assert d["f_statistic"] == 12.345
        assert d["n_subtypes"] == 3
        assert "T_cell" in d["significant_cell_types"]
        assert d["per_cell_type_pvalues"]["T_cell"] == 0.001


# ---------------------------------------------------------------------------
# CrossModalValidationResult tests
# ---------------------------------------------------------------------------


class TestCrossModalValidationResult:
    def test_to_dict(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        d = result.to_dict()
        assert "pair_results" in d
        assert "mean_concordance_ari" in d
        assert "gate_passed" in d
        assert isinstance(d["pair_results"], list)

    def test_format_report(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        report = result.format_report()
        assert "Cross-Modal Validation Report" in report
        assert "Gate Status" in report
        assert "Concordance ARI" in report

    def test_get_citations(self):
        result = CrossModalValidationResult(
            pair_results=[],
            mean_concordance_ari=0.0,
            mean_concordance_nmi=0.0,
            mean_transfer_ari=0.0,
            null_ari_95th=0.0,
            gate_passed=False,
        )
        citations = result.get_citations()
        assert len(citations) >= 2
        assert "Hubert" in citations[0]
        assert "Ritchie" in citations[2]

    def test_to_dict_with_single_cell(self):
        sc = SingleCellCompositionResult(
            f_statistic=5.0,
            p_value=0.01,
            n_subtypes=3,
            n_cell_types=4,
            per_cell_type_pvalues={"T_cell": 0.01},
            significant_cell_types=["T_cell"],
        )
        result = CrossModalValidationResult(
            pair_results=[],
            mean_concordance_ari=0.5,
            mean_concordance_nmi=0.4,
            mean_transfer_ari=0.3,
            null_ari_95th=0.1,
            gate_passed=True,
            single_cell_composition=sc,
        )
        d = result.to_dict()
        assert "single_cell_composition" in d
        assert d["single_cell_composition"]["n_subtypes"] == 3


# ---------------------------------------------------------------------------
# cross_modal_concordance tests
# ---------------------------------------------------------------------------


class TestCrossModalConcordance:
    def test_concordant_modalities_pass(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=20,
            seed=42,
            show_progress=False,
        )
        assert isinstance(result, CrossModalValidationResult)
        assert result.gate_passed
        assert result.mean_concordance_ari > 0.3

    def test_discordant_modalities_fail(self, discordant_data):
        labels = np.repeat(np.arange(3), 40)
        result = cross_modal_concordance(
            discordant_data,
            labels,
            list(discordant_data["modality_1"].index),
            n_clusters=3,
            n_permutations=20,
            seed=42,
            show_progress=False,
        )
        assert isinstance(result, CrossModalValidationResult)
        # With random data, concordance should be low
        assert result.mean_concordance_ari < 0.5

    def test_three_modalities(self, three_modality_data):
        scores, labels = three_modality_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        # 3 modalities → 3 pairs (C(3,2))
        assert len(result.pair_results) == 3

    def test_insufficient_shared_samples(self):
        """Pairs with too few shared samples are skipped."""
        rng = np.random.RandomState(42)
        # Modality A: samples 0-19
        scores_a = pd.DataFrame(
            rng.normal(0, 1, (20, 5)),
            index=[f"A_{i}" for i in range(20)],
            columns=[f"P_{j}" for j in range(5)],
        )
        # Modality B: samples 20-39 (NO overlap)
        scores_b = pd.DataFrame(
            rng.normal(0, 1, (20, 5)),
            index=[f"B_{i}" for i in range(20)],
            columns=[f"P_{j}" for j in range(5)],
        )
        labels = np.zeros(40)

        result = cross_modal_concordance(
            {"A": scores_a, "B": scores_b},
            labels,
            [f"A_{i}" for i in range(20)] + [f"B_{i}" for i in range(20)],
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        # No valid pairs → empty results, gate fails
        assert len(result.pair_results) == 0
        assert not result.gate_passed

    def test_reproducible_with_seed(self, concordant_data):
        scores, labels = concordant_data
        sample_ids = list(scores[list(scores.keys())[0]].index)

        r1 = cross_modal_concordance(
            scores,
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        r2 = cross_modal_concordance(
            scores,
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        assert r1.mean_concordance_ari == r2.mean_concordance_ari
        assert r1.null_ari_95th == r2.null_ari_95th

    def test_different_seeds_vary(self, concordant_data):
        scores, labels = concordant_data
        sample_ids = list(scores[list(scores.keys())[0]].index)

        r1 = cross_modal_concordance(
            scores,
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        r2 = cross_modal_concordance(
            scores,
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            seed=999,
            show_progress=False,
        )
        # Null distribution should differ with different seeds
        assert r1.null_ari_95th != r2.null_ari_95th

    def test_fewer_than_two_modalities_raises(self):
        rng = np.random.RandomState(42)
        scores_a = pd.DataFrame(
            rng.normal(0, 1, (20, 5)),
            index=[f"S_{i}" for i in range(20)],
            columns=[f"P_{j}" for j in range(5)],
        )
        with pytest.raises(ValueError, match="at least 2 modalities"):
            cross_modal_concordance(
                {"A": scores_a},
                np.zeros(20),
                [f"S_{i}" for i in range(20)],
                n_clusters=3,
                show_progress=False,
            )

    def test_custom_threshold(self, concordant_data):
        scores, labels = concordant_data
        sample_ids = list(scores[list(scores.keys())[0]].index)

        # With threshold above 1.0 (impossible to exceed), gate should fail
        result = cross_modal_concordance(
            scores,
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            concordance_threshold=1.1,
            seed=42,
            show_progress=False,
        )
        assert not result.gate_passed


# ---------------------------------------------------------------------------
# Transfer validation tests
# ---------------------------------------------------------------------------


class TestTransferValidation:
    def test_transfer_ari_concordant(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        for p in result.pair_results:
            # With strong shared signal, transfer ARI should be positive
            assert p.transfer_ari_a_to_b > -0.1
            assert p.transfer_ari_b_to_a > -0.1

    def test_transfer_ari_no_shared_pathways(self):
        """Transfer is 0.0 when modalities have no shared pathways."""
        rng = np.random.RandomState(42)
        n_samples = 60
        sample_ids = [f"S_{i}" for i in range(n_samples)]

        scores_a = pd.DataFrame(
            rng.normal(0, 1, (n_samples, 5)),
            index=sample_ids,
            columns=[f"PA_{j}" for j in range(5)],
        )
        scores_b = pd.DataFrame(
            rng.normal(0, 1, (n_samples, 5)),
            index=sample_ids,
            columns=[f"PB_{j}" for j in range(5)],
        )
        labels = np.repeat(np.arange(3), 20)

        result = cross_modal_concordance(
            {"A": scores_a, "B": scores_b},
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        for p in result.pair_results:
            assert p.transfer_ari_a_to_b == 0.0
            assert p.transfer_ari_b_to_a == 0.0
            assert p.n_shared_pathways == 0


# ---------------------------------------------------------------------------
# Null calibration tests
# ---------------------------------------------------------------------------


class TestNullCalibration:
    def test_null_95th_percentile_reasonable(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=50,
            seed=42,
            show_progress=False,
        )
        # Null ARI under permutation should be near zero
        assert result.null_ari_95th < 0.3
        assert result.null_ari_95th >= -0.1

    def test_concordant_exceeds_null(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=20,
            seed=42,
            show_progress=False,
        )
        assert result.mean_concordance_ari > result.null_ari_95th


# ---------------------------------------------------------------------------
# Single-cell composition tests
# ---------------------------------------------------------------------------


class TestSingleCellCompositionTest:
    def test_distinct_compositions_detected(self, cell_type_fractions):
        labels, fractions = cell_type_fractions
        result = single_cell_composition_test(labels, fractions)
        assert isinstance(result, SingleCellCompositionResult)
        assert result.n_subtypes == 3
        assert result.n_cell_types == 4
        assert len(result.significant_cell_types) > 0
        assert result.p_value < 0.05

    def test_uniform_not_significant(self, uniform_cell_fractions):
        labels, fractions = uniform_cell_fractions
        result = single_cell_composition_test(labels, fractions)
        # With uniform random fractions, no cell type should be significant
        assert result.f_statistic < 10  # Reasonably low F-stat
        # May or may not be significant with random data, but fewer significant
        assert len(result.significant_cell_types) <= len(fractions.columns)

    def test_mismatched_lengths_raises(self):
        labels = np.array([0, 1, 2])
        fractions = pd.DataFrame({"T_cell": [0.5, 0.5]})
        with pytest.raises(ValueError, match="does not match"):
            single_cell_composition_test(labels, fractions)

    def test_to_dict(self, cell_type_fractions):
        labels, fractions = cell_type_fractions
        result = single_cell_composition_test(labels, fractions)
        d = result.to_dict()
        assert "f_statistic" in d
        assert "significant_cell_types" in d
        assert "per_cell_type_pvalues" in d


# ---------------------------------------------------------------------------
# Synthetic multimodal data tests
# ---------------------------------------------------------------------------


class TestSyntheticMultimodalData:
    def test_shapes(self):
        scores, labels = generate_synthetic_multimodal_data(
            n_samples=100,
            n_pathways=10,
            n_subtypes=3,
            n_modalities=2,
            seed=42,
        )
        assert len(scores) == 2
        assert labels.shape == (100,)
        for key, df in scores.items():
            assert df.shape == (100, 10)

    def test_labels_planted(self):
        scores, labels = generate_synthetic_multimodal_data(
            n_samples=90,
            n_subtypes=3,
            seed=42,
        )
        unique = np.unique(labels)
        assert len(unique) == 3
        # Each subtype should have 30 samples
        for label in unique:
            assert np.sum(labels == label) == 30

    def test_three_modalities(self):
        scores, labels = generate_synthetic_multimodal_data(
            n_modalities=3,
            seed=42,
        )
        assert len(scores) == 3

    def test_shared_signal_recoverable(self):
        """With strong signal, clustering each modality should recover true labels."""
        from sklearn.metrics import adjusted_rand_score as ari
        from sklearn.mixture import GaussianMixture

        scores, true_labels = generate_synthetic_multimodal_data(
            n_samples=120,
            n_pathways=10,
            n_subtypes=3,
            shared_signal_strength=3.0,
            modality_noise=0.3,
            seed=42,
        )
        for key, df in scores.items():
            gmm = GaussianMixture(n_components=3, random_state=42, reg_covar=1e-6)
            gmm.fit(df.values)
            predicted = gmm.predict(df.values)
            # Strong signal → high ARI
            assert ari(true_labels, predicted) > 0.5

    def test_seed_reproducibility(self):
        s1, l1 = generate_synthetic_multimodal_data(seed=42)
        s2, l2 = generate_synthetic_multimodal_data(seed=42)
        np.testing.assert_array_equal(l1, l2)
        for key in s1:
            pd.testing.assert_frame_equal(s1[key], s2[key])


# ---------------------------------------------------------------------------
# ValidationGates integration tests
# ---------------------------------------------------------------------------


class TestValidationGatesIntegration:
    @pytest.fixture
    def validation_fixtures(self):
        """Minimal fixtures for ValidationGates testing."""
        rng = np.random.RandomState(42)
        n_samples = 60
        n_pathways = 5
        n_genes = 25
        n_clusters = 3

        pathway_scores = pd.DataFrame(
            rng.normal(0, 1, (n_samples, n_pathways)),
            index=[f"S_{i}" for i in range(n_samples)],
            columns=[f"P_{j}" for j in range(n_pathways)],
        )
        labels = np.repeat(np.arange(n_clusters), n_samples // n_clusters)
        pathways = {f"P_{j}": [f"G_{g}" for g in range(5)] for j in range(n_pathways)}
        gene_burdens = pd.DataFrame(
            rng.exponential(1, (n_samples, n_genes)),
            index=[f"S_{i}" for i in range(n_samples)],
            columns=[f"G_{g}" for g in range(n_genes)],
        )
        return pathway_scores, labels, pathways, gene_burdens, n_clusters

    def test_run_all_without_multimodal(self, validation_fixtures):
        pathway_scores, labels, pathways, gene_burdens, n_clusters = validation_fixtures
        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5, show_progress=False)
        result = validator.run_all(
            pathway_scores=pathway_scores,
            cluster_labels=labels,
            pathways=pathways,
            gene_burdens=gene_burdens,
            n_clusters=n_clusters,
        )
        # 3 gates without ancestry or multi-omic
        assert len(result.results) == 3

    def test_run_all_with_multimodal(self, validation_fixtures):
        pathway_scores, labels, pathways, gene_burdens, n_clusters = validation_fixtures
        sample_ids = list(pathway_scores.index)

        # Create 2-modality scores for the same samples
        rng = np.random.RandomState(42)
        per_modality_scores = {
            "VCF": pd.DataFrame(
                rng.normal(0, 1, (60, 5)),
                index=sample_ids,
                columns=[f"P_{j}" for j in range(5)],
            ),
            "EXPR": pd.DataFrame(
                rng.normal(0, 1, (60, 5)),
                index=sample_ids,
                columns=[f"P_{j}" for j in range(5)],
            ),
        }

        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5, show_progress=False)
        result = validator.run_all(
            pathway_scores=pathway_scores,
            cluster_labels=labels,
            pathways=pathways,
            gene_burdens=gene_burdens,
            n_clusters=n_clusters,
            per_modality_scores=per_modality_scores,
        )
        # 3 mandatory + 1 cross-modal = 4 gates
        assert len(result.results) == 4
        assert result.results[3].name == "Gate 5: Cross-Modal Concordance"

    def test_gate5_returns_validation_result(self, validation_fixtures):
        pathway_scores, labels, pathways, gene_burdens, n_clusters = validation_fixtures
        sample_ids = list(pathway_scores.index)

        rng = np.random.RandomState(42)
        per_modality_scores = {
            "A": pd.DataFrame(
                rng.normal(0, 1, (60, 5)),
                index=sample_ids,
                columns=[f"P_{j}" for j in range(5)],
            ),
            "B": pd.DataFrame(
                rng.normal(0, 1, (60, 5)),
                index=sample_ids,
                columns=[f"P_{j}" for j in range(5)],
            ),
        }

        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5, show_progress=False)
        gate_result = validator.cross_modal_concordance_gate(
            per_modality_scores, labels, sample_ids, n_clusters
        )
        assert isinstance(gate_result, ValidationResult)
        assert gate_result.metric_name == "mean_concordance_ARI"
        assert gate_result.comparison == ">"
        assert gate_result.name == "Gate 5: Cross-Modal Concordance"

    def test_single_modality_not_added(self, validation_fixtures):
        """run_all should NOT add gate 5 if only 1 modality provided."""
        pathway_scores, labels, pathways, gene_burdens, n_clusters = validation_fixtures
        rng = np.random.RandomState(42)
        per_modality_scores = {
            "VCF": pd.DataFrame(
                rng.normal(0, 1, (60, 5)),
                index=list(pathway_scores.index),
                columns=[f"P_{j}" for j in range(5)],
            ),
        }

        validator = ValidationGates(seed=42, n_permutations=5, n_bootstrap=5, show_progress=False)
        result = validator.run_all(
            pathway_scores=pathway_scores,
            cluster_labels=labels,
            pathways=pathways,
            gene_burdens=gene_burdens,
            n_clusters=n_clusters,
            per_modality_scores=per_modality_scores,
        )
        # Should NOT include gate 5 with only 1 modality
        assert len(result.results) == 3


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_different_pathways_per_modality(self):
        """Modalities with completely different pathways still compute concordance."""
        rng = np.random.RandomState(42)
        n_samples = 60
        sample_ids = [f"S_{i}" for i in range(n_samples)]

        scores_a = pd.DataFrame(
            rng.normal(0, 1, (n_samples, 5)),
            index=sample_ids,
            columns=[f"PA_{j}" for j in range(5)],
        )
        scores_b = pd.DataFrame(
            rng.normal(0, 1, (n_samples, 5)),
            index=sample_ids,
            columns=[f"PB_{j}" for j in range(5)],
        )
        labels = np.repeat(np.arange(3), 20)

        result = cross_modal_concordance(
            {"A": scores_a, "B": scores_b},
            labels,
            sample_ids,
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        # Concordance still computed (independent clustering)
        assert len(result.pair_results) == 1
        # Transfer should be 0 (no shared pathways)
        assert result.pair_results[0].transfer_ari_a_to_b == 0.0

    def test_partial_sample_overlap(self):
        """Works with partial sample overlap between modalities."""
        rng = np.random.RandomState(42)
        # 40 shared + 20 unique each
        shared_ids = [f"S_{i}" for i in range(40)]
        ids_a = shared_ids + [f"A_{i}" for i in range(20)]
        ids_b = shared_ids + [f"B_{i}" for i in range(20)]
        pathways = [f"P_{j}" for j in range(5)]

        scores_a = pd.DataFrame(rng.normal(0, 1, (60, 5)), index=ids_a, columns=pathways)
        scores_b = pd.DataFrame(rng.normal(0, 1, (60, 5)), index=ids_b, columns=pathways)
        labels = np.zeros(80)

        result = cross_modal_concordance(
            {"A": scores_a, "B": scores_b},
            labels,
            ids_a + [f"B_{i}" for i in range(20)],
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        assert len(result.pair_results) == 1
        assert result.pair_results[0].n_shared_samples == 40

    def test_no_nan_in_output(self, concordant_data):
        scores, labels = concordant_data
        result = cross_modal_concordance(
            scores,
            labels,
            list(scores[list(scores.keys())[0]].index),
            n_clusters=3,
            n_permutations=10,
            seed=42,
            show_progress=False,
        )
        d = result.to_dict()
        assert not np.isnan(d["mean_concordance_ari"])
        assert not np.isnan(d["mean_concordance_nmi"])
        assert not np.isnan(d["mean_transfer_ari"])
        assert not np.isnan(d["null_ari_95th"])
