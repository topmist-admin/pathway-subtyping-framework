"""
Tests for v0.6 F9 sequence-level off-target scoring.

Covers:
    - DNA validation and Hamming similarity primitives
    - SimilarityBackend: baseline scoring on synthetic candidates
    - SimulatedEvo2Backend: functional score surfaces seed-mismatch cases
      that the Hamming baseline cannot distinguish
    - Evo2OffTargetScorer (official stub): informative error without deps
    - Roadmap acceptance proxy: AUROC improvement >= 0.03 vs baseline on
      a synthetic CRISPR panel designed around seed-region conservation
"""

from __future__ import annotations

import numpy as np
import pytest

from pathway_subtyping.qc.offtarget_sequence import (
    Evo2OffTargetScorer,
    OffTargetScore,
    SimilarityBackend,
    SimulatedEvo2Backend,
    auroc,
    compare_backends,
)

# --------------------------------------------------------------------------- #
# AUROC primitive
# --------------------------------------------------------------------------- #


class TestAuroc:

    def test_perfect_separation(self):
        scores = [0.1, 0.2, 0.8, 0.9]
        labels = [0, 0, 1, 1]
        assert auroc(scores, labels) == pytest.approx(1.0)

    def test_random_scoring(self):
        rng = np.random.default_rng(0)
        scores = rng.uniform(size=1000).tolist()
        labels = rng.integers(0, 2, size=1000).tolist()
        # Should be close to 0.5 for random data
        assert 0.4 < auroc(scores, labels) < 0.6

    def test_constant_score_nan_or_half(self):
        scores = [0.5, 0.5, 0.5, 0.5]
        labels = [0, 0, 1, 1]
        # Tied scores -> U = 4 (n_pos * n_neg / 2), AUROC = 0.5
        assert auroc(scores, labels) == pytest.approx(0.5)

    def test_all_one_class_returns_nan(self):
        assert np.isnan(auroc([0.1, 0.2], [1, 1]))


# --------------------------------------------------------------------------- #
# Similarity backend
# --------------------------------------------------------------------------- #


class TestSimilarityBackend:

    def test_exact_match_scores_1(self):
        be = SimilarityBackend()
        rows = be.score("AAACCCGGG", {"site1": "AAACCCGGG"})
        assert rows[0].similarity_score == pytest.approx(1.0)

    def test_full_mismatch_scores_0(self):
        be = SimilarityBackend()
        rows = be.score("AAAAAAAAA", {"site1": "TTTTTTTTT"})
        assert rows[0].similarity_score == pytest.approx(0.0)

    def test_rejects_non_dna(self):
        be = SimilarityBackend()
        with pytest.raises(ValueError, match="non-DNA"):
            be.score("AAAXCCCGGG", {"s": "AAACCCGGG"})

    def test_combined_score_falls_back_to_similarity(self):
        be = SimilarityBackend()
        rows = be.score("AAACCCGGG", {"site1": "AAACCCGGG"})
        # No functional_score was set -> combined falls back to similarity
        assert rows[0].combined_score == pytest.approx(rows[0].similarity_score)


class TestSimulatedEvo2:

    def test_rejects_bad_seed_length(self):
        with pytest.raises(ValueError, match="seed_length"):
            SimulatedEvo2Backend(seed_length=0)

    def test_seed_matched_scores_higher_than_seed_mismatched(self):
        """Crucial F9 property: when a candidate preserves the seed
        region, the simulated functional score should surface it even
        if the overall similarity is modest."""
        be = SimulatedEvo2Backend(seed_length=8)
        # Guide and two candidates at same overall Hamming similarity,
        # one with seed match and one with seed mismatch.
        guide = "AAAAAAAACCCCCCCC"  # 16 nt, seed = 'AAAAAAAA'
        seed_match = "AAAAAAAATTTTTTTT"  # seed matches, tail mismatches
        seed_mismatch = "TTTTTTTTCCCCCCCC"  # seed mismatches, tail matches
        rows = be.score(guide, {"A": seed_match, "B": seed_mismatch})
        row_a = next(r for r in rows if r.site_id == "A")
        row_b = next(r for r in rows if r.site_id == "B")
        # Overall similarity equal by construction
        assert row_a.similarity_score == row_b.similarity_score
        # But functional score favours seed-match
        assert row_a.functional_score > row_b.functional_score


# --------------------------------------------------------------------------- #
# Official Evo 2 stub
# --------------------------------------------------------------------------- #


class TestEvo2Stub:

    def test_informative_error(self):
        scorer = Evo2OffTargetScorer()
        try:
            scorer._lazy_load()
        except (ImportError, NotImplementedError) as exc:
            msg = str(exc).lower()
            assert "qc-sequence" in msg or "evo2" in msg or "checkpoint" in msg
        else:
            pytest.skip("Evo 2 installed — skipping stub error test")


# --------------------------------------------------------------------------- #
# Roadmap acceptance: AUROC improvement >= 0.03
# --------------------------------------------------------------------------- #


def _synthetic_offtarget_panel(seed_len: int = 8, n_per_class: int = 25, rng_seed: int = 0):
    """Synthetic CRISPR panel where seed-region conservation is the
    true determinant of functional off-target activity."""
    rng = np.random.default_rng(rng_seed)
    guide_length = 20
    guide = "".join(rng.choice(list("ACGT"), size=guide_length))
    guide_seed = guide[:seed_len]
    guide_tail = guide[seed_len:]

    candidates = {}
    labels = {}

    # Positives: preserve seed, variable tail
    for i in range(n_per_class):
        tail_mutated = list(guide_tail)
        n_muts = rng.integers(0, len(tail_mutated) // 2 + 1)
        positions = rng.choice(len(tail_mutated), size=n_muts, replace=False)
        for p in positions:
            orig = tail_mutated[p]
            options = [b for b in "ACGT" if b != orig]
            tail_mutated[p] = rng.choice(options)
        seq = guide_seed + "".join(tail_mutated)
        site_id = f"POS_{i}"
        candidates[site_id] = seq
        labels[site_id] = 1

    # Negatives: mutate seed region (broken seed -> no off-target activity),
    # possibly preserve tail for similar Hamming distance.
    for i in range(n_per_class):
        seed_mut = list(guide_seed)
        n_muts = rng.integers(2, len(seed_mut) + 1)
        positions = rng.choice(len(seed_mut), size=n_muts, replace=False)
        for p in positions:
            orig = seed_mut[p]
            options = [b for b in "ACGT" if b != orig]
            seed_mut[p] = rng.choice(options)
        seq = "".join(seed_mut) + guide_tail
        site_id = f"NEG_{i}"
        candidates[site_id] = seq
        labels[site_id] = 0

    return guide, candidates, labels


class TestRoadmapAcceptance:

    def test_simulated_evo2_beats_similarity_by_3pct(self):
        guide, candidates, labels = _synthetic_offtarget_panel()
        res = compare_backends(
            baseline=SimilarityBackend(),
            contender=SimulatedEvo2Backend(seed_length=8),
            guide=guide,
            candidates=candidates,
            true_labels=labels,
        )
        assert res["delta"] >= 0.03, (
            f"contender AUROC {res['contender_auroc']:.3f} - "
            f"baseline AUROC {res['baseline_auroc']:.3f} = "
            f"delta {res['delta']:+.3f} (roadmap requires >= +0.03)"
        )

    def test_compare_backends_returns_both_aurocs(self):
        guide, candidates, labels = _synthetic_offtarget_panel(n_per_class=10)
        res = compare_backends(
            baseline=SimilarityBackend(),
            contender=SimulatedEvo2Backend(seed_length=8),
            guide=guide,
            candidates=candidates,
            true_labels=labels,
        )
        assert 0.0 <= res["baseline_auroc"] <= 1.0
        assert 0.0 <= res["contender_auroc"] <= 1.0
        assert res["n_sites"] == len(candidates)


class TestScoreToDict:

    def test_score_fields_present(self):
        be = SimulatedEvo2Backend()
        rows = be.score("AAAAAAAA", {"s": "AAAAAAAA"})
        d = rows[0].to_dict()
        assert {
            "guide",
            "site_id",
            "candidate_sequence",
            "similarity_score",
            "functional_score",
            "combined_score",
        }.issubset(d.keys())
