"""
Tests for Autism-Specific Interpretation Layer (#26, #27).

Tests cover:
    - Biological rules (R1-R7)
    - Condition evaluation
    - Rule engine (autism-only enforcement)
    - Explanation generation
    - Neurosymbolic combiner
    - Evidence scoring with safety flags
    - Drug mapping
    - Therapeutic hypothesis ranking
"""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import pytest

from pathway_subtyping.autism.rules.biological_rules import (
    BiologicalRules,
    Conclusion,
    ConclusionType,
    Rule,
)
from pathway_subtyping.autism.rules.conditions import Condition, ConditionEvaluator
from pathway_subtyping.autism.rules.engine import (
    AutismRuleEngine,
    BiologicalContext,
    FiredRule,
    IndividualData,
)
from pathway_subtyping.autism.rules.explanation import ExplanationGenerator, ReasoningChain


# ── Mock variant for testing ──────────────────────────────────────────

@dataclass
class MockVariant:
    gene: str = ""
    consequence: str = ""
    is_lof: bool = False
    is_damaging: bool = False


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def context():
    return BiologicalContext(
        gene_constraints={"CHD8": 0.99, "SHANK3": 0.95, "NRXN1": 0.92, "SCN2A": 0.88},
        sfari_genes={"CHD8": 1, "SHANK3": 1, "NRXN1": 1, "SCN2A": 1, "DYRK1A": 2},
        prenatal_expression={"CHD8": 3.5, "SHANK3": 2.0, "NRXN1": 2.5},
        cell_type_expression={
            "SHANK3": {"excitatory_neuron": 4.0, "inhibitory_neuron": 1.0},
            "NRXN1": {"excitatory_neuron": 3.0},
        },
        pathway_genes={
            "CHROMATIN": {"CHD8", "KDM6A", "ARID1B"},
            "SYNAPTIC": {"SHANK3", "NRXN1", "NLGN3"},
        },
        paralog_map={"SHANK3": ["SHANK1", "SHANK2"]},
        chd8_targets={"KDM6A", "ARID1B", "TBR1"},
        syngo_genes={"SHANK3", "NRXN1", "NLGN3", "SYN1"},
    )


@pytest.fixture
def individual_chd8():
    return IndividualData(
        sample_id="ASD_001",
        variants=[
            MockVariant(gene="CHD8", consequence="frameshift", is_lof=True, is_damaging=True),
        ],
        pathway_scores={"CHROMATIN": -2.5, "SYNAPTIC": 0.3},
    )


@pytest.fixture
def individual_shank3():
    return IndividualData(
        sample_id="ASD_002",
        variants=[
            MockVariant(gene="SHANK3", consequence="stop_gained", is_lof=True, is_damaging=True),
        ],
        pathway_scores={"SYNAPTIC": -3.0},
    )


# ══════════════════════════════════════════════════════════════════════
# BIOLOGICAL RULES
# ══════════════════════════════════════════════════════════════════════

class TestBiologicalRules:

    def test_get_all_rules(self):
        rules = BiologicalRules.get_all_rules()
        assert len(rules) == 8
        ids = {r.id for r in rules}
        assert ids == {"R1", "R2", "R3", "R3b", "R4", "R5", "R6", "R7"}

    def test_get_core_rules(self):
        rules = BiologicalRules.get_core_rules()
        assert len(rules) == 7
        assert "R7" not in {r.id for r in rules}

    def test_get_rule_by_id(self):
        r1 = BiologicalRules.get_rule_by_id("R1")
        assert r1 is not None
        assert r1.base_confidence == 0.90

    def test_nonexistent_rule(self):
        assert BiologicalRules.get_rule_by_id("R99") is None

    def test_get_rules_by_conclusion_type(self):
        disruption = BiologicalRules.get_rules_by_conclusion_type(
            ConclusionType.PATHWAY_DISRUPTION
        )
        assert len(disruption) >= 2

    def test_rule_to_dict(self):
        r1 = BiologicalRules.get_rule_by_id("R1")
        d = r1.to_dict()
        assert d["id"] == "R1"
        assert "conditions" in d
        assert "conclusion" in d

    def test_rule_priorities(self):
        rules = BiologicalRules.get_all_rules()
        # R1 should be highest priority (10)
        r1 = next(r for r in rules if r.id == "R1")
        assert r1.priority == 10


# ══════════════════════════════════════════════════════════════════════
# RULE ENGINE
# ══════════════════════════════════════════════════════════════════════

class TestAutismRuleEngine:

    def test_autism_only_enforcement(self, context):
        rules = BiologicalRules.get_all_rules()
        with pytest.raises(ValueError, match="autism-only"):
            AutismRuleEngine(rules, context, disease="schizophrenia")

    def test_valid_disease_names(self, context):
        rules = BiologicalRules.get_all_rules()
        for disease in ("autism", "asd", "autism_spectrum_disorder"):
            engine = AutismRuleEngine(rules, context, disease=disease)
            assert engine is not None

    def test_evaluate_chd8(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        assert len(fired) > 0
        rule_ids = {fr.rule.id for fr in fired}
        # CHD8 LoF should fire R1 (constrained + prenatal) and R7 (SFARI high-conf)
        assert "R1" in rule_ids or "R7" in rule_ids

    def test_evaluate_shank3(self, context, individual_shank3):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_shank3)
        rule_ids = {fr.rule.id for fr in fired}
        # SHANK3 is synaptic + constrained + SFARI 1 + has paralogs
        assert len(fired) > 0

    def test_evaluate_batch(self, context, individual_chd8, individual_shank3):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        results = engine.evaluate_batch([individual_chd8, individual_shank3])
        assert "ASD_001" in results
        assert "ASD_002" in results

    def test_evaluate_subset_rules(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8, rule_ids=["R1"])
        for fr in fired:
            assert fr.rule.id == "R1"

    def test_fired_rule_to_dict(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        if fired:
            d = fired[0].to_dict()
            assert "rule_id" in d
            assert "confidence" in d

    def test_get_summary(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        summary = engine.get_summary(fired)
        assert "n_fired" in summary
        assert "genes_affected" in summary

    def test_no_variants_no_rules(self, context):
        empty = IndividualData(sample_id="CTRL_001")
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(empty)
        assert len(fired) == 0


# ══════════════════════════════════════════════════════════════════════
# EXPLANATION GENERATOR
# ══════════════════════════════════════════════════════════════════════

class TestExplanationGenerator:

    def test_generate_chain(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        gen = ExplanationGenerator()
        chain = gen.generate_reasoning_chain("ASD_001", fired)
        assert isinstance(chain, ReasoningChain)
        assert chain.individual_id == "ASD_001"
        assert chain.n_rules_fired == len(fired)

    def test_chain_to_dict(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        gen = ExplanationGenerator()
        chain = gen.generate_reasoning_chain("ASD_001", fired)
        d = chain.to_dict()
        assert "individual_id" in d
        assert "steps" in d

    def test_clinical_summary(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        gen = ExplanationGenerator()
        chain = gen.generate_reasoning_chain("ASD_001", fired)
        summary = gen.generate_clinical_summary(chain)
        assert "ASD_001" in summary
        assert "Research hypothesis" in summary or "not a clinical" in summary

    def test_disclaimer_included(self, context, individual_chd8):
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        gen = ExplanationGenerator(include_disclaimers=True)
        chain = gen.generate_reasoning_chain("ASD_001", fired)
        assert "DISCLAIMER" in chain.explanation_text

    def test_empty_fired_rules(self):
        gen = ExplanationGenerator()
        chain = gen.generate_reasoning_chain("CTRL_001", [])
        assert chain.n_rules_fired == 0


# ══════════════════════════════════════════════════════════════════════
# NEUROSYMBOLIC COMBINER
# ══════════════════════════════════════════════════════════════════════

class TestNeurosymbolicCombiner:

    def test_weighted_sum(self):
        from pathway_subtyping.autism.neurosymbolic.combiner import (
            CombinerConfig,
            NeurosymbolicCombiner,
        )
        combiner = NeurosymbolicCombiner()
        result = combiner.combine(
            neural_scores={"CHD8": 0.9, "SHANK3": 0.7},
            symbolic_scores={"CHD8": 0.92, "NRXN1": 0.85},
        )
        assert "CHD8" in result.combined_scores
        assert "SHANK3" in result.combined_scores
        assert "NRXN1" in result.combined_scores

    def test_max_combination(self):
        from pathway_subtyping.autism.neurosymbolic.combiner import (
            CombinationMethod,
            CombinerConfig,
            NeurosymbolicCombiner,
        )
        config = CombinerConfig(method=CombinationMethod.MAX)
        combiner = NeurosymbolicCombiner(config)
        result = combiner.combine(
            neural_scores={"A": 0.3},
            symbolic_scores={"A": 0.8},
        )
        # Max of 0.3, 0.8 = 0.8 (normalized to 1.0 since only gene)
        assert result.combined_scores["A"] > 0

    def test_create_symbolic_score_vector(self, context, individual_chd8):
        from pathway_subtyping.autism.neurosymbolic.combiner import (
            create_symbolic_score_vector,
        )
        engine = AutismRuleEngine(BiologicalRules.get_all_rules(), context)
        fired = engine.evaluate(individual_chd8)
        genes = ["CHD8", "SHANK3", "NRXN1"]
        vec = create_symbolic_score_vector(fired, genes)
        assert len(vec) == 3

    def test_to_dict(self):
        from pathway_subtyping.autism.neurosymbolic.combiner import NeurosymbolicCombiner
        combiner = NeurosymbolicCombiner()
        result = combiner.combine({"A": 0.5}, {"B": 0.5})
        d = result.to_dict()
        assert "method" in d
        assert "top_genes" in d


# ══════════════════════════════════════════════════════════════════════
# EVIDENCE SCORING (#27)
# ══════════════════════════════════════════════════════════════════════

class TestEvidenceScorer:

    def test_basic_scoring(self):
        from pathway_subtyping.autism.therapeutic.evidence import EvidenceScorer
        scorer = EvidenceScorer()
        result = scorer.score(
            drug_id="SORAFENIB",
            target_pathway="MAPK",
            target_genes=["RAF"],
            mechanism="inhibitor",
        )
        assert result.overall >= 0
        assert result.level is not None

    def test_safety_flags(self):
        from pathway_subtyping.autism.therapeutic.evidence import EvidenceScorer, SafetyFlag
        scorer = EvidenceScorer(
            safety_db={
                "DRUG_X": [SafetyFlag.CONTRAINDICATED_PEDIATRIC.value, SafetyFlag.TERATOGENIC.value],
            }
        )
        result = scorer.score(
            drug_id="DRUG_X", target_pathway="PW", target_genes=["G"],
            mechanism="inhibitor",
        )
        assert result.has_critical_safety_flags
        assert len(result.safety_flags) >= 2

    def test_literature_boost(self):
        from pathway_subtyping.autism.therapeutic.evidence import EvidenceScorer
        scorer = EvidenceScorer(literature_db={"DRUG_A": 0.8})
        result = scorer.score(
            drug_id="DRUG_A", target_pathway="PW",
            target_genes=["G"], mechanism="inhibitor",
        )
        assert result.literature_support == 0.8

    def test_to_dict(self):
        from pathway_subtyping.autism.therapeutic.evidence import EvidenceScorer
        scorer = EvidenceScorer()
        result = scorer.score(
            drug_id="D", target_pathway="P", target_genes=["G"], mechanism="mod",
        )
        d = result.to_dict()
        assert "overall" in d
        assert "safety_flags" in d
        assert "has_critical_safety" in d


# ══════════════════════════════════════════════════════════════════════
# DRUG MAPPING (#27)
# ══════════════════════════════════════════════════════════════════════

class TestDrugMapping:

    @pytest.fixture
    def drug_db(self):
        from pathway_subtyping.autism.therapeutic.drug_mapping import (
            DrugCandidate,
            DrugMechanism,
            DrugStatus,
            DrugTargetDatabase,
        )
        db = DrugTargetDatabase()
        db.add_drug(DrugCandidate(
            drug_id="SORAFENIB", drug_name="Sorafenib",
            target_genes=["RAF", "VEGFR"],
            mechanism="inhibitor", mechanism_type=DrugMechanism.INHIBITOR,
            status=DrugStatus.APPROVED, asd_relevance_score=0.3,
        ))
        db.add_drug(DrugCandidate(
            drug_id="MEMANTINE", drug_name="Memantine",
            target_genes=["GRIN1", "GRIN2A"],
            mechanism="antagonist", mechanism_type=DrugMechanism.ANTAGONIST,
            status=DrugStatus.APPROVED, asd_relevance_score=0.7,
        ))
        return db

    def test_add_and_get_drug(self, drug_db):
        drug = drug_db.get_drug("SORAFENIB")
        assert drug is not None
        assert drug.drug_name == "Sorafenib"

    def test_get_drugs_for_gene(self, drug_db):
        drugs = drug_db.get_drugs_for_gene("RAF")
        assert len(drugs) == 1
        assert drugs[0].drug_id == "SORAFENIB"

    def test_search_by_status(self, drug_db):
        from pathway_subtyping.autism.therapeutic.drug_mapping import DrugStatus
        approved = drug_db.search_drugs(status=DrugStatus.APPROVED)
        assert len(approved) == 2

    def test_pathway_mapper(self, drug_db):
        from pathway_subtyping.autism.therapeutic.drug_mapping import PathwayDrugMapper
        mapper = PathwayDrugMapper(drug_db)
        candidates = mapper.map_pathway("MAPK", pathway_genes=["RAF", "MEK"])
        assert len(candidates) >= 1

    def test_statistics(self, drug_db):
        stats = drug_db.get_statistics()
        assert stats["n_drugs"] == 2


# ══════════════════════════════════════════════════════════════════════
# THERAPEUTIC RANKING (#27)
# ══════════════════════════════════════════════════════════════════════

class TestTherapeuticRanking:

    @pytest.fixture
    def ranker(self):
        from pathway_subtyping.autism.therapeutic.drug_mapping import (
            DrugCandidate,
            DrugMechanism,
            DrugStatus,
            DrugTargetDatabase,
            PathwayDrugMapper,
        )
        from pathway_subtyping.autism.therapeutic.evidence import EvidenceScorer
        from pathway_subtyping.autism.therapeutic.ranking import HypothesisRanker

        db = DrugTargetDatabase()
        db.add_drug(DrugCandidate(
            drug_id="DRUG_A", drug_name="Drug A",
            target_genes=["GENE_1", "GENE_2"],
            mechanism="inhibitor", mechanism_type=DrugMechanism.INHIBITOR,
            status=DrugStatus.APPROVED, asd_relevance_score=0.6,
        ))
        db.add_drug(DrugCandidate(
            drug_id="DRUG_B", drug_name="Drug B",
            target_genes=["GENE_3"],
            mechanism="agonist", mechanism_type=DrugMechanism.AGONIST,
            status=DrugStatus.APPROVED, asd_relevance_score=0.4,
        ))

        mapper = PathwayDrugMapper(db)
        scorer = EvidenceScorer()
        return HypothesisRanker(drug_mapper=mapper, evidence_scorer=scorer)

    def test_rank_produces_hypotheses(self, ranker):
        result = ranker.rank(
            pathway_scores={"PW1": 2.5, "PW2": -1.8},
            pathway_genes={"PW1": ["GENE_1", "GENE_2"], "PW2": ["GENE_3"]},
        )
        assert len(result.hypotheses) > 0

    def test_requires_validation_enforced(self, ranker):
        result = ranker.rank(
            pathway_scores={"PW1": 3.0},
            pathway_genes={"PW1": ["GENE_1"]},
        )
        for h in result.hypotheses:
            assert h.requires_validation is True

    def test_ranking_order(self, ranker):
        result = ranker.rank(
            pathway_scores={"PW1": 3.0},
            pathway_genes={"PW1": ["GENE_1", "GENE_2"]},
        )
        if len(result.hypotheses) >= 2:
            assert result.hypotheses[0].score >= result.hypotheses[1].score

    def test_diversity_caps(self, ranker):
        from pathway_subtyping.autism.therapeutic.ranking import RankingConfig
        ranker.config = RankingConfig(max_per_pathway=1)
        result = ranker.rank(
            pathway_scores={"PW1": 3.0},
            pathway_genes={"PW1": ["GENE_1", "GENE_2", "GENE_3"]},
        )
        pw1_count = sum(1 for h in result.hypotheses if h.target_pathway == "PW1")
        assert pw1_count <= 1

    def test_result_summary(self, ranker):
        result = ranker.rank(
            pathway_scores={"PW1": 2.5},
            pathway_genes={"PW1": ["GENE_1"]},
        )
        summary = result.summary()
        assert "hypotheses" in summary.lower()
        assert "not clinical" in summary.lower() or "Computational" in summary

    def test_to_dict(self, ranker):
        result = ranker.rank(
            pathway_scores={"PW1": 2.5},
            pathway_genes={"PW1": ["GENE_1"]},
        )
        d = result.to_dict()
        assert "n_hypotheses" in d
        assert "hypotheses" in d

    def test_below_zscore_threshold(self, ranker):
        result = ranker.rank(
            pathway_scores={"PW1": 0.5},  # Below 1.5 threshold
            pathway_genes={"PW1": ["GENE_1"]},
        )
        assert len(result.hypotheses) == 0
