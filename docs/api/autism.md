# Autism-Specific Interpretation Layer API

> **Module**: `pathway_subtyping.autism`
> **Install**: `pip install pathway-subtyping[autism]`
> **Scope**: AUTISM-ONLY — refuses to run on non-autism datasets

Optional interpretation layer providing autism-specific biological rules, neuro-symbolic reasoning, and therapeutic hypothesis generation on top of the disease-agnostic PSF pipeline. Implements Mohit's standing rule: tool is disease-agnostic, biological interpretation is autism-only.

---

## Quick Example

```python
from pathway_subtyping.autism import (
    BiologicalRules, AutismRuleEngine, BiologicalContext,
    IndividualData, ExplanationGenerator,
)

# Set up biological context (reference databases)
context = BiologicalContext(
    gene_constraints={"CHD8": 0.99, "SHANK3": 0.95},
    sfari_genes={"CHD8": 1, "SHANK3": 1},
    prenatal_expression={"CHD8": 3.5, "SHANK3": 2.0},
    syngo_genes={"SHANK3", "NRXN1"},
    chd8_targets={"KDM6A", "ARID1B", "TBR1"},
)

# Create engine (autism-only enforced)
engine = AutismRuleEngine(
    rules=BiologicalRules.get_all_rules(),
    context=context,
    disease="autism",  # "schizophrenia" would raise ValueError
)

# Evaluate an individual
individual = IndividualData(
    sample_id="ASD_001",
    variants=[variant_with_gene_CHD8],
    pathway_scores={"CHROMATIN": -2.5},
)
fired = engine.evaluate(individual)

# Generate explanation
gen = ExplanationGenerator(include_disclaimers=True)
chain = gen.generate_reasoning_chain("ASD_001", fired)
print(chain.explanation_text)
```

---

## Biological Rules (#26)

### `BiologicalRules`

Factory class providing 8 curated autism-specific rules.

| Rule | Name | Confidence | Evidence |
|------|------|-----------|----------|
| R1 | Constrained LoF in Developing Cortex | 0.90 | Samocha 2014, Willsey 2013 |
| R2 | Pathway Convergence | 0.85 | Voineagu 2011, Sanders 2015 |
| R3 | CHD8 Chromatin Cascade | 0.92 | Cotney 2015, Sugathan 2014 |
| R3b | CHD8 Target Disruption | 0.75 | Cotney 2015 |
| R4 | Synaptic Excitatory Disruption | 0.82 | SynGO, Rubenstein 2003 |
| R5 | Paralog Compensation | 0.70 | Diss & Bhatt 2020 |
| R6 | Therapeutic Pathway Targeting | 0.60 | — |
| R7 | SFARI High-Confidence Gene | 0.88 | SFARI Gene Database |

```python
rules = BiologicalRules.get_all_rules()         # All 8 rules
rules = BiologicalRules.get_core_rules()         # R1-R6 (without R7)
rule = BiologicalRules.get_rule_by_id("R1")      # Specific rule
rules = BiologicalRules.get_rules_by_conclusion_type(ConclusionType.PATHWAY_DISRUPTION)
```

**Conclusion types:** `PATHWAY_DISRUPTION` | `PATHWAY_CONVERGENCE` | `SUBTYPE_INDICATOR` | `EFFECT_MODIFIER` | `THERAPEUTIC_HYPOTHESIS`

### Condition Predicates (18 types)

| Category | Predicates |
|----------|-----------|
| Variant | `has_variant`, `has_lof_variant`, `has_missense_variant`, `has_damaging_variant`, `has_multiple_hits` |
| Gene | `is_constrained`, `is_sfari_gene`, `is_high_confidence_sfari`, `has_paralog`, `is_chd8_target`, `is_synaptic_gene` |
| Expression | `prenatally_expressed`, `cell_type_enriched` |
| Pathway | `gene_in_pathway`, `pathway_disrupted`, `hits_are_independent` |
| Drug | `drug_targets`, `mechanism_alignment` |

---

## Rule Engine

### `AutismRuleEngine`

Evaluates rules against individual variant/expression data. **Autism-only enforcement** — raises `ValueError` for non-autism disease labels.

```python
engine = AutismRuleEngine(
    rules: List[Rule],
    context: BiologicalContext,
    disease: str = "autism",  # Also accepts "asd", "autism_spectrum_disorder"
)

fired = engine.evaluate(individual_data, rule_ids=["R1", "R3"])  -> List[FiredRule]
batch = engine.evaluate_batch(cohort)  -> Dict[str, List[FiredRule]]
summary = engine.get_summary(fired)    -> Dict with n_fired, genes_affected, etc.
```

### `BiologicalContext`

Reference biological data for condition evaluation. All fields are optional — rules that need missing data simply don't fire.

```python
context = BiologicalContext(
    gene_constraints: Dict[str, float],          # gene -> pLI score
    sfari_genes: Dict[str, int],                 # gene -> SFARI score (1-3)
    prenatal_expression: Dict[str, float],       # gene -> prenatal expression level
    cell_type_expression: Dict[str, Dict[str, float]],  # gene -> {cell_type: level}
    pathway_genes: Dict[str, Set[str]],          # pathway -> gene set
    paralog_map: Dict[str, List[str]],           # gene -> paralog list
    chd8_targets: Set[str],                      # Genes regulated by CHD8
    syngo_genes: Set[str],                       # SynGO synaptic genes
    drug_targets: Dict[str, Set[str]],           # drug -> target genes
)
```

### `IndividualData`

Variant and expression data for one individual.

```python
individual = IndividualData(
    sample_id: str,
    variants: List[Any],                    # Objects with .gene, .is_lof, .is_damaging
    gene_burdens: Dict[str, float],         # Per-gene burden scores
    pathway_scores: Dict[str, float],       # Per-pathway z-scores
)
individual.get_genes_with_variants()         # Set of genes with any variant
individual.get_genes_with_damaging_variants()  # Set with LoF/damaging
```

---

## Explanation Generator

### `ExplanationGenerator`

Produces human-readable reasoning chains with mandatory disclaimers.

```python
gen = ExplanationGenerator(
    include_evidence=True,
    include_disclaimers=True,  # DISCLAIMER always included
)

chain = gen.generate_reasoning_chain("ASD_001", fired_rules) -> ReasoningChain
chain.n_rules_fired
chain.average_confidence
chain.genes_affected           # Set of affected genes
chain.pathway_conclusions      # Dict[pathway, confidence]
chain.subtype_indicators       # e.g., ["chromatin_remodeling", "e_i_imbalance"]
chain.explanation_text         # Full human-readable text

summary = gen.generate_clinical_summary(chain)  # Brief clinical-style summary
```

---

## Neurosymbolic Combiner (#26)

### `NeurosymbolicCombiner`

Merges GNN predictions with symbolic rule conclusions. Pure Python — no PyTorch dependency.

```python
from pathway_subtyping.autism.neurosymbolic.combiner import (
    NeurosymbolicCombiner, CombinerConfig, CombinationMethod,
)

combiner = NeurosymbolicCombiner(CombinerConfig(
    method=CombinationMethod.WEIGHTED_SUM,
    neural_weight=0.6,
    symbolic_weight=0.4,
))

result = combiner.combine(
    neural_scores={"CHD8": 0.9, "SHANK3": 0.7},
    symbolic_scores={"CHD8": 0.92, "NRXN1": 0.85},
)
result.combined_scores  # Merged gene -> score mapping
```

| Method | Description |
|--------|-------------|
| `WEIGHTED_SUM` | `w_n * neural + w_s * symbolic` (default) |
| `MAX` | `max(neural, symbolic)` per gene |
| `PRODUCT` | `neural * symbolic` (both must agree) |
| `RULE_GUIDED` | Boost neural where rules fire; penalize where absent |

---

## Therapeutic Hypotheses (#27)

### Evidence Scoring

Multi-criteria scoring with pediatric safety flags.

```python
from pathway_subtyping.autism.therapeutic.evidence import (
    EvidenceScorer, EvidenceLevel, SafetyFlag,
)

scorer = EvidenceScorer(
    literature_db={"MEMANTINE": 0.6},
    safety_db={"DRUG_X": [SafetyFlag.CONTRAINDICATED_PEDIATRIC.value]},
)

score = scorer.score(
    drug_id="MEMANTINE", target_pathway="GLUTAMATE",
    target_genes=["GRIN1"], mechanism="antagonist",
)
score.overall       # 0-1 weighted composite
score.level         # EvidenceLevel.HIGH / MODERATE / LOW / INSUFFICIENT
score.safety_flags  # List of SafetyFlag values
score.has_critical_safety_flags  # True for pediatric contraindications
```

**Safety flags (11 types):** `CNS_EFFECTS`, `DEVELOPMENTAL_CONCERNS`, `DRUG_INTERACTIONS`, `CONTRAINDICATED_PEDIATRIC`, `BLACK_BOX_WARNING`, `OFF_LABEL_USE`, `IMMUNOSUPPRESSION`, `HEPATOTOXICITY`, `CARDIOTOXICITY`, `TERATOGENIC`, `WITHDRAWAL_RISK`

### Drug Mapping

```python
from pathway_subtyping.autism.therapeutic.drug_mapping import (
    DrugTargetDatabase, DrugCandidate, PathwayDrugMapper,
    DrugMechanism, DrugStatus,
)

db = DrugTargetDatabase()
db.add_drug(DrugCandidate(
    drug_id="MEMANTINE", drug_name="Memantine",
    target_genes=["GRIN1", "GRIN2A"],
    mechanism="antagonist", mechanism_type=DrugMechanism.ANTAGONIST,
    status=DrugStatus.APPROVED, asd_relevance_score=0.7,
))

mapper = PathwayDrugMapper(db)
candidates = mapper.map_pathway("GLUTAMATE", pathway_genes=["GRIN1", "GRIN2A"])
```

### Hypothesis Ranking

```python
from pathway_subtyping.autism.therapeutic.ranking import HypothesisRanker, RankingConfig

ranker = HypothesisRanker(
    drug_mapper=mapper, evidence_scorer=scorer,
    config=RankingConfig(max_hypotheses=20, max_per_pathway=5),
)
result = ranker.rank(
    pathway_scores={"GLUTAMATE": 2.5, "GABA": -1.8},
    pathway_genes={"GLUTAMATE": ["GRIN1", "GRIN2A"]},
)

for h in result.hypotheses:
    print(h.summary())
    assert h.requires_validation  # ALWAYS True — cannot be overridden
```

**`TherapeuticHypothesis.requires_validation` is always `True`** — enforced in `__post_init__`. These are computational hypotheses, NOT clinical recommendations.

Ranking uses diversity constraints (`max_per_pathway`, `max_per_mechanism`) to ensure coverage across pathways.

---

**DISCLAIMER:** This layer generates computational research hypotheses only. All findings require independent validation by qualified researchers and clinicians. Not for clinical diagnosis or treatment decisions.
