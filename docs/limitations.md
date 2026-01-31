# Limitations

## Overview

This document explicitly acknowledges the limitations of the Pathway Subtyping Framework. Understanding these constraints is essential for appropriate interpretation and use.

---

## Fundamental Limitations

### 1. Hypothesis-Generating Only

This framework is designed to **generate hypotheses**, not to:
- Establish causality
- Make clinical predictions
- Guide treatment decisions
- Replace experimental validation

All findings require independent validation through:
- Functional experiments
- Animal models
- Clinical studies
- Replication in diverse cohorts

### 2. Correlation, Not Causation

Pathway associations indicate **statistical enrichment**, not causal relationships. Possible explanations for associations include:
- True biological effect
- Confounding by ancestry, batch, or technical factors
- Chance (despite statistical correction)
- Pathway annotation artifacts

### 3. No Individual-Level Prediction

The framework identifies **group-level patterns**, not individual diagnoses or prognoses. It cannot:
- Predict who will develop disease
- Determine disease severity for an individual
- Guide personalized treatment
- Provide genetic counseling

---

## Biological Limitations

### 1. Incomplete Pathway Definitions

Curated pathways have significant gaps:
- **Study bias**: Well-funded diseases are over-represented
- **Gene bias**: Some genes are studied more than others
- **Update lag**: Databases do not capture recent discoveries
- **Arbitrary boundaries**: Pathway membership is often fuzzy

### 2. Static Models

The framework treats biology as static, ignoring:

| Dynamic Factor | What is Missed |
|----------------|----------------|
| Developmental timing | Effects may be stage-specific |
| Cell-type specificity | Different cells have different pathway activity |
| Compensatory mechanisms | Organisms adapt to perturbations |
| Pleiotropy | Genes have multiple functions |
| Epistasis | Gene-gene interactions |

### 3. Missing Biology

The framework cannot capture:
- **Non-coding variation**: Most genetic variation is non-coding
- **Epigenetics**: Methylation, chromatin state, etc.
- **Environmental factors**: Diet, toxins, social environment
- **Gene-environment interactions**: GxE effects
- **Somatic mosaicism**: Post-zygotic mutations

---

## Technical Limitations

### 1. Data Quality Dependence

Outputs are only as good as inputs:
- Sequencing errors propagate
- Variant calling is imperfect
- Annotation databases have errors
- Phenotype data may be inconsistent

### 2. Sequencing Technology Artifacts

Different technologies produce different results:
- WES vs. WGS coverage differences
- Capture kit biases
- Short-read vs. long-read
- Batch effects from processing pipelines

### 3. Population Structure

Genetic analyses are confounded by:
- Ancestry differences between cases and controls
- Population stratification
- Cryptic relatedness
- Founder effects

### 4. Sample Size Constraints

| Analysis Type | Typical Minimum | Limitation |
|---------------|-----------------|------------|
| Pathway enrichment | 100-300 samples | Low power for weak effects |
| Subtype discovery | 300-1000 samples | Unstable clusters if too small |
| Rare variant association | 1000+ samples | Most variants too rare |
| Cross-cohort replication | Multiple cohorts | Availability |

### 5. Multiple Testing

With many pathways and genes:
- Many false positives expected
- Stringent correction reduces power
- Winner's curse inflates effect sizes

---

## Analytical Limitations

### 1. Parameter Sensitivity

Results depend on arbitrary choices:
- Variant weighting schemes
- Pathway databases used
- Clustering algorithm and K
- Significance thresholds
- Normalization methods

### 2. Clustering Limitations

Unsupervised clustering has inherent issues:
- **Cluster number**: Often unclear how many subtypes exist
- **Boundary artifacts**: Continuous variation forced into discrete groups
- **Instability**: Small changes can shift assignments
- **Interpretation**: Clusters may not correspond to biology

### 3. Overfitting Risk

Complex models can:
- Fit noise rather than signal
- Fail to generalize to new data
- Appear successful on training data only

---

## Interpretational Limitations

### 1. Biological Plausibility is Not Truth

A pathway association that "makes biological sense" is not necessarily true. Confirmation bias leads us to find stories for any result.

### 2. Replication Challenges

True biological effects may not replicate due to:
- Different populations
- Different sequencing technologies
- Different phenotype definitions
- Insufficient power

Failed replication does not prove null (may be underpowered).
Successful replication does not prove causation.

### 3. Publication Bias

Published literature is biased toward:
- Positive results
- Novel findings
- Large effect sizes

This affects both our prior expectations and pathway annotations.

---

## Ethical and Social Limitations

### 1. Representation Bias

Most genetic studies over-represent:
- European ancestry individuals
- Families with resources to participate
- Certain phenotypic presentations

Findings may not generalize to underrepresented groups.

### 2. Reductionism Risk

Focusing on genetics may:
- Underemphasize environmental factors
- Reinforce medical model of disease
- Miss protective factors
- Ignore patient perspectives

### 3. Dual Use Concerns

Genetic subtyping could potentially be misused for:
- Discrimination
- Prenatal selection
- Insurance decisions
- Educational gatekeeping

---

## What This Framework Cannot Do

| Cannot Do | Why |
|-----------|-----|
| Diagnose disease | Not validated for clinical use |
| Predict individual outcomes | Group-level patterns only |
| Determine causality | Observational data only |
| Replace genetic counseling | Requires professional interpretation |
| Guide treatment | No clinical validation |
| Capture all biology | Static, incomplete models |

---

## Appropriate Use

Given these limitations, the framework is appropriate for:

:white_check_mark: Generating research hypotheses
:white_check_mark: Prioritizing pathways for experimental follow-up
:white_check_mark: Informing study design
:white_check_mark: Educational purposes
:white_check_mark: Exploratory data analysis

And inappropriate for:

:x: Clinical diagnosis
:x: Treatment decisions
:x: Individual predictions
:x: Prenatal screening
:x: Insurance/employment decisions

---

## Mitigation Strategies

| Limitation | Mitigation |
|------------|------------|
| Incomplete pathways | Use multiple pathway databases; acknowledge gaps |
| Static models | Integrate developmental data when available |
| Overfitting | Cross-validation; replication; validation gates |
| Multiple testing | FDR correction; pre-registration |
| Confounding | Careful covariate adjustment; negative controls |
| Interpretation | Require experimental validation; temper claims |

---

## Conclusion

Acknowledging limitations is not weakness—it is scientific integrity. This framework offers a principled approach to studying genetic heterogeneity, but its outputs are starting points for further research, not conclusions.

Users should:
1. Read and understand these limitations
2. Communicate uncertainty in all findings
3. Seek independent validation
4. Avoid overinterpretation
5. Consider ethical implications

---

## See Also

- [DISCLAIMER.md](../DISCLAIMER.md) - Full research-only disclaimer
- [Outputs Dictionary](outputs_dictionary.md) - Interpretation guidance
- [Validation Gates Guide](guides/validation-gates.md) - Statistical validation
