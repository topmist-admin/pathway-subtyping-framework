# Multi-Omic Fusion API

> **Module**: `pathway_subtyping.multi_omic`

Fuses pathway scores from multiple data modalities (VCF variant burden, bulk RNA-seq, single-cell scRNA-seq, deconvolution) into a combined feature matrix for unified subtype discovery.

---

## Quick Example

```python
from pathway_subtyping import (
    ModalityType, FusionStrategy, MissingStrategy,
    prepare_modality, fuse_modalities,
)

# Prepare modalities (each has pathway_scores DataFrame)
vcf_mod = prepare_modality(ModalityType.VCF, vcf_scores, label="WES")
expr_mod = prepare_modality(ModalityType.EXPRESSION, expr_scores, label="RNA-seq")

# Fuse
result = fuse_modalities(
    [vcf_mod, expr_mod],
    strategy=FusionStrategy.CONCATENATE,
    missing_strategy=MissingStrategy.IMPUTE_ZERO,
    seed=42,
)

# result.fused_pathway_scores: unified samples x pathways DataFrame
# result.per_modality_scores: dict for cross-modal validation
print(result.format_report())
```

---

## Enums

### `ModalityType`

| Value | Description |
|-------|-------------|
| `VCF` | Variant burden pathway scores |
| `EXPRESSION` | Bulk RNA-seq expression pathway scores |
| `SINGLE_CELL` | Single-cell scRNA-seq pathway scores |
| `DECONVOLUTION` | Deconvolution-derived cell-type proportions |

### `FusionStrategy`

| Value | Description | Output Dimensions |
|-------|-------------|-------------------|
| `CONCATENATE` | Column-bind with modality prefixes (e.g., `WES:PATHWAY_1`) | All samples x all pathways |
| `WEIGHTED_AVERAGE` | Weighted mean of shared pathway scores | Shared samples x shared pathways |
| `INTERSECTION_ONLY` | Restrict to shared samples and pathways only | Shared samples x shared pathways |

### `MissingStrategy`

| Value | Description |
|-------|-------------|
| `IMPUTE_ZERO` | Fill missing values with 0 (default) |
| `IMPUTE_MEAN` | Fill with column mean |
| `DROP` | Drop samples missing from any modality |

---

## Functions

### `prepare_modality()`

```python
def prepare_modality(
    modality_type: ModalityType,
    pathway_scores: pd.DataFrame,
    gene_data: Optional[pd.DataFrame] = None,
    quality_report: Optional[Any] = None,
    label: Optional[str] = None,
) -> ModalityInput
```

Wrap a modality's pathway scores into a `ModalityInput` for fusion.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `modality_type` | ModalityType | required | Which omic type this represents |
| `pathway_scores` | pd.DataFrame | required | Z-normalized scores (samples x pathways) |
| `gene_data` | pd.DataFrame or None | None | Gene-level data for characterization |
| `quality_report` | Any or None | None | Quality report from the scoring step |
| `label` | str or None | None | Human-readable label (e.g., `"WES"`, `"RNA-seq"`) |

**Returns:** `ModalityInput` ready for fusion.

**Raises:** `ValueError` if data is empty or has fewer than 2 samples.

---

### `fuse_modalities()`

```python
def fuse_modalities(
    modalities: List[ModalityInput],
    strategy: FusionStrategy = FusionStrategy.CONCATENATE,
    missing_strategy: MissingStrategy = MissingStrategy.IMPUTE_ZERO,
    weights: Optional[Dict[str, float]] = None,
    renormalize: bool = True,
    seed: Optional[int] = None,
) -> MultiOmicFusionResult
```

Fuse multiple modalities into a single feature matrix.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `modalities` | List[ModalityInput] | required | Two or more prepared modalities |
| `strategy` | FusionStrategy | `CONCATENATE` | Fusion strategy |
| `missing_strategy` | MissingStrategy | `IMPUTE_ZERO` | How to handle missing samples |
| `weights` | Dict[str, float] or None | None | Per-modality weights (for `WEIGHTED_AVERAGE`; must sum to 1.0) |
| `renormalize` | bool | True | Z-normalize output columns |
| `seed` | int or None | None | Random seed |

**Returns:** `MultiOmicFusionResult` with fused scores and quality report.

**Raises:** `ValueError` if fewer than 2 modalities, no shared samples, or invalid weights.

---

### `correlation_analysis()`

```python
def correlation_analysis(
    modalities: List[ModalityInput],
    method: str = "pearson",
) -> List[Dict[str, Any]]
```

Compute pairwise pathway-level correlations between modalities.

**Returns:** List of dicts with correlation stats per modality pair.

---

### `compute_sample_overlap()`

```python
def compute_sample_overlap(
    modalities: List[ModalityInput],
) -> SampleOverlapStats
```

Compute sample overlap statistics across modalities.

---

## Dataclasses

### `MultiOmicFusionResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `fused_pathway_scores` | pd.DataFrame | Combined feature matrix (samples x pathways) |
| `fused_gene_data` | pd.DataFrame or None | Combined gene-level data |
| `per_modality_scores` | Dict[str, pd.DataFrame] | Per-modality scores for cross-modal validation |
| `strategy` | FusionStrategy | Strategy used |
| `quality_report` | MultiOmicQualityReport | Quality metrics |
| `n_modalities` | int | Number of modalities fused |
| `effective_weights` | Dict[str, float] | Final weights applied |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `ModalityInput`

| Attribute | Type | Description |
|-----------|------|-------------|
| `modality_type` | ModalityType | Omic type |
| `pathway_scores` | pd.DataFrame | Z-normalized scores |
| `gene_data` | pd.DataFrame or None | Gene-level data |
| `label` | str or None | Human-readable label |
| `sample_ids` | List[str] | Property: sample IDs from index |
| `pathway_names` | List[str] | Property: pathway names from columns |
| `n_samples` | int | Property: number of samples |
| `n_pathways` | int | Property: number of pathways |

### `MultiOmicQualityReport`

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_modalities` | int | Number of modalities |
| `total_samples` | int | Union of all sample IDs |
| `shared_samples` | int | Intersection of all sample IDs |
| `sample_overlap_fraction` | float | shared / total |
| `pathway_overlap_count` | int | Number of shared pathways |
| `warnings` | List[str] | Quality warnings (e.g., low overlap) |
| `is_usable` | bool | Whether fusion result is viable |

**Methods:** `to_dict()`

---

## Fusion Strategy Details

### Concatenate (Default)

Column-binds all modality scores with label prefixes:

```
WES:PATHWAY_1  WES:PATHWAY_2  RNA-seq:PATHWAY_1  RNA-seq:PATHWAY_2
```

- Preserves all information from all modalities
- Handles partial sample overlap (missing values filled per `missing_strategy`)
- Output has more columns than any single modality

### Weighted Average

Computes weighted mean across modalities for shared pathways only:

```
fused[pathway] = w_vcf * vcf[pathway] + w_expr * expr[pathway]
```

- Requires shared pathways across modalities
- Weights default to uniform (1/n_modalities) if not specified
- Output has same dimensionality as single-modality

### Intersection Only

Restricts to samples and pathways present in all modalities. No imputation needed. Conservative but avoids any missing data assumptions.

---

## End-to-End Multi-Omic Workflow

```python
from pathway_subtyping import (
    load_expression_matrix, score_pathways_from_expression,
    ModalityType, FusionStrategy, prepare_modality, fuse_modalities,
    run_clustering, ValidationGates,
)

# 1. Score each modality
expr, _ = load_expression_matrix("expression.csv")
expr_result = score_pathways_from_expression(expr, pathways, seed=42)

# 2. Prepare modalities
vcf_mod = prepare_modality(ModalityType.VCF, vcf_scores, label="WES")
expr_mod = prepare_modality(ModalityType.EXPRESSION, expr_result.pathway_scores, label="RNA-seq")

# 3. Fuse
fusion = fuse_modalities([vcf_mod, expr_mod], strategy=FusionStrategy.CONCATENATE, seed=42)

# 4. Cluster fused scores
clustering = run_clustering(fusion.fused_pathway_scores.values, n_clusters=3, seed=42)

# 5. Validate (including cross-modal gate)
gates = ValidationGates(seed=42)
result = gates.run_all(
    pathway_scores=fusion.fused_pathway_scores,
    cluster_labels=clustering.labels,
    pathways=pathways,
    gene_burdens=fusion.fused_gene_data,
    n_clusters=3,
    per_modality_scores=fusion.per_modality_scores,  # enables Gate 5
)
```

---

## YAML Configuration

```yaml
data:
  input_type: multi_omic
  pathway_db: data/pathways/autism_pathways.gmt

multi_omic:
  modalities:
    - type: vcf
      path: data/variants.vcf
      label: WES
    - type: expression
      path: data/expression.csv
      label: RNA-seq
      expression_input_type: tpm
      expression_scoring_method: ssgsea
  fusion_strategy: concatenate
  missing_strategy: impute_zero
```

---

## See Also

- [Expression Scoring](expression.md) — bulk RNA-seq pathway scoring
- [Single-Cell Scoring](single_cell.md) — scRNA-seq pathway scoring
- [Deconvolution](deconvolution.md) — cell-type proportion estimation
- [Cross-Modal Validation](cross_modal_validation.md) — Gate 5 validation
