# Cross-Modal Validation API

> **Module**: `pathway_subtyping.cross_modal_validation`

Tests whether molecular subtypes are consistent across different data modalities (e.g., VCF variant burden vs bulk RNA-seq expression). If subtypes discovered from one data type also appear when clustering a completely different data type, they likely reflect real biology rather than modality-specific artifacts.

This implements **Validation Gate 5** for multi-omic analyses.

---

## Quick Example

```python
from pathway_subtyping import cross_modal_concordance

# per_modality_scores: dict of modality_label -> pathway_scores DataFrame
result = cross_modal_concordance(
    per_modality_scores={"WES": vcf_scores, "RNA-seq": expr_scores},
    cluster_labels=fused_labels,
    fused_sample_ids=list(fused_scores.index),
    n_clusters=3,
    n_permutations=100,
    seed=42,
)

print(f"Gate passed: {result.gate_passed}")
print(f"Mean concordance ARI: {result.mean_concordance_ari:.3f}")
print(f"Null 95th percentile: {result.null_ari_95th:.3f}")
print(result.format_report())
```

---

## Functions

### `cross_modal_concordance()`

```python
def cross_modal_concordance(
    per_modality_scores: Dict[str, pd.DataFrame],
    cluster_labels: np.ndarray,
    fused_sample_ids: List[str],
    n_clusters: int,
    n_permutations: int = 100,
    concordance_threshold: Optional[float] = None,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> CrossModalValidationResult
```

Main entry point for cross-modal validation.

**How it works:**
1. For each pair of modalities, cluster each independently using GMM
2. Measure agreement via ARI and NMI (concordance)
3. Test transfer: train on modality A, predict on modality B (and vice versa)
4. Build null distribution by permuting labels across `n_permutations` runs
5. Gate passes if observed ARI > 95th percentile of null distribution

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `per_modality_scores` | Dict[str, pd.DataFrame] | required | Modality label to pathway scores (2+ entries) |
| `cluster_labels` | np.ndarray | required | Cluster labels from fused analysis |
| `fused_sample_ids` | List[str] | required | Sample IDs used in fused analysis |
| `n_clusters` | int | required | Number of clusters |
| `n_permutations` | int | 100 | Permutations for null distribution |
| `concordance_threshold` | float or None | None | Override gate threshold (default: null 95th percentile) |
| `seed` | int or None | None | Random seed |
| `show_progress` | bool | True | Show tqdm progress bar |

**Returns:** `CrossModalValidationResult`

**Raises:** `ValueError` if fewer than 2 modalities provided.

---

### `single_cell_composition_test()`

```python
def single_cell_composition_test(
    cell_type_fractions: pd.DataFrame,
    bulk_labels: np.ndarray,
    significance_threshold: float = 0.05,
) -> SingleCellCompositionResult
```

Test whether bulk subtypes have distinct cell-type compositions using one-way ANOVA (Kruskal-Wallis) per cell type with Bonferroni correction.

---

### `generate_synthetic_multimodal_data()`

```python
def generate_synthetic_multimodal_data(
    n_samples: int = 60,
    n_pathways: int = 10,
    n_subtypes: int = 3,
    n_modalities: int = 2,
    signal_strength: float = 2.0,
    seed: Optional[int] = None,
) -> Tuple[Dict[str, pd.DataFrame], np.ndarray]
```

Generate synthetic multi-modal data with shared planted subtypes for testing.

**Returns:** `(per_modality_scores, true_labels)` tuple.

---

## Dataclasses

### `CrossModalValidationResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `pair_results` | List[CrossModalPairResult] | Per-pair concordance metrics |
| `mean_concordance_ari` | float | Average ARI across all modality pairs |
| `mean_concordance_nmi` | float | Average NMI across all modality pairs |
| `mean_transfer_ari` | float | Average transfer ARI |
| `null_ari_95th` | float | 95th percentile of null ARI distribution |
| `gate_passed` | bool | Whether observed ARI exceeds null threshold |
| `single_cell_composition` | SingleCellCompositionResult or None | Cell-type composition test result |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

### `CrossModalPairResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `modality_a` | str | First modality label |
| `modality_b` | str | Second modality label |
| `concordance_ari` | float | ARI between independent clusterings |
| `concordance_nmi` | float | NMI between independent clusterings |
| `transfer_ari_a_to_b` | float | ARI when training on A, predicting on B |
| `transfer_ari_b_to_a` | float | ARI when training on B, predicting on A |
| `n_shared_samples` | int | Samples present in both modalities |
| `n_shared_pathways` | int | Pathways present in both modalities |

**Methods:** `to_dict()`

### `SingleCellCompositionResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `f_statistic` | float | Overall F-statistic |
| `p_value` | float | Overall p-value |
| `n_subtypes` | int | Number of subtypes tested |
| `n_cell_types` | int | Number of cell types |
| `per_cell_type_pvalues` | Dict[str, float] | Per-cell-type p-values |
| `significant_cell_types` | List[str] | Cell types with p < Bonferroni threshold |

**Methods:** `to_dict()`

---

## Using Gate 5 via ValidationGates

Gate 5 integrates automatically when you pass `per_modality_scores` to `ValidationGates.run_all()`:

```python
from pathway_subtyping import ValidationGates

gates = ValidationGates(seed=42, n_permutations=100)
result = gates.run_all(
    pathway_scores=fused_scores,
    cluster_labels=labels,
    pathways=pathways,
    gene_burdens=gene_data,
    n_clusters=3,
    per_modality_scores={"WES": vcf_scores, "RNA-seq": expr_scores},
)
# Gate 5 is automatically included when >= 2 modalities provided
```

---

## Interpretation

| Metric | Good | Concerning |
|--------|------|------------|
| `mean_concordance_ari` | > 0.3 | < 0.1 (subtypes may be modality-specific) |
| `gate_passed` | True | False (subtypes not replicated across modalities) |
| `transfer_ari` | > 0.2 | 0.0 (no shared feature space or no transferability) |
| `null_ari_95th` | < 0.15 | > 0.3 (data may have structure artifacts) |

---

## See Also

- [Multi-Omic Fusion](multi_omic.md) — fusing modalities before validation
- [Validation Gates](validation.md) — all validation gates
- [Deconvolution](deconvolution.md) — cell-type proportion estimation
