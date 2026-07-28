# Public Dataset Validation API

> **Module**: `pathway_subtyping.validation_datasets`

Downloads open-access reference data (ClinVar, Reactome) and validates that a project's pathway definitions correspond to real biology — that the genes in a custom GMT are actually disease-relevant, and that its pathways overlap curated equivalents — rather than being an arbitrary partition of gene space.

References:
- Landrum MJ et al. ClinVar: improvements to accessing data. *Nucleic Acids Res.* 2020;48(D1):D835-D844.
- Gillespie M et al. The Reactome Pathway Knowledgebase 2022.

---

## Quick Example

```python
from pathway_subtyping.validation_datasets import run_full_validation

report = run_full_validation(
    gmt_path="data/pathways/autism_pathways.gmt",
    disease_name="autism",
    n_samples=100,
    n_subtypes=3,
    seed=42,
)

print(report.overall_pass)
for w in report.warnings:
    print("WARN:", w)
```

Network access is required unless `skip_download=True`. Downloads are cached; see [Caching and offline use](#caching-and-offline-use).

---

## Loading reference data

### `download_dataset`

```python
download_dataset(dataset_key, cache_dir=None, force=False, timeout=60) -> Path
```

Fetches a registered dataset to the cache and returns its path. Raises `urllib.error.URLError` on failure.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dataset_key` | `str` | — | Registered dataset name (e.g. `"reactome_pathways"`) |
| `cache_dir` | `Optional[Path]` | `None` | Cache location; a default is used when omitted |
| `force` | `bool` | `False` | Re-download even if cached |
| `timeout` | `int` | `60` | Per-request timeout, seconds |

### `load_clinvar_gene_summary`

```python
load_clinvar_gene_summary(file_path=None, cache_dir=None) -> Dict[str, ClinVarGeneSummary]
```

Returns per-gene pathogenic/benign variant counts keyed by gene symbol.

### `load_reactome_pathways`

```python
load_reactome_pathways(file_path=None, cache_dir=None, species="Homo sapiens") -> Dict[str, List[str]]
```

Returns `{pathway_name: [gene_symbol, ...]}` for the requested species.

---

## Validating pathway definitions

### `validate_pathway_coverage`

```python
validate_pathway_coverage(pathways, clinvar_genes) -> List[PathwayCoverageResult]
```

For each pathway, how many of its genes appear in ClinVar at all, and how many carry pathogenic variants.

### `validate_pathway_against_reactome`

```python
validate_pathway_against_reactome(
    pathways, reactome_pathways, min_overlap=0.1
) -> Dict[str, Dict[str, Any]]
```

Cross-references custom pathway definitions against Reactome, reporting the best-matching curated pathway per input set. `min_overlap` is the Jaccard-style floor below which a match is not reported.

### `run_biological_plausibility_check`

```python
run_biological_plausibility_check(
    pathway_scores, cluster_labels, pathways, clinvar_genes,
    disease_name="autism", seed=42,
) -> BiologicalPlausibilityResult
```

Asks whether discovered subtypes are enriched for pathways whose genes carry known pathogenic variation — i.e. whether the partition is biologically interpretable, not merely statistically stable.

> This is a **plausibility** check, not a validation gate. It does not test whether the subtypes are discrete or confound-free. Use it alongside [`../discreteness_gate.md`](../discreteness_gate.md) (Gate A) and the confound gate in [`validation.md`](validation.md), not instead of them.

---

## Synthetic data with realistic structure

### `generate_disease_realistic_synthetic`

```python
generate_disease_realistic_synthetic(
    pathways, clinvar_genes, n_samples=100, n_subtypes=3,
    effect_size=1.0, noise_level=1.0, seed=42,
) -> SimulatedData
```

Generates synthetic cohorts whose signal is placed in pathways that ClinVar says actually carry pathogenic variation, rather than in arbitrary gene sets. Returns the same `SimulatedData` type as [`simulation.md`](simulation.md).

---

## End-to-end

### `run_full_validation`

```python
run_full_validation(
    gmt_path=None, disease_name="autism", n_samples=100, n_subtypes=3,
    effect_size=1.0, cache_dir=None, seed=42, skip_download=False,
) -> ValidationReport
```

Runs coverage, Reactome cross-reference, synthetic validation and the plausibility check, returning one `ValidationReport`. Pass `skip_download=True` to run only the parts that need no network.

---

## Result types

### `ValidationReport`

| Field | Type | Description |
|---|---|---|
| `timestamp` | `str` | Run timestamp |
| `data_sources` | `List[DatasetInfo]` | Which reference datasets were used |
| `pathway_coverage` | `List[PathwayCoverageResult]` | Per-pathway ClinVar coverage |
| `reactome_cross_ref` | `Dict[str, Dict[str, Any]]` | Per-pathway Reactome match |
| `synthetic_validation` | `Dict[str, Any]` | Synthetic-cohort recovery results |
| `biological_plausibility` | `Optional[BiologicalPlausibilityResult]` | Plausibility outcome |
| `overall_pass` | `bool` | Aggregate verdict |
| `warnings` | `List[str]` | Non-fatal issues — **read these even when `overall_pass` is True** |

### `PathwayCoverageResult`

`pathway_name`, `total_genes`, `genes_in_clinvar`, `genes_with_pathogenic`, `coverage_fraction`, `pathogenic_fraction`, `gene_details`.

### `ClinVarGeneSummary`

`symbol`, `gene_id`, `n_pathogenic`, `n_likely_pathogenic`, `n_uncertain`, `n_benign`, `n_likely_benign`.

### `BiologicalPlausibilityResult`

`disease_name`, `n_subtypes`, `n_enriched_pathways`, `pathway_gene_clinvar_overlap`, `subtypes_biologically_distinct`, `details`.

### `DatasetInfo`

`name`, `url`, `description`, `file_format`, `license_note`.

---

## Caching and offline use

Downloads are cached under `cache_dir`. Pass `skip_download=True` to `run_full_validation` for a network-free run.

> **Known flake.** `tests/test_validation_datasets.py::TestNetworkIntegration::test_reactome_download_and_parse` depends on a live Reactome endpoint and intermittently fails with `HTTP 403 Forbidden`. That is an upstream availability issue, not a defect in this module. Deselect it for offline or CI-constrained runs:
>
> ```bash
> pytest --deselect tests/test_validation_datasets.py::TestNetworkIntegration::test_reactome_download_and_parse
> ```

---

## See Also

- [`simulation.md`](simulation.md) — the `SimulatedData` type returned by the generator
- [`validation.md`](validation.md) — the validation gates proper
- [`../discreteness_gate.md`](../discreteness_gate.md) — Gate A
