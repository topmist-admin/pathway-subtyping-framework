# Genetics API

> **Module**: `pathway_subtyping.genetics`

Genetic-anchoring primitives behind **Gate 7**: do the discovered subtypes line up with real genetic strata (GWAS risk-gene enrichment at the feature level, or somatic driver strata at the sample level) rather than with technical structure?

The gates themselves live in [`validation.md`](validation.md) (`genetic_anchoring_gate`, `somatic_anchoring_gate`); this module holds the statistics they call.

---

## Quick Example

```python
from pathway_subtyping.genetics import (
    feature_level_anchoring,
    hypergeometric_enrichment,
    somatic_alignment,
)

# Feature level — are risk genes concentrated in particular pathways?
results = feature_level_anchoring(
    gene_sets=pathways,            # {pathway_name: [gene, ...]}
    risk_genes=sfari_genes,
    universe=all_measured_genes,
)
for name, r in results.items():
    print(name, r.fold, r.p_value)

# Sample level — does the partition track somatic driver strata?
alignment = somatic_alignment(
    cluster_labels,
    somatic_strata={"BRAF_V600E": braf_calls, "MSI": msi_calls},
    cramers_v_min=0.30,
)
```

---

## Feature-level anchoring

### `hypergeometric_enrichment`

```python
hypergeometric_enrichment(
    test_set, risk_genes, universe, label="", null=""
) -> EnrichmentResult
```

One-sided hypergeometric test for over-representation of `risk_genes` within `test_set`, against `universe`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `test_set` | `Iterable[str]` | — | Genes under test (e.g. one pathway's members) |
| `risk_genes` | `Iterable[str]` | — | Reference risk-gene list |
| `universe` | `Iterable[str]` | — | Background gene universe |
| `label` | `str` | `""` | Free-text label carried into the result |
| `null` | `str` | `""` | Free-text description of the null, carried into the result |

### `feature_level_anchoring`

```python
feature_level_anchoring(
    gene_sets, risk_genes, universe, reference_universe=None
) -> Dict[str, EnrichmentResult]
```

Runs `hypergeometric_enrichment` across every entry in `gene_sets` and returns one `EnrichmentResult` per set.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gene_sets` | `Dict[Any, Iterable[str]]` | — | `{name: genes}` to test |
| `risk_genes` | `Iterable[str]` | — | Reference risk-gene list |
| `universe` | `Iterable[str]` | — | Background universe |
| `reference_universe` | `Optional[Iterable[str]]` | `None` | Alternative universe for the null, when it should differ from the test universe |

### `EnrichmentResult`

| Field | Type | Description |
|---|---|---|
| `label` | `str` | Identifier passed in |
| `null` | `str` | Description of the null used |
| `universe_n` | `int` | Universe size |
| `risk_in_universe` | `int` | Risk genes present in the universe |
| `testset_n` | `int` | Size of the tested set |
| `risk_hits` | `int` | Risk genes in the tested set |
| `expected` | `float` | Expected hits under the null |
| `fold` | `float` | Observed / expected |
| `p_value` | `float` | Hypergeometric p |

> **Read `fold` next to `risk_in_universe`.** A large fold change over a handful of risk genes in the universe is not strong evidence, and the result object reports both so the denominator cannot be lost.

---

## Sample-level somatic anchoring

### `somatic_alignment`

```python
somatic_alignment(
    cluster_labels, somatic_strata, cramers_v_min=0.30, alpha=0.05
) -> Dict[str, Any]
```

Tests the partition against each named somatic stratum by chi-square with a Cramér's V effect size.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `cluster_labels` | array-like | — | Cluster assignment per sample |
| `somatic_strata` | `Dict[str, Any]` | — | `{stratum_name: per-sample categorical values}` |
| `cramers_v_min` | `float` | `0.30` | Effect-size floor for a stratum to count as aligned |
| `alpha` | `float` | `0.05` | Significance level |

### `StratumAlignment`

| Field | Type | Description |
|---|---|---|
| `stratum` | `str` | Stratum name |
| `chi2` | — | Chi-square statistic |
| `p_value` | — | Chi-square p |
| `cramers_v` | `float` | Effect size (0 = independent, 1 = perfect association) |
| `dof` | `int` | Degrees of freedom |
| `n` | `int` | Samples contributing |

> **Statistical significance is not enough on its own here.** At large *n* a trivial association reaches p < 0.05, which is why the alignment requires both `p < alpha` **and** `cramers_v >= cramers_v_min`. This is the same two-part criterion the confound gate uses, for the same reason.

---

## See Also

- [`validation.md`](validation.md) — `genetic_anchoring_gate` and `somatic_anchoring_gate`, which consume these
- [`../discreteness_gate.md`](../discreteness_gate.md) — the gate that decides whether there is a partition worth anchoring in the first place
