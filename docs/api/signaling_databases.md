# Signaling Databases API

> **Module**: `pathway_subtyping.signaling_databases`

Loads cell-cell signaling interaction databases (CellPhoneDB, CellChatDB) and converts them into pathway gene sets (`Dict[str, List[str]]`) compatible with the existing pathway scoring and clustering infrastructure.

---

## Quick Example

```python
from pathway_subtyping import (
    load_cellphonedb,
    load_cellchatdb,
    merge_signaling_databases,
    SignalingDatabase,
)

# CellPhoneDB — auto-download from GitHub
cpdb = load_cellphonedb(min_genes_per_pathway=3)
print(cpdb.format_report())
# cpdb.pathway_gene_sets: Dict[str, List[str]]

# CellChatDB — user-exported CSV from R
ccdb = load_cellchatdb("cellchatdb_human.csv")

# Merge both databases
merged = merge_signaling_databases(cpdb, ccdb)
# merged: Dict[str, List[str]] — ready for scoring and clustering
```

---

## Enums

### `SignalingDatabase`

| Value | Description |
|-------|-------------|
| `CELLPHONEDB` | CellPhoneDB ligand-receptor interactions (Efremova et al., 2020). CSV data auto-downloaded from GitHub. |
| `CELLCHATDB` | CellChatDB ligand-receptor interactions (Jin et al., 2021). Requires user-exported CSV from R. |

---

## Functions

### `load_cellphonedb()`

```python
def load_cellphonedb(
    interactions_path: Optional[Path] = None,
    gene_path: Optional[Path] = None,
    complex_path: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
    force_download: bool = False,
    timeout: int = 60,
    min_genes_per_pathway: int = 2,
) -> SignalingDatabaseResult
```

Load CellPhoneDB interactions and convert to pathway gene sets. Downloads three CSV files from the CellPhoneDB GitHub repository, resolves multi-subunit complexes to gene symbols, and groups interactions by their `classification` field.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `interactions_path` | Path or None | None | Explicit path to interaction_input.csv (skips download) |
| `gene_path` | Path or None | None | Explicit path to gene_input.csv (skips download) |
| `complex_path` | Path or None | None | Explicit path to complex_input.csv (skips download) |
| `cache_dir` | Path or None | `data/validation_cache/` | Cache directory for downloads |
| `force_download` | bool | False | Re-download even if cached |
| `timeout` | int | 60 | Download timeout in seconds |
| `min_genes_per_pathway` | int | 2 | Minimum genes per pathway to include |

**Returns:** `SignalingDatabaseResult` with pathway gene sets grouped by signaling classification.

**Raises:** `FileNotFoundError` if explicit paths don't exist. `ValueError` if CSVs are empty or missing required columns. `urllib.error.URLError` if download fails.

---

### `load_cellchatdb()`

```python
def load_cellchatdb(
    file_path: Path,
    species: str = "human",
    min_genes_per_pathway: int = 2,
) -> SignalingDatabaseResult
```

Load CellChatDB interactions from a user-exported CSV file. CellChatDB data is distributed as R binary `.rda` files. Export from R with:

```r
library(CellChat)
write.csv(CellChatDB.human$interaction, "cellchatdb_human.csv", row.names = FALSE)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file_path` | Path | required | Path to exported CellChatDB CSV |
| `species` | str | "human" | Species label for reporting |
| `min_genes_per_pathway` | int | 2 | Minimum genes per pathway to include |

**Returns:** `SignalingDatabaseResult` with pathway gene sets grouped by `pathway_name`.

**Raises:** `FileNotFoundError` if file doesn't exist. `ValueError` if required columns (`pathway_name`, `ligand`, `receptor`) are missing.

---

### `convert_interactions_to_pathways()`

```python
def convert_interactions_to_pathways(
    interactions: List[SignalingInteraction],
    group_by: str = "pathway_name",
    min_genes: int = 2,
    prefix: str = "Signaling",
) -> Dict[str, List[str]]
```

Convert a list of signaling interactions to pathway gene sets by grouping and collecting the union of all gene symbols per group.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `interactions` | List[SignalingInteraction] | required | Interactions to convert |
| `group_by` | str | "pathway_name" | Field to group by: `"pathway_name"` or `"annotation"` |
| `min_genes` | int | 2 | Minimum genes per group to include |
| `prefix` | str | "Signaling" | Prefix for pathway names (empty string to disable) |

**Returns:** `Dict[str, List[str]]` — pathway name to sorted list of unique gene symbols.

---

### `merge_signaling_databases()`

```python
def merge_signaling_databases(
    *results: SignalingDatabaseResult,
) -> Dict[str, List[str]]
```

Merge pathway gene sets from multiple signaling database results. Overlapping pathway names have their gene sets unioned.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `*results` | SignalingDatabaseResult | One or more results to merge |

**Returns:** `Dict[str, List[str]]` — merged pathway name to sorted list of unique gene symbols.

---

## Dataclasses

### `SignalingInteraction`

| Attribute | Type | Description |
|-----------|------|-------------|
| `interaction_id` | str | Unique identifier for the interaction |
| `partner_a` | str | First partner (gene symbol, UniProt ID, or complex name) |
| `partner_b` | str | Second partner |
| `genes_a` | List[str] | Resolved gene symbols for partner A |
| `genes_b` | List[str] | Resolved gene symbols for partner B |
| `pathway_name` | str | Signaling pathway classification |
| `source_database` | SignalingDatabase | Which database this came from |
| `annotation` | str | Additional annotation |

**Properties:** `all_genes` — all unique gene symbols (sorted).

**Methods:** `to_dict()`

### `SignalingDatabaseResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `database` | SignalingDatabase | Which database was loaded |
| `pathway_gene_sets` | Dict[str, List[str]] | Pathway name → gene list (the main output) |
| `n_interactions` | int | Total interactions parsed |
| `n_pathways` | int | Number of unique pathway gene sets |
| `n_unique_genes` | int | Total unique genes across all pathways |
| `interactions` | List[SignalingInteraction] | Full list of parsed interactions |
| `warnings` | List[str] | Any warnings during parsing |

**Methods:** `to_dict()`, `format_report()`, `get_citations()`

---

## Workflow

### Using Signaling Pathways for Subtype Discovery

Signaling pathway gene sets are `Dict[str, List[str]]` — the same format used throughout the framework. They can be passed directly to pathway scoring:

```python
from pathway_subtyping import (
    load_cellphonedb,
    load_expression_matrix,
    score_pathways_from_expression,
    run_clustering,
    ExpressionScoringMethod,
    ExpressionInputType,
)

# 1. Load signaling pathways
cpdb = load_cellphonedb(min_genes_per_pathway=5)
signaling_pathways = cpdb.pathway_gene_sets

# 2. Score expression data against signaling pathways
expr, qc = load_expression_matrix("expression.csv", input_type=ExpressionInputType.TPM)
result = score_pathways_from_expression(
    expr, signaling_pathways,
    method=ExpressionScoringMethod.SSGSEA, seed=42,
)

# 3. Cluster to discover signaling-based subtypes
clustering = run_clustering(result.pathway_scores.values, n_clusters=3, seed=42)
```

### Combining with Reactome Pathways

```python
from pathway_subtyping import load_cellphonedb, load_reactome_pathways

# Load both sources
reactome = load_reactome_pathways(species="Homo sapiens")
cpdb = load_cellphonedb(min_genes_per_pathway=3)

# Merge into one pathway dictionary
combined = {**reactome, **cpdb.pathway_gene_sets}
# Use combined for scoring
```

---

## Complex Resolution

CellPhoneDB interactions often involve multi-subunit receptor or ligand complexes (e.g., `integrin_a2b1_complex`). The loader automatically resolves these:

1. Downloads `gene_input.csv` (UniProt → gene symbol mapping)
2. Downloads `complex_input.csv` (complex → UniProt subunit mapping)
3. For each interaction partner:
   - If it's a complex name → resolve to all constituent gene symbols
   - If it's a UniProt accession → look up gene symbol
   - Otherwise → treat as a gene name directly

This ensures that all genes contributing to a signaling pathway are captured in the gene set.

---

## Interpretation

| Metric | Typical | Notes |
|--------|---------|-------|
| CellPhoneDB pathways | 50-100 | Grouped by `classification` |
| CellPhoneDB unique genes | 500-1000 | After complex resolution |
| CellChatDB pathways | 100-250 | Grouped by `pathway_name` |
| Genes per signaling pathway | 2-50 | Varies by pathway |

Pathway names are prefixed with `"Signaling: "` to distinguish from metabolic/disease pathways when merged.

---

## References

- Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. *Nat Protoc*. 2020;15(4):1484-1506.
- Jin S, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun*. 2021;12(1):1088.

---

## See Also

- [Expression Scoring](expression.md) — score pathways from bulk RNA-seq
- [Validation Datasets](../api/index.md) — Reactome pathway loading
- [Multi-Omic Fusion](multi_omic.md) — combine multiple pathway sources
- [Clustering](clustering.md) — discover subtypes from pathway scores
