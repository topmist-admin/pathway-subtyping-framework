# VCF Data Quality API

> **Module**: `pathway_subtyping.data_quality`

Loads and validates real-world VCF files before analysis: multi-allelic variant expansion, genotype parsing, annotation coverage reporting (GENE / CONSEQUENCE / CADD), and error messages that carry actionable fix suggestions.

---

## Quick Example

```python
from pathway_subtyping import (
    load_vcf_with_quality_check,
    validate_vcf_for_pipeline,
    VCFDataQualityError,
)

# Pre-flight check: prints a summary and returns fix suggestions
is_valid, report, fix_suggestions = validate_vcf_for_pipeline("cohort.vcf", verbose=True)

# Full load
variants_df, genotypes_df, samples, report = load_vcf_with_quality_check(
    "cohort.vcf",
    strict=False,
)
print(report.gene_coverage, report.is_usable)
print(report.summary())

# strict=True turns an unusable file into an exception carrying the report
try:
    load_vcf_with_quality_check("unannotated.vcf", strict=True)
except VCFDataQualityError as e:
    print(e.fix_suggestions)
```

`DataQualityReport`, `VCFDataQualityError`, `load_vcf_with_quality_check`, and `validate_vcf_for_pipeline` are re-exported from the top-level `pathway_subtyping` package. The lower-level parsers below are imported from `pathway_subtyping.data_quality` directly.

---

## Functions

### `load_vcf_with_quality_check()`

```python
def load_vcf_with_quality_check(
    vcf_path: str,
    strict: bool = False,
    min_gene_coverage: float = 50.0,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], DataQualityReport]
```

Load a VCF (plain or `.gz`, detected by file extension) with comprehensive quality checking.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vcf_path` | str | required | Path to VCF file |
| `strict` | bool | False | If True, raise `VCFDataQualityError` when `report.is_usable` is False |
| `min_gene_coverage` | float | 50.0 | Minimum percentage of variants with gene annotation |

**Returns:** `(variants_df, genotypes_df, samples, quality_report)`.

- `variants_df` columns: `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `filter`, `gene`, `consequence`, `cadd`. Missing `qual` and `cadd` values become `0.0`; missing `gene` / `consequence` become `""`.
- `genotypes_df` has one row per parsed (expanded) variant and one column per sample, holding the count of the relevant alternate allele.
- `samples` is the list of sample IDs from the `#CHROM` header line.

**Key behaviors:**
- Multi-allelic records are expanded into one bi-allelic record per ALT, with IDs suffixed `_1`, `_2`, … and `original_alts` retained internally during expansion.
- Variant lines with fewer than 9 columns are counted as skipped and recorded as a warning.
- Coverage warnings are appended to the report below 50% GENE coverage, below 50% CONSEQUENCE coverage, and below 30% CADD coverage.

**Raises:**
- `FileNotFoundError` if the file does not exist (message includes a path/permissions checklist).
- `VCFDataQualityError` on a corrupt gzip (`gzip.BadGzipFile`), on a `UnicodeDecodeError`, or when `strict=True` and the data is not usable.

> Note: `min_gene_coverage` is accepted for API compatibility but is not read by the function body. The usability decision comes from `DataQualityReport.is_usable`, which applies a fixed 50% GENE-coverage threshold, and the coverage warnings use fixed thresholds as listed above.

---

### `validate_vcf_for_pipeline()`

```python
def validate_vcf_for_pipeline(
    vcf_path: str,
    verbose: bool = True,
) -> Tuple[bool, DataQualityReport, List[str]]
```

Pre-flight validation. Calls `load_vcf_with_quality_check(vcf_path, strict=False)` and adds cohort-shape warnings.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vcf_path` | str | required | Path to VCF file |
| `verbose` | bool | True | Print `report.summary()` and any suggestions to stdout |

**Returns:** `(is_valid, report, fix_suggestions)` where `is_valid` is `report.is_usable and not report.errors`.

**Additional warnings raised by this function:**
- Fewer than 10 samples — clustering may not be meaningful.
- Fewer than 100 parsed variants — analysis may be underpowered.

**Failure modes (returns rather than raises):** a missing file returns `(False, report, ["Verify the file path and permissions"])`; a `VCFDataQualityError` returns `(False, e.report, e.fix_suggestions)`.

---

### `parse_genotype()`

```python
def parse_genotype(gt_string: str, target_allele: int = 1) -> Tuple[int, bool]
```

Parse a VCF genotype string and return the count of a specific allele.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gt_string` | str | required | Genotype string, e.g. `"0/1"`, `"1\|2"`; only the first colon-separated field is used |
| `target_allele` | int | 1 | Which alternate allele to count (1 = first ALT, 2 = second, …) |

**Returns:** `(allele_count, is_valid)`. Counts occurrences of the specific allele number, not any non-reference allele. Missing (`.`, `./.`, `.|.`), empty, non-diploid, or unparseable genotypes return `(0, False)`.

| Input | `target_allele` | Result |
|-------|-----------------|--------|
| `"0/0"` | 1 | `(0, True)` |
| `"0/1"` | 1 | `(1, True)` |
| `"1/1"` | 1 | `(2, True)` |
| `"1\|2"` | 2 | `(1, True)` |
| `"./."` | 1 | `(0, False)` |

---

### `parse_info_field()`

```python
def parse_info_field(info_string: str) -> Dict[str, Any]
```

Parse a VCF INFO field into a dictionary. `Key=Value` pairs become string values; comma-separated values become a list of strings; `Key=.` becomes `None`; valueless flag fields become `True`. An empty string or `"."` returns an empty dict. Keys keep their original case.

---

### `expand_multi_allelic()`

```python
def expand_multi_allelic(
    chrom: str,
    pos: int,
    vid: str,
    ref: str,
    alts: List[str],
    info_dict: Dict[str, Any],
    sample_genotypes: List[str],
    samples: List[str],
) -> List[Dict[str, Any]]
```

Expand one multi-allelic variant into a list of bi-allelic records, one per entry in `alts`.

Each returned record has keys `chrom`, `pos`, `id`, `ref`, `alt`, `original_alts`, `genotypes`, plus lowercased INFO keys. When an INFO value is a list whose length matches `len(alts)`, the allele-matching element is used; other list values fall back to their first element. `genotypes` maps each sample ID to the count of that specific ALT allele (missing or malformed genotypes count as `0`). IDs are suffixed `_1`, `_2`, … when there is more than one ALT, and `original_alts` holds the comma-joined original ALT string in that case (otherwise `None`).

---

## Dataclasses

### `DataQualityReport`

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `total_variants` | int | 0 | VCF variant lines seen |
| `parsed_variants` | int | 0 | Records produced after multi-allelic expansion |
| `skipped_variants` | int | 0 | Lines skipped as malformed |
| `variants_with_gene` | int | 0 | Records with a usable GENE annotation |
| `variants_with_consequence` | int | 0 | Records with a usable CONSEQUENCE annotation |
| `variants_with_cadd` | int | 0 | Records with a float-parseable CADD score |
| `multi_allelic_variants` | int | 0 | Multi-allelic source records |
| `multi_allelic_expanded` | int | 0 | Bi-allelic records produced from them |
| `missing_genotypes` | int | 0 | Genotype calls recorded as missing |
| `malformed_genotypes` | int | 0 | Genotype calls that failed to parse |
| `warnings` | List[str] | `[]` | Accumulated warnings |
| `errors` | List[str] | `[]` | Accumulated errors |

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `gene_coverage` | float | Percentage of parsed variants with a gene annotation (0.0 when none parsed) |
| `consequence_coverage` | float | Percentage with a consequence annotation |
| `cadd_coverage` | float | Percentage with a CADD score |
| `is_usable` | bool | `gene_coverage >= 50.0` and `parsed_variants > 0` |

**Methods:**

- `add_warning(message)` / `add_error(message)` — append (de-duplicated) and log at the corresponding level.
- `to_dict()` — nested dictionary with `annotation_coverage`, `multi_allelic`, `genotype_issues`, `warnings`, `errors`, and `is_usable`. Coverage values are pre-formatted percentage strings such as `"100.0%"`.
- `summary()` — human-readable report ending in `Data Quality Status: PASS` or `FAIL`.

---

## Exceptions

### `VCFDataQualityError`

```python
VCFDataQualityError(message: str, report: DataQualityReport, fix_suggestions: List[str])
```

Raised when VCF data quality is insufficient. Attributes: `message`, `report` (the `DataQualityReport` at the point of failure), and `fix_suggestions` (list of remediation strings). The exception text is composed from all three — it lists the gene-coverage shortfall when below 50%, the first three warnings, and the numbered fix suggestions.

---

## Notes

- Internal helpers `_safe_float()` and `_clean_annotation()` are implementation details and are not documented here.
- Annotation lookup is INFO-key based: `GENE`, `CONSEQUENCE`, and `CADD` are read from the INFO field (lowercased during parsing). VCFs annotated only via a packed `CSQ`/`ANN` field will report low coverage — run VEP or ANNOVAR so that these keys are present.

---

## See Also

- [Variant QC](variant_qc.md) — quality filters applied before burden computation
- [Limitations](../limitations.md) — where data quality bounds interpretation
- [Troubleshooting](../troubleshooting.md) — data issue walkthroughs
