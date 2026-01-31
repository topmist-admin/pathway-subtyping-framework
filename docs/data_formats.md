# Data Formats Specification

> This document specifies the input and output data formats used in the Pathway Subtyping Framework.

## Input Formats

### VCF (Variant Call Format)

Standard VCF format (v4.2+) with required annotation fields.

#### Basic Structure

```
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD phred score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	12345	rs123	A	G	50	PASS	GENE=SHANK3;CONSEQUENCE=missense_variant;CADD=25.3	GT	0/1	0/0
```

#### Required Fields

| Field | Location | Description |
|-------|----------|-------------|
| CHROM | Column 1 | Chromosome (chr1-22, chrX, chrY) |
| POS | Column 2 | Position (1-based) |
| REF | Column 4 | Reference allele |
| ALT | Column 5 | Alternate allele(s) |
| GT | FORMAT | Genotype (0/0, 0/1, 1/1) |

#### Required INFO Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| GENE | String | Gene symbol (HGNC) | `GENE=SHANK3` |
| CONSEQUENCE | String | Variant consequence | `CONSEQUENCE=missense_variant` |

#### Recommended INFO Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| CADD | Float | CADD phred score | `CADD=25.3` |
| AF | Float | Allele frequency | `AF=0.0001` |
| gnomAD_AF | Float | gnomAD allele frequency | `gnomAD_AF=0.00005` |

#### Supported Consequence Terms

The framework recognizes VEP/Ensembl consequence terms:

**High Impact (LoF)**
- `frameshift_variant`
- `stop_gained`
- `splice_donor_variant`
- `splice_acceptor_variant`
- `start_lost`
- `stop_lost`

**Moderate Impact**
- `missense_variant`
- `inframe_insertion`
- `inframe_deletion`

**Low Impact**
- `synonymous_variant`
- `splice_region_variant`

#### Compression

- Uncompressed: `.vcf`
- Compressed: `.vcf.gz` (bgzip recommended)
- Index: `.vcf.gz.tbi` (tabix) for random access

#### Example VCF

```
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_001	SAMPLE_002	SAMPLE_003
chr22	51135990	.	G	A	99	PASS	GENE=SHANK3;CONSEQUENCE=missense_variant;CADD=28.5	GT	0/1	0/0	0/0
chr14	21853913	.	C	T	99	PASS	GENE=CHD8;CONSEQUENCE=frameshift_variant;CADD=35.2	GT	0/0	0/1	0/0
chr2	50145678	.	A	G	99	PASS	GENE=NRXN1;CONSEQUENCE=splice_donor_variant;CADD=32.1	GT	0/0	0/0	0/1
```

---

### GMT (Gene Matrix Transposed)

Pathway definition format with one pathway per line.

#### Format

```
PATHWAY_ID<TAB>DESCRIPTION<TAB>GENE1<TAB>GENE2<TAB>GENE3<TAB>...
```

#### Fields

| Position | Field | Description |
|----------|-------|-------------|
| 1 | Pathway ID | Unique identifier |
| 2 | Description | Human-readable name or URL |
| 3+ | Genes | Gene symbols (HGNC) |

#### Example GMT

```
SYNAPTIC_TRANSMISSION	Genes involved in synaptic signaling	SHANK3	SHANK2	NRXN1	NLGN1	GRIN2B	SCN2A
CHROMATIN_REMODELING	Chromatin modification genes	CHD8	CHD2	ARID1B	KMT2A	HDAC4	SETD5
ION_CHANNEL	Ion channel complex genes	SCN1A	SCN2A	KCNQ2	CACNA1A	GRIN2A	GRIN2B
IMMUNE_RESPONSE	Immune signaling pathway	IL6	TNF	IFNG	IL1B	NFKB1	STAT3
```

#### Validation Rules

- Minimum 3 genes per pathway
- Gene symbols should be HGNC-approved
- No duplicate pathway IDs
- Tab-separated (not spaces)

#### Sources for Pathway Files

| Source | URL | Format |
|--------|-----|--------|
| MSigDB | https://www.gsea-msigdb.org/ | GMT |
| Gene Ontology | http://geneontology.org/ | GAF → GMT |
| Reactome | https://reactome.org/ | GMT |
| KEGG | https://www.genome.jp/kegg/ | Custom → GMT |

---

### Phenotypes CSV

Sample metadata file (optional but recommended).

#### Format

```csv
sample_id,sex,age,cohort,planted_subtype
SAMPLE_001,M,25,discovery,synaptic
SAMPLE_002,F,32,discovery,chromatin
SAMPLE_003,M,28,replication,synaptic
```

#### Required Columns

| Column | Type | Description |
|--------|------|-------------|
| sample_id | string | Must match VCF sample IDs |

#### Optional Columns

| Column | Type | Description |
|--------|------|-------------|
| sex | string | M/F or Male/Female |
| age | numeric | Age at assessment |
| cohort | string | Cohort identifier |
| planted_subtype | string | Ground truth label (for validation) |

#### Notes

- Sample IDs must exactly match VCF column headers
- Missing values can be empty or `NA`
- Additional columns are preserved but not used

---

### Configuration YAML

Pipeline configuration file.

#### Full Schema

```yaml
# Run identification
run_name: my_analysis          # Required: unique run identifier

# Input files
input:
  vcf_path: path/to/variants.vcf        # Required: VCF file
  phenotypes_path: path/to/pheno.csv    # Optional: phenotypes
  pathways_path: path/to/pathways.gmt   # Required: pathway definitions

# Output settings
output:
  output_dir: outputs/my_analysis       # Output directory

# Pipeline settings
pipeline:
  seed: 42                    # Random seed (null = random)
  min_samples_per_cluster: 10 # Minimum cluster size

# Clustering settings
clustering:
  n_clusters_range: [2, 8]    # Range of K to try
  covariance_type: full       # GMM covariance type

# Validation settings
validation:
  run_validation: true        # Run validation gates
  n_permutations: 100         # Label shuffle iterations
  n_bootstrap: 50             # Bootstrap iterations
  ari_threshold: 0.7          # Minimum stability ARI
```

#### Minimal Configuration

```yaml
run_name: quick_test

input:
  vcf_path: data/variants.vcf
  pathways_path: data/pathways.gmt

pipeline:
  seed: 42
```

---

## Output Formats

### pathway_scores.csv

Sample-by-pathway score matrix.

#### Format

```csv
,PATHWAY_A,PATHWAY_B,PATHWAY_C,PATHWAY_D
SAMPLE_001,1.23,-0.45,0.87,0.12
SAMPLE_002,-0.67,1.89,-0.12,0.45
SAMPLE_003,0.34,0.23,1.56,-0.89
```

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| (index) | string | Sample ID |
| [pathway_name] | float | Z-score normalized pathway score |

### subtype_assignments.csv

Cluster assignments with confidence scores.

#### Format

```csv
sample_id,cluster_id,cluster_label,confidence,planted_subtype
SAMPLE_001,0,synaptic,0.92,synaptic
SAMPLE_002,1,chromatin,0.78,chromatin
SAMPLE_003,0,synaptic,0.85,synaptic
```

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| sample_id | string | Sample identifier |
| cluster_id | int | Numeric cluster (0-indexed) |
| cluster_label | string | Dominant pathway label |
| confidence | float | GMM posterior probability [0-1] |
| planted_subtype | string | Ground truth (if provided) |

### report.json

Machine-readable results summary.

#### Schema

```json
{
  "run_name": "string",
  "timestamp": "ISO8601 datetime",
  "seed": "integer",
  "summary": {
    "n_samples": "integer",
    "n_pathways": "integer",
    "n_clusters": "integer"
  },
  "clusters": {
    "0": {"count": "integer", "label": "string"}
  },
  "validation": {
    "all_passed": "boolean",
    "label_shuffle": {"status": "string", "value": "float"},
    "random_genes": {"status": "string", "value": "float"},
    "bootstrap_stability": {"status": "string", "value": "float"}
  }
}
```

### report.md

Human-readable Markdown report.

Contains sections for:
- Run metadata
- Input summary
- Clustering results
- Validation gate results
- Disclaimer

---

## File Size Guidelines

| Dataset Size | Expected VCF Size | Memory Needed |
|--------------|-------------------|---------------|
| 100 samples | ~10-50 MB | 2 GB |
| 1,000 samples | ~100-500 MB | 4 GB |
| 10,000 samples | ~1-5 GB | 16 GB |

---

## Validation Tools

### Validate VCF

```bash
# Check VCF structure
bcftools view -h your_file.vcf

# Check required INFO fields
bcftools query -f '%INFO/GENE\t%INFO/CONSEQUENCE\n' your_file.vcf | head
```

### Validate GMT

```python
from pathway_subtyping.io import load_pathways

# Will raise error if invalid
pathways = load_pathways("pathways.gmt")
print(f"Loaded {len(pathways)} pathways")
```

### Validate Config

```python
import yaml
from pathway_subtyping import PipelineConfig

# Load and validate
with open("config.yaml") as f:
    raw = yaml.safe_load(f)
config = PipelineConfig(**raw)
```

---

## See Also

- [Quickstart Guide](quickstart.md) - Getting started
- [Troubleshooting](troubleshooting.md) - Common format issues
- [API Reference](api/index.md) - Loading functions
