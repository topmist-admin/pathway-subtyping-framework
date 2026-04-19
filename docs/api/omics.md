# `omics` — Multi-omics pathway-score fusion (F10)

ATAC-seq and proteomics pathway scoring, weighted fusion with RNA, and explicit discordance detection for pathways where modalities disagree.

**Import path:** `pathway_subtyping.omics`
**Source:** [`src/pathway_subtyping/omics/`](../../src/pathway_subtyping/omics/)
**Guide:** [multi-omics.md](../guides/multi-omics.md) · [Notebook 28](../../examples/notebooks/28_multi_omics.ipynb)
**Install:** core (no extra needed)

---

## Public surface

```python
from pathway_subtyping.omics import (
    # Per-modality scorers
    ATACScorer, score_atac_pathways,
    ProteomicsScorer, score_proteomics_pathways,
    # Fusion
    FusionWeights, FusionResult, MultiOmicsFusion,
    # Discordance
    DiscordanceReport, flag_discordant_pathways,
)
```

---

## `ATACScorer`

```python
ATACScorer(
    peak_to_gene: Optional[Mapping[str, Sequence[str]]] = None,
    aggregation: str = "sum",          # or "mean"
)
```

- `peak_to_gene`: mapping `{peak_id: [gene_1, gene_2, ...]}`. Bring your own peak-to-gene annotation (bedtools intersect, ArchR, etc.).
- If `peak_to_gene=None`, assumes the ATAC matrix columns are already gene symbols.

**Method:**

```python
scores = scorer.score(
    accessibility: pd.DataFrame,         # samples × peaks (or genes if peak_to_gene=None)
    pathways: Mapping[str, Iterable[str]],
    min_peaks_per_pathway: int = 2,
)
# scores: DataFrame(samples × pathways), z-normalised per pathway
```

The convenience function `score_atac_pathways(...)` is a one-shot wrapper that builds the scorer internally.

---

## `ProteomicsScorer`

```python
ProteomicsScorer(
    gene_to_protein: Optional[Mapping[str, str]] = None,
)
```

- `gene_to_protein`: mapping from gene symbol → protein accession (UniProt ID etc.). When supplied, pathway membership (expressed as gene symbols) is translated before lookup.
- If `gene_to_protein=None`, pathway members are assumed to match protein-matrix column names directly.

**Method:**

```python
scores = scorer.score(
    abundance: pd.DataFrame,              # samples × proteins (log-normalised)
    pathways: Mapping[str, Iterable[str]],
    min_proteins_per_pathway: int = 2,
)
# scores: DataFrame(samples × pathways), z-normalised per pathway
```

`score_proteomics_pathways(...)` is the one-shot convenience wrapper.

---

## `FusionWeights`

```python
@dataclass
class FusionWeights:
    rna: float = 1.0
    atac: float = 0.0
    protein: float = 0.0
```

- Auto-normalises to sum to 1.0 (`FusionWeights(rna=1, atac=1)` → `rna=0.5, atac=0.5`).
- Zero weights are legal (useful to drop a modality without removing the input).
- Negative weights raise.

`.as_dict()` returns `{"rna": ..., "atac": ..., "protein": ...}`.

---

## `MultiOmicsFusion`

```python
fusion = MultiOmicsFusion()
```

**Methods:**

### `.fuse(rna=None, atac=None, protein=None, weights=None) -> FusionResult`

Weighted fusion on the intersection of sample indices and pathway columns. At least one modality must be supplied. `weights=None` means equal non-zero weights across supplied modalities.

### `.learn_weights(labels, rna=None, atac=None, protein=None, grid_step=0.1) -> FusionWeights`

Simplex grid-search optimising leave-one-out 1-NN cell-type classification accuracy. `grid_step=0.1` → 11 × 11 candidate triples (fast; ~seconds on cohorts < 1k). Smaller step = finer search, more compute.

---

## `FusionResult`

| Attribute | Type | Description |
|---|---|---|
| `.fused` | `DataFrame(samples × pathways)` | The fused score matrix |
| `.weights` | `FusionWeights` | Weights actually used |
| `.per_modality` | `Dict[str, DataFrame]` | Aligned per-modality matrices (post-intersection) |
| `.union_pathways` | `List[str]` | Pathway columns present in `.fused` |
| `.to_dict()` | — | JSON-serialisable summary |

---

## `flag_discordant_pathways`

```python
report = flag_discordant_pathways(
    rna: pd.DataFrame,           # samples × pathways
    protein: pd.DataFrame,       # samples × pathways
    magnitude_threshold: float = 0.5,
    direction_threshold: float = 0.0,
)
```

Per-pathway classification into:

- **concordant** — RNA and protein agree in sign and within `magnitude_threshold`.
- **magnitude-discordant** — same sign, difference > `magnitude_threshold`.
- **direction-discordant** — opposite signs (after `direction_threshold`).

**`DiscordanceReport`:**

| Attribute / method | Type | Description |
|---|---|---|
| `.concordant_pathways` | `List[str]` | — |
| `.discordant_pathways` | `List[str]` | magnitude- or direction-discordant |
| `.direction_discordant` | `List[str]` | direction-discordant only (the biologically "interesting" ones) |
| `.per_pathway` | `DataFrame` | Per-pathway classification + deltas |
| `.summary()` | `str` | — |
| `.to_dict()` | `dict` | — |

---

## Roadmap acceptance (F10)

- Fused − RNA-only 1-NN classification accuracy uplift ≥ 3 pp, with 95% bootstrap CI lower bound > 0.
- Packaged run: 10x `pbmc_1k_protein_v3` CITE-seq (630 gated cells). **56.5% → 79.5%**, uplift +23.0 pp, CI +18.1..+27.6 pp. JSON: [`results/f10_validation/fusion_uplift.json`](../../results/f10_validation/fusion_uplift.json).
- Reproduce: `python scripts/validate_f10_real_data.py`.

---

## See also

- [Guide: multi-omics](../guides/multi-omics.md) — beginner-facing walkthrough
- [`cross_modal_validation`](cross_modal_validation.md) — Validation Gate 5 on the fused scores
- [`embed`](embed.md) — `embed_joint` for Nicheformer-style cross-modality joint embedding upstream of fusion
