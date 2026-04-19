# Multi-Omics Fusion (F10)

> v0.6 Phase 3 F10 — pathway-score fusion across RNA + ATAC + proteomics, plus an explicit discordance flag for pathways where the modalities disagree.

**No extra required** — the `omics` layer is part of the core install.

**When to use:** you have pathway scores from ≥2 modalities (e.g., RNA-seq + CITE-seq ADT, or bulk RNA + ATAC) and want either (a) a single fused score per pathway for downstream clustering, or (b) an explicit signal showing *where* RNA and protein disagree (which is often the biologically interesting answer, not noise).

**API cheat sheet:**

| Class / function | Role |
|---|---|
| `ATACScorer` | Score pathways from a samples × peaks ATAC matrix. |
| `ProteomicsScorer` | Score pathways from a samples × proteins abundance matrix. |
| `MultiOmicsFusion.fuse(...)` | Weighted fusion across RNA + ATAC + protein. Auto-intersects on samples + pathways. |
| `MultiOmicsFusion.learn_weights(...)` | Grid-search fusion weights that maximise 1-NN accuracy against a labelled reference. |
| `FusionWeights` | Per-modality weights; auto-normalises to sum to 1. |
| `flag_discordant_pathways(...)` | Per-pathway RNA-vs-protein concordance check; returns a `DiscordanceReport`. |

---

## 5-minute example

```python
import pandas as pd
from pathway_subtyping.omics import (
    ATACScorer, ProteomicsScorer, MultiOmicsFusion, FusionWeights,
    flag_discordant_pathways,
)

# Inputs (samples × features per modality — your usual normalised matrices)
rna = pd.read_csv("rna_expression.csv", index_col=0)
atac = pd.read_csv("atac_accessibility.csv", index_col=0)
protein = pd.read_csv("protein_abundance.csv", index_col=0)

pathways = {
    "TCELL_ACTIVATION": ["CD3D", "CD3E", "CD4", "IL2RA"],
    "BCELL": ["CD19", "MS4A1", "CD79A"],
    # ...
}

# 1) Score each modality — same API shape, one row per sample, one column per pathway
rna_scores = score_rna_pathways(rna, pathways)          # existing v0.5 scorer
atac_scores = ATACScorer().score(atac, pathways)
protein_scores = ProteomicsScorer().score(protein, pathways)

# 2) Fuse at equal weight
fusion = MultiOmicsFusion()
result = fusion.fuse(
    rna=rna_scores,
    atac=atac_scores,
    protein=protein_scores,
    weights=FusionWeights(rna=0.5, atac=0.25, protein=0.25),
)
fused = result.fused   # samples × pathways — plug into GMM clustering / MSV

# 3) Flag pathways where RNA and protein disagree (often the interesting ones)
disc = flag_discordant_pathways(rna=rna_scores, protein=protein_scores)
print(disc.summary())
print(disc.discordant_pathways)  # list of pathway names
```

---

## Learning fusion weights against a labelled reference

If you have a paired reference with cell-type labels (e.g., a CITE-seq dataset), let the fusion learn weights that maximise 1-NN accuracy:

```python
labels = pd.Series(cell_types, index=rna_scores.index)

weights = fusion.learn_weights(
    labels=labels,
    rna=rna_scores,
    atac=atac_scores,
    protein=protein_scores,
    grid_step=0.1,            # grid over the simplex; 0.1 -> 11 x 11 candidates
)
print(weights)  # FusionWeights(rna=0.6, atac=0.0, protein=0.4) — auto-normalised

result = fusion.fuse(rna=rna_scores, atac=atac_scores, protein=protein_scores, weights=weights)
```

`learn_weights` does a leave-one-out 1-NN search on the simplex, so small cohorts (≤100 labelled cells) run in seconds.

---

## Real-data acceptance run (roadmap F10)

The shipped validation script runs this flow against the 10x Genomics `pbmc_1k_protein_v3` CITE-seq reference (713 cells, 17-antibody panel):

```bash
python scripts/validate_f10_real_data.py
```

Outcome on the packaged data (see [`results/f10_validation/fusion_uplift.json`](../../results/f10_validation/fusion_uplift.json) after running):

| | RNA-only | Fused (RNA + protein) |
|---|---|---|
| 5-fold 1-NN cell-type accuracy | 56.5% | **79.5%** |
| Bootstrap 95% CI on uplift | — | +18.1 .. +27.6 pp |

Gate: uplift ≥ 3 pp with CI lower bound > 0. Passes by a wide margin.

---

## Interpreting discordance

`flag_discordant_pathways` returns one of three statuses per pathway: **concordant** (RNA and protein agree in magnitude + sign), **magnitude-discordant** (same direction, different strength), or **direction-discordant** (opposite sign — classic post-transcriptional regulation, translation control, protein stability, etc.).

In practice, **direction-discordant pathways are rarely noise**. They are a short list of biological candidates for follow-up: the pathway is transcribed but not translated, or vice versa. Use `.discordant_pathways` as a shortlist, then inspect the per-pathway per-sample scores.

---

## Known gotchas

- **Modality inputs must share sample index and pathway columns.** `MultiOmicsFusion.fuse` auto-restricts to the intersection of samples AND pathways across all supplied modalities. If the intersection is empty, it raises — typically a sign that your pathway definitions differ across modalities.
- **Protein pathways use gene-symbol membership by default.** If your protein matrix is keyed by UniProt IDs, pass a `gene_to_protein` mapping to `ProteomicsScorer`.
- **ATAC scoring assumes peaks are annotated to genes.** `ATACScorer` reduces peaks → gene by summing/averaging over annotated peaks, then mean-Zs the gene-level scores against each pathway. Bring your own peak-to-gene annotation (bedtools intersect or similar); no peak-calling is done inside the module.
- **`FusionWeights` auto-normalises.** Passing `FusionWeights(rna=1.0, atac=1.0)` is equivalent to `FusionWeights(rna=0.5, atac=0.5)`. Zero weights are legal (useful to drop a modality without deleting the input).

---

## See also

- [API reference: `omics`](../api/omics.md) — full `ATACScorer` / `ProteomicsScorer` / `MultiOmicsFusion` / `flag_discordant_pathways` signatures
- [Notebook 28 — multi_omics](../../examples/notebooks/28_multi_omics.ipynb) — paired CITE-seq walkthrough
- [docs/api/cross_modal_validation.md](../api/cross_modal_validation.md) — Validation Gate 5 for fused scores
- [CHANGELOG F10 validation entry](../../CHANGELOG.md) — real-cohort acceptance numbers
