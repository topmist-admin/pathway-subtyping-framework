# `harmonize` — Cross-platform pathway-score alignment (F2)

Embedding-anchored per-platform alignment of pathway scores across 10x / Smart-seq2 / bulk RNA-seq / spatial transcriptomics.

**Import path:** `pathway_subtyping.harmonize`
**Source:** [`src/pathway_subtyping/harmonize/`](../../src/pathway_subtyping/harmonize/)
**Guide:** [cross-platform.md](../guides/cross-platform.md) · [Notebook 22](../../examples/notebooks/22_cross_platform.ipynb)
**Install:** `pip install pathway-subtyping[harmonize]` (optional — `FallbackEmbedder` works without it)

---

## Public surface

```python
from pathway_subtyping.harmonize import (
    # Embedders
    UCEEmbedder, FallbackEmbedder,
    # Aligner
    CrossPlatformAligner, AlignmentResult,
    # Reports + benchmarking
    HarmonizationReport, CrossPlatformBenchmark, simulate_platform_distortion,
)
```

---

## `UCEEmbedder`

Wraps the Universal Cell Embeddings model (Rosen et al. 2024) to produce platform-invariant 1280-dim per-cell embeddings.

```python
UCEEmbedder(
    checkpoint: str = "uce-33l-v2",
    device: str = "cpu",
    batch_size: int = 32,
)
```

- Requires `[harmonize]` extra + locally-cached UCE checkpoint. `_lazy_load()` raises a clear `ImportError` with install instructions if either is missing.
- CPU inference works but is slow; GPU recommended for >10k cells.

**Method:** `.embed(expression: pd.DataFrame) -> np.ndarray` of shape `(n_cells, 1280)`.

---

## `FallbackEmbedder`

Deterministic PCA-based embedder. Not a substitute for UCE on biological claims — it's an API harness so tests, CI, and offline development run uniformly without the checkpoint.

```python
FallbackEmbedder(embedding_dim: int = 64, seed: int = 0)
```

Fits lazily on the first `.embed()` call. Calls with the same DataFrame + same seed produce identical outputs (useful for reproducibility tests).

```python
emb = FallbackEmbedder(embedding_dim=32, seed=42)
emb.fit(expression_df)        # optional; auto-fits on first embed
X_emb = emb.embed(expression_df)   # (n_cells, 32)
```

---

## `CrossPlatformAligner`

Fits a per-platform linear correction on an embedding anchor, then maps cells from any platform into the reference platform's frame.

```python
CrossPlatformAligner(
    ridge_alpha: float = 1e-3,
    reference_platform: Optional[str] = None,   # defaults to the most-populated platform
)
```

**Workflow:**

```python
aligner = CrossPlatformAligner(reference_platform="10x")
aligner.fit(scores, platforms, embeddings)       # learn per-platform betas
result = aligner.transform(scores, platforms, embeddings)
# result.aligned_scores  — same shape as scores, but in reference frame
# result.per_cell_shift  — per-cell magnitude of correction
# result.per_platform_drift  — dict of per-platform offsets (for reporting)
```

Or one-shot:

```python
result = aligner.fit_transform(scores, platforms, embeddings)
```

**`AlignmentResult`:**

| Attribute | Type | Description |
|---|---|---|
| `.aligned_scores` | `DataFrame(n_cells, n_pathways)` | Drift-removed pathway scores in reference frame |
| `.platform` | `Series(n_cells)` | Platform label per cell (for downstream grouping) |
| `.per_cell_shift` | `Series(n_cells)` | Mean absolute shift per cell |
| `.per_platform_drift` | `Dict[str, Dict[str, float]]` | Per-platform per-pathway intercept offset |
| `.to_dict()` | — | JSON-serialisable summary |

**Accessors after fit:**
- `.platforms` — list of platforms seen at fit time
- `.pathway_names` — pathway columns from fit
- `.platform_offset(platform)` — `Series` of per-pathway offsets vs reference

**Unseen platforms at transform time** are left unchanged with a warning logged — useful for new cells from a platform that wasn't present at fit.

---

## `HarmonizationReport`

Per-cell confidence of the harmonisation (useful as a filter for downstream analysis).

```python
report = HarmonizationReport.from_alignment(
    result,
    pre_scores=raw_scores,     # optional: for pre-vs-post uplift
)

report.per_cell_confidence        # Series: 1.0 = high confidence, 0.0 = noisy
report.flagged_cells              # cells below the default threshold
print(report.summary())
```

---

## Benchmarking

`CrossPlatformBenchmark` runs a synthetic-distortion stress test end-to-end:

```python
from pathway_subtyping.harmonize import simulate_platform_distortion, CrossPlatformBenchmark

distorted = simulate_platform_distortion(
    reference_scores, platforms=["A", "B"], drift_strength=0.5, seed=0,
)
bench = CrossPlatformBenchmark(aligner=CrossPlatformAligner())
report = bench.run(reference=reference_scores, distorted=distorted, platforms=["A", "B"])
# report.pre_rho / post_rho / uplift across pathways
```

---

## Real-data acceptance (roadmap F2)

- Pathway-mean Spearman rho uplift ≥ 0.10 post-alignment.
- Aspirational target: post-alignment rho ≥ 0.75 on paired-cell matched cortex (requires paired-cell data; not enforced in CI).
- Packaged run: GSE28521 × GSE80655 — pre rho = −0.02, post rho = +0.52, uplift +0.55. JSON: [`results/f2_validation/harmonize_spearman.json`](../../results/f2_validation/harmonize_spearman.json).
- Reproduce: `python scripts/validate_f2_real_data.py`.

---

## See also

- [Guide: cross-platform harmonization](../guides/cross-platform.md) — beginner-facing walkthrough
- [`embed`](embed.md) — scGPT / Nicheformer embedders that can feed the aligner's embedding anchor
- [`uncertainty`](uncertainty.md) — use `HarmonizationReport.per_cell_confidence` as an uncertainty signal into F12 active selection
