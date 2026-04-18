# Cross-Platform Harmonization Guide

*PSF v0.6 Phase 1 — Rigor Layer, F2*

PSF pathway scores computed on different single-cell platforms (10x
Chromium, Smart-seq2, bulk RNA-seq, spatial transcriptomics) are
directionally comparable after standard preprocessing but drift in
magnitude: platform-specific detection biases couple with cell biology
in ways that z-scoring alone cannot remove. Users combining public
datasets across platforms have long reported a 0.55–0.65 per-pathway
correlation where ≥ 0.80 is expected.

The `pathway_subtyping.harmonize` module closes this gap with a
UCE-anchored per-platform linear aligner and a diagnostic suite that
flags low-confidence cells.

See [examples/notebooks/22_cross_platform.ipynb](../../examples/notebooks/22_cross_platform.ipynb)
for a runnable walkthrough.

---

## When to reach for the harmonize layer

| Scenario | Use |
|---|---|
| Combining public datasets spanning several platforms | Always — drift is real and measurable |
| Single-platform cohort, multiple batches | Use `pathway_subtyping.batch_correction` instead |
| Cross-species meta-analysis | Real UCE (not the fallback) via `[harmonize]` extra |
| Quick sanity check without UCE installed | `FallbackEmbedder` — deterministic, no heavyweight deps |

---

## Mental model

```
scores_platform_p  ≈  biology(cell)  +  platform_offset_p(cell)  +  noise
```

Two design choices follow:

1. **Get a platform-invariant estimate of `biology(cell)`.** In production
   this is `UCEEmbedder` (Rosen et al. 2024 — 1280-dim embeddings from a
   model pretrained on 36M cells across 33 species). For local testing,
   `FallbackEmbedder` provides a deterministic PCA substitute.
2. **Learn the per-platform linear map** from embedding → pathway scores,
   then express every cell in a chosen reference-platform frame. That map
   is what `CrossPlatformAligner` fits.

---

## API walkthrough

```python
from pathway_subtyping.harmonize import (
    CrossPlatformAligner,
    FallbackEmbedder,
    HarmonizationReport,
)

# 1. Get cell embeddings (real UCE in prod; fallback for testing)
embedder = FallbackEmbedder(embedding_dim=64).fit(reference_scores)
embeddings = embedder.embed(pooled_scores)

# 2. Fit + transform pathway scores
aligner = CrossPlatformAligner()  # optional: reference_platform="10x"
result = aligner.fit_transform(
    pooled_scores,           # DataFrame (n_cells x n_pathways)
    platforms,               # sequence[str] length n_cells
    embeddings,              # (n_cells, embedding_dim) array
)

# 3. Diagnostics
report = HarmonizationReport.from_alignment(result)
print(report.summary())
# HarmonizationReport(n=800, mean_confidence=0.497, ...)
```

### Choosing the reference platform

By default, `CrossPlatformAligner` picks the platform with the most cells
as the reference. Override with `CrossPlatformAligner(reference_platform="10x")`
if you have a canonical frame — e.g., the platform that most of your
downstream analysis was historically done on.

### Platforms unseen at fit time

Cells whose platform was not in the fit set pass through unshifted, with
a single warning in the log. In practice, re-fit the aligner any time you
add a new platform.

### Diagnostics to report

- `report.mean_platform_drift` — one number per platform summarising how
  far that platform's raw scores deviated from the reference on average.
- `report.confidence` — per-cell `1 / (1 + z-scored shift)`; values near 1
  mean the cell's raw scores needed little correction.
- `report.correlate_with_quality(quality_series)` — sanity check that
  low-confidence cells are not random; pass a per-cell quality covariate
  (read depth, mitochondrial fraction, etc.) and read the Spearman rho.

---

## Acceptance criteria (v0.6 roadmap)

The layer is considered complete when:

- [x] Harmonized pathway-level rho across platform pairs exceeds 0.75 on
      the synthetic cross-platform benchmark (pre-harmonization baseline
      0.55–0.65 typical).
- [x] Harmonization confidence correlates with a per-cell quality covariate
      (tested via shift-magnitude proxy in `test_harmonize.py`).
- [x] No degradation of within-platform performance — single-platform
      inputs keep drift ≈ 0 after `fit_transform`.
- [x] Cross-platform benchmark runs on 4 platforms × 1 biological system
      (synthetic cortex-shaped) end-to-end in `CrossPlatformBenchmark`.

### Synthetic acceptance harness

`CrossPlatformBenchmark` applies deterministic platform distortions to a
reference pathway-score matrix, fits the aligner, and reports pre- and
post-harmonization Spearman rho:

```python
from pathway_subtyping.harmonize import CrossPlatformBenchmark

bench = CrossPlatformBenchmark(
    reference_scores=my_reference_scores,
    platforms=["10x", "smartseq2", "bulk_rnaseq", "spatial"],
)
summary = bench.run_many(n_seeds=5)
# summary["post_rho_mean"] > 0.75 on a cohort with 3+ biological clusters
```

The benchmark simulates the ideal "UCE-provides-platform-invariant-embedding"
case by freezing an embedder fit on the undistorted reference. Real UCE
approximates this because it is pretrained across dozens of platforms
and does not leak platform signal into the embedding space. The
`FallbackEmbedder` substitute is PCA-only; on data where platform drift
exceeds biological variance, its embeddings will contain platform signal
and alignment quality drops — install the `[harmonize]` extra to use the
real UCE model when coupling is severe.

### Real-data evaluation (next iteration)

Matched cross-platform cohorts for 10x × Smart-seq2 cortex comparison are
not shipped with PSF. Users with matched data (e.g., Trujillo 2019 cortex
organoid subset × Gouwens 2020 Patch-seq) can reuse `CrossPlatformBenchmark`
by passing their own pathway-score DataFrame; the protocol is identical.

---

## Gotchas

- **Do not fit the embedder on pooled distorted data.** If your embedder
  is PCA-style (the fallback), pooled distorted data will contain platform
  variance that the embedder then surfaces as "biology" — alignment
  quality collapses. Always use a reference-fit (or UCE-pretrained)
  embedder and project distorted data through it.
- **Minimum 3 cells per platform at fit time.** Platforms below that
  threshold are skipped with a warning; transforms on those platforms
  fall back to the reference frame.
- **Check confidence distribution.** If the mean harmonization confidence
  is below ~0.3, either your platform distortion is too extreme for
  linear alignment (consider per-platform re-preprocessing) or your
  embedding does not capture the relevant biology.

---

## See also

- [examples/notebooks/22_cross_platform.ipynb](../../examples/notebooks/22_cross_platform.ipynb)
  — runnable walkthrough
- [docs/roadmap-v06-codeberg.md](../roadmap-v06-codeberg.md) — Phase 1 F2
- [docs/guides/uncertainty.md](uncertainty.md) — the F1 sibling module
- Rosen J et al. (2024). *Universal Cell Embeddings: A Foundation Model
  for Cell Biology.* Nature.
