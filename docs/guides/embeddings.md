# Cell Embeddings Guide

*PSF v0.6 Phase 2 — Foundation-Model Interface, F6*

PSF's primary output is **pathway scores** — interpretable,
z-normalised, version-pinned per-pathway numbers. For some downstream
analyses (rare cell detection, trajectory inference, cross-platform
harmonization) learned cell embeddings carry complementary signal.
The `pathway_subtyping.embed` layer provides thin, stable wrappers
around foundation models — starting with scGPT — so those analyses
don't need ad-hoc integrations.

See [examples/notebooks/24_embeddings.ipynb](../../examples/notebooks/24_embeddings.ipynb)
for a runnable walkthrough.

---

## Architecture

```
expression_df  →  scGPTEmbedder  →  EmbeddingResult
                       │                  │
                       │                  ▼
                       │         (n_cells, embedding_dim)
                       ▼
              OfficialSCGPTBackend (scGPT checkpoint)
                       OR
              FallbackSCGPTEmbedder (PCA substitute)
```

All backends implement the same `Embedder` interface — if a downstream
function accepts one, it accepts all. This is why the F2 cross-platform
aligner can swap UCE for scGPT (or vice versa) without glue code.

---

## Choosing a backend

| Backend | When to use | Install |
|---|---|---|
| `OfficialSCGPTBackend` | Production. Biology-grade embeddings. | `pip install pathway-subtyping[embed]` + cached scGPT checkpoint |
| `FallbackSCGPTEmbedder` | CI, tests, offline experimentation. | Always available — pure numpy |

The fallback is a deterministic PCA projection. It preserves the API
and is great for testing but is not a substitute for scGPT when making
biological claims.

---

## Minimal usage

```python
from pathway_subtyping.embed import scGPTEmbedder

embedder = scGPTEmbedder()                       # fallback by default
result = embedder.embed(expression_df)
# result.embeddings is (n_cells, embedding_dim)
# result.as_dataframe() aligns with expression_df.index
```

## Caching across reruns

Foundation-model inference is the slowest step in the PSF pipeline.
The content-hashed cache at `pathway_subtyping.embed.EmbeddingCache`
skips the expensive re-inference when nothing has changed.

```python
from pathway_subtyping.embed import EmbeddingCache, cache_key_for

cache = EmbeddingCache("~/.cache/pathway-subtyping/embed")
key = cache_key_for(backend=result.backend, expression=expression_df)

if cache.has(key):
    result = cache.get(key)
else:
    result = embedder.embed(expression_df)
    cache.put(key, result)
```

The key covers the backend identifier (which includes checkpoint
version) and the content of the expression matrix. Changing either
produces a fresh key; there is no manual cache invalidation to manage.

---

## Integrating with other PSF layers

### F2 cross-platform harmonization

`CrossPlatformAligner` accepts an arbitrary embedding array. Use any
`Embedder` backend:

```python
from pathway_subtyping.embed import scGPTEmbedder
from pathway_subtyping.harmonize import CrossPlatformAligner

embeddings = scGPTEmbedder().embed(expression_df).embeddings
aligned = CrossPlatformAligner().fit_transform(
    pathway_scores, platform_labels, embeddings,
)
```

### F5 perturbation screens

`MSVFromEmbedding` likewise accepts any embedding source. See the
[perturbation guide](perturbation.md) for a full pipeline example.

---

## Acceptance criteria (v0.6 roadmap)

- [x] Stable `Embedder` interface; scGPT is the first concrete backend.
- [x] Deterministic output on the fallback backend (same seed → same
      embedding across reruns).
- [x] Content-hashed cache invalidates on backend change or expression
      change; cross-run stable.
- [x] Integrates with F2 harmonize without glue code (tested in
      `test_embed.py::TestHarmonizeIntegration`).
- [ ] Official-backend bit parity against frozen fixtures (gated on
      scGPT checkpoint distribution).

---

## See also

- [examples/notebooks/24_embeddings.ipynb](../../examples/notebooks/24_embeddings.ipynb)
  — runnable walkthrough
- [docs/guides/cross-platform.md](cross-platform.md) — F2 integration
- [docs/guides/perturbation.md](perturbation.md) — F5 integration
- [docs/roadmap-v06-codeberg.md](../roadmap-v06-codeberg.md) — F6 spec
- Cui H et al. (2024). *scGPT: toward building a foundation model for
  single-cell multi-omics using generative AI.* Nature Methods.
