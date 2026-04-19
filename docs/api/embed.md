# `embed` — Single-cell embedding wrappers (F6, F8)

Backend-agnostic cell embedders (scGPT, Nicheformer) with a deterministic PCA fallback and a content-hashed on-disk cache. The embedding is the anchor used by F2 harmonization, F5 perturbation, F11 causal, and more.

**Import path:** `pathway_subtyping.embed`
**Source:** [`src/pathway_subtyping/embed/`](../../src/pathway_subtyping/embed/)
**Guide:** [embeddings.md](../guides/embeddings.md) · [Notebook 24 (scGPT)](../../examples/notebooks/24_embeddings.ipynb) · [Notebook 26 (Nicheformer)](../../examples/notebooks/26_spatial_joint.ipynb)
**Install:** `pip install pathway-subtyping[embed]` (optional — fallbacks work without it)

---

## Public surface

```python
from pathway_subtyping.embed import (
    # ABC + result type
    Embedder, EmbeddingResult,
    # scGPT — F6
    scGPTEmbedder, OfficialSCGPTBackend, FallbackSCGPTEmbedder,
    # Nicheformer — F8 (dissociated + spatial joint)
    NicheformerEmbedder, OfficialNicheformerBackend, FallbackNicheformerEmbedder,
    # Cache
    EmbeddingCache, cache_key_for,
)
```

---

## `Embedder` (ABC)

```python
class Embedder:
    embedding_dim: int
    backend_id: str
    def embed(self, expression: pd.DataFrame) -> EmbeddingResult: ...
```

Every backend exposes `.embed(DataFrame) -> EmbeddingResult`. Stable public contract; custom backends that pass the parity harness (see `tests/test_embed*.py`) can drop in.

---

## `EmbeddingResult`

| Attribute | Type | Description |
|---|---|---|
| `.embeddings` | `(n_cells, embedding_dim) ndarray` | The per-cell embedding |
| `.cell_index` | `pd.Index` | Original cell IDs preserved from input DataFrame |
| `.backend` | `str` | e.g. `'scgpt:official/v1.2.0'` or `'scgpt:fallback'` |
| `.meta` | `Dict[str, Any]` | Backend-specific metadata (tokenizer version, checkpoint path, etc.) |
| `.n_cells` / `.embedding_dim` | `int` | Convenience properties |

---

## `scGPTEmbedder` (F6)

```python
scGPTEmbedder(
    checkpoint: str = "scgpt-whole-human-v1.2.0",
    device: str = "cpu",
    batch_size: int = 32,
    fallback_when_unavailable: bool = True,
)
```

Attempts to load the real scGPT model; if the `[embed]` extra or the checkpoint is missing and `fallback_when_unavailable=True`, transparently falls back to `FallbackSCGPTEmbedder`. Set `fallback_when_unavailable=False` to surface `ImportError` / `FileNotFoundError` explicitly.

```python
embedder = scGPTEmbedder()
result = embedder.embed(expression)          # EmbeddingResult
```

---

## `FallbackSCGPTEmbedder`

Deterministic PCA on the input. Same `.embed(...) -> EmbeddingResult` contract; `.backend_id = "scgpt:fallback"`. Used by default in CI and tests.

```python
FallbackSCGPTEmbedder(embedding_dim: int = 64, seed: int = 0)
```

---

## `NicheformerEmbedder` (F8)

For paired dissociated + spatial transcriptomics. Produces a joint embedding space so cells from both modalities can be compared without a harmonisation step afterward.

```python
NicheformerEmbedder(
    checkpoint: str = "nicheformer-v1",
    device: str = "cpu",
    fallback_when_unavailable: bool = True,
)

# Single-modality embed, or...
result = embedder.embed(dissociated_expression)

# ...paired embed (recommended for downstream F10 cross-modal use):
from pathway_subtyping.embed import embed_joint
joint = embed_joint(
    embedder,
    dissociated=dissociated_df,
    spatial=spatial_df,
)
# joint.embeddings: (n_dissociated + n_spatial, d); joint.meta['modality']: array of 'dissociated' / 'spatial'
```

---

## `EmbeddingCache`

Filesystem cache for `EmbeddingResult` objects. The same cache powers the F5 `OfficialBackend(cache_dir=...)` pattern under the hood.

```python
from pathway_subtyping.embed import EmbeddingCache, cache_key_for

cache = EmbeddingCache(cache_dir="~/.cache/pathway-subtyping/embed")

key = cache_key_for(
    backend="scgpt:official/v1.2.0",   # include checkpoint version — upgrades invalidate automatically
    expression=expression_df,
    extra={"batch_size": 32},          # optional extras bake into the key
)

if cache.has(key):
    result = cache.get(key)
else:
    result = embedder.embed(expression_df)
    cache.put(key, result)
```

**Key stability:** the SHA-256 covers `backend + column order + row order + float64 bytes + extras JSON`. Same expression, same backend, same extras → same key. Any change anywhere → cache miss → fresh inference.

**Storage:** one `.npy` + one `.json` sidecar per entry. Delete per-entry with `cache.delete(key)` or `cache.clear()`.

---

## `cache_key_for`

```python
cache_key_for(
    backend: str,
    expression: pd.DataFrame,
    extra: Optional[Dict[str, Any]] = None,
) -> str
```

Exposed as a stand-alone helper so you can use the cache key scheme outside `EmbeddingCache` (e.g. for downstream keyed artefacts).

---

## Backend identifier convention

`{model_name}:{backend_kind}/{version}`:

- `scgpt:official/v1.2.0`
- `scgpt:fallback`
- `nicheformer:official/v1`
- `geneformer:official:Geneformer-V2-104M:cls:maxlen4096` (F5 `OfficialBackend._backend_id()` — includes config)

Changing any component produces a different cache key automatically. No manual invalidation needed across upgrades.

---

## See also

- [Guide: embeddings](../guides/embeddings.md) — beginner-facing walkthrough
- [`harmonize`](harmonize.md) — uses `Embedder` output as the alignment anchor
- [`perturb`](perturb.md) — uses `OfficialBackend` (Geneformer) + `MSVFromEmbedding` on top of embeddings
- [`omics`](omics.md) — joint embedding via `embed_joint` feeds cross-modal fusion
