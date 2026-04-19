# `perturb` — In-silico gene perturbation (F5)

Geneformer-backed perturbation of PSF pathway-score predictions via a Molecular-State-Vector (MSV) head and a panel-level screen driver.

**Import path:** `pathway_subtyping.perturb`
**Source:** [`src/pathway_subtyping/perturb/`](../../src/pathway_subtyping/perturb/)
**Guide:** [perturbation.md](../guides/perturbation.md) · [Notebook 23](../../examples/notebooks/23_perturbation.ipynb)
**Install:** `pip install pathway-subtyping[perturb]` (optional — `FallbackPerturber` works without it)

---

## Public surface

```python
from pathway_subtyping.perturb import (
    # High-level wrapper
    GeneformerPerturber, PerturbationMode, PerturbationResult,
    # Backends
    GeneformerBackend, OfficialBackend, FallbackPerturber,
    # MSV head + screen
    MSVFromEmbedding, PerturbationScreen, ScreenResult, PerturbationReport,
)
```

---

## `GeneformerPerturber` (high-level)

```python
GeneformerPerturber(backend: Optional[GeneformerBackend] = None)
```

Default backend is `FallbackPerturber` (PCA-based, no deps). For real Geneformer inference, pass `OfficialBackend(...)`.

**Methods:**

| Method | Returns | Notes |
|---|---|---|
| `.embed(expression)` | `np.ndarray (n_cells, d)` | Baseline embeddings |
| `.perturb(expression, gene, mode="knockout")` | `PerturbationResult` | In-silico knockout or overexpress |

**`PerturbationResult`:**

| Attribute | Type | Description |
|---|---|---|
| `.gene` | `str` | The gene that was perturbed |
| `.mode` | `PerturbationMode` | `KNOCKOUT` or `OVEREXPRESS` |
| `.baseline_embedding` | `(n, d)` | Unperturbed embedding |
| `.perturbed_embedding` | `(n, d)` | Embedding after perturbation |
| `.delta_embedding` | `(n, d)` | perturbed − baseline |
| `.l2_impact_per_cell` | `(n,)` | Per-cell L2 of delta |
| `.mean_l2_impact` | `float` | Mean L2 across cells |

---

## `OfficialBackend` (production)

Real Geneformer V2 104M wrapper — rank tokenisation + CLS-token embeddings + direct knockout-by-zero-count.

```python
OfficialBackend(
    model_directory: Optional[str] = None,    # path to Geneformer-V2-104M dir
    device: str = "cpu",
    forward_batch_size: int = 8,
    max_input_len: int = 4096,
    emb_mode: str = "cls",                    # or "mean"
    cache_dir: Optional[str] = None,          # on-disk content-hashed embedding cache
)
```

**Checkpoint discovery:** clone `https://huggingface.co/ctheodoris/Geneformer` locally and point `model_directory` at the `Geneformer-V2-104M` subdir, or set `GENEFORMER_MODEL_DIR` env var before running [`scripts/validate_f5_real_data.py`](../../scripts/validate_f5_real_data.py).

**Cache:** with `cache_dir=` set, CLS embeddings are persisted to disk keyed by SHA-256 over `(checkpoint_basename + emb_mode + max_input_len + expression bytes)`. Reruns on unchanged input return in sub-millisecond time. Default cache path in the validation script: `~/.cache/pathway-subtyping/geneformer`. Verified 12 000× speedup (2.49 s cold → 0.0002 s warm on a 5-cell toy).

**Embedding dim:** `768` (inferred from checkpoint `hidden_size`; auto-updates on load).

---

## `FallbackPerturber` (fallback)

Deterministic PCA-based substitute. API-compatible with `OfficialBackend` — tests, CI, and offline runs use it without the heavy deps.

```python
FallbackPerturber(embedding_dim: int = 64, overexpression_factor: float = 3.0, seed: int = 0)
```

- Knockout: zero the gene's column before projecting into the PCA basis.
- Overexpress: multiply by `overexpression_factor`.
- Directional behaviour on well-separated markers is correct; absolute magnitudes are not biology-calibrated.

---

## `MSVFromEmbedding`

Ridge-regression head mapping embeddings → per-pathway MSV scores.

```python
head = MSVFromEmbedding(ridge_alpha: float = 1e-2)
head.fit(embeddings, pathway_scores)       # pathway_scores = samples × pathways DataFrame
head.transform(new_embeddings)             # returns DataFrame(samples × pathways)
head.fit_transform(embeddings, pathway_scores)
head.delta(baseline_emb, perturbed_emb)    # returns per-pathway ΔMSV
```

**Floor:** needs ≥ 5 samples to fit. On small cohorts (e.g. n=6 in GSE123753 Neu samples), fit on a superset (WT+KO pooled) and predict on the subset; the real-data validator demonstrates this pattern.

---

## `PerturbationScreen`

Run a panel of perturbations and rank genes by MSV impact.

```python
screen = PerturbationScreen(perturber, msv_head, pathway_target="HALLMARK_MYC_TARGETS_V1")
result = screen.run(expression, gene_list=["MYC", "TP53", "E2F1", "CDK1", "CCNE1"])
# ScreenResult: per-gene ranked by |ΔMSV| on the target pathway
```

**`ScreenResult`:**
- `.ranking` — DataFrame with columns `gene`, `delta_msv`, `abs_impact`, sorted by absolute impact
- `.per_gene_details` — per-gene raw `PerturbationResult` for inspection
- `.to_dict()` — JSON-serialisable

---

## `PerturbationReport`

Narrative wrapper around a `ScreenResult` + a directional-signature check (does KO of a known master regulator move its target pathway in the expected direction?).

```python
report = PerturbationReport.from_screen(result, directional_expectations={"MYC": -1})
print(report.summary())      # markdown summary
print(report.directional_agreement_rate)
```

---

## Validation scripts

- `scripts/validate_f5_real_data.py --backend {fallback, geneformer}` — runs a TCGA-COAD in-silico directional test + a GSE123753 (Boxer 2020 MECP2-KO iPSC neurons) real WT-vs-KO comparison, writes `perturbation_directional.json` (Fallback) or `perturbation_directional_geneformer.json` (Geneformer).
- Geneformer real-data result: **50/50 pathways directionally agree** on GSE123753, Spearman predicted-vs-observed ΔMSV rho = **+0.85**.
- Tests: [`tests/test_perturb_real_data.py`](../../tests/test_perturb_real_data.py).

---

## See also

- [Guide: perturbation](../guides/perturbation.md) — beginner-facing walkthrough
- [`uncertainty`](uncertainty.md) — wrap `MSVFromEmbedding.transform` with `ConformalPathwayPredictor` for calibrated intervals on perturbed MSV predictions
- [CHANGELOG](../../CHANGELOG.md) — v0.6.3 release notes and acceptance numbers
