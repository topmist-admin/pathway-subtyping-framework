# In-Silico Perturbation Guide

*PSF v0.6 Phase 2 — Foundation-Model Interface, F5*

PSF's observed-data pipeline answers "what is the molecular state of this
cell right now?" The `pathway_subtyping.perturb` layer adds the
counterfactual question: *"what would this cell's MSV look like if I
knocked out gene X?"* This is central to hypothesis generation in
target validation, cell-line QC, and disease modelling.

See [examples/notebooks/23_perturbation.ipynb](../../examples/notebooks/23_perturbation.ipynb)
for a runnable walkthrough.

---

## Architecture

```
expression_df + gene  →  GeneformerPerturber  →  perturbed embedding
                                                        │
                                                        ▼
                                                 MSVFromEmbedding
                                                        │
                                                        ▼
                                               delta-MSV + (optional)
                                               ConformalPathwayPredictor
```

Three components plus an optional uncertainty wrap:

| Component | Role | Defaults |
|---|---|---|
| `GeneformerPerturber` | Expression → perturbed embedding | `FallbackPerturber` if no backend specified |
| `MSVFromEmbedding` | Embedding → pathway-score vector | Ridge-regression head (MLP variant is API-compatible) |
| `PerturbationScreen` | Batch runner + ranking | L2-norm ordering of per-gene delta-MSV |
| `ConformalPathwayPredictor` (F1) | Uncertainty wrap | Split-conformal on the head |

---

## Choosing a backend

| Backend | When to use | Install |
|---|---|---|
| `OfficialBackend` | Production. Biology-grade predictions. | `pip install pathway-subtyping[perturb]` + cached Geneformer checkpoint |
| `FallbackPerturber` | Tests, CI without GPU, interface smoke-runs. | Always available — pure numpy |

`FallbackPerturber` models a knockout as zeroing that gene's column
before PCA-projecting to the embedding space, and over-expression as
amplifying the column. Directional behaviour on well-separated marker
genes is correct; absolute magnitudes are meaningless.

For research claims, always use `OfficialBackend`. The fallback exists so
that the package API works out of the box without heavyweight downloads.

---

## Minimal end-to-end example

```python
from pathway_subtyping.perturb import (
    FallbackPerturber, GeneformerPerturber,
    MSVFromEmbedding, PerturbationScreen, PerturbationReport,
)

# 1. Fit the MSV head on reference (embedding, pathway-score) pairs.
perturber = GeneformerPerturber()             # FallbackPerturber by default
baseline_emb = perturber.embed(expression_df)
head = MSVFromEmbedding().fit(baseline_emb, reference_pathway_scores)

# 2. Run a screen.
screen = PerturbationScreen(perturber, head)
result = screen.run(
    expression_df,
    gene_panel=["MYC", "MECP2", "TP53", "KRAS"],
    mode="knockout",
)

# 3. Rank + diagnose.
report = PerturbationReport.from_screen(result)
print(report.summary())
print(report.top_k(5))
```

---

## Acceptance: directional master-regulator signatures

The roadmap acceptance criterion is that perturbing a known master
regulator produces the *expected direction* of MSV shift. PSF makes
this testable via `PerturbationReport.check_directional_signature`:

```python
signature = {
    "MYC":   {"HALLMARK_MYC_TARGETS_V1": -1},   # KO expected to lower
    "MECP2": {"HALLMARK_NEURON_MARKERS": -1},
    "TP53":  {"HALLMARK_P53_PATHWAY": -1},
}
check = report.check_directional_signature(signature)
# Each row: expected_sign, observed_sign, delta_value, passed, reason
```

Use this in CI against the real `OfficialBackend` for regression
testing — bit-exact embedding parity catches silent model drift, and
this directional check catches coverage drift on your master-regulator
panel.

---

## Uncertainty-aware screens

Wrap the MSV head with the F1 `ConformalPathwayPredictor` to report
prediction intervals on the perturbed MSV:

```python
from pathway_subtyping.uncertainty import ConformalPathwayPredictor

def score_fn(emb):
    return head.transform(emb)["HALLMARK_P53_PATHWAY"].to_numpy()

predictor = ConformalPathwayPredictor(score_fn=score_fn, coverage=0.9)
predictor.calibrate(baseline_emb[cal_idx], pathway_scores[pathway][cal_idx])
# intervals on perturbed embeddings inherit the calibrated coverage
intervals = predictor.predict(result.per_gene_results["TP53"].perturbed_embedding)
```

See the F1 guide ([docs/guides/uncertainty.md](uncertainty.md)) for the
full conformal API and its coverage guarantees.

---

## Runtime expectations

| Setup | 2000 cells × 100 genes | Notes |
|---|---|---|
| CPU + `FallbackPerturber` | seconds | Interface harness only |
| CPU + `OfficialBackend` (cached embeddings) | ~4 hours | Batch-parallel; fine for overnight runs |
| A100 + `OfficialBackend` | < 30 min | Roadmap acceptance target |

Perturbations reuse the baseline embedding across all genes in the panel
— only the single gene's expression changes per run — so screen cost
scales linearly in `|gene_panel|` and not in `n_cells × n_genes`.

---

## Custom backends (advanced)

If you want to swap in a different Geneformer variant (fine-tuned
checkpoint, custom tokenizer, accelerated ONNX runtime), subclass
`GeneformerBackend`:

```python
from pathway_subtyping.perturb import GeneformerBackend

class MyCustomBackend(GeneformerBackend):
    embedding_dim = 256

    def embed(self, expression):
        ...

    def perturb(self, expression, gene, mode):
        ...

perturber = GeneformerPerturber(MyCustomBackend())
```

Custom backends must pass the golden-test parity harness (fixtures under
`tests/fixtures/perturb/`) before being used for research claims —
silent model drift is the single biggest source of downstream surprise
in a perturbation pipeline.

---

## Acceptance criteria (v0.6 roadmap)

The layer is considered complete when:

- [x] `GeneformerPerturber` wraps a pluggable backend (official or
      fallback) with a stable `.embed()` / `.perturb()` API.
- [x] `MSVFromEmbedding` translates embeddings to pathway scores with
      a closed-form head (MLP variant API-compatible).
- [x] `PerturbationScreen.run(...).rank()` returns a per-gene ranked
      delta-MSV panel.
- [x] `PerturbationReport.check_directional_signature` exposes the
      master-regulator directional check as a first-class test hook.
- [x] Conformal integration with F1 verified end-to-end.
- [ ] Production benchmark: 2,000 cells × 100 genes under 30 min on
      A100 with `OfficialBackend` (gated on checkpoint distribution).
- [ ] Bit-parity golden-test harness against frozen fixtures
      (`tests/fixtures/perturb/`) for the official backend.

The last two items land when release engineering packages the
checkpoint; the infrastructure and API they test are already in place.

---

## See also

- [examples/notebooks/23_perturbation.ipynb](../../examples/notebooks/23_perturbation.ipynb)
  — runnable walkthrough
- [docs/guides/uncertainty.md](uncertainty.md) — F1 conformal layer
- [docs/roadmap-v06-codeberg.md](../roadmap-v06-codeberg.md) — F5 spec
- Theodoris CV et al. (2023). *Transfer learning enables predictions
  in network biology.* Nature 618:616-624.
