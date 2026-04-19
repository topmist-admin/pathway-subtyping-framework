# Gene-Set Expansion (F7)

> v0.6 Phase 2 F7 — suggest candidate additions to a hand-curated seed gene set using either a DNA-sequence-predicted regulatory model (Borzoi) or a co-expression fallback.

**Install:** `pip install pathway-subtyping[genesets]`

**When to use:** you have a small hand-curated gene set (5–30 genes) and want literature-agnostic candidate additions. The module does not replace curation — it ranks candidates by similarity to the seed and flags membership in known databases so you can curate informed additions.

**API cheat sheet:**

| Class | Role |
|---|---|
| `RegulatoryGeneSetExpander` | User-facing entry point — wraps a backend and runs the expansion. |
| `CoexpressionBackend` | Fallback backend. Scores similarity as Spearman rho on a user-supplied expression matrix. Always available. |
| `BorzoiBackend` | Opt-in backend. Delegates similarity to predicted regulatory activity. Requires a Borzoi / Enformer checkpoint. |
| `ExpansionResult` | Container returned by `.expand()`. Has `.as_dataframe()`, `.top(n)`, `.recall_against(ground_truth)`. |

---

## 5-minute example (co-expression fallback)

```python
import pandas as pd
from pathway_subtyping.genesets import (
    CoexpressionBackend, RegulatoryGeneSetExpander,
)

# You already have a samples × genes expression matrix. Any scale works —
# the backend computes per-gene Spearman rho over samples internally.
expr = pd.read_csv("your_expression_matrix.csv", index_col=0)

# Seed: a MYC-related curated set
seed = ["MYC", "MYCN", "MAX", "MXI1"]

backend = CoexpressionBackend(expression=expr)
expander = RegulatoryGeneSetExpander(backend=backend)

result = expander.expand(
    seed_genes=seed,
    candidate_genes=expr.columns.tolist(),  # pool to rank
    top_n=20,
)

print(result.as_dataframe().head())
# columns: gene, score, in_curated_db, source_1, source_2, source_3
```

The result's `source_1..source_k` columns name the seed genes contributing most to each candidate's score — a simple "why was this suggested" signal.

---

## Production backend (Borzoi)

```python
from pathway_subtyping.genesets import BorzoiBackend, RegulatoryGeneSetExpander

backend = BorzoiBackend(checkpoint="borzoi-human-v1.0")
expander = RegulatoryGeneSetExpander(backend=backend)
result = expander.expand(seed_genes=seed, candidate_genes=candidates, top_n=20)
```

Borzoi requires the `[genesets]` extra plus a locally-cached checkpoint. The `CoexpressionBackend` is intentionally kept as the deterministic fallback so CI runs and offline development don't depend on a heavyweight model.

---

## Flagging curated-database membership

If a candidate already appears in a known curated database, you usually want a weaker "gentle reminder" rather than a strong recommendation. Pass `known_members` to flag these:

```python
known = {
    "MSigDB_MYC_TARGETS": {"BRCA1", "CCNB1", "CDK4", ...},
    "KEGG_CELL_CYCLE": {"CCND1", "CDK4", "RB1", ...},
}

result = expander.expand(
    seed_genes=seed,
    candidate_genes=candidates,
    known_members=known,
    top_n=20,
)

df = result.as_dataframe()
# df['in_curated_db'] == True for candidates present in any listed DB
novel = df[~df["in_curated_db"]]
```

A common pattern is: take the top-20 by score, split into `in_curated_db == True` (confirmatory picks) vs. `False` (novel suggestions), and write them both into your curated-set spreadsheet with a provenance note.

---

## Validating against a ground truth

If you have a held-out "gold" set (e.g., an orthogonal pathway like MSigDB MYC_TARGETS_V1 that your seed is designed to approximate), compute top-N recall:

```python
recall_at_20 = result.recall_against(ground_truth=gold_set)
```

The roadmap target is ≥30% top-20 recall on literature-grounded gold sets. On synthetic benchmarks the co-expression fallback achieves ~100% when the correlation signal is clean.

---

## Known gotchas

- **Seed genes must be present in the similarity matrix.** For `CoexpressionBackend`, that means every seed gene must be a column of the expression matrix. If some are missing, the expander raises rather than silently down-weighting the seed.
- **`candidate_genes` must be a superset of `seed_genes`** — the expander filters seed genes out of the final candidate list automatically, but requires the seed to be present in the pool to compute similarity.
- **Co-expression similarity is symmetric, not causal.** If you need directional regulator → target claims, use the Borzoi backend or F11 (causal) for invariant prediction.

---

## See also

- [API reference: `genesets`](../api/genesets.md) — full class + method signatures
- [Notebook 25 — gene_set_expansion](../../examples/notebooks/25_gene_set_expansion.ipynb) — cell-by-cell walkthrough
- [F11 causal inference](causal.md) — for directional parent-of-target claims instead of similarity
