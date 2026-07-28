# `genesets` — Regulatory gene-set expansion (F7)

Suggest candidate additions to a curated seed gene set using either DNA-sequence-predicted regulatory similarity (Borzoi) or expression co-expression (fallback).

**Import path:** `pathway_subtyping.genesets`
**Source:** [`src/pathway_subtyping/genesets/`](../../src/pathway_subtyping/genesets/)
**Guide:** [genesets.md](../guides/genesets.md) · [Notebook 25](../../examples/notebooks/25_gene_set_expansion.ipynb)
**Install:** `pip install pathway-subtyping[genesets]` (optional — `CoexpressionBackend` works without it)

---

## Public surface

```python
from pathway_subtyping.genesets import (
    RegulatoryGeneSetExpander,
    # Backends
    RegulatoryBackend, BorzoiBackend, CoexpressionBackend,
    # Types
    ExpansionCandidate, ExpansionResult,
)
```

---

## `RegulatoryGeneSetExpander`

```python
RegulatoryGeneSetExpander(
    backend: RegulatoryBackend,
    source_k: int = 3,              # report top-k source seed genes per candidate
)
```

**Method:**

```python
result = expander.expand(
    seed_genes: Sequence[str],
    candidate_genes: Sequence[str],
    top_n: int = 20,
    known_members: Optional[Mapping[str, Iterable[str]]] = None,
    min_score: float = 0.0,
) -> ExpansionResult
```

- `seed_genes` must be a subset of `candidate_genes` (the expander removes them from the final candidate list automatically but requires them in the similarity matrix).
- `known_members`: dict `{db_name: iterable_of_genes}` — candidates in any listed DB are flagged `in_curated_db=True` but still ranked.
- `min_score`: candidates below this are excluded.

**Scoring:** mean similarity of each candidate to the seed set (higher = more similar).

---

## `ExpansionResult`

| Attribute / method | Returns | Description |
|---|---|---|
| `.seed_genes` | `List[str]` | Echo of seed |
| `.candidates` | `List[ExpansionCandidate]` | Ranked descending by score |
| `.backend` | `str` | Backend identifier |
| `.as_dataframe()` | `DataFrame` | Columns: `gene`, `score`, `in_curated_db`, `source_1..source_k` |
| `.recall_against(ground_truth)` | `float` | Top-N recall vs a held-out gold set |
| `.top(n)` | `List[ExpansionCandidate]` | Prefix of the ranking |
| `.to_dict()` | `dict` | JSON-serialisable |

---

## `ExpansionCandidate`

```python
@dataclass
class ExpansionCandidate:
    gene: str
    score: float
    in_curated_db: bool
    sources: List[str]       # top-k seed genes contributing to the score
```

The `sources` field is the "why was this suggested" explanation: the seed genes most similar to this candidate, ranked by contribution.

---

## `RegulatoryBackend` (ABC)

```python
class RegulatoryBackend:
    backend_id: str
    def similarity_matrix(self, genes: Sequence[str]) -> pd.DataFrame: ...
```

Must return a square `(n_genes, n_genes)` DataFrame indexed by gene symbol. Symmetric matrices are expected but not required.

---

## `CoexpressionBackend` (fallback)

```python
CoexpressionBackend(expression: pd.DataFrame)
```

- `expression`: samples × genes DataFrame (log-TPM or any normalised scale).
- Similarity is Spearman rho over samples. Always available; used by default in tests and when the `[genesets]` extra isn't installed.

---

## `BorzoiBackend` (production)

```python
BorzoiBackend(checkpoint: str = "borzoi-human-v1.0")
```

- Lazy-loads the Borzoi DNA-sequence regulatory model via the `[genesets]` extra.
- Raises `ImportError` with install instructions if the extra is missing; raises `FileNotFoundError` if the checkpoint isn't cached.
- Similarity derives from predicted regulatory activity profiles — captures signal that co-expression misses (e.g. regulator-target pairs that aren't always co-expressed across environments).

---

## Roadmap acceptance (F7)

- Top-20 recall ≥ 30% against a literature-grounded gold set.
- On synthetic benchmarks with a clean correlation signal, `CoexpressionBackend` hits **100%** top-20 recall (tests: [`tests/test_genesets_regulatory.py`](../../tests/test_genesets_regulatory.py)). Borzoi-backed gold-set acceptance is a post-release follow-up.

---

## See also

- [Guide: gene-set expansion](../guides/genesets.md) — beginner-facing walkthrough
- [`causal`](causal.md) — complementary tool when you need directional parent-of-target claims instead of symmetric similarity
- [`knowledge_graph`](../guides/cross-cohort-validation.md) — for an alternative KG-based expansion path (pinned OmniPath / SIGNOR / Reactome snapshots)
