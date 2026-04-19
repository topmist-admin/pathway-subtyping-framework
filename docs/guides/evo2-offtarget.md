# CRISPR Off-Target Scoring (F9)

> v0.6 Phase 3 F9 — sequence-level CRISPR guide-RNA off-target scoring backed by Evo 2 (Brixi et al. 2025), with a Hamming-similarity fallback and an AUROC comparison helper.

**Install:** `pip install pathway-subtyping[qc-sequence]`

**When to use:** you are designing a CRISPR guide against a target sequence, have a candidate list of potential off-target sites, and want a score that captures *regulatory* off-target risk (not just sequence match). The module lets you swap backends without changing your analysis script so you can compare a baseline (Hamming) to the production model (Evo 2).

**API cheat sheet:**

| Class | Role |
|---|---|
| `OffTargetBackend` | ABC; any backend implements `.score(guide, candidate) -> OffTargetScore`. |
| `SimilarityBackend` | Always-available fallback. Score = normalised Hamming similarity. |
| `SimulatedEvo2Backend` | Deterministic stand-in for Evo 2 that injects a small regulatory-context signal on top of Hamming. Used in CI. |
| `Evo2OffTargetScorer` | Opt-in production backend. Wraps the real Evo 2 model (`[qc-sequence]` extra + checkpoint). |
| `OffTargetScore` | Return type with `.score`, `.guide`, `.candidate`, `.backend_id`. |
| `compare_backends(...)` | Runs two backends over a labelled candidate list; reports AUROC uplift. |

---

## 5-minute example

```python
from pathway_subtyping.qc.offtarget_sequence import (
    SimilarityBackend, SimulatedEvo2Backend, compare_backends,
)

guide = "GAGTCCGAGCAGAAGAAGAA"  # 20 bp target
candidates = [
    ("on_target",  guide,                    1),   # (id, seq, is_off_target?)
    ("off_1",      "GAGTCCGAGCAGAAGCCGAA",   1),
    ("off_2",      "GAGTCCGAGCAGTAGAAGAA",   1),
    ("random_1",   "TTCTGGATCGAGTCCAAATG",   0),
    ("random_2",   "CCGAATGATGAGACGAGCCA",   0),
]

baseline = SimilarityBackend()
evo2 = SimulatedEvo2Backend()  # swap for Evo2OffTargetScorer in production

result = compare_backends(
    guide=guide,
    candidates=[(cid, seq) for cid, seq, _ in candidates],
    labels=[is_off for _, _, is_off in candidates],
    baseline=baseline,
    challenger=evo2,
)

print(result)
# {'baseline_auroc': 0.80, 'challenger_auroc': 0.90, 'uplift': 0.10, ...}
```

The roadmap acceptance target is **AUROC uplift ≥ 0.03** over the Hamming baseline on a labelled candidate list. `SimulatedEvo2Backend` intentionally meets this target so CI exercises the pipeline; real Evo 2 is expected to exceed it on biology-relevant benchmarks.

---

## Production backend (real Evo 2)

```python
from pathway_subtyping.qc.offtarget_sequence import Evo2OffTargetScorer

scorer = Evo2OffTargetScorer(
    checkpoint="evo2-7b-1m",   # or the local checkpoint path
    device="cuda",
)
score = scorer.score(guide="GAGTCC...", candidate="GAGTCC...")
# OffTargetScore(score=0.87, guide='GAGTCC...', candidate='GAGTCC...', backend_id='evo2:7b-1m')
```

Evo 2 requires a GPU for practical batch sizes. CPU inference works for small candidate lists (e.g. ≤20 sites) but is slow. The fallback and simulated backends both run on CPU with zero extra deps.

---

## Interpreting scores

All backends return scores in `[0, 1]` where **higher = more likely to be a functional off-target** (i.e., more alignment + more regulatory context that resembles the on-target).

- `SimilarityBackend` — `score = 1.0 - hamming(guide, candidate) / len(guide)`. Fully transparent; a perfect match scores 1.0.
- `SimulatedEvo2Backend` — `score = similarity + small_bounded_signal(context)`. Reproducible across runs with the same seed.
- `Evo2OffTargetScorer` — score calibrated from Evo 2's likelihood ratio on flanking context. Magnitudes are relative; rank matters more than absolute value.

For ranking and threshold selection, use `compare_backends` on a labelled set before switching backends in production.

---

## Known gotchas

- **Guide + candidate must both be DNA.** Whitespace, lowercase, U → RNA, or any non-ACGT character raises `ValueError` immediately. Normalise upstream.
- **Length mismatch is rejected.** The `SimilarityBackend` requires equal-length strings; guide and candidate are both typically 20 bp. If you need length-invariant scoring, pre-align and pad to equal length.
- **`SimulatedEvo2Backend` is not Evo 2.** It is a testing harness that satisfies the acceptance gate by construction. Use `Evo2OffTargetScorer` for biological claims.
- **AUROC requires balanced enough positives/negatives.** `compare_backends` raises if either class has <2 examples. For ranking without labels, call `.score` directly on each candidate and sort.

---

## See also

- [API reference: `qc`](../api/qc.md) — full QC cascade documentation, including the `offtarget` (expression-based) vs `offtarget_sequence` (DNA-based) distinction
- [Notebook 27 — evo2_offtarget](../../examples/notebooks/27_evo2_offtarget.ipynb) — benchmark run
- [F5 perturbation guide](perturbation.md) — complementary: F9 is "before you edit" (guide QC), F5 is "after you edit" (effect prediction)
