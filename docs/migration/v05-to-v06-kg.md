# KG Migration Guide — v0.5 → v0.6

*PSF v0.6 Phase 1, Feature F3 (Knowledge Graph Refresh)*

v0.6 bumps three upstream data sources that feed the PSF knowledge
graph:

| Source     | v0.5 release | v0.6 release | Notes |
|------------|-------------|-------------|-------|
| OmniPath   | 2024        | 2025        | +~15 % curated interactions; edge-direction fixes in FGF/Notch/WNT subnetworks |
| SIGNOR     | 2.4         | 3.0         | Adds neuronal-signalling subnetwork (~2,000 new edges) |
| Reactome   | 2025        | 2026        | Refactored cell-cycle hierarchy; expect minor pathway-ID migrations |

Three of the twelve molecular-QC features are KG-sensitive:
`F1_cascade`, `F3_tension`, and `F9_crosstalk`. If you rely on any of
these, read this guide before upgrading — otherwise the upgrade is
transparent.

---

## What changes in your code

Nothing changes in the pathway-scoring API. What changes is the content
of `pathway_subtyping.knowledge_graph.KG_SOURCES` — the pinned manifest
that identifies which upstream releases PSF builds its KG from.

Inspect the v0.6 manifest:

```python
from pathway_subtyping.knowledge_graph import KG_SOURCES, manifest_digest

for src in KG_SOURCES.values():
    print(src.source_id, src.release, src.release_date)
print("manifest digest:", manifest_digest())
```

Two KGs are bitwise-reproducible iff their manifest digests match. Pin
`manifest_digest()` in any downstream artifact that needs KG provenance
(publication supplementary materials, Zenodo bundles, CLI output).

---

## What to do if you've built downstream artifacts on v0.5

### 1. Diff the two KG snapshots

```python
from pathway_subtyping.knowledge_graph import diff_kgs

diff = diff_kgs(kg_v05, kg_v06)
print(diff.summary)
# -> {'nodes_added_total': 142, 'nodes_removed_total': 7,
#     'edges_added_total': 2104, 'edges_removed_total': 63,
#     'edges_flipped_total': 18}

print(diff.edges_flipped[:5])
# -> [('FGF2', 'FGFR1', 'gene_regulates_gene'), ...]
```

`edges_flipped` is the most important field. Direction reversals in
signalling edges can flip the sign of cascade-score contributions —
they are the single biggest source of silent-but-meaningful score
changes in a KG upgrade.

### 2. Run the regression harness

If you have a benchmark cohort (pathway-score matrix, TCGA-COAD
expression matrix, autism validation cohort), pass it through
`run_kg_regression` with your score functions to quantify per-score
drift:

```python
from pathway_subtyping.knowledge_graph import run_kg_regression

def cascade_score(kg, input_df):
    # Your existing cascade scoring function
    return {"F1_cascade_score": compute_cascade(kg, input_df)}

report = run_kg_regression(
    v1=kg_v05,
    v2=kg_v06,
    score_fns=[cascade_score],
    benchmark_inputs=[tcga_coad_expression],
    tolerance=0.05,   # roadmap default: flag >5% relative change
)

print(report.summary())
print("flagged:", [d.score_name for d in report.flagged_scores])
```

`report.passed` is True iff no score moved by more than the tolerance.
The roadmap acceptance target is 5 % — any score that moves more than
that should be investigated before v0.6 lands in production.

### 3. If regression reports changes

Three useful follow-ups:

1. **Inspect the flipped edges within the pathway**: most regressions
   trace to a handful of edges whose direction reversed. Fetch them
   from `diff.edges_flipped` and check the subnetwork.
2. **Check the citation/notes for the release**: `KG_SOURCES[src].notes`
   contains a short human summary of what changed — e.g. "Refactored
   cell-cycle hierarchy" gives you a head-start on what to expect.
3. **Pin your artifact to v0.5**: if you need stability before you can
   absorb the v0.6 changes, swap `KG_SOURCES` for `KG_SOURCES_V05` at
   runtime. The registry is a public symbol for exactly this reason.

```python
from pathway_subtyping.knowledge_graph import KG_SOURCES_V05, KnowledgeGraphBuilder
# Build v0.5 KG on v0.6 PSF — stable during migration window.
kg = KnowledgeGraphBuilder().build_from_sources(KG_SOURCES_V05)
```

---

## Reproducibility contract

Pinned hashes (`KGSource.sha256`) make the v0.6 KG reproducible for at
least 18 months — the upstream sources guarantee that their archive
URLs remain stable for that window. PSF additionally mirrors every
pinned bundle into its own archive and the `KGSource.verify_archive()`
method validates the integrity of local bundles before use.

> The `sha256` values shipped in `KG_SOURCES_V06` at the v0.6 release
> are the authoritative hashes. If you see the placeholder hash (64
> zeros), you are running a development snapshot and must not rely on
> the KG for any reproducibility-critical output.

---

## Release engineering checklist (maintainers)

When bumping a source release:

1. Download the bundle locally.
2. Compute SHA-256: `shasum -a 256 bundle.tar.gz`.
3. Update the relevant entry in `src/pathway_subtyping/knowledge_graph/sources.py`.
4. Regenerate `manifest_digest()` and pin it in the release notes.
5. Run `tests/test_kg_refresh.py` to confirm the manifest is well-formed.
6. Build a v0.5 and v0.6 KG; run `diff_kgs()` and record the diff
   summary in the changelog.
7. Run `run_kg_regression()` against the TCGA-COAD and autism benchmark
   cohorts; flag any score delta > 5 % for manual review.

---

## See also

- [docs/roadmap-v06-codeberg.md](../roadmap-v06-codeberg.md) — v0.6
  roadmap, Feature F3.
- [src/pathway_subtyping/knowledge_graph/sources.py](../../src/pathway_subtyping/knowledge_graph/sources.py)
  — pinned manifest.
- [tests/test_kg_refresh.py](../../tests/test_kg_refresh.py) —
  acceptance test suite.
