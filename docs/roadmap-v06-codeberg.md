# PSF v0.6 Public Roadmap (Codeberg)

**Date:** April 16, 2026
**Author:** Rohit Chauhan
**Status:** Proposed
**Repo:** https://codeberg.org/pathways/pathway-subtyping-framework
**Target release:** v0.6.0
**Version string (main):** `0.6.0`
**Concept DOI:** 10.5281/zenodo.18442426
**Base state:** post-cleanup `main` (commit `fcfda1a` equivalent, molecular QC layer implemented)

---

## Strategic Rationale

PSF v0.4.0 established pathway-level subtyping on bulk and single-cell transcriptomics. v0.5 added the molecular QC layer — twelve features that distinguish pathway *activation* from pathway *resolution* (cascade, temporal, tension, gates, drift, offtarget, heterogeneity, dosage, crosstalk, feedback, stress, atlas).

v0.6 answers three limitations users and reviewers have raised against v0.5:

1. **Claims without confidence.** MSV components return point estimates. Researchers making clinical or experimental decisions need calibrated uncertainty — not a single number, but a distribution.
2. **Cross-platform generalization is asserted, not measured.** PSF scores from 10x, Smart-seq2, bulk RNA-seq, and spatial platforms are treated as comparable without a principled harmonization layer.
3. **Gene-set–only views miss context.** Pathway scores collapse cell state into a per-pathway number. Cells with identical pathway scores can have different underlying programs; cells with different pathway scores can share program motifs. Contextual embeddings from genomic foundation models complement — but do not replace — interpretable pathway scores.

v0.6 adds a **rigor layer** (uncertainty, cross-platform harmonization, knowledge graph refresh) and a **foundation-model interface** (Geneformer, scGPT, Enformer/Borzoi, UCE, Nicheformer, AlphaMissense, Evo 2). Pathway scores stay the primary output; foundation-model embeddings are complementary, not a replacement.

---

## Design Principles

1. **Interpretability first.** Pathway scores remain the primary output. Foundation-model embeddings are an optional enrichment channel, not a substitute.
2. **Additive, not breaking.** v0.5 users should upgrade to v0.6 with zero code changes. All new capabilities are opt-in.
3. **No mandatory GPU.** Every v0.6 feature must have a CPU path, even if slow. Foundation-model inference can download precomputed embeddings for benchmark datasets.
4. **Reproducible by default.** Every new feature ships with a pinned model version, a frozen test dataset, and a deterministic seed.
5. **License-clean.** Every integrated model must be Apache-2.0, MIT, or BSD-licensed, or have a documented commercial-research waiver. No GPL, no research-only-no-redistribution.
6. **No silent model downloads.** Models are fetched only after explicit user consent via `psf.models.download("geneformer")`.

---

## Scope — 12 enhancements, 3 phases

| # | Enhancement | Module | Phase | Est. LOC |
|---|---|---|---|---|
| 1 | Uncertainty quantification | `uncertainty/` | 1 | ~900 |
| 2 | UCE cross-platform harmonization | `harmonize/uce.py` | 1 | ~700 |
| 3 | Knowledge graph refresh (OmniPath 2025, SIGNOR 3.0, Reactome 2026) | `kg/` | 1 | ~400 |
| 4 | AlphaMissense in cascade scoring | `qc/cascade.py` (extension) | 1 | ~250 |
| 5 | Geneformer in-silico perturbation API | `perturb/` | 2 | ~1,200 |
| 6 | scGPT embedding API | `embed/scgpt.py` | 2 | ~600 |
| 7 | Enformer/Borzoi gene-set augmentation | `genesets/regulatory.py` | 2 | ~500 |
| 8 | Nicheformer spatial upgrade | `spatial.py` (extension) | 3 | ~450 |
| 9 | Evo 2 off-target sequence scoring | `qc/offtarget.py` (extension) | 3 | ~350 |
| 10 | Multi-omics fusion (ATAC, proteomics) | `omics/` | 3 | ~1,400 |
| 11 | Causal inference layer | `causal/` | 3 | ~800 |
| 12 | Active learning framework | `active/` | 3 | ~600 |
| — | **Total** | — | — | **~8,150** |

Tests target ~30% of implementation LOC (~2,450 additional lines). Expected net addition: ~10,600 LOC across ~45 new files and extensions to ~15 existing modules.

---

## Phase 1 — Rigor Layer

Phase 1 is the gating phase. Without it, the foundation-model features in Phase 2 inherit poorly calibrated uncertainty and untested cross-platform assumptions.

### Feature 1: Uncertainty Quantification

**Module:** `pathway_subtyping.uncertainty`
**Priority:** HIGH — foundational for reviewer and clinician trust
**Depends on:** existing MSV modules, pathway scoring

#### Problem

Every MSV output today (`maturation_score`, `stress_score`, `drift_score`, etc.) is a single scalar. A researcher reading `stress_score = 0.72` cannot distinguish "this cell is confidently in a stress state" from "the model is uncertain and the true value could be anywhere from 0.4 to 0.9."

#### Design

Three complementary uncertainty channels:

```
uncertainty/
├── bootstrap.py          # Non-parametric bootstrap over cells for MSV scores
├── conformal.py          # Split-conformal prediction intervals for pathway scores
├── bayesian_pathway.py   # Bayesian pathway_gmm variant with posterior over component assignment
└── calibration.py        # Reliability diagrams, ECE, Brier score
```

- **Bootstrap:** resample cells with replacement, recompute MSV, report 2.5–97.5 percentile intervals
- **Conformal:** split-conformal wrapper for any pathway score function, returns distribution-free prediction intervals with user-specified coverage
- **Bayesian GMM:** replace point estimates in `pathway_gmm` with posterior samples over component weights and means
- **Calibration:** tools to assess whether reported uncertainty matches held-out performance

#### Deliverables

- [ ] `BootstrapMSV` class with parallelized resampling
- [ ] `ConformalPathwayPredictor` — wraps any scoring function
- [ ] `BayesianPathwayGMM` — drop-in replacement for `pathway_gmm` with `.sample()` method
- [ ] `CalibrationReport` with reliability diagrams, ECE, Brier score
- [ ] Notebook 21: "Uncertainty in PSF outputs"
- [ ] Documentation: `docs/guides/uncertainty.md`
- [ ] Tests covering coverage guarantees on synthetic ground truth

#### Acceptance Criteria

- Conformal intervals achieve nominal coverage (e.g., 90%) on held-out TCGA-COAD and autism benchmark datasets within ±2%
- Bootstrap intervals on MSV scores remain stable (< 5% width variance) across 3 independent resampling seeds
- Bayesian GMM reproduces point-estimate results within 1% when posterior is collapsed to mode
- No breaking change to existing scoring APIs

---

### Feature 2: UCE Cross-Platform Harmonization

**Module:** `pathway_subtyping.harmonize.uce`
**Priority:** HIGH — resolves a long-standing user request and a claim in the PSF manuscript
**Depends on:** UCE checkpoint (Rosen et al. 2024), scanpy, anndata

#### Problem

PSF currently assumes pathway scores from 10x Chromium, Smart-seq2, bulk RNA-seq, and spatial transcriptomics are comparable after standard preprocessing. Empirically this holds within a platform but the cross-platform rho claim in the manuscript rests on a small number of paired samples. Users combining public datasets across platforms have reported score drift.

#### Design

Universal Cell Embeddings (UCE) provides platform-invariant cell embeddings trained on 36M cells across 33 species. Use UCE embeddings to:

1. Produce a platform-invariant "reference embedding" for every cell
2. Fit a linear map from platform-specific pathway scores → reference-space pathway scores
3. Report a **harmonization confidence** that quantifies how much a cell's raw scores needed to shift to align with the reference

```
harmonize/
├── uce.py              # UCE embedding extraction wrapper
├── align.py            # Linear alignment from platform-specific to reference-space scores
├── confidence.py       # Harmonization confidence scoring
└── benchmark.py        # Cross-platform evaluation suite
```

#### Deliverables

- [ ] `UCEEmbedder` — extracts UCE embeddings with CPU fallback
- [ ] `CrossPlatformAligner` — fits platform-specific → reference alignment
- [ ] `HarmonizationReport` — per-cell confidence + per-platform score drift summary
- [ ] Cross-platform benchmark: 4 platforms × 2 biological systems (cortex + colon)
- [ ] Notebook 22: "Combining datasets across platforms"
- [ ] Documentation: `docs/guides/cross-platform.md`

#### Acceptance Criteria

- Harmonized pathway-level rho across 10x vs Smart-seq2 on matched cortex samples exceeds 0.75 (pre-harmonization baseline: 0.55–0.65 typical)
- Harmonization confidence correlates with known quality issues (e.g., low-read-depth cells flagged as low-confidence)
- No degradation of within-platform performance on TCGA-COAD benchmark

---

### Feature 3: Knowledge Graph Refresh

**Module:** `pathway_subtyping.kg`
**Priority:** MEDIUM — low effort, immediately benefits all KG-dependent features
**Depends on:** OmniPath 2025 release, SIGNOR 3.0, Reactome 2026 release

#### Problem

The current KG was assembled against OmniPath 2024, SIGNOR 2.4, and Reactome 2025. Three of the twelve molecular QC features (cascade, tension, crosstalk) are sensitive to KG coverage. Updated releases add ~15% more curated interactions and fix edge-direction errors in several pathways.

#### Design

- Update KG loader to target OmniPath 2025 (April 2026 release) as the default upstream
- Pin exact OmniPath/SIGNOR/Reactome snapshot hashes in `pyproject.toml` metadata
- Add KG diff tooling to report what changed vs. the v0.5 KG
- Regression test: cascade and tension scores on TCGA-COAD must change by < 5% or flag for review

#### Deliverables

- [ ] Updated KG loader with new source versions
- [ ] `kg.diff(v1, v2)` utility
- [ ] Regression test suite against TCGA-COAD and autism benchmarks
- [ ] Migration guide: `docs/migration/v05-to-v06-kg.md`

#### Acceptance Criteria

- All v0.5 KG-dependent tests pass with the v0.6 KG (either with identical results or with flagged changes documented)
- Pinned hashes make KG reproducible for at least 18 months

---

### Feature 4: AlphaMissense in Cascade Scoring

**Module:** `pathway_subtyping.qc.cascade` (extension)
**Priority:** MEDIUM — small extension, large interpretive payoff
**Depends on:** AlphaMissense precomputed missense scores (DeepMind, public)

#### Problem

Cascade scoring (F1) currently treats all genes in a pathway as either expressed or not. For variants carriers (e.g., autism families with rare missense variants), whether a missense variant is likely pathogenic should modulate how much weight that gene contributes to cascade completion.

#### Design

- Ingest the AlphaMissense precomputed score table (one score per protein-coding missense variant, public)
- Extend `CascadeAnalyzer` to accept optional `variant_table` input
- Down-weight a gene's cascade contribution when a carrier has a high-AM-score variant in that gene

#### Deliverables

- [ ] `AlphaMissenseScorer` loader and lookup
- [ ] `CascadeAnalyzer.score(..., variants=...)` extension
- [ ] Unit tests with synthetic variant tables
- [ ] Example notebook: autism trio with AlphaMissense-modulated cascade score

#### Acceptance Criteria

- Existing cascade tests unchanged when `variants=None`
- On SPARK/SFARI autism trios, AM-modulated cascade score differs from unmodulated in carriers of curated ClinVar pathogenic variants

---

## Phase 2 — Foundation-Model Interface

### Feature 5: Geneformer In-Silico Perturbation API

**Module:** `pathway_subtyping.perturb`
**Priority:** HIGH — single highest-impact new capability
**Depends on:** Geneformer checkpoint (Theodoris et al. 2023), Phase 1 uncertainty layer

#### Problem

Users can compute PSF scores on observed data but cannot ask "what would the MSV look like if I knocked out gene X?" This counterfactual question is central to hypothesis generation in cell line QC, disease modeling, and target validation.

#### Design

Geneformer's in-silico perturbation produces a perturbed cell embedding given a gene deletion or over-expression. We wrap this to return a perturbed MSV, not just a perturbed embedding.

```
perturb/
├── geneformer_wrapper.py    # Load Geneformer, run perturbation
├── msv_from_embedding.py    # Map perturbed embedding → MSV scores
├── screen.py                # Batch perturbation screens
└── report.py                # Perturbation impact reports
```

Pipeline:
1. User supplies observed expression matrix + gene of interest
2. Geneformer produces perturbed embedding for each cell
3. A light head (fit on Phase 1 reference data) maps embeddings → MSV scores
4. Output: perturbed MSV with uncertainty intervals (via Phase 1 conformal wrapper)

#### Deliverables

- [ ] `GeneformerPerturber` class
- [ ] `MSVFromEmbedding` head (fit during package build, cached)
- [ ] `PerturbationScreen` — runs a panel of gene KOs and ranks by MSV impact
- [ ] Notebook 23: "In-silico perturbation screens with PSF"
- [ ] Documentation: `docs/guides/perturbation.md`
- [ ] GPU-optional; CPU fallback with cached embeddings for benchmark datasets

#### Acceptance Criteria

- Perturbing known master regulators (e.g., MYC in cancer, MECP2 in neurons) produces directionally expected MSV shifts
- Perturbation screens complete on benchmark dataset (2,000 cells × 100 genes) in under 30 minutes on single A100, under 4 hours on CPU with cached embeddings
- Conformal intervals on perturbed MSV remain calibrated

---

### Feature 6: scGPT Embedding API

**Module:** `pathway_subtyping.embed.scgpt`
**Priority:** MEDIUM — foundation for future integrations
**Depends on:** scGPT checkpoint (Cui et al. 2024)

#### Problem

Some downstream analyses (rare cell detection, cell-state trajectory inference) benefit from learned embeddings rather than pathway scores. Users currently compute these with external packages and pass them back to PSF — lossy and error-prone.

#### Design

Thin, official scGPT wrapper with a stable interface. Minimal logic; the value is (a) stable API, (b) cached benchmark embeddings, (c) integration with Phase 1 uncertainty.

```
embed/
├── scgpt.py          # scGPT embedding wrapper
├── base.py           # Abstract Embedder interface (for future models)
└── cache.py          # Embedding cache (content-hashed)
```

#### Deliverables

- [ ] `scGPTEmbedder` with stable interface
- [ ] Content-hashed embedding cache (safe across package upgrades)
- [ ] Integration: pass scGPT embeddings into `harmonize/align.py` as an alternative to UCE
- [ ] Notebook 24: "When to use embeddings alongside pathway scores"

#### Acceptance Criteria

- Embedding output is byte-identical across reruns with frozen seed
- Cache invalidation correctly triggered by scGPT checkpoint changes

---

### Feature 7: Enformer/Borzoi Gene-Set Augmentation

**Module:** `pathway_subtyping.genesets.regulatory`
**Priority:** MEDIUM — strengthens gene-set curation for users building custom sets
**Depends on:** Borzoi checkpoint, human genome sequence

#### Problem

PSF ships with curated gene sets but users defining custom sets rely on keyword matching or existing databases. Borzoi (Linder et al. 2024) predicts expression from DNA sequence and can identify genes under shared regulatory control — a complementary signal to curated pathway membership.

#### Design

- Given a user-supplied seed gene set, Borzoi-predicted regulatory similarity identifies candidate additions
- Report candidates with Borzoi-predicted co-regulation score, flagged by whether they are already in curated databases
- Human-in-the-loop: this is a *suggestion* tool, not an auto-expansion tool

#### Deliverables

- [ ] `RegulatoryGeneSetExpander`
- [ ] Notebook 25: "Expanding a custom gene set with regulatory evidence"

#### Acceptance Criteria

- On held-out Reactome pathways, expander's top-20 suggestions contain at least 30% true Reactome members (recall proxy)

---

## Phase 3 — Extensions

### Feature 8: Nicheformer Spatial Upgrade

**Module:** `pathway_subtyping.spatial` (extension)
**Priority:** MEDIUM
**Depends on:** Nicheformer checkpoint (Theis lab 2024)

Extend the existing 933-LOC spatial module with Nicheformer embeddings for joint dissociated + spatial analysis. Key capability: users who have both 10x dissociated scRNA-seq and Visium spatial of the same tissue can now score pathways in a common embedding space.

**Deliverables:** `NicheformerEmbedder`, joint-space alignment API, notebook 26.

**Acceptance:** Dissociated-vs-spatial pathway score rho exceeds 0.7 on paired cortex reference.

---

### Feature 9: Evo 2 Off-Target Sequence Scoring

**Module:** `pathway_subtyping.qc.offtarget` (extension)
**Priority:** LOW — niche, but valuable to CRISPR-editing workflows
**Depends on:** Evo 2 checkpoint (Arc Institute 2025)

Augment F6 off-target scoring with Evo 2's long-context sequence predictions. Current F6 uses sequence similarity; Evo 2 adds functional prediction for ambiguous off-target sites.

**Deliverables:** `Evo2OffTargetScorer`, comparison against current F6, notebook 27.

**Acceptance:** On a held-out CRISPR off-target benchmark, Evo 2–augmented F6 improves AUROC by at least 0.03 over current F6.

---

### Feature 10: Multi-Omics Fusion (ATAC, Proteomics)

**Module:** `pathway_subtyping.omics`
**Priority:** MEDIUM — long-standing limitation
**Depends on:** nothing new; uses standard formats

PSF v0.5 is RNA-centric. Pathway activity inferred from transcript abundance misses post-transcriptional and post-translational regulation. Fuse chromatin accessibility (ATAC-seq) and proteomics where available.

```
omics/
├── atac.py         # ATAC-seq pathway scoring
├── proteomics.py   # Proteomics pathway scoring
├── fusion.py       # Weighted fusion of RNA + ATAC + proteomics
└── discordance.py  # Flag pathways where RNA and protein disagree
```

Discordance between RNA and protein-level pathway scores is informative signal, not noise — flag these pathways explicitly.

**Deliverables:** ATAC + proteomics scorers, fusion weights learned from paired RNA+ATAC+protein reference data, notebook 28.

**Acceptance:** Fused score on paired CITE-seq reference exceeds RNA-only score on downstream cell-type classification by at least 3% accuracy.

---

### Feature 11: Causal Inference Layer

**Module:** `pathway_subtyping.causal`
**Priority:** LOW — advanced, long-tail
**Depends on:** `dowhy` or custom invariant-prediction implementation

Cascade scoring (F1) is currently correlational: "upstream and downstream are co-activated." A causal statement — "upstream activation causes downstream activation" — requires intervention or invariant prediction across environments.

Implement invariant causal prediction for pathway-level relationships when users provide environment labels (e.g., cell line, batch, donor).

**Deliverables:** `InvariantPathwayPredictor`, notebook 29.

**Acceptance:** On simulated data with ground truth causal structure, identifies causal parents at recall ≥ 0.7.

---

### Feature 12: Active Learning Framework

**Module:** `pathway_subtyping.active`
**Priority:** LOW — research-facing, not core
**Depends on:** Phase 1 uncertainty

Given a pool of unlabeled samples and a fixed label budget, select samples whose labels would most reduce PSF subtype uncertainty. Valuable for cohort expansion decisions and targeted sample selection.

**Deliverables:** `ActiveSampleSelector` with uncertainty-based, diversity-based, and hybrid strategies; notebook 30.

**Acceptance:** On a held-out autism cohort simulation, active learning reaches 90% of full-cohort accuracy using 40% of labels.

---

## Dependency Graph

```
Phase 1 (Rigor):
    F1 Uncertainty ────────┐
    F2 UCE Harmonization ──┤
    F3 KG Refresh ─────────┼─── Phase 2 (Foundation Models):
    F4 AlphaMissense ──────┘        F5 Geneformer Perturbation (needs F1)
                                     F6 scGPT Embeddings (needs F2)
                                     F7 Enformer/Borzoi Gene Sets

                                 ───── Phase 3 (Extensions):
                                        F8 Nicheformer (needs F6)
                                        F9 Evo 2 Off-Target
                                        F10 Multi-Omics Fusion
                                        F11 Causal Layer (needs F1)
                                        F12 Active Learning (needs F1)
```

Phase 1 is hard-blocking for Phase 2. Phase 3 can proceed in parallel once Phase 1 is complete.

---

## Timeline

| Phase | Target window | Deliverables |
|---|---|---|
| Phase 1 (Rigor) | June 2026 | F1–F4, v0.6.0-alpha |
| Phase 2 (Foundation models) | July 2026 | F5–F7, v0.6.0-beta |
| Phase 3 (Extensions) | August 2026 | F8–F12, v0.6.0-rc1 |
| Final release | September 2026 | v0.6.0 stable on all 6 platforms |

Rationale: grant deadlines through end of May 2026 dominate availability. Development begins in earnest in June.

---

## Testing Strategy

| Test level | Coverage target | v0.6 additions |
|---|---|---|
| L1 unit | per-function | ~1,500 LOC new unit tests |
| L2 integration | per-module | ~600 LOC |
| L3 cross-module | MSV + uncertainty end-to-end | ~200 LOC |
| L4 reference benchmark | TCGA-COAD, autism, cortex | ~150 LOC new asserts |

Every foundation-model feature ships with **cached embeddings** for at least one public benchmark dataset, so CI can run without GPU. Model downloads are gated behind opt-in.

New reliability tests:
- Conformal coverage holds on TCGA-COAD (F1)
- Cross-platform rho improvement on paired reference (F2)
- KG-dependent score drift within tolerance (F3)
- Geneformer perturbation directionality on known regulators (F5)

---

## Backward Compatibility

- All v0.5 APIs unchanged. `CascadeAnalyzer.score(expr)` works identically; `CascadeAnalyzer.score(expr, variants=...)` is the new opt-in path.
- MSV output adds optional `.uncertainty` attribute; default access patterns unchanged.
- Default KG stays v0.5 KG for any user who pins; v0.6 KG is opt-in via `psf.kg.load(version="2026-q2")` or becomes default on next semver bump.
- `pathway_gmm` unchanged; `BayesianPathwayGMM` is a new class.

Any user upgrading with zero configuration changes gets identical v0.5 behavior.

---

## Release Channels

- **PyPI:** `pip install pathway-subtyping==0.6.0`
- **Codeberg:** tagged release + GitHub Releases–style notes
- **Zenodo:** new version DOI under concept 10.5281/zenodo.18442426
- **Docker Hub:** `pathways/psf:0.6.0`
- **bio.tools:** updated entry with new capabilities
- **RRID:** SCR_028051 updated descriptor

Release notes must include: migration guide, new benchmarks, model checkpoint provenance table, and an explicit "what has NOT changed" section to reassure existing users.

---

## Success Criteria

v0.6.0 is ready to release when:

- [ ] All Phase 1 acceptance criteria pass on three benchmark datasets
- [ ] At least three Phase 2 features pass acceptance
- [ ] Zero regressions on v0.5 test suite (1,363 tests)
- [ ] New test suite (~2,450 LOC) passes in under 20 minutes on CI
- [ ] Documentation updated: new guides for uncertainty, harmonization, perturbation, embeddings
- [ ] Six notebooks added (21–30 where applicable)
- [ ] Migration guide published
- [ ] Concept-DOI citation updated on README

---

## Out of Scope for v0.6

Explicit non-goals — defer to v0.7 or later:

- Real-time / streaming inference
- Online learning / continuous adaptation
- Federated training across institutions
- Protein structure modeling (ESM-3, AlphaFold 3 integration)
- Training our own foundation model
- Wet-lab integration APIs (LIMS, robotics)
- Multi-tenant deployment tooling

---

## Appendix: Foundation-Model Provenance Table

| Model | Version | License | Size | Paper | Checkpoint URL |
|---|---|---|---|---|---|
| UCE | 4-layer, 33M cells | MIT | 0.8B params | Rosen et al. 2024 | HuggingFace |
| Geneformer | V2 | Apache-2.0 | 10M params | Theodoris et al. 2023, Nature | HuggingFace |
| scGPT | whole-human | MIT | 51M params | Cui et al. 2024, Nat. Methods | GitHub |
| Nicheformer | public | Apache-2.0 | TBD | Schaar et al. 2024 | HuggingFace |
| AlphaMissense | precomputed scores | CC BY-NC-SA-4.0 (research) | — | Cheng et al. 2023, Science | EBI |
| Borzoi | public | Apache-2.0 | 300M params | Linder et al. 2024 | Google Storage |
| Evo 2 | 7B public weights | Apache-2.0 | 7B params | Arc Institute 2025 | HuggingFace |

Provenance is pinned in `pyproject.toml` metadata and verified at model download time via checksum.

---

**End of v0.6 Public Roadmap**
