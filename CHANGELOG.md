# Changelog

All notable changes to the Pathway Subtyping Framework will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

_Nothing yet._

---

## [0.9.0] - 2026-07-29

**Minor, not patch:** this release adds a new public module (`kg_sensitivity`,
Gate K) alongside the fixes, so it is a feature release under semver.

Artefacts (fill in at release time):
- **PyPI:** https://pypi.org/project/pathway-subtyping/0.9.0/ — `pip install pathway-subtyping==0.9.0`
- **Source:** tag `v0.9.0` on GitHub (`topmist-admin/pathway-subtyping-framework`) and Codeberg (`pathways/pathway-subtyping-framework`) · RRID:SCR_028051
- **Zenodo:** _versioned DOI pending_ — deposit under concept DOI `10.5281/zenodo.18638048`

> **No manuscript number changes in this release.** All four offline
> reproduction packages (`benchmark_audit`, `ablation_honest`,
> `flagship_donor_level`, `flagship_stability`) were re-run against the fixed
> code on 2026-07-29 and are **byte-identical** to the deposited artifacts. The
> corrected code paths are not reachable from the reproduction bundle: it calls
> `select_n_clusters` without `min_cluster_fraction`, its somatic tables are
> well-powered (n=430, χ²≈187, far from the sparse regime), and it does not
> import network propagation, the bootstrap CI, conformal prediction, threshold
> calibration, or Gate K at all.

### ⚠️ Breaking changes — do not take this release blindly

This is a **correctness** release. The headline is the six statistical fixes
below, not the new module. It is a MINOR bump rather than a patch because code
written against 0.8.0 can break, in three ways — one of which is silent.

**1. Calls that used to succeed now raise.**
- `statistical_rigor.sensitivity_analysis_weights(...)` raises `NotImplementedError`.
  It previously returned hard-coded zeros for every input, and
  `SensitivityResult.is_robust()` reported `True` on them — "your results are
  robust to weight scheme", confidently, for any input. Raising is the fix.
- `clustering.select_n_clusters(..., method=<unrecognised>)` raises `ValueError`
  instead of silently falling through to silhouette.

**2. A return type changed sign.**
- `statistical_rigor.bootstrap_effect_size_ci(...)` now returns a **signed**
  interval and can be negative (e.g. `(-3.308, -2.290)`). It previously returned
  `max(|d|)` bounds, always ≥ 0. Any caller that assumed non-negativity, or that
  compared the lower bound to 0 as a significance test, must be revisited — that
  comparison was meaningless before (the bound excluded 0 on 40/40 null datasets)
  and is meaningful now.

**3. Numbers change, silently — the one to actually worry about.**
- `network_propagation` random-walk-with-restart and pagerank return **different
  values**. The previous implementation did not conserve probability mass and
  suppressed hubs. Nothing errors; results simply differ, and the old ones were
  wrong. **Re-run any analysis that used network propagation.**

**Also:** `CrosstalkDetector`, `CrosstalkResult` and `InterferenceEdge` are no
longer in `pathway_subtyping.qc.__all__`, so `from pathway_subtyping.qc import *`
no longer provides them. Direct imports still work — see the F9 deprecation entry.

**Not affected:** pathway scoring (ssGSEA / GSVA / mean-Z), the discreteness gate,
clustering results under the default configuration, and every number in the
cautionary-manuscript reproduction bundle (re-verified byte-identical above).

### Fixed
- **Split-conformal under-covered when the calibration set was too small.** The
  valid quantile is the `ceil((n+1)(1-alpha))`-th order statistic of the
  calibration scores; when that rank exceeds `n` no such statistic exists and the
  correct quantile is **+∞** (the all-labels set). Clamping to the max observed
  score instead silently missed the coverage guarantee — measured **0.9035 at
  n=10** against a 0.95 target, 0.9167 at n=12, 0.9373 at n=15. Coverage is the
  one guarantee conformal prediction offers, so quietly missing it is worse than
  returning a wide set. Now returns the infinite interval with a warning.
  `ReframedMembershipGate` sets `n_cal = max(10, round(0.30n))` against a 0.95
  target, so cohorts of roughly n ≲ 63 were in the affected regime.
- **`calibrate_thresholds` produced impossible thresholds for non-default alpha.**
  The linear rescale was unclamped, so `alpha=0.5` gave a stability floor of
  **−0.5000** (a bar every partition clears) and a null ceiling of 0.58 (a bar
  almost nothing exceeds) — both silently disabling the gate they configure. Now
  clamped to [0, 1], and documented as a heuristic rescale rather than a
  recalibration: a threshold for a different alpha is a different percentile of
  the null, not a scalar multiple.

### Documented (behaviour unchanged, but the claim was wrong)
- **The calibration lookup tables do not reproduce from their stated recipe.** The
  comment claimed `generate_calibration_table(n_simulations=500, seed=42)`.
  Re-running that exact recipe gives systematically smaller 95th percentiles —
  (50,2) 0.0920→0.0544, (100,3) 0.0580→0.0257, (100,5) 0.0830→0.0244, (30,4)
  0.2150→0.0944 — so the embedded values are **1.7×–3.4× more permissive** than
  the simulation they were said to come from, in the anti-conservative direction
  for a gate ceiling. **Deliberately not regenerated:** doing so would tighten
  every threshold 2–3× and change pass/fail for existing users, which is a
  decision, not a bug fix. Scope: **not used anywhere in the manuscript
  reproduction bundle**, so no published result depends on it; it does feed
  `pipeline.py`'s null-control gate.
- **ssGSEA is sensitive to gene column order when values are tied**, and the
  reproduction pipeline creates ties: `log2(CPM + 1)` maps every zero count to
  exactly 0. Measured on a realistic bulk matrix: **29% exact-zero entries**, and
  permuting gene order shifted scores by up to 1.68 with sample-ordering Spearman
  **0.70**. On dense continuous data ssGSEA is exactly invariant (max |Δ| = 0);
  mean-Z is invariant either way. `consolidation-cautionary/RUNME.md` now states
  that gene column order must be preserved. **Not changed in code** — altering
  tie handling would move deposited numbers mid-submission.
- **`run_statistical_analysis(normalization=...)` is recorded, never applied.** It
  is copied to `result.normalization_method` and nothing else; no aggregation
  happens in that function and `aggregate_pathway_scores` has no production call
  site. A result reporting `normalization_method: size_normalized` is a label, not
  a fact — annotated at the parameter so it stops implying work that did not happen.
- **`_score_mechanistic_alignment` ignores its `pathway` argument.** It scores the
  mechanism string alone against a fixed unsourced lookup, so "alignment" between
  a mechanism and a pathway is never computed. Annotated as a mechanism-type prior.

### Fixed
- **Random-walk-with-restart applied a row-stochastic matrix without transposing.**
  `W = D⁻¹A` is row-stochastic, so the walk on a column probability vector must
  apply `Wᵀ`. Using `W·p` made each node divide its **incoming** mass by its own
  degree instead of the sender's: probability was not conserved (a 5-node test
  summed to **0.719** instead of 1.0) and **hubs were suppressed rather than
  boosted** — the opposite of what propagation is for. Node B (degree 3) received
  exactly `0.30288/3`. Now matches both the closed form `α·(I−(1−α)Wᵀ)⁻¹p₀` and
  `networkx.pagerank` exactly. Affects `random_walk` and `pagerank`; heat
  diffusion and insulated diffusion were already correct. An un-converged run now
  logs a warning instead of returning silently.
- **`bootstrap_effect_size_ci` was one-sided by construction.** It bootstrapped
  `max(|Cohen's d|)` over all cluster pairs — a maximum of absolute values,
  bounded below by zero and positively biased — so the percentile interval
  **excluded zero on 40 of 40 pure-null datasets**, e.g. `(0.039, 0.606)` where no
  effect exists. It could never signal "no effect". Now pre-specifies one cluster
  pair on the full data and bootstraps the **signed** effect for that fixed pair.
  Null false-positive rate **40/40 → 3/40** (near the nominal 5%), with power
  retained (20/20 large effects still detected). The residual pair-selection
  optimism is documented in the docstring.
- **`somatic_alignment` ran chi-square on sparse tables, where it is invalid and
  anticonservative.** With min expected count 3.0 it returned p=0.0079 against
  Fisher exact p=0.0202 — a 2.6× overstatement, on the *positive-evidence* gate
  for cancer, where a rare driver in a small cohort could fabricate an anchor. Now
  reports `min_expected_count` and `sparse_table` always, falls back to **Fisher's
  exact** for 2×2 (valid there), and refuses to anchor on an r×c table whose
  p came from an invalid chi-square.
- **`select_n_clusters` silently returned a k it had just rejected.** When every k
  failed `min_cluster_fraction`, it fell back to `k_range[0]` with only a log line,
  and no caller checked. The fallback is now loud and machine-detectable via
  `optimal_k in result.rejected_k`. Also: `np.bincount(labels)` could not see an
  **empty trailing cluster**, so a degenerate k=3 fit producing only labels {0,1}
  passed the guard — now `minlength=k`. An unrecognised `method` string used to
  fall through to silhouette silently and now raises.
- **Gate K: the `degree`/`module` null perturbed at 2× the requested size.** A
  Maslov–Sneppen double-edge swap removes **two** edges and adds **two**, but each
  swap was counted as consuming one removal and one addition. Measured on a
  200-node graph, requesting 5/5 changes actually changed 10/10 — exactly 2.00×
  (`uniform` was unaffected at 1.00×). This defeats the size-matching the module's
  entire argument rests on: an over-perturbed null pushes null agreement down,
  which raises `P(null ≤ observed)` and therefore **systematically suppresses the
  `kg-sensitive` verdict**. Now `min(n_remove, n_add) // 2` swaps, each consuming
  two units per side. Only an **even** balanced budget is fully swap-able; an odd
  remainder falls to the degree-weighted residual, which is documented.
- **Gate K: a null that could not perturb the graph returned a confident
  `kg-sensitive`.** `rewire_kg` skips work it cannot do — no schema-valid exemplar
  for an edge type, graph too dense — and logs at DEBUG, so a draw can come back
  identical to the baseline. When a release *introduces* an edge type the baseline
  lacks, every null agreement was 1.0, making the observed value look extreme:
  reproduced as `kg-sensitive`, `p=0.0196`, `testable=True` on a null that never
  perturbed anything. Gate K now measures achieved perturbation per draw and
  **abstains** when the median is zero, warns below half the requested change, and
  reports `null_perturbation_requested` / `null_perturbation_median` on the result
  and in `to_dict()`.
  - Both found by adversarial review of code written the same week; both pinned by
    regression tests (`test_degree_mode_respects_the_requested_budget`,
    `test_abstains_when_rewiring_cannot_perturb_the_graph`). No published result is
    affected — Gate K has never been run outside its own tests.
- **`statistical_rigor.sensitivity_analysis_weights` now raises `NotImplementedError`
  instead of returning fabricated results.** It ignored all four of its inputs
  (`gene_burdens`, `pathways`, `cluster_labels`, `seed`) and appended a literal `0.0`
  / `0` per weight scheme, so two entirely different cohorts returned byte-identical
  output. Because `SensitivityResult.is_robust()` takes the `mean_val == 0` branch
  and returns `np.std(values) < threshold`, that all-zero vector reported
  **`is_robust(...) → True`** — i.e. "your results are robust to weight scheme",
  confidently, for any input. The signature and docstring are retained because the
  intent is worth preserving; the docstring now records what a real implementation
  requires, including that the already-weighted `gene_burdens` frame is insufficient
  input (the raw variant table is needed, which is the likely reason it was stubbed).
  - **No published result affected:** the function had zero call sites and was not
    exported from the package root.
  - `SensitivityResult.is_robust()` now documents its degenerate-vector trap:
    `[0.85, 0.85, 0.85]` and `[0.0, 0.0, 0.0]` both return `True`, but only the
    first is good news. Same shape as the degenerate-ground-truth ARI artifact
    behind the 2026-07 correction. Behaviour unchanged (zero callers); documented
    so whoever implements the function does not walk into it.
  - Pinned by `TestSensitivityAnalysisWeightsNotImplemented`.

### Deprecated
- **F9 (molecular QC layer) `CrosstalkDetector` soft-deprecated — its `competition_score` does not measure
  what it claims.** `_compute_competition` subtracts the shared gene set from both
  pathways and returns `mean(A-exclusive) - mean(B-exclusive)`, so the shared genes
  never enter the arithmetic and serve only as a presence gate. Demonstrated: varying
  shared-node expression across four orders of magnitude (0.1 → 500) leaves the score
  bit-identical, while holding shared nodes perfectly balanced and varying only the
  exclusive genes swings it from **+8.0 to −8.0** and flips `dominant`.
  - `competition_score`, `dominant` and the PASS/FAIL summary are **uninterpretable**.
    `n_significant` thresholds `abs(score)` at 0.3, which on log2-scale expression
    would fail most real batches for an unrelated reason.
  - Constructing the detector now emits a **`FutureWarning`** (not
    `DeprecationWarning`, which Python hides by default outside `__main__`).
  - The three names are withheld from `pathway_subtyping.qc.__all__` but remain
    **importable**, so existing code does not break — only `import *` and generated
    API docs are affected.
  - **Not a regression:** the module shipped this way in the v0.5 QC layer. A
    `shared_mean` local was computed and never used from the first commit, and was
    removed as dead code in the 2026-07-28 lint cleanup — which is what surfaced it.
  - **No published result is affected.** `CrosstalkDetector` is invoked only in its
    own unit tests; no deposited artifact contains `competition_score`. The existing
    tests could not have caught this: all three use i.i.d. `normal(5.0, 0.5)`
    expression, under which the score is ~0 for any formula, and none asserts on the
    score at all.
  - Unrelated and unaffected: `KnowledgeGraph.get_pathway_crosstalk()`,
    `get_shared_genes()`, and the separately-numbered **v0.6 F9**
    (`qc.offtarget_sequence`, Evo 2 off-target scoring) — "F9" is overloaded.
  - Fixing requires deciding what the score should mean; candidate formulas are in
    the module docstring. `TestCrosstalkSoftDeprecation` pins the current defect, so
    a genuine fix will make it fail — that is the signal to re-export F9.

### Added
- **`pathway_subtyping.kg_sensitivity` — Gate K, time-sliced knowledge-graph
  sensitivity.** Tests whether a partition survives being recomputed against a
  different KG version. Every existing gate holds the knowledge graph fixed and
  varies the data; none checks whether a finding is an artifact of the database
  release it was derived under.
  - `kg_timeslice_sensitivity(v1, v2, partition_fn, cohort, ...)` returns one of
    **`robust`**, **`kg-sensitive`**, **`generically-fragile`**, or a
    `not-testable (...)` abstention, with a `testable` flag so abstentions can be
    excluded from failure-rate denominators.
  - `rewire_kg()` builds the reference: a **size-matched random perturbation**
    that changes the same number of edges of the same types as the observed diff,
    chosen at random. Without it the raw agreement number cannot be read — the
    test suite constructs two cases with an identical observed ARI of −0.026 that
    receive opposite verdicts, separated only by the null (median 1.000 vs −0.026).
  - **Three null models via `rewiring=`**, weakest to strongest:
    `uniform` (node set + edge-type counts only), **`degree` (default)** —
    Maslov–Sneppen double-edge swaps preserving every node's exact in/out-degree,
    and `module` — degree-preserving swaps additionally matched to the observed
    diff's within-/cross-module split, with modules derived from
    `GENE_IN_PATHWAY` membership by `module_map_from_pathways()` /
    `within_module_fractions()`.
  - **The null choice can change the verdict**, so it is not cosmetic and it is
    recorded on the result (`KGSensitivityResult.rewiring`). On identical data
    with an identical observed ARI, `uniform` returns `generically-fragile`
    (p=0.066) where `degree` returns `kg-sensitive` (p=0.016) — pinned by
    `test_null_choice_can_change_the_verdict`. `uniform` destroys the degree
    sequence and so misreads changes at both ends of it: hub edges look specially
    targeted when losing one is ordinary, and peripheral edges look ordinary when
    losing one is specific. It is retained as a comparison baseline only.
  - `module` **abstains** when the graph carries no `GENE_IN_PATHWAY` edges rather
    than silently falling back to a weaker null.
  - Reuses `diff_kgs()` to size the perturbation and optionally
    `run_kg_regression()` for a scalar view. **Neither decides the verdict**; a
    regression test asserts a flagged scalar score leaves the partition verdict
    unchanged.
  - Decision rule, in full: `robust` if `observed_ari >= ari_min`, else
    `kg-sensitive` if `empirical_p < alpha`, else `generically-fragile`. Both
    terms are live — stated explicitly because the v0.8.0 gate shipped documenting
    three criteria of which only one decided.
  - Docs: `docs/kg_sensitivity_gate.md`. Tests: `tests/test_kg_sensitivity.py`.

### Documentation
- **Nine previously undocumented modules now have API reference pages**, closing the
  gap where public top-level exports had no reference at all: `kg_sensitivity`,
  `discreteness`, `clustering_dl`, `network_propagation`, `data_quality`,
  `validation_datasets`, `genetics`, `utils`. All are linked from `docs/api/index.md`,
  which gains a **v0.8 Layers** section and a "Proposed, not implemented" section so a
  design spec can never be mistaken for a capability.
- **Accuracy fixes to existing API pages**, each verified against source:
  `validation.md` (`run_all()` documented 6 of 11 parameters; five gates —
  ancestry-independence, cross-modal, confound-association, genetic and somatic
  anchoring — plus `cramers_v` were entirely undocumented, and the sample output
  implied only three gates exist); `pipeline.md` (`PipelineConfig` documented 22 of 39
  fields, so every expression and multi-omic entry point was undiscoverable; `run()`
  described a VCF-only flow); `config.md` (`ancestry` and `multi_omic` validation
  sections were missing). Two dead cross-links repaired in `genesets.md` and `qc.md`.
- `docs/METHODS.md` carries a **coverage notice**: its body describes v0.2/v0.3 and
  omits the discreteness gate, Gate K, Gates 5/6/7, conformal prediction,
  harmonization, perturbation, embeddings and causal inference. The Validation
  Framework section is updated to enumerate the current battery and to flag that
  bootstrap stability was demoted. The remaining sections are unchanged and the
  document should not be cited as a complete v0.8.0 methodology.
- `docs/terminology.md`: gate battery updated from three gates to the current eight,
  plus new entries for **abstention (`not-testable`)**, **size-matched null** and
  **degree-preserving rewiring**.
- `README.md` documents Gate K in the gate battery and feature table.
- Cross-links added between the three gate docs (A, K, T) so the shared null-design
  argument is discoverable from any of them.
- Fixed three dead relative links (`zenodo_setup.md` ×2, `contributor-kit/06`) and a
  stale line citation in `discreteness_gate.md`. **All 124 tracked markdown files now
  have zero dead relative links.**
- `docs/roadmap-trajectory-gate.md` — proposed Gate T spec (trajectory validation).
  **Design only; no implementation exists.**

## [0.8.0] - 2026-07-25

Artefacts:
- **PyPI:** https://pypi.org/project/pathway-subtyping/0.8.0/ — `pip install pathway-subtyping==0.8.0` (published via GitHub Actions Trusted Publishing).
- **Source:** tag `v0.8.0` on GitHub (`topmist-admin/pathway-subtyping-framework`) and Codeberg (`pathways/pathway-subtyping-framework`) · RRID:SCR_028051.
- **Zenodo:** `10.5281/zenodo.21566406` (https://doi.org/10.5281/zenodo.21566406) — v0.8.0, published 2026-07-25 under concept DOI `10.5281/zenodo.18638048`.

Corrects the stability gate's null. Adversarial methods review established that the
bootstrap-stability null tested the wrong hypothesis at small n, and this release
fixes it. See `docs/discreteness_gate.md` for the full rationale.

### Corrected
- **The stability null tested independence, not discreteness.** The bootstrap-stability
  null permuted each pathway column independently, which destroys inter-pathway
  correlation while preserving marginals — a null of *"the pathways are mutually
  independent,"* **not** *"there are discrete clusters."* At small n these diverge:
  a single correlated Gaussian blob or a 1-D continuous gradient (tumor purity,
  immune infiltration, proliferation) has no discrete clusters, yet a mixture model
  reproducibly bisects it every bootstrap (high observed ARI) while the feature-permuted
  null collapses (low ARI) — so a **continuum is falsely certified as a "reproducible
  subtype."** The independence null is **retained but demoted** to a confound / marginal
  control; it no longer decides discreteness.

### Added
- **`pathway_subtyping.discreteness` subpackage** with a discreteness-aware Gate A
  (`DiscretenessGateA`) that keeps the identical observed statistic (mean bootstrap-ARI
  at the fixed selected k) and replaces the reference with:
  - a **single-Gaussian (SigClust) reference** — primary (Liu, Hayes, Nobel & Marron,
    JASA 2008, `doi:10.1198/016214508000000454`);
  - the **gap statistic** — complementary (Tibshirani, Walther & Hastie, JRSS B 2001,
    `doi:10.1111/1467-9868.00293`);
  - **Hartigan's dip test** — optional corroborating unimodality flag (Ann. Statist.
    1985, `doi:10.1214/aos/1176346577`), via the new `diptest` extra.

  **Only the single-Gaussian reference decides the verdict** (`passed = obs > sg_p95
  and sg_p < alpha`). Gap and dip are computed and reported but never enter the
  decision rule — thresholding the SigClust *p* alone reproduces the ablation's TPR
  and FPR exactly. The gate has three outcomes (certify / reject / `not-testable`);
  `not-testable` is an abstention and must not be counted as a rejection.
- **Small-n hardening** (identical for observed data and every reference draw): PCA to
  d = min(10, ⌊n/10⌋) ≪ n (a 50-dim full-covariance mixture is singular at n in the
  tens); k fixed across observed and reference; silhouette-based k-stability routing to
  `not-testable` when no reproducible modal k exists.
- **`ReframedMembershipGate`** — the conformal membership gate reframed as *assignment
  sharpness conditional on a Gate-A/B-certified partition* (not a standalone guarantee),
  at a single 0.90 operating point with a bootstrap coverage CI and a minimum-n rule.
- **`pathway_subtyping.clustering_dl`** — deep-learning clustering baselines, DEC
  (Xie, Girshick & Farhadi, ICML 2016) and VAE-GMM / VaDE (Jiang et al., IJCAI 2017),
  exposed as `run_dec` / `run_vae_gmm`. These let the validation gates be run over
  DL-produced partitions (the gate is clusterer-agnostic — it tests the data, not the
  algorithm) and provide the method-comparison baselines for the cancer worked example.
  Torch is an optional dependency; the baselines are skipped gracefully when absent.
- **Synthetic-control tests** (`tests/test_discreteness_gate.py`): the new gate must fail
  a single Gaussian and a 1-D continuum (which the old null wrongly passed) and pass two
  separated clusters — with ground truth known by construction.
- **`pathway_subtyping.genetics` subpackage + Gate 7 — Genetic Anchoring** (feature-level).
  `ValidationGates.genetic_anchoring_gate()` tests whether a subtype's *defining genes* are
  over-represented for disease genetic-risk genes (hypergeometric over-representation,
  BH-adjusted across subtypes) against a **background-matched** null. Unlike every other
  gate this is **positive evidence**: a germline-variant enrichment cannot be manufactured
  by any postmortem/technical confound (PMI, RIN, dissection, batch), so a pass is
  confound-immune evidence a subtype axis is genetically implicated — necessary, not
  sufficient, and a specific low-power test (a null is weak evidence of absence). The null
  is background-matched because brain-expressed genes are already enriched for brain-disease
  risk; a genome-wide reference is reported per subtype for contrast but does not decide the
  gate. Reusable core in `genetics/gwas_enrichment.py` (`feature_level_anchoring`,
  `hypergeometric_enrichment`, `EnrichmentResult`); wired into `run_all()` behind an optional
  `genetic_anchoring=` argument (backward compatible). Scope is feature-level only;
  subject-level anchoring (rare-variant burden / PRS on same-donor genotypes) is a later,
  data-use-agreement-gated addition.
- **Positive-control tests** (`tests/test_genetic_anchoring_gate.py`): reconstruct the
  Voineagu discrimination (a neuronal/synaptic set enriches for autism risk under a
  brain-expressed-matched null while a glial/immune set does not) and verify the gate passes
  on the anchored axis. The refactored code reproduces the v0.7.0 reproduction-bundle
  headline exactly on the real deposited inputs (neuronal fold 16.5×, p≈1.0e-23; glial 2.6×,
  n.s.).
- **Gate 7 somatic mode — `ValidationGates.somatic_anchoring_gate()`** (+
  `genetics/somatic_anchoring.py`: `somatic_alignment`, `StratumAlignment`). The cancer
  counterpart to feature-level anchoring: instead of testing whether a subtype's *genes* are
  enriched for germline risk, it tests whether the *tumors* in a subtype carry a **somatic
  stratum** (driver mutation like BRAF-V600E/KRAS/MSI, CNA, or mutational-signature class)
  more than the others — a subject-level chi-square + Bergsma Cramér's V association,
  BH-adjusted across strata. It is the confound gate's statistic with **inverted polarity**:
  a strong nuisance association fails Gate 6, a strong somatic-driver association passes
  Gate 7. High-power (somatic signal is strong per-tumor) and matched to oncology, where the
  psychiatry-appropriate germline mode is weak. Positive-evidence; wired into `run_all()`
  behind an optional `somatic_anchoring=` argument (backward compatible). Docstrings carry
  the confound caveat (a somatic driver can co-vary with tissue-of-origin / tumor purity —
  run Gate 6 first). Germline *subject-level* anchoring (rare-variant burden / PRS on
  same-donor genotypes) remains the later, data-use-agreement-gated addition.
- **Somatic-anchoring tests** (`tests/test_somatic_anchoring_gate.py`): a synthetic
  colorectal-style scenario — a subtype aligned with a BRAF-V600E stratum anchors (Cramér's
  V high, significant) while a random stratum does not; per-stratum missing-value handling;
  `run_all` integration.

### Dependencies
- Added `joblib` as an explicit core dependency (previously relied on scikit-learn's
  transitive install; the discreteness gate imports it directly).
- Added a `discreteness` optional extra pinning `diptest` (**GPLv2+** — intentionally NOT
  a core dependency so the default MIT install and Docker image stay MIT-clean; the dip
  flag is skipped gracefully when absent).

### Attribution
- The discreteness gate and reframed membership gate were developed by a Topmist
  computational analyst under work-for-hire; © Topmist LLC, MIT-licensed with the rest
  of the framework.

## [0.7.0] - 2026-07-09

Post-correction hardening (2026-07-09). Two source fixes that close the failure
modes behind the 2026-07 benchmark/model correction, plus a frozen, verified
dependency set.

Artefact DOIs:
- **Zenodo:** `10.5281/zenodo.21279842` (https://doi.org/10.5281/zenodo.21279842) — v0.7.0, published 2026-07-09 under concept DOI `10.5281/zenodo.18638048`.
- **PyPI:** https://pypi.org/project/pathway-subtyping/0.7.0/ — `pip install pathway-subtyping==0.7.0` (published via GitHub Actions Trusted Publishing).

### Corrected / retracted (carried since 0.6.3)
- The **adaptive bootstrap-threshold model** (`threshold_model_real47.json`,
  reported at R²=0.889) does **not reproduce** and is **retracted** (see
  `src/pathway_subtyping/RETRACTED_threshold_model_real47.md`). The pipeline
  never consumed it. The **47-dataset benchmark** is **corrected** (invalid rows
  flagged; see `CORRECTION_2026-07/` + `ERRATUM_2026-07-08.md`). Corrected data
  is on Zenodo v2.0 (`10.5281/zenodo.21262112`). Correction notices ship in
  `README.md` and `KNOWN-ISSUES.md`.

### Added
- **Confound Association Gate (validation Gate 6)** — `ValidationGates.confound_association_gate()`
  in `validation.py`, with a `cramers_v()` helper (Bergsma-corrected Cramér's V).
  Tests each partition against named confounds (e.g. brain region, batch) via
  chi-square + Cramér's V, BH-adjusted across confounds; a partition **fails**
  if any *nuisance* confound is both significant and non-trivial (V ≥ 0.30).
  Diagnosis (keys `diagnosis`/`dx`/`disease`/`condition`/`phenotype`) is treated
  as biology-of-interest and never fails the gate. Wired into `run_all()` behind
  a new optional `confounds=` argument (backward compatible; the gate only runs
  when confounds are supplied). This is the gate whose absence let a brain-region
  artifact pass the full battery on GSE80655 (bootstrap ARI ≈ 0.92, yet Cramér's
  V ≈ 0.67 vs region and independent of diagnosis, p ≈ 0.41).
- **Guarded ARI metrics** — new `pathway_subtyping.utils.metrics` module:
  `safe_adjusted_rand_score`, `ari_with_validity`, `ari_degenerate_reason`.
  These return `NaN` + a reason (not a misleading "perfect" 1.0) on degenerate
  ground truth. The guard keys on ground-truth **structure**
  (`n_true_clusters < 2`, too few samples, or empty), so it also catches the
  degenerate case that returns ARI = 0.0 — which a value-based (`== 1.0`) guard
  would miss.
- **Regression tests** — `tests/test_metrics_ari_guard.py` (incl. the GSE136196
  ARI = 0.0 fixture) and `tests/test_confound_gate.py` (19 tests total, all pass).

### Fixed
- **Empty-input ARI artifact fixed at the source.** The `adjusted_rand_score`
  degeneracy that produced the spurious "perfect" scores in the 47-dataset
  calibration benchmark is now guarded framework-wide and wired into
  `benchmark.py` (degenerate ground truth → `ari=None`, never a fake score) and
  `pipeline.py` (`_generate_reports` planted-subtype ARI).
- **`KNOWN-ISSUES.md`** — corrected the count to **14 degenerate rows
  (`n_true_clusters=1`), 13 of which carried a spurious ARI = 1.0** (the 14th,
  GSE136196, returned 0.0), and added the "fixed at source" note.
- **Public erratum + Zenodo v2.0 description (2026-07-09)** — propagated the same
  count correction to `CORRECTION_2026-07/ERRATUM_2026-07-08.md` and
  `CORRECTION_2026-07/ZENODO_v2_description.md`, which previously stated all 14
  degenerate rows carried ARI = 1.0. NOTE: the Zenodo deposit
  `10.5281/zenodo.21262112` must be re-uploaded with the corrected erratum to
  make this live on the DOI record (manual step).

### Changed
- **Dependencies (frozen).** `requirements.txt` now pins exact versions verified
  against the test suite on 2026-07-09 (299 passed / 1 skipped; 19/19 new tests).
  Added `requests` (a declared core dependency that was missing from this file).
- **`pandas` cap raised `<3.0.0` → `<4.0.0`** in `pyproject.toml` and
  `requirements.txt`: the suite passes on pandas 3.0.2, and the previous cap
  excluded the tested version.
- **Test/plotting dependencies confirmed installed** in the working environment:
  `pytest`, `scikit-learn`, `matplotlib`, `seaborn` (the last three are already
  declared core/dev deps in `pyproject.toml`; two figure-generating tests
  require `matplotlib`+`seaborn` and are otherwise skipped-red in a minimal env).

## [0.6.3] - 2026-04-18

Artefact DOIs:
- **PyPI:** https://pypi.org/project/pathway-subtyping/0.6.3/
- **Zenodo:** `10.5281/zenodo.19648024` (https://doi.org/10.5281/zenodo.19648024)
  — released under concept DOI `10.5281/zenodo.18442426` (retired — now resolves to nothing; current concept DOI `10.5281/zenodo.18638048`)
- **Docker Hub:** `rohitdataops/pathway-subtyping:0.6.3-runtime`,
  `:0.6.3-jupyter`, `:latest`

### Validation

Real-data acceptance runs added for F2 (harmonize), F5 (perturb), and F10
(multi-omics fusion), plus a production Geneformer-backed F5 run on a real
WT vs MECP2-KO cohort. All land as new skip-on-absent test suites mirroring
F1's blueprint; CI remains deterministic when cohort artefacts are not
present locally.

- **F2** — `scripts/validate_f2_real_data.py` +
  `tests/test_harmonize_real_data.py`. Cohorts: GSE28521 (Affymetrix U133,
  post-mortem frontal cortex, n=79) × GSE80655 (Illumina HiSeq 2000
  RNA-seq, DLPFC, n=281). Pathway-mean Spearman rho lifts from
  -0.024 (95% CI -0.30..+0.25) to +0.52 (95% CI +0.24..+0.73) after
  alignment — uplift +0.55, passing the +0.10 uplift gate. The stricter
  roadmap 0.75 post-rho target requires paired-cell data and is tracked
  as aspirational in the JSON artefact.
- **F5** — `scripts/validate_f5_real_data.py` +
  `tests/test_perturb_real_data.py`. Cohort: TCGA-COAD (n=57, log1p TPM).
  FallbackPerturber backend + MSVFromEmbedding head produce directional
  agreement 13/14 = **92.9%** across curated (gene, pathway) edges
  anchored in MYC / TP53 / E2F1 / CCNE1 / CDK1 literature — clears the
  70% gate. Perturbed MSV conformal oracle deviation -0.0012 at 90%
  target — preserves the F1 calibration guarantee through the
  perturbation wrapper.
- **F10** — `scripts/validate_f10_real_data.py` +
  `tests/test_omics_real_data.py` +
  `data/omics/cite_adt_to_pathway.yaml`. Cohort: 10x
  `pbmc_1k_protein_v3` CITE-seq (713 cells, 17-antibody panel; 630
  after ADT-gating into 5 PBMC types). 1-NN cell-type classification
  accuracy rises from 56.5% (RNA-only) to **79.5%** (fused) — uplift
  +23.0 pp (95% CI +18.1..+27.6 pp), passing both the 3% uplift and
  strictly-positive CI-lower-bound gates.

### Features

- **Real Geneformer-backed F5 perturbation.** ``OfficialBackend`` is
  now a working Geneformer V2 104M wrapper (CLS-token embeddings + rank-
  tokenization + direct knockout-by-zero-count). Requires the optional
  ``[perturb]`` extra plus a locally cloned ``Geneformer-V2-104M``
  checkpoint — configured via the ``GENEFORMER_MODEL_DIR`` env var or the
  ``--geneformer-model-dir`` CLI flag on ``scripts/validate_f5_real_data.py``.
- ``scripts/validate_f5_real_data.py --backend {fallback,geneformer}``
  selects the perturbation backend; the Geneformer path also runs a new
  **real WT vs MECP2-KO** comparison on GSE123753 (Rodrigues et al. 2020,
  isogenic iPSC-derived cortical neurons with MECP2 deletion). Predicted
  in-silico MECP2-KO ΔMSV is compared to observed ΔMSV from the real RTT
  vs WT cohort on 50 hallmark pathways.

### Validation (Geneformer-backed F5)

- **GSE123753 WT vs MECP2-KO** (neurons; 3 WT + 3 KO, MSV head fit on
  all 11 GSE123753 samples): **50/50 pathways directionally agree**
  between in-silico KO and observed KO (gate: ≥70%); Spearman
  predicted-vs-observed ΔMSV rho = **+0.85**. Gates enforced in new
  ``test_wt_vs_ko_*`` cases in ``tests/test_perturb_real_data.py``.

Public test count: 1,612 → **1,634** (+13 real-data tests, +4
Geneformer WT-vs-KO tests, +5 Geneformer cache tests). 3 skipped are
the wt_vs_ko subtests that only run against the Geneformer artefact.

### Features (continued)

- **Content-hashed Geneformer embedding cache.** ``OfficialBackend``
  now accepts ``cache_dir=`` and, when set, transparently caches CLS-
  token embeddings keyed by (checkpoint + emb_mode + max_input_len +
  expression bytes). Reruns on the same cohort return in sub-
  millisecond time instead of the ~40-minute CPU forward pass.
  ``scripts/validate_f5_real_data.py`` wires a ``--geneformer-cache-dir``
  flag with a default at ``~/.cache/pathway-subtyping/geneformer`` and
  a ``GENEFORMER_CACHE_DIR`` env override. Disable with an empty string.

---

## [0.6.2] - 2026-04-18

Artefact DOIs:
- **PyPI:** https://pypi.org/project/pathway-subtyping/0.6.2/
- **Zenodo:** `10.5281/zenodo.19646697` (https://doi.org/10.5281/zenodo.19646697)
  — released under concept DOI `10.5281/zenodo.18442426` (retired — now resolves to nothing; current concept DOI `10.5281/zenodo.18638048`)
- **Docker Hub:** `rohitdataops/pathway-subtyping:0.6.2-runtime`,
  `:0.6.2-jupyter`, `:latest`

### Fixed

- Sync `__version__` in `src/pathway_subtyping/__init__.py` and the `version`
  field in `CITATION.cff` with `pyproject.toml`. In v0.6.1 these were left at
  `0.5.0`, so `import pathway_subtyping; pathway_subtyping.__version__`
  reported the wrong string on the v0.6.1 wheel. Publication-facing
  metadata (CITATION.cff) also now reflects v0.6.2 and the 2026-04-18
  release date. Still a packaging-only patch — tested behaviour is
  identical to v0.6.0.

---

## [0.6.1] - 2026-04-18

### Fixed

- `pyproject.toml`: declare the five foundation-model extras that the v0.6.0
  README advertised but didn't actually ship — `[harmonize]`, `[perturb]`,
  `[embed]`, `[genesets]`, `[qc-sequence]`. Each extra installs the PyTorch
  substrate (`torch>=2.0.0`; `perturb` additionally pulls `transformers>=4.35`)
  so the corresponding `Official*Backend` can lazy-load; the upstream
  model-specific package (geneformer, scgpt, nicheformer, borzoi, evo2, uce)
  must be installed separately until those packages ship stable PyPI wheels.
  The deterministic fallback implementations work without any of these
  extras. `[all]` is updated to include the new extras.

Known issue — `__version__` in `__init__.py` and the `version` field in
`CITATION.cff` were not synced to 0.6.1. Fixed in v0.6.2.

No code changes — v0.6.0 tested behaviour is preserved. This is a
packaging-only patch.

---

## [0.6.0] - 2026-04-18

v0.6 adds a **Rigor layer** (uncertainty, cross-platform harmonization,
KG refresh, AlphaMissense cascade) and a **Foundation-Model Interface**
(Geneformer perturbation, scGPT embeddings, Borzoi gene-set expansion,
Nicheformer spatial join, Evo 2 off-target, multi-omics fusion, causal
inference, active learning). All twelve roadmap features ship with
test-asserted acceptance criteria; every foundation-model wrapper has
an opt-in production backend (gated on the relevant extra + checkpoint)
and a deterministic PCA-based fallback so CI runs without heavyweight
model downloads.

Public-edition test count: 1,362 → 1,612 (+250 tests, all synthetic or
real-data acceptance).

### Added

#### Phase 1 — Rigor Layer

- **F1 Uncertainty quantification (`pathway_subtyping.uncertainty`)** —
  `ConformalPathwayPredictor` (split-conformal prediction intervals),
  `BootstrapMSV` (non-parametric bootstrap with per-cell and aggregate
  modes), `BayesianPathwayGMM` (drop-in Bayesian replacement for the
  point-estimate GMM with posterior sampling), `CalibrationReport`
  (ECE + Brier + reliability diagrams). Real-data acceptance on
  TCGA-COAD (n=57) and GSE28521 autism cortex (n=79): oracle-adjusted
  conformal coverage within ±1% of target. See
  `docs/guides/uncertainty.md` and `examples/notebooks/21_uncertainty.ipynb`.
- **F2 Cross-platform harmonization (`pathway_subtyping.harmonize`)** —
  `UCEEmbedder` (opt-in, `[harmonize]` extra) / `FallbackEmbedder`,
  `CrossPlatformAligner` (per-platform biology-regression), and
  `HarmonizationReport`. Synthetic 4-platform harmonized rho > 0.75
  (baseline 0.3–0.5). See `docs/guides/cross-platform.md` and
  `examples/notebooks/22_cross_platform.ipynb`.
- **F3 Knowledge-graph refresh (`pathway_subtyping.knowledge_graph.{sources,diff,regression}`)** —
  pinned v0.5 + v0.6 source manifests (OmniPath 2025, SIGNOR 3.0,
  Reactome 2026) with SHA-256 verification, `diff_kgs` utility
  reporting node/edge/direction changes, `run_kg_regression` with
  configurable threshold flagging, `manifest_digest` for whole-
  manifest reproducibility hashing. Migration guide at
  `docs/migration/v05-to-v06-kg.md`.
- **F4 AlphaMissense-modulated cascade (`pathway_subtyping.qc.alphamissense`)** —
  `AlphaMissenseScorer` loads per-variant pathogenicity scores and
  produces per-cell per-gene weight matrices. `CascadeAnalyzer` gains
  an optional `gene_weights` parameter; `gene_weights=None` is
  bit-identical to the variant-naive v0.5 baseline.

#### Phase 2 — Foundation-Model Interface

- **F5 Geneformer in-silico perturbation (`pathway_subtyping.perturb`)** —
  `GeneformerPerturber` with pluggable backend (`OfficialBackend` gated
  on `[perturb]` extra, `FallbackPerturber` for tests), `MSVFromEmbedding`
  ridge-regression head, `PerturbationScreen` batch runner,
  `PerturbationReport` with directional-signature check. Synthetic
  master-regulator tests pass for three simulated cluster markers.
  See `docs/guides/perturbation.md` and `examples/notebooks/23_perturbation.ipynb`.
- **F6 scGPT embeddings + cache (`pathway_subtyping.embed`)** — abstract
  `Embedder` interface, `scGPTEmbedder` with `OfficialSCGPTBackend` /
  `FallbackSCGPTEmbedder`, content-hashed `EmbeddingCache` safe across
  package upgrades. See `docs/guides/embeddings.md` and
  `examples/notebooks/24_embeddings.ipynb`.
- **F7 Regulatory gene-set expansion (`pathway_subtyping.genesets.regulatory`)** —
  `RegulatoryGeneSetExpander` with `BorzoiBackend` (opt-in) and
  `CoexpressionBackend` (fallback). On the synthetic Reactome-like
  fixture, top-20 recovers 100% of held-out pathway members (roadmap
  ≥30%). See `examples/notebooks/25_gene_set_expansion.ipynb`.

#### Phase 3 — Extensions

- **F8 Nicheformer spatial-joint (`pathway_subtyping.embed.nicheformer`)** —
  `NicheformerEmbedder` with shared-basis `embed_joint(dissociated,
  spatial)`. Dissociated-vs-spatial pathway-score rho > 0.7 on paired
  cortex synthetic reference. See `examples/notebooks/26_spatial_joint.ipynb`.
- **F9 Evo 2 sequence-level off-target (`pathway_subtyping.qc.offtarget_sequence`)** —
  `Evo2OffTargetScorer` (opt-in, `[qc-sequence]` extra), `SimulatedEvo2Backend`
  for tests, `SimilarityBackend` baseline, and `compare_backends` /
  `auroc` primitives. On a seed-conservation CRISPR panel, the
  simulated contender beats the similarity baseline by AUROC ≥ 0.03.
  See `examples/notebooks/27_evo2_offtarget.ipynb`.
- **F10 Multi-omics fusion (`pathway_subtyping.omics`)** — `ATACScorer`,
  `ProteomicsScorer`, `MultiOmicsFusion` with `learn_weights` grid
  search, and `flag_discordant_pathways` for RNA/protein disagreement.
  CITE-seq-style 1-NN cell-type accuracy uplift ≥ 3% over RNA-only.
  See `examples/notebooks/28_multi_omics.ipynb`.
- **F11 Invariant causal prediction (`pathway_subtyping.causal`)** —
  `InvariantPathwayPredictor` with combined mean+variance invariance
  test. On a 2-environment synthetic cohort, identifiable causal
  parents recalled with 1.0 recall (roadmap ≥ 0.7). See
  `examples/notebooks/29_causal_inference.ipynb`.
- **F12 Active-learning selection (`pathway_subtyping.active`)** —
  `ActiveSampleSelector` with uncertainty / diversity / hybrid
  strategies. 40% label budget reaches ≥ 90% of full-cohort 1-NN
  accuracy on an autism-style synthetic cohort. See
  `examples/notebooks/30_active_learning.ipynb`.

### Added — Scripts, tests, docs

- `scripts/validate_f1_real_data.py` — reproducible conformal coverage
  validation on TCGA-COAD + GSE28521.
- `tests/test_uncertainty*.py`, `tests/test_harmonize.py`,
  `tests/test_kg_refresh.py`, `tests/test_alphamissense.py`,
  `tests/test_perturb.py`, `tests/test_embed.py`,
  `tests/test_genesets_regulatory.py`, `tests/test_nicheformer.py`,
  `tests/test_offtarget_sequence.py`, `tests/test_omics_fusion.py`,
  `tests/test_causal.py`, `tests/test_active.py` — one test module
  per v0.6 feature.
- New guides: `docs/guides/{uncertainty,cross-platform,embeddings,perturbation}.md`.
- Migration guide: `docs/migration/v05-to-v06-kg.md`.

### Deprecated

- None. v0.6 is additive. No public v0.5 API changes.

### Removed

- None.

### Fixed

- None (feature release, not a bugfix release).

---

## [0.5.0] - 2026-04-13

### Added

#### Molecular QC Layer (`[qc]` extra)
- **12-feature molecular QC** for manufactured and engineered cells (`pathway_subtyping/qc/`)
- **F1 CascadeAnalyzer**: Topology-aware incomplete cascade detection using KG directed edges
- **F2 TemporalTracker**: Trajectory classification (resolving/stalled/reversing/oscillating)
- **F3 TensionScorer**: Molecular tension from open signaling loops
- **F4 ResolutionGate**: Unified RELEASE/HOLD/REJECT integrating all QC signals
- **F5 DriftDetector**: Cumulative pathway drift from baseline across passages
- **F6 OffTargetDetector**: INTENDED/TOLERATED/OFF_TARGET/EXCLUDED_VIOLATION classification
- **F7 HeterogeneityProfiler**: Batch uniformity with DBSCAN subpopulation detection
- **F8 DosageAnalyzer**: UNDER/IN_RANGE/OVER with therapeutic windows and stoichiometry
- **F9 CrosstalkDetector**: Shared node competition between pathways via KG
- **F10 FeedbackMonitor**: Activator-inhibitor correlation (intact/decoupled/inverted)
- **F11 StressFingerprinter**: Matches pathway patterns to 6 stressor signatures with remediation
- **F12 AtlasComparator**: Distance from reference atlas with nearest-type mapping
- **ManufacturingSimulator**: 9 injectable defect types, 12x12 orthogonality matrix, severity titration
- **3 scenario tests**: CAR T manufacturing, neural organoid, 20-passage stability

#### GNN & Graph Embeddings (`[gnn]` extra) — Experimental
- **TransEModel**: Translational KG embeddings (pure numpy, no PyTorch required)
- **RotatEModel**: Relational rotation embeddings in complex space
- **OntologyAwareGNN**: Heterogeneous GNN with edge-type-aware message passing (requires PyTorch)
- **BiologicalAttention**: Multi-head attention with biological prior bias injection (pLI, expression, SFARI)
- **PathwayCoAttention**: Bi-directional gene-pathway cross-attention
- **EmbeddingFusion**: Multi-source fusion (concat, weighted sum, PCA)
- **GNNTrainer**: Training loop with AdamW optimizer and evaluation metrics

#### Autism Interpretation Layer (`[autism]` extra) — Autism-Only
- **BiologicalRules**: 8 curated rules (R1-R7 + R3b) with literature citations
- **AutismRuleEngine**: Priority-sorted rule evaluation with autism-only enforcement
- **ConditionEvaluator**: 18 predicate types (variant, gene, expression, pathway, drug)
- **ExplanationGenerator**: Human-readable reasoning chains with mandatory disclaimers
- **NeurosymbolicCombiner**: 4 combination methods (weighted_sum, max, product, rule_guided)
- **EvidenceScorer**: Multi-criteria scoring with 11 pediatric safety flags
- **DrugTargetDatabase**: In-memory drug-gene mapping with mechanism classification
- **HypothesisRanker**: Diversity-constrained ranking (requires_validation always True)

#### Knowledge Graph Enhancements
- **GENE_REGULATES** edge type for directed signaling/regulatory relationships
- **Topology-aware methods**: `partition_pathway_genes()`, `topological_sort_pathway()`, `find_cascade_paths()`
- **Centrality scoring**: `compute_centrality()` (degree, betweenness, closeness, PageRank)
- **Topology-weighted scoring**: `topology_weighted_pathway_score()` — hub genes contribute more
- **Hierarchical queries**: `get_pathway_hierarchy()`, `get_all_descendant_genes()`
- **Cross-omics resolution**: `resolve_entity_chain()`, `get_drug_targets_in_pathway()`
- **Crosstalk quantification**: `get_pathway_crosstalk()`, `get_shared_genes()`
- **Builder methods**: `add_signaling_edges()`, `add_signaling_edges_from_dict()`

#### Documentation
- **docs/how-it-works.md**: Plain-language conceptual guide (pathway scoring, validation gates, 5-layer architecture)
- **docs/api/qc.md**: Full API reference for 12 QC features
- **docs/api/gnn.md**: Full API reference for GNN embeddings, model, attention
- **docs/api/autism.md**: Full API reference for rules engine, therapeutic ranking

#### Infrastructure
- **`[qc]`, `[gnn]`, `[autism]` optional extras** in pyproject.toml
- **Codeberg issues #30-#41**: 12 QC features with labels (qc, manufacturing, safety, v0.5)
- **Warning filters**: Suppressed 92K+ benign sklearn/umap warnings in pytest config
- **Test count**: 1,054 → 1,363 tests (309 new)

### Fixed
- **20 broken GitHub links** replaced with Codeberg URLs (account suspended since Feb 22)
- **3 broken internal cross-references** (data/pathways/README.md path, one-pager notebook link)
- **Python version badge**: 3.9+ → 3.8+ (matches pyproject.toml)
- **.gitignore**: Added audit reports, task files, QC roadmap docs

## [0.4.0] - 2026-03-04

### Added

#### Notebooks
- **Notebook 13** (`13_geo_blood_ados.ipynb`): GSE111175 blood transcriptomics with ADOS clinical correlation
  - 141 blood leukocyte samples (28 ASD, 113 control), Illumina BeadChip
  - Pathway subtypes detectable in peripheral tissue
  - Synaptic transmission correlates with ADOS Social Affect (rho=-0.52, FDR p=0.032)
  - ADOS severity stratification (Section 10b): chi-square, Fisher exact, ARI, contingency heatmap
- **Notebook 14** (`14_geo_blood_large_cohort.ipynb`): GSE18123 large blood cohort validation
  - 285 blood samples (72 ASD, 213 control) across two Affymetrix platforms
  - Cross-cohort projection from GSE111175 achieves ARI=0.374 (exceeds 0.3 threshold)
  - Cross-tissue pathway correlation with brain subtypes (rho=0.371)
- **Notebook 15** (`15_geo_scz_replication.ipynb`): GSE53987 SCZ replication on Affymetrix microarray
  - 205 samples across 4 diagnoses (SCZ, BD, MDD, CTL), 3 brain regions
  - Cross-platform projection from GSE80655 RNA-seq: ARI=0.319 (PASS)
  - Cross-disease ARI=0.792 (SCZ/ASD pathway convergence)
  - Multi-diagnosis pooled clustering: k=5, silhouette=0.450
- **Notebook 16** (`16_knowledge_graph_analysis.ipynb`): Cross-disease knowledge graph & drug repurposing
  - 336 genes from 21 unique pathways (15 ASD + 14 SCZ, 8 shared)
  - STRING PPI network: 4,378 edges, 97.9% gene coverage
  - 20 hub genes (11 cross-disease bridges) ranked by betweenness centrality
  - DGIdb drug repurposing: 1,546 unique drugs targeting 44 hub/core genes
  - 6 Louvain communities (all cross-disease), convergence subnetwork (260 genes)
  - Manuscript-ready network figure, drug overlay, pathway crosstalk heatmap
- **Notebook 17** (`17_tcga_cancer_validation.ipynb`): TCGA-COAD colorectal cancer subtyping
  - 452 primary tumor samples from NCI GDC REST API (STAR Counts, TPM)
  - 50 MSigDB Hallmark cancer pathway gene sets scored via ssGSEA
  - k=3 subtypes: Stromal/EMT (19%), Immune-cold (29%), Proliferative/Metabolic (52%)
  - CMS external validation via pure Python NTP classifier (482 CMScaller marker genes)
    - Subtype 0 → CMS4 at 76% (Fisher OR=16.7, p=1.4e-25)
    - 436/452 classified (96.5%, FDR ≤ 0.05), global ARI=0.10 (k=3 vs k=4 mismatch)
  - Survival analysis: Kaplan-Meier curves, log-rank test, Cox PH regression
  - k=4 sensitivity analysis: global ARI worsened to 0.082, k=3 confirmed as primary
  - Benchmark: pathway_gmm wins, 2/3 validation gates pass
  - Demonstrates disease-agnostic framework applicability to cancer data
  - Standalone notebook (no prior notebooks required)
- **Notebook 18** (`18_geo_clinical_phenotype.ipynb`): GSE15402 clinical phenotype validation
  - 116 lymphoblastoid cell lines (87 ASD, 29 controls), TIGR 40K platform
  - ADI-R severity subgroup association: chi-square p=0.001
  - S5→L subgroup OR=13.2 (raw p=0.008)
- **Notebook 19** (`19_scz_blood_multi_cohort.ipynb`): SCZ blood multi-cohort Hertzberg replication
  - 407 total samples (177 SCZ) across 5 GEO datasets, 5 microarray platforms
  - Per-dataset subtyping + merged analysis: k=7, silhouette=0.088, 2/3 gates
  - Cross-cohort projection: mean ARI=0.205 (GSE38484→GSE18312 ARI=0.469 PASS)
  - Hertzberg concordance: 4 Immune-like (A) + 3 Neuro-like (B) subtypes mapped

#### Scripts
- **Cytoscape figure generator** (`scripts/generate_cytoscape_figures.py`): Publication-ready network figures via py4cytoscape
  - Hub gene subnetwork (20 hubs + top neighbors, community-colored)
  - Drug-hub target map (bipartite gene-drug network with interaction types)
  - Community overview (336 nodes with community-aware layout)
  - Requires `[graph]` extra + Cytoscape desktop app

#### Dependencies
- New `[graph]` optional extra: `networkx>=3.0`, `py4cytoscape>=1.0.0`
- Updated `requirements.txt` to match `pyproject.toml` (removed stale `pysam` from core)

#### Infrastructure
- **Docker multi-stage build** (`Dockerfile`): 4 stages (builder → runtime → development → jupyter)
  - Runtime image: minimal production CLI (374 MB compressed)
  - Jupyter image: notebook server with all analysis dependencies (617 MB compressed)
  - 1003 tests pass in container; NB00 quick demo executes successfully
- **Docker Hub:** Published as `rohitdataops/pathway-subtyping` (`:0.4.0-runtime`, `:0.4.0-jupyter`, `:latest`)
- **Docker Compose** (`docker-compose.yml`): pipeline, dev, test, and jupyter services
- **SciCrunch RRID:** `RRID:SCR_028051` — added to CITATION.cff, README.md, pyproject.toml, Dockerfile
- **bio.tools registration:** `biotools:pathway-subtyping` — added to pyproject.toml

#### Documentation
- Notebook execution guide with dependency diagram (`docs/notebook-guide.md`)
- Notebook execution registry (`research-results/NOTEBOOK-EXECUTION-REGISTRY.md`)
- Updated JOSS paper with all 6 validation datasets (1,075 samples)

### Fixed
- **Notebook 17 MSI ARI bug:** `adjusted_rand_score([], [])` returned 1.0 when MSI column existed but had 0 valid samples; added `valid_msi.sum() >= 2` guard before ARI computation and `ari_vs_msi is not None` guard on visualization panel
- **Notebook 17 Cox PH ZeroDivisionError:** `tumor_stage` entirely NaN in GDC TCGA-COAD; added `dropna(axis=1, how='all')` and sample count guard
- **Notebook 13 missing `import json`:** ADOS severity cell 10b.3 lacked json import for results_summary update

## [0.3.1] - 2026-02-23

### Changed
- Migrated all repository URLs from GitHub to Codeberg (GitHub account suspended Feb 22)
- Updated PyPI metadata, CITATION.cff, README, CONTRIBUTING, CHANGELOG, and all docs
- Replaced Colab notebook links with local notebook paths (Colab requires GitHub hosting)

## [0.3.0] - 2026-02-16

### Added

#### Signaling Pathway Databases (#32)
- **Signaling databases module** (`signaling_databases.py`): Load cell-cell signaling databases as pathway gene sets
  - `SignalingDatabase`: Enum for CellPhoneDB and CellChatDB
  - `SignalingInteraction`: Ligand-receptor interaction dataclass with gene resolution
  - `SignalingDatabaseResult`: Result with pathway gene sets, citations (`.to_dict()`, `.format_report()`, `.get_citations()`)
  - `load_cellphonedb()`: Auto-download from GitHub, resolve complexes to genes, group by signaling classification
  - `load_cellchatdb()`: Load user-exported CellChatDB CSV (R binary not Python-readable)
  - `convert_interactions_to_pathways()`: Group interactions into pathway gene sets
  - `merge_signaling_databases()`: Merge gene sets from multiple databases
- Signaling pathway gene sets directly compatible with `score_pathways_from_expression()`, `run_clustering()`, and validation gates

#### Bulk Deconvolution Integration (#31)
- **Deconvolution module** (`deconvolution.py`): Estimate cell-type proportions from bulk RNA-seq
  - `DeconvolutionMethod`: Enum for NNLS deconvolution
  - `DeconvolutionQualityReport`: Reference coverage, proportion validation, gene overlap
  - `DeconvolutionResult`: Cell-type proportion matrix with `.to_dict()`, `.format_report()`, `.get_citations()`
  - `build_reference_profile()`: Aggregate single-cell reference to cell-type mean expression profiles
  - `deconvolve_bulk()`: Main entry point — NNLS deconvolution with quality checks
  - `combine_features()`: Merge pathway scores + cell-type proportions into unified feature matrix
  - `generate_synthetic_bulk()`: Create synthetic bulk from known proportions for testing
- **Multi-omic integration**: Added `ModalityType.DECONVOLUTION` for deconvolution-derived features

#### Cross-Modal Validation Gate (#33)
- **Cross-modal validation module** (`cross_modal_validation.py`): Gate 5 — tests whether subtypes replicate across data modalities
  - `CrossModalPairResult`: Per-pair concordance metrics (ARI, NMI, bidirectional transfer ARI)
  - `CrossModalValidationResult`: Top-level result with gate pass/fail, null calibration, format_report/get_citations
  - `SingleCellCompositionResult`: ANOVA-based test for distinct cell-type compositions across subtypes
  - `cross_modal_concordance()`: Main entry point — independent clustering per modality, pairwise ARI/NMI, transfer validation, permutation null calibration
  - `single_cell_composition_test()`: One-way ANOVA per cell type with Bonferroni correction
  - `generate_synthetic_multimodal_data()`: Planted shared subtype structure for testing and calibration
- **Validation integration**: Gate 5 added to `ValidationGates.run_all()` when `per_modality_scores` is provided (>= 2 modalities)
- **Pipeline integration**: `run_validation_gates()` passes `per_modality_scores` from multi-omic fusion result

#### Multi-Omic Pipeline Integration (#35)
- **Multi-omic fusion module** (`multi_omic.py`): Fuse pathway scores from VCF, bulk RNA-seq, and scRNA-seq modalities
  - `ModalityType`: Enum for VCF, expression, and single-cell modalities
  - `FusionStrategy`: Enum for concatenate, weighted_average, intersection_only fusion strategies
  - `MissingStrategy`: Enum for handling samples missing from some modalities (impute_zero, impute_mean, drop)
  - `ModalityInput`: Dataclass wrapping a single modality's scored output with validation
  - `SampleOverlapStats`: Statistics about sample overlap across modalities
  - `MultiOmicFusionResult`: Result container with fused scores, per-modality scores, overlap stats (`.to_dict()`, `.format_report()`, `.get_citations()`)
  - `MultiOmicQualityReport`: Per-modality quality reports, sample/pathway overlap statistics, warnings
  - `prepare_modality()`: Validate and wrap modality output for fusion
  - `fuse_modalities()`: Main entry point — concatenation (default), weighted average, or intersection-only fusion
  - `compute_sample_overlap()`: Analyze sample overlap across modalities with pairwise statistics
  - `correlation_analysis()`: Cross-modality pathway correlation (Pearson/Spearman)
- **Pipeline integration**: `input_type: "multi_omic"` in YAML config with `multi_omic` section for modality list, fusion strategy, weights
- **Config validation**: New `multi_omic` section validation with modality-type, path, weight, and strategy checks

#### Single-Cell Pathway Scoring (#30)
- **Single-cell module** (`single_cell.py`): Per-cell and pseudobulk pathway scoring from scRNA-seq data
  - `SingleCellScoringMethod`: Enum for per-cell mean_z and pseudobulk (mean_z, ssGSEA, GSVA) methods
  - `SingleCellInputType`: Enum for raw counts, log-normalized, and h5ad input types
  - `SingleCellQualityReport`: QC report with sparsity, cell counts, gene detection, pathway coverage
  - `SingleCellScoringResult`: Result container with cell-type-level and optional per-cell scores (`.to_dict()`, `.format_report()`, `.get_citations()`)
  - `load_single_cell_data()`: Load and validate h5ad or CSV/TSV single-cell data with auto-normalization
  - `score_single_cell_pathways()`: Main scoring entry point; pseudobulk reuses expression.py internals
  - Memory-efficient chunked per-cell scoring for datasets up to 50K cells
  - Sparse matrix support (scipy CSR/CSC) throughout
- **Optional `[sc]` dependency group**: `anndata>=0.9.0`

#### Advanced Visualization (#9)
- **Visualization module** (`visualization.py`): Interactive and publication-quality visualizations
  - `DimReductionMethod`: Enum for PCA, t-SNE, UMAP dimensionality reduction
  - `FigureFormat`: Enum for PNG, SVG, PDF, HTML export formats
  - `ReportConfig`: Configuration dataclass for interactive report generation
  - `VisualizationResult`: Result dataclass with output paths and metadata (with `.to_dict()`)
  - `compute_dim_reduction()`: PCA, t-SNE, or UMAP embedding with metadata (explained variance, KL divergence, etc.)
  - `plot_interactive_scatter()`: Plotly scatter plot of samples in reduced space (PCA/t-SNE/UMAP)
  - `plot_interactive_heatmap()`: Plotly heatmap of mean pathway Z-scores per subtype
  - `plot_cluster_distribution()`: Plotly bar chart of cluster sizes
  - `plot_subtype_trajectories()`: Plotly radar chart of subtype pathway profiles
  - `plot_static_scatter()`: Matplotlib fallback scatter plot (works without Plotly)
  - `export_figure()`: Multi-format export (PNG, SVG, PDF, HTML) for both Plotly and matplotlib figures
  - `create_interactive_report()`: Self-contained HTML report with all charts, summary statistics, and styling
  - `generate_all_figures()`: Convenience function to generate all visualizations at once
- **Optional `[viz]` dependency group**: `plotly>=5.15.0`, `umap-learn>=0.5.0`, `kaleido>=0.2.0`
- **Pipeline integration**: `generate_interactive_report` and `interactive_dim_reduction` config options in `PipelineConfig`
  - YAML config: `output.generate_interactive_report: true` and `output.interactive_dim_reduction: "umap"`
  - Graceful fallback when Plotly not installed (warning + skip)
- 53 new tests covering all visualization functions, dimensionality reduction, export formats, edge cases
- All interactive features degrade gracefully to static matplotlib when `[viz]` extra not installed

#### Pipeline Input Validation
- **Phenotype file validation**: Checks for empty files, missing `sample_id` column, duplicate sample IDs
- **Minimum sample size guard**: Error if samples < 2× max clusters, warning if < 5× max clusters
- **Sample overlap reporting**: Logs overlap between data samples and phenotype samples
- 6 new tests for pipeline input validation

#### CI and Code Quality
- Fixed all black/isort/flake8 lint issues across 16 files
- Added `-m "not network"` to CI test commands to skip flaky network tests
- Removed unused imports in expression.py, threshold_calibration.py, validation_datasets.py
- Closed GitHub issues #7 (cross-cohort), #8 (performance), #26 (variant QC)
- Created `py.typed` PEP 561 marker file

#### CI and Testing
- Added `requires_plotly` and `requires_umap` skip markers in `test_visualization.py` so visualization tests skip gracefully when `[viz]` extras are not installed (fixes CI failures)
- 49 new single-cell tests, 768 total tests passing

### Changed
- Total test count: 1,003 (up from 971)

## [0.2.3] - 2026-02-14

### Added

#### Cross-Cohort Validation Enhancements (#7)
- **CohortResult.to_dict()**: Serialize cohort results to dictionary
- **CrossCohortResult methods**: Added `.to_dict()`, `.format_report()`, `.get_citations()` for publication-ready output
- **`generate_synthetic_cohort_pair()`**: Convenience function for creating matched synthetic cohort pairs for testing and demos
- **Cross-cohort example script** (`scripts/cross_cohort_example.py`): CLI example with argparse
- **Cross-cohort user guide** (`docs/guides/cross-cohort-validation.md`): Interpretation tables, real-world workflow, common pitfalls
- **Cross-cohort API reference** (`docs/api/cross_cohort.md`): Full API documentation
- 12 new tests for cross-cohort enhancements

#### Performance Optimization (#8)
- **tqdm progress bars**: Added to validation gates (label shuffle, random gene sets, bootstrap stability), simulation analysis (Type I error, power analysis, sample size analysis), expression scoring (ssGSEA, GSVA), and sensitivity analysis (feature LOO)
- **`show_progress` parameter**: Added to `ValidationGates`, `estimate_type_i_error()`, `run_power_analysis()`, `run_sample_size_analysis()`, `score_pathways_from_expression()`, `vary_feature_subset()`, `run_sensitivity_analysis()`
- **Chunked processing**: `PipelineConfig` now accepts `use_chunked_processing` and `chunk_size` options, delegating to existing `compute_gene_burdens_chunked()` for large VCF files
- **Benchmark script** (`scripts/benchmark_performance.py`): Measures wall-clock time and peak memory via `tracemalloc`, reports pass/fail against targets (30 min, 8 GB for 10K samples)
- **Hardware guide** (`docs/guides/performance-and-hardware.md`): Recommended hardware table, memory estimation, chunked processing guide, Colab constraints, performance tips
- 9 new tests for performance parameters

#### Variant Quality Control (#26)
- **Variant QC module** (`variant_qc.py`): Standard genetic variant quality control filters applied before burden computation
  - `VariantQCConfig`: Dataclass with min_qual, min_call_rate, hwe_p_threshold, max_maf, min_gq, min_dp (with `.to_dict()`)
  - `VariantQCResult`: Dataclass with removal counts, per-variant metrics, retention rate (with `.to_dict()`, `.format_report()`, `.get_citations()`)
  - `compute_call_rate()`: Per-variant genotype call rate
  - `compute_maf()`: Minor allele frequency computation for diploid samples
  - `check_hwe()`: Hardy-Weinberg equilibrium chi-squared test per variant
  - `apply_genotype_filters()`: Per-genotype GQ/DP masking
  - `filter_variants()`: Applies QUAL, call rate, HWE, MAF filters in sequence with per-filter removal tracking
- **Pipeline integration**: `variant_qc_enabled` config option runs QC between data loading and burden computation
  - `PipelineConfig`: Added `variant_qc_enabled`, `variant_qc_min_qual`, `variant_qc_min_call_rate`, `variant_qc_hwe_p_threshold`, `variant_qc_max_maf` fields
  - `from_yaml()`: Now parses `variant_qc:` section
  - QC results included in JSON and Markdown reports
- **Config validation** (`config.py`): `_validate_variant_qc_section()` validates all QC parameter ranges
- 40 new tests covering all QC functions, config validation, and package imports

#### Data-Driven Validation Threshold Calibration (#19)
- **Threshold calibration module** (`threshold_calibration.py`): Replaces hard-coded validation thresholds with data-driven values that adjust for sample size and number of clusters
  - `CalibratedThresholds`: Dataclass with null ARI threshold, stability threshold, calibration method, interpolation flag (with `.to_dict()`, `.format_report()`, `.get_citations()`)
  - `CalibrationSimulationResult`: Dataclass for simulation distributions
  - `calibrate_thresholds(n_samples, n_clusters, ...)`: Auto-calibrate thresholds via lookup table, interpolation, or simulation fallback
  - `get_default_thresholds()`: Returns legacy 0.15/0.8 values for backward compatibility
  - `generate_calibration_table()`: Regenerate lookup tables from simulations
  - Pre-computed lookup tables: 56-entry grid (8 sample sizes × 7 cluster counts) with empirically-derived 95th percentile null ARI and 5th percentile stability ARI
  - Bilinear interpolation for intermediate configurations
  - On-the-fly simulation fallback for out-of-range configurations
- **Pipeline integration**: Auto-calibration in `run_validation_gates()` when thresholds are `null` in config
  - `PipelineConfig`: Added `validation_calibrate`, `validation_stability_threshold`, `validation_null_ari_max`, `validation_alpha`, `validation_n_permutations`, `validation_n_bootstrap` fields
  - `from_yaml()`: Now parses `validation:` section (previously ignored)
  - Calibration info included in JSON and Markdown reports
- **Config validation** (`config.py`): `_validate_validation_section()` validates threshold ranges, alpha, iteration counts
- **Threshold calibration tests** (`tests/test_threshold_calibration.py`): 46 tests covering lookup tables, interpolation, simulation, calibration modes, reproducibility
- **Table generation script** (`scripts/generate_calibration_table.py`): CLI script to regenerate lookup tables

### Fixed

#### ClinVar and Reactome Parser Updates
- **ClinVar parser** (`validation_datasets.py`): Handle NCBI's updated `gene_specific_summary.txt` column format
  - New format uses `Alleles_reported_Pathogenic_Likely_pathogenic` (combined column) instead of separate `Number_Pathogenic`/`Number_Likely_Pathogenic` columns
  - Parser auto-detects format and handles both old and new column names
  - Handles `Number_uncertain` (new) vs `Number_Uncertain_Significance` (old) column naming
- **Reactome parser** (`validation_datasets.py`): Handle Reactome's updated GMT layout
  - New format: `Pathway Name\tR-HSA-ID\tGenes` (R-HSA-ID moved to column 1)
  - Old format: `R-HSA-ID\tHomo sapiens: Description\tGenes`
  - Parser now checks species prefix in description field alongside existing name checks

## [0.2.0] - 2026-02-09

### Added

#### Ancestry / Population Stratification Correction
- **Ancestry correction module** (`ancestry.py`): PCA-based population stratification detection and correction
  - `AncestryMethod`: Enum for correction approaches (REGRESS_OUT, COVARIATE_AWARE, STRATIFIED)
  - `compute_ancestry_pcs()`: PCA on genotype matrix with monomorphic variant handling
  - `adjust_pathway_scores()`: Regression residualization to remove ancestry-correlated variance
  - `check_ancestry_independence()`: Kruskal-Wallis test with Bonferroni correction for cluster-ancestry independence
  - `stratified_analysis()`: Per-ancestry-group clustering with cross-group concordance
  - `compute_ancestry_correlation()`: Pearson correlation matrix between pathways and ancestry PCs
  - Dataclasses: `AncestryPCs`, `AncestryAdjustmentResult`, `AncestryStratificationReport` (all with `.to_dict()`, `.format_report()`, `.get_citations()`)
- **Ancestry validation gate** (`validation.py`): 4th validation gate tests cluster-ancestry independence
- **Ancestry simulation** (`simulation.py`): Configurable ancestry confounding for synthetic data
  - `n_ancestry_groups`, `ancestry_effect_size`, `ancestry_confounding` parameters in `SimulationConfig`
  - Simulated ancestry PCs and group labels in `SimulatedData`
- **Pipeline integration** (`pipeline.py`): Optional ancestry correction between pathway scoring and clustering
  - `ancestry_pcs_path`, `ancestry_correction`, `ancestry_n_pcs` in `PipelineConfig`
  - Ancestry section in JSON and Markdown reports
- **Config validation** (`config.py`): Ancestry section validation (method, PCs file, n_pcs)
- **Ancestry test suite** (`tests/test_ancestry.py`): 44 tests covering all functions, edge cases, and end-to-end validation

#### Batch Correction & Sensitivity Analysis
- **Batch correction module** (`batch_correction.py`): ComBat-style batch effect detection and correction
  - `BatchCorrectionMethod`: Enum for correction approaches (COMBAT, MEAN_CENTER, STANDARDIZE)
  - `detect_batch_effects()`: ANOVA-based batch effect detection with eta-squared variance explained
  - `correct_batch_effects()`: ComBat empirical Bayes correction, mean centering, or standardization
  - `validate_batch_correction()`: Post-correction validation of variance reduction and signal preservation
  - Dataclasses: `BatchEffectReport`, `BatchCorrectionResult` (all with `.to_dict()`, `.format_report()`, `.get_citations()`)
- **Sensitivity analysis module** (`sensitivity.py`): Systematic parameter robustness testing
  - `SensitivityParameter`: Enum for parameter axes (CLUSTERING_ALGORITHM, N_CLUSTERS, NORMALIZATION, FEATURE_SUBSET)
  - `vary_clustering_algorithm()`: Compare GMM, K-means, Hierarchical across algorithms
  - `vary_n_clusters()`: Sweep cluster count range with pairwise ARI
  - `vary_feature_subset()`: Leave-one-out pathway sensitivity
  - `vary_normalization()`: Compare z-score, min-max, robust, rank normalization
  - `run_sensitivity_analysis()`: Full sensitivity analysis with robustness scoring
  - Dataclasses: `ParameterVariationResult`, `SensitivityAnalysisResult` (all with `.to_dict()`, `.format_report()`, `.get_citations()`)
- **Batch correction tests** (`tests/test_batch_correction.py`): 34 tests covering detection, correction, validation, edge cases
- **Sensitivity analysis tests** (`tests/test_sensitivity.py`): 27 tests covering all parameter axes, reproducibility, dataclasses

#### Scientific Rigor Modules (Publication Readiness)
- **Statistical rigor module** (`statistical_rigor.py`): Publication-quality statistics
  - `benjamini_hochberg()`: FDR correction for multiple testing
  - `BurdenWeightScheme`: Literature-based variant weighting (DEFAULT, GNOMAD_CONSTRAINT, ACMG_INSPIRED, UNIFORM)
  - `PathwayNormalization`: Multiple aggregation methods (MEAN, MEDIAN, SIZE_NORMALIZED, PCA)
  - `compute_pathway_effect_sizes()`: Cohen's d with bootstrap confidence intervals
  - `compute_pathway_pvalues()`: Permutation-based p-value computation
  - `run_statistical_analysis()`: Comprehensive statistical analysis pipeline
- **Multiple clustering algorithms** (`clustering.py`): Algorithm comparison framework
  - `ClusteringAlgorithm`: GMM, K-means, Hierarchical, Spectral clustering
  - `run_clustering()`: Unified interface for all algorithms
  - `select_n_clusters()`: BIC or silhouette-based model selection
  - `cross_validate_clustering()`: K-fold cross-validation for stability
  - `compare_algorithms()`: Pairwise ARI comparison, consensus labels
- **Simulation framework** (`simulation.py`): Ground truth validation
  - `SimulationConfig`: Configurable synthetic data generation
  - `generate_synthetic_data()`: Planted subtype structure with effect size control
  - `estimate_type_i_error()`: False positive rate estimation under null
  - `run_power_analysis()`: Power curves across effect sizes
  - `run_sample_size_analysis()`: Sample size recommendations for target power
  - `validate_framework()`: Comprehensive framework validation
- **Formal methods documentation** (`docs/METHODS.md`): Statistical methodology for publications

#### Additional Disease Pathways (Week 5)
- **Parkinson's Disease pathways** (`parkinsons_pathways.gmt`): 14 pathways, ~280 genes
  - Alpha-synuclein aggregation, mitochondrial function, autophagy-lysosomal pathway
  - Dopamine metabolism, endolysosomal trafficking, immune/inflammation
  - Sources: Nalls et al. 2019 (Lancet Neurol), Blauwendraat et al. 2020, IPDGC
- **Bipolar Disorder pathways** (`bipolar_pathways.gmt`): 14 pathways, ~290 genes
  - Calcium signaling, circadian rhythm, WNT/GSK3 signaling
  - Glutamate/GABA signaling, HPA stress response, neuroplasticity
  - Sources: Mullins et al. 2021 (Nat Genet), Stahl et al. 2019, BDgene
- **Literature citations** added to autism pathway file header
- Updated pathway documentation with new disease recommendations

#### Real-World Data Support (Week 4)
- **Multi-allelic variant support**: Automatically expands multi-allelic variants (e.g., A→G,T) into separate bi-allelic records
- **Data quality module** (`data_quality.py`): Comprehensive VCF parsing with quality checks
  - `DataQualityReport`: Reports annotation coverage, multi-allelic handling, and data usability
  - `VCFDataQualityError`: User-friendly exceptions with fix suggestions
  - `load_vcf_with_quality_check()`: Robust VCF loading with quality validation
  - `validate_vcf_for_pipeline()`: Pre-flight validation function
- **Graceful handling of missing annotations**: Pipeline continues with warnings instead of failing
- **Enhanced annotation helper** (`scripts/annotate_vcf.py`):
  - Verbose mode with detailed statistics
  - Validation-only mode to check existing VCF
  - Better VEP/ANNOVAR format detection
  - Comprehensive error messages with fix suggestions

#### Performance Module Enhancements (`utils/performance.py`)
- **Gzip support**: `chunked_vcf_reader()` now handles `.vcf.gz` files automatically
- **Multi-allelic expansion**: Chunked reader expands multi-allelic variants with allele-specific genotype counting
- **Consistent genotype parsing**: Uses `parse_genotype()` with `target_allele` parameter for accurate multi-allelic handling
- **CADD defaults in chunked processing**: `compute_gene_burdens_chunked()` now applies consequence-based CADD defaults (35/20/10)
- **Zero-variance pathway filtering**: `parallel_pathway_scores()` removes zero-variance pathways before Z-score normalization with clear error messages

#### Configuration Validation (`config.py`)
- **ConfigValidationError class**: Custom exception with field tracking and actionable fix suggestions
- **Enhanced `validate_config()`**: Now accepts `check_files` parameter to skip file existence checks during testing
- **`validate_gmt_file()` function**: Validates GMT pathway files with detailed error reporting
  - Checks for minimum 3 tab-separated fields per line
  - Validates minimum 2 genes per pathway
  - Reports duplicate pathway names and parsing errors

#### Analytical Reliability Improvements
- **GMM convergence checking**: All GMM fits now verify convergence and log warnings if not converged
- **GMM covariance regularization**: Added `reg_covar=1e-6` to all GMM calls for numerical stability
- **Zero-variance pathway handling**: Automatically detects and removes pathways with zero variance before normalization
- **CADD missing value handling**: Uses consequence-based defaults (35/20/10) instead of silent zeros
- **Consistent genotype parsing**: Unified allele-specific counting between bi-allelic and multi-allelic variants
- **Empty ARI array handling**: Validation gates now handle edge cases where no GMM fits converge

#### Test Suite Expansion
- **CLI test suite** (`tests/test_cli.py`): 20+ tests for command-line interface
  - Version and help display tests
  - Config loading and validation tests
  - Command-line override tests (--output, --seed, --quiet)
  - Error handling and exit code tests
- **Performance module tests** (`tests/test_performance.py`): 25+ tests for performance utilities
  - Chunked VCF reader tests (plain and gzipped)
  - Multi-allelic expansion and genotype parsing tests
  - Gene burden computation tests
  - Parallel pathway scoring tests with zero-variance handling
  - Memory estimation and downsampling tests
  - Progress tracking tests
- **Updated config tests** (`tests/test_config.py`): Tests for new validation functions
- **Updated data quality tests**: Tests for multi-allelic `parse_genotype()` with `target_allele`
- **Statistical rigor tests** (`tests/test_statistical_rigor.py`): 32 tests for FDR, burden weights, effect sizes
- **Clustering tests** (`tests/test_clustering.py`): 26 tests for algorithms, CV, comparison
- **Simulation tests** (`tests/test_simulation.py`): 24 tests for synthetic data, power analysis
- **Ancestry tests** (`tests/test_ancestry.py`): 44 tests for PCA, correction, independence, stratified analysis
- **Batch correction tests** (`tests/test_batch_correction.py`): 34 tests for detection, correction, validation
- **Sensitivity analysis tests** (`tests/test_sensitivity.py`): 27 tests for parameter variation, robustness
- **Total test count**: 347 tests (up from 242)

### Changed
- Pipeline now uses `data_quality` module for VCF loading
- Pipeline reports include data quality section
- Version bumped to 0.2.0-dev
- `parse_genotype()` now takes `target_allele` parameter for consistent multi-allelic handling
- `validate_config()` raises `ConfigValidationError` instead of generic `ValueError`

### Documentation
- Updated troubleshooting guide with comprehensive real-world data section
- Added VCF validation instructions
- Added multi-allelic variant handling explanation
- Added CADD score coverage guidance
- Updated API documentation for config and validation modules

### Other
- PyPI package publishing preparation

## [0.1.0] - 2026-01-29

### Added

#### Core Pipeline
- Complete pathway subtyping pipeline with VCF → clustering → report workflow
- Gene burden computation with LoF/missense weighting and CADD score normalization
- Pathway score aggregation using GMT file definitions
- GMM clustering with automatic cluster selection via BIC
- Cluster labeling based on top contributing pathways

#### Validation Gates
- Negative Control 1: Label shuffle test (ARI should be < 0.15)
- Negative Control 2: Random gene sets test (ARI should be < 0.15)
- Stability Test: Bootstrap resampling (ARI should be >= 0.8)
- Comprehensive validation reporting in JSON and Markdown

#### Disease Support
- Autism pathway definitions (15 pathways, ~200 genes) - Validated
- Schizophrenia pathway template
- Epilepsy pathway template
- Intellectual disability pathway template
- Guide for adapting to new diseases

#### Infrastructure
- `pyproject.toml` for modern Python packaging
- CLI entry points (`psf`, `pathway-subtyping`)
- YAML-based configuration system
- Reproducibility features (seed control, metadata logging)

#### Testing
- 64 unit and integration tests
- Test fixtures for synthetic data
- CI/CD with GitHub Actions
- Multi-OS (Ubuntu, macOS) and multi-Python (3.9-3.12) testing

#### Documentation
- Comprehensive README with quick start guide
- Disease adaptation guide
- Pathway curation guide
- Validation gates explanation
- Getting-started Jupyter notebook tutorial

#### Containerization
- Multi-stage Dockerfile (runtime, dev, jupyter targets)
- docker-compose.yml for easy orchestration
- Pre-configured development environment

#### Sample Data
- Synthetic VCF with 60 samples and 30 variants
- Synthetic phenotypes with 4 planted subtypes
- Test configuration ready to run out of the box

### Technical Details
- Python 3.8+ support
- Dependencies: numpy, pandas, scikit-learn, scipy, pysam, matplotlib, seaborn
- MIT License

## [0.0.1] - 2026-01-29

### Added
- Initial project structure
- Core module scaffolding
- Basic documentation

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| 0.5.0 | 2026-04-13 | Molecular QC layer, GNN embeddings, autism interpretation, KG topology |
| 0.4.0 | 2026-03-04 | Clinical validation (NB13-19), network biology, Docker, RRID/bio.tools |
| 0.3.1 | 2026-02-23 | Repository URL migration from GitHub to Codeberg |
| 0.3.0 | 2026-02-16 | Multi-omic integration: single-cell, deconvolution, signaling databases, cross-modal validation, fusion |
| 0.2.3 | 2026-02-14 | Cross-cohort validation, performance optimization, threshold calibration, expression scoring |
| 0.2.0 | 2026-02-09 | Scientific rigor, ancestry/batch correction, benchmarks, sensitivity analysis |
| 0.1.0 | 2026-01-29 | First public release with full pipeline |
| 0.0.1 | 2026-01-29 | Initial project setup |

[Unreleased]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.5.0...HEAD
[0.5.0]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.4.0...v0.5.0
[0.4.0]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.3.1...v0.4.0
[0.3.1]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.3.0...v0.3.1
[0.3.0]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.2.3...v0.3.0
[0.2.3]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.2.0...v0.2.3
[0.2.0]: https://codeberg.org/pathways/pathway-subtyping-framework/compare/v0.1.0...v0.2.0
[0.1.0]: https://codeberg.org/pathways/pathway-subtyping-framework/releases/tag/v0.1.0
[0.0.1]: https://codeberg.org/pathways/pathway-subtyping-framework/releases/tag/v0.0.1
