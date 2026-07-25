# RUNME — third-party reproduction index (cautionary-framework paper)

**Purpose.** One entry point from which an independent reviewer can regenerate
every number in the cautionary-validation paper. Each package below is
self-contained (script + deposited reference output + README explaining what the
numbers do and do not support), uses **public data only**, and is deterministic
at **seed 42**.

**Scope note.** This file indexes the *new* cautionary-framework evidence
(Results 1–4 below). The sibling [`README.md`](README.md) is the older, separate
reproduction package for the **"Stable but confounded"** paper (GSE28521 /
GSE64018 / GSE80655) and is pinned to the public **v0.7.0** release. The two are
not interchangeable — see the version blocker immediately below.

---

## ⚠️ BLOCKER — these packages need v0.8.0, which is not publicly released

Every cross-domain package here imports `pathway_subtyping.discreteness`
(Gate A) and, for the DL baselines, `pathway_subtyping.clustering_dl`. **Neither
module exists in v0.7.0**, which is the newest public release on PyPI / Codeberg
/ Zenodo. Verified:

```
git ls-tree -r --name-only v0.7.0 -- src/pathway_subtyping/discreteness   # empty
git ls-tree -r --name-only v0.7.0 -- src/pathway_subtyping/clustering_dl.py   # empty
```

Consequences a reviewer will hit:

- [`requirements.txt`](requirements.txt) pins `pathway-subtyping==0.7.0`. Following
  it and then running any cross-domain script fails at import.
- The **linchpin** gate-calibration result, the discreteness verdicts in the
  cancer and GTEx packages, and the entire ablation are unreproducible from a
  public install.

**Until v0.8.0 is tagged and published, the claim "third-party reproducible" holds
only for the v0.7.0-based `README.md` package, not for the cautionary-framework
evidence.** Publishing v0.8.0 (currently written, tested, untagged, and held) is
therefore a hard prerequisite for submission, not an independent release decision.

Interim instructions for a reviewer with repo access: install from source at the
commit that produced these outputs rather than from PyPI —

```bash
git clone <repo> && cd pathway-subtyping-framework
git checkout 39cd5a2          # head of the rebuild series (d13bec9..39cd5a2)
pip install -e .
python -c "import pathway_subtyping as p; print(p.__version__)"   # 0.8.0 (unreleased)
```

---

## Packages → paper Results

| Paper section | Package | Deposited output | Network |
|---|---|---|---|
| **Result 1** ablation (honest three-way + head-to-head) | [`cross-domain/gate_ablation/`](cross-domain/gate_ablation/) | **`ablation_honest.json`** (authoritative), `separation_sweep.json`, `gate_ablation_raw.csv`, figure | none (synthetic) |
| **Result 1** real-data calibration (**within-study**) | [`cross-domain/gate_calibration/`](cross-domain/gate_calibration/) | **`gate_calibration_within_study.json`** (authoritative; the pooled `gate_calibration.json` is a withdrawn batch artifact) | cBioPortal (public, no auth) |
| **Result 2** benchmark audit | [`cross-domain/benchmark_audit/`](cross-domain/benchmark_audit/) | `benchmark_audit.json` (incl. column-validity diagnostic) | none (reads deposited CSV) |
| **Result 4** flagship donor-level stats | [`cross-domain/flagship_stats/`](cross-domain/flagship_stats/) | `flagship_donor_level.json` | none (reads deposited labels) |
| **Result 3** cancer worked example | [`cross-domain/cancer_r38/`](cross-domain/cancer_r38/) | `brca_pam50_validation.json`, `cptac_brca_multiomic.json` | cBioPortal |
| **Result 4** psychiatry flagship | [`README.md`](README.md) (v0.7.0 package) + [`genetic-anchoring/`](genetic-anchoring/) | see that README | GEO |
| Large-N calibration point | [`cross-domain/gtex_brain/`](cross-domain/gtex_brain/) | `gtex_brain_region_confound.json` | recount3 (R) |
| Scoping (negative result) | [`cross-domain/psychiatric_meta/`](cross-domain/psychiatric_meta/) | `track_a_recount3.tsv` | NCBI E-utilities |
| Gate-6 domain remap | [`cross-domain/`](cross-domain/) | `results/confound_remap_results.json` | none (seeded) |

---

## Headline numbers, traced to source

Every figure below is read directly from the deposited artifact named in its row.
None is restated from prose.

**Gate ablation (authoritative: `ablation_honest.json`)** — three-way accounting
- Negatives (n=30): stability-only false-certifies **11** (FPR 0.367, Wilson 95%
  CI [0.22, 0.54]); the recalibrated gate certifies **0**. Paired exact McNemar
  b=11 c=0 **p=0.001** — real reduction.
- ⚠️ The gate reaches that by **abstaining on 28/30 (93%)** negatives; FPR excluding
  abstentions is 0/2, CI [0.00, 0.66] — nearly uninformative. Do NOT quote "FPR 0.000".
- No detectable TPR cost: 0.967 vs 1.000, McNemar **p=1.0**.
- Head-to-head: the SigClust p-value **alone** reproduces the composite gate exactly
  → the contribution is a null recalibration, not a new instrument.
- Separation sweep (`separation_sweep.json`, 20 reps/step): gate resolves (certify 0.00
  through δ≤2.0 → 0.15 at δ=2.5 → 0.55 at δ=3.0; median p 0.78→0.23→0.04) but is
  **conservative** — the transition is confined to δ=2.5–3.0. Script default `--reps 20`.
- The old `gate_ablation_results.json` "FPR 0.000 / TPR held" framing is SUPERSEDED.

**Gate calibration (authoritative: `gate_calibration_within_study.json`)**
- Discrete positive control, **within one study**: IDH-mut vs IDH-wt low-grade glioma
  (TCGA-LGG, n=507); recovery ARI **0.418** (within-study z-scores); Gate A **certified**.
- Continuum negative control: LUAD immune infiltration (n=510), dip p 0.98; **not certified**.
- ⚠️ The old pooled `gate_calibration.json` (ARI 0.921) is a **withdrawn batch artifact**
  (3 studies, 1 tumor type each; within-study z-score collapsed it to 0.05). Do not cite it.

**Cancer worked example** — `cross-domain/cancer_r38/results/brca_pam50_validation.json`
- TCGA-BRCA n=1082 (981 PAM50-labelled), k=5: PSF pathway-GMM recovery
  **ARI 0.218** vs **LumA single-subtype enrichment 87.6%** (OR 11.64, p 7.8e-47)
  — *the same partition*. Metric choice, not method quality, drives the apparent
  result. Single-subtype enrichment (what the field and the original manuscript
  reported) flatters subtyping; k-way ARI is the honest test.
- Do **not** claim PSF beats deep learning: VAE-GMM edges PSF on enrichment
  (89.1% vs 87.6%) while PSF leads on ARI. "Competitive" is the defensible claim.
- CPTAC-BRCA (n=122, `cptac_brca_multiomic.json`): protein modality ARI 0.172 /
  LumA enrichment 87.0%; mRNA ARI 0.189 / Basal 88.2%. Discreteness gate certified
  in **both** modalities; bootstrap stability **fails** the 0.80 bar in both
  (mRNA 0.457, protein 0.369). Expression↔protein concordance ARI 0.166.

**GTEx brain large-N** — `cross-domain/gtex_brain/results/gtex_brain_region_confound.json`
- n=2931, 13 regions. At k=13: region recovery **ARI 0.151** (weak), driven by
  cerebellum (enrichment 0.589, OR 22.14). At BIC k=3: ARI 0.074, Gate A **not
  certified** (continuum).
- Honest reading: this is **not** a clean "region dominates at large N" rebuttal.
  Cerebellum separates; cortical regions behave more like a continuum. Retained as
  a large-N calibration point and a third independent recurrence of the
  metric-dependence finding.

**Flagship donor-level (Result 4)** — `cross-domain/flagship_stats/results/flagship_donor_level.json`
- GSE80655: 141 samples = **48 donors**. Region (sample-level, valid): Bergsma V **0.660**.
  Diagnosis at donor level: permutation p **0.234** (indistinguishable from chance), obs
  V 0.258, null 95th 0.258. Region-not-diagnosis holds under correct inference.

**Psychiatric meta-cohort scoping** — `cross-domain/psychiatric_meta/results/track_a_recount3.tsv`
- 57 candidates → 30 bulk-RNA-seq keepers → only **6 studies / 243 samples** in
  recount3; **~595** with the GSE80655 anchor (n=352) added. Not large-N.

**Benchmark audit (Result 2)** — `cross-domain/benchmark_audit/results/benchmark_audit.json`
- Retracted adaptive-threshold model refit: R² **0.111** (all rows) / **0.015**
  (valid) / **0.001** (stricter screen), slope reversing sign between screens.
  Published claim was R² 0.889, slope +0.914. Reproduces the erratum exactly.
- Reproducibility across valid cohorts: median **−0.002**, max **0.391**,
  **0 of 22 reach 0.5** (0 of 15 under the stricter screen).
- ⚠️ The correction was incomplete: **7 of 22 "valid" rows** have a ground-truth
  label count ≥ half the sample count (GSE2109 and GSE5204 at exactly one label per
  sample), which the erratum's `n_true_clusters > n_samples` rule missed. Every
  headline is reported under both screens and none move.
- ⚠️ The statistic is a **5th percentile** (conservative lower bound). Claims must be
  about the gate criterion, not about biology. See that package's caveat section.

---

## What this bundle does NOT contain

- **CommonMind (n=986)** and any other controlled-access psychiatric cohort. Access
  is gated on the usual controlled-access approval. When that analysis runs, controlled
  data cannot be redistributed — the package will ship the **code plus access
  instructions**, not the data, following the pattern of the packages above.
- Any claim that the framework *finds* strong validated subtypes. The paper's
  contribution is the calibrated instrument and the honest audit; see
  `cross-domain/README.md` and each package's honesty notes.

## Open items before this bundle is submission-ready

1. **Publish v0.8.0** — hard prerequisite (see blocker above). Currently untagged
   and held.
2. **Add a cautionary-bundle pin file.** `requirements.txt` covers only the
   v0.7.0 package; the cross-domain packages have no pinned environment.
3. **Zenodo deposit** of this bundle, cited from the paper.
4. ~~**Result 2 writeup**~~ **DONE 2026-07-23** — `cross-domain/benchmark_audit/`
   plus the Result 2 section of the rebuild draft.
5. ~~**Rewrite the abstract**~~ **DRAFTED 2026-07-23** in
   `AIForAutismOutReach/.../REBUILD-DRAFT-2026-07-23-abstract-result2-result4.md`,
   carrying none of the retracted figures. Still needs the 3 hostile-review rounds
   before submission.
6. **Re-issue the benchmark correction.** The audit found the 2026-07-08 erratum's
   ground-truth rule (`n_true_clusters > n_samples`) is too weak — 7 rows with
   near-singleton labels remain marked valid. The conclusions are unaffected, but
   the deposited Zenodo v2.0 artifact should be superseded by a v2.1 applying the
   ratio screen, or the manuscript must disclose the gap. **PI decision.**
