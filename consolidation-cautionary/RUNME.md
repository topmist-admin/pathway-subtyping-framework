# RUNME — third-party reproduction index (cautionary-framework paper)

**Purpose.** One entry point from which an independent reviewer can regenerate
every number in the cautionary-validation paper. Each package below is
self-contained (script + deposited reference output + README explaining what the
numbers do and do not support), uses **public data only**, and is deterministic
at **seed 42**.

> ⚠️ **Preserve gene column order.** Determinism at seed 42 holds *for a given
> gene ordering*. ssGSEA ranks genes within each sample, and where expression
> values are exactly **tied** the ranking falls back to column position — so
> permuting the gene columns can change scores even with the seed fixed. This
> bites here because the pipeline computes `log2(CPM + 1)` and `log2(0 + 1) = 0`,
> turning every zero-count gene into an exact tie. On a realistic bulk matrix
> processed that way, **29% of entries are exact zeros**, and permuting gene order
> moved scores by up to 1.68 with sample-ordering Spearman **0.70** (measured
> 2026-07-29). On dense continuous data with no ties, ssGSEA is exactly invariant
> (max |Δ| = 0.0000); mean-Z is invariant either way.
>
> **What this means for you.** Reading the deposited matrices as shipped
> reproduces the deposited numbers exactly — that is why the independent
> clean-room run was byte-identical. But re-parsing from GEO, sorting genes
> alphabetically, or joining a different annotation can reorder columns and shift
> the scores. If your numbers differ, check gene order before suspecting anything
> else. This is a property of ssGSEA tie-handling, not of the data or the seed.

**Scope note.** This file indexes the *new* cautionary-framework evidence
(Results 1–4 below). The sibling [`README.md`](README.md) is the older, separate
reproduction package for the **"Stable but confounded"** paper (GSE28521 /
GSE64018 / GSE80655). The two are not interchangeable.

> ✅ **v0.9.0 is published.** Released 2026-07-29 — `pip install
> pathway-subtyping==0.9.0` resolves from PyPI, and the pins below are live rather
> than staged. Verified from a clean venv against the published wheel.

**One version, one install: `pathway-subtyping==0.9.0`.** Everything in this
bundle — both this file's packages and the sibling README's — installs and runs
against v0.9.0, which is a superset of v0.7.0. Where the older scripts and their
READMEs mention v0.7.0, that is **provenance** (the release under which their
deposited reference outputs were generated), *not* an instruction to install it.
Do not downgrade: v0.7.0 lacks `pathway_subtyping.discreteness` and
`pathway_subtyping.clustering_dl`, so Results 1–3 cannot run on it at all.

As the sibling README already notes for its Layer B, exact partition
reproduction is version-sensitive — those deposited numbers were produced under
earlier releases, and the invariant to check is the *conclusion*, not a
bit-identical partition.

> 📄 **Re-running everything from source?** Read
> [`DATA-ACQUISITION.md`](DATA-ACQUISITION.md) first. It lists what actually needs
> downloading (almost nothing — six packages reproduce from deposited inputs), the
> measured runtime of every job, and **the acceptance criterion per result**. That last
> point matters: byte-identity is the right bar for the six deposited packages and the
> *wrong* bar for anything re-fetched from a live API, where the invariant is the
> conclusion. Applying the wrong bar manufactures disagreements.

---

## Install (public — everything runs from PyPI)

The framework is on PyPI as **`pathway-subtyping==0.9.0`** (the version these packages
need — it contains `pathway_subtyping.discreteness` (Gate A) and
`pathway_subtyping.clustering_dl` (the DL baselines), which the v0.7.0 line did not).

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt      # pins pathway-subtyping==0.9.0 + numpy/pandas/sklearn/scipy/statsmodels/requests
# torch is pinned in requirements.txt (==2.11.0) and needed ONLY by cancer_r38's DL
# baselines. They are not deterministic across torch releases; the version actually
# used is recorded in each result under provenance.environment.torch.
python -c "import pathway_subtyping as p; print(p.__version__)"   # 0.9.0
```

**Verified reviewer reproduction (2026-07-25):** from a clean virtual environment with
only the PyPI wheel installed (no local source tree), every no-network package
reproduces its deposited headline numbers exactly — benchmark audit (22 valid rows),
honest ablation (McNemar p=0.001, SigClust-p-alone == composite gate), donor-level
flagship (region V 0.660, diagnosis permutation p 0.234), and the framework-gate
flagship stability (0.921). The network packages (calibration, cancer, sweep) fetch
public cBioPortal/GEO/recount3 data with no authentication.

- **PyPI:** https://pypi.org/project/pathway-subtyping/0.9.0/
- **Source (tag v0.9.0):** GitHub `topmist-admin/pathway-subtyping-framework` ·
  Codeberg `pathways/pathway-subtyping-framework` · RRID:SCR_028051
- **This reproduction bundle, citable:** **`10.5281/zenodo.18638048`** — the **concept
  DOI**, which always resolves to the current version. Cite this, not a version DOI
- **Corrected 47-dataset benchmark:** `10.5281/zenodo.19323753` — the **concept DOI**,
  which always resolves to the newest version of the record. Cite this, not a version
  DOI. The analysis in this bundle read version 2.0 (`10.5281/zenodo.21262112`), and
  the deposited `benchmark_audit.json` records that version DOI deliberately, as
  provenance of the exact file it read

---

## Packages → paper Results

| Paper section | Package | Deposited output | Network |
|---|---|---|---|
| **Result 1** ablation (honest three-way + head-to-head) | [`cross-domain/gate_ablation/`](cross-domain/gate_ablation/) | **`ablation_honest.json`** (authoritative), `separation_sweep.json`, `gate_ablation_raw.csv`, figure | none (synthetic) |
| **Result 1** real-data calibration (**within-study**) | [`cross-domain/gate_calibration/`](cross-domain/gate_calibration/) | **`gate_calibration_within_study.json`** (authoritative; the pooled `gate_calibration.json` is a withdrawn batch artifact) | cBioPortal — or **none** with `--cache-in`, which replays the deposited inputs |
| **Result 2** benchmark audit | [`cross-domain/benchmark_audit/`](cross-domain/benchmark_audit/) | `benchmark_audit.json` (incl. column-validity diagnostic) | none (reads deposited CSV) |
| **Result 4** flagship donor-level stats | [`cross-domain/flagship_stats/`](cross-domain/flagship_stats/) | `flagship_donor_level.json` | none (reads deposited labels) |
| **Result 3** cancer worked example | [`cross-domain/cancer_r38/`](cross-domain/cancer_r38/) | `brca_pam50_validation.json`, `cptac_brca_multiomic.json` | cBioPortal |
| **Result 4** autism subgroup (GSE28521) | [`results/autism_subgroup/`](results/autism_subgroup/) | `autism_subgroup_result.json` + `autism_subgroup_genes.csv.gz` | **none** — re-runs offline via `--genes` |
| **Result 4** psychiatry flagship | [`README.md`](README.md) (deposited outputs produced under v0.7.0; runs on v0.9.0) + [`genetic-anchoring/`](genetic-anchoring/) | see that README | GEO |
| Large-N calibration point | [`cross-domain/gtex_brain/`](cross-domain/gtex_brain/) | `gtex_brain_region_confound.json` | **none** — the analysis reads the deposited `gtex_brain_pathway_scores.tsv`; recount3 (R) is only needed to *regenerate* that input |
| Scoping (negative result) | [`cross-domain/psychiatric_meta/`](cross-domain/psychiatric_meta/) | `track_a_recount3.tsv` | NCBI E-utilities |
| Gate-6 domain remap | [`cross-domain/`](cross-domain/) | `results/confound_remap_results.json` | none (seeded) |
| Gate-7 somatic anchoring (real TCGA-CRC positive control: BRAF-V600E / KRAS / MSI) | [`cross-domain/tcga_crc/`](cross-domain/tcga_crc/) | `tcga_crc_somatic_result.json` | cBioPortal |

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
- TPR moves **1.000 → 0.967** (one `discrete_k3` replicate lost); exact McNemar on
  positives **p=1.0**, so the cost is not distinguishable from zero and the study
  cannot exclude one. Report both figures; write neither "TPR held" nor a settled
  "3% cost".
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

1. ~~**Publish the framework**~~ **DONE 2026-07-29** — `pathway-subtyping==0.9.0` is on
   PyPI and tagged `v0.9.0` on GitHub, Codeberg and Bitbucket. v0.9.0 is a correctness
   release (six statistical fixes); **no manuscript number changes** — all offline
   packages are byte-identical under it. Re-verified from a clean PyPI install.
2. ~~**Environment pin**~~ **DONE** — `requirements.txt` now pins
   `pathway-subtyping==0.9.0` and the analysis deps; it covers the whole bundle.
3. ~~**Zenodo deposit**~~ **DONE** — deposited under concept DOI
   **`10.5281/zenodo.18638048`**, which always resolves to the current version. The
   manuscript cites the **concept** DOI, not a version DOI, so the citation survives
   future deposits without an edit.
4. ~~**Result 2 writeup**~~ **DONE 2026-07-23** — `cross-domain/benchmark_audit/`
   plus the Result 2 section of the rebuild draft.
5. ~~**Rewrite the abstract**~~ **DRAFTED 2026-07-23** in the manuscript working copy
   (maintained outside this repository), carrying none of the retracted figures.
   Still needs the hostile-review rounds before submission.
6. **Re-issue the benchmark correction.** *(Decided 2026-07-29 — description-only
   v2.1; data unchanged.)* Two separate problems were bundled under this item:

   - **The record's description carried a retracted claim.** The v2.0 description
     ended by asserting that reproducibility is "low on most independent cohorts" —
     the same claim Result 2 withdraws. After the repo README was corrected, the
     Zenodo record was the *only* public artifact still making it, and the one a
     reader reaches by following the DOI. **Fixed by v2.1, published 2026-07-29**
     (version DOI `10.5281/zenodo.21694795`; the concept DOI resolves to it). Text
     archived at [`CORRECTION_2026-07/ZENODO_v2.1_description.md`](../CORRECTION_2026-07/ZENODO_v2.1_description.md).
   - **The erratum's ground-truth rule is too weak** — 7 near-singleton rows remain
     marked valid. **Disclosed, not fixed in the data.** Re-issuing the CSV would
     change `benchmark_audit`'s input and break the byte-identical reproduction
     certified for this bundle, for no gain: the audit already applies the ratio
     screen in code and the retraction conclusion holds under it (R²≈0.001). The
     manuscript discloses the gap in Availability.

   Citations now use the **concept DOI** `10.5281/zenodo.19323753` throughout, so
   they survive v2.1 and any later re-issue without edits. Version DOIs are retained
   only where they record provenance — notably `benchmark_audit.json`, which names
   the exact version its input file came from.

   Zenodo DOIs cannot be deleted, only superseded; v1.0, v1.1 and v2.0 remain
   permanently resolvable and nothing cites them.
