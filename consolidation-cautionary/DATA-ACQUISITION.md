# Data acquisition, runtime budget, and what counts as a successful reproduction

One page for a reviewer who intends to re-derive everything from source. It answers
three questions the per-package READMEs answer only in pieces: **what do I download**,
**how long will it take**, and **what result means I succeeded**.

---

## 1. What you actually need to download

**Nothing, for most of it.** Six packages reproduce byte-identically from inputs already
inside the bundle. Downloads are needed only to *re-derive* what those inputs were built
from, and for the three cBioPortal packages.

### 1a. Already deposited — no download

| Input | Path | Size | Used by |
|---|---|---|---|
| Autism gene matrix (32 × 9,147) | `results/autism_subgroup/autism_subgroup_genes.csv.gz` | 1.4 MB | Result 4 autism arm |
| GTEx brain pathway scores (n=2,931) | `cross-domain/gtex_brain/results/gtex_brain_pathway_scores.tsv` | 2.7 MB | large-N calibration |
| Flagship partition + metadata | `data/partition/sample_metadata_with_subtypes.csv` | 16 KB | Result 4 flagship |
| Corrected benchmark | `../CORRECTION_2026-07/corrected_benchmark_47datasets_v2.csv` | 5 KB | Result 2 |
| Ablation raw grids | `cross-domain/gate_ablation/results/*.csv` | 12 KB | Result 1 |
| recount3 candidate table | `cross-domain/psychiatric_meta/results/track_a_recount3.tsv` | 4 KB | scoping result |
| Hallmark panel (50 sets) | `panels/hallmark_200genes.gmt` | 48 KB | autism arm, Result 3 |
| Schizophrenia panel (14 sets) | `panels/schizophrenia_pathways.gmt` | 4 KB | Result 4 flagship |

### 1b. Live downloads — needed only for the network packages

| Source | Identifier | Needed by | Auth |
|---|---|---|---|
| cBioPortal | `lgg_tcga_pan_can_atlas_2018` | `gate_calibration` (discrete control) | none |
| cBioPortal | `luad_tcga_pan_can_atlas_2018` | `gate_calibration` (continuum control) | none |
| cBioPortal | `brca_tcga_pan_can_atlas_2018` | `cancer_r38` (Result 3) | none |
| cBioPortal | `brca_cptac_2020` | `cancer_r38` multi-omic | none |
| cBioPortal | `coadread_tcga_pan_can_atlas_2018` | `tcga_crc` (Gate 7) | none |
| NCBI E-utilities | psychiatric study search | `psychiatric_meta` candidate list only | none |
| NCBI GEO | GSE28521, GSE80655 | only to *rebuild* deposited inputs | none |
| recount3 (R) | GTEx brain | only to *rebuild* the deposited score matrix | none |

**GSE64018 is not in this table on purpose.** It appears in bundle material inherited
from the withdrawn methodology paper and supports no result reported in the manuscript.

> ⚠️ **Do not re-fetch GEO expecting a match.** GEO revises series matrices and platform
> annotation files **in place, with no version bump**. The autism arm records the SHA-256
> of both source files in `GEO_SOURCE_PROVENANCE.json` precisely so drift is detectable —
> but the reproducible route is the deposited matrix, not a re-fetch.

> ⚠️ **cBioPortal is not a versioned API.** A study can gain or lose samples between
> fetches. Every network result now records endpoint, study, UTC fetch date, post-filter
> shape, and a SHA-256 of the assembled matrix under `provenance.fetch`. **If your numbers
> disagree, compare `matrix_sha256` first** — a different hash means the upstream data
> moved, not that the analysis is wrong.

---

## 2. Runtime budget

Measured on a 2023-class laptop. Nothing below hangs; three jobs are simply long.

| Job | Wall clock | Network |
|---|---|---|
| `benchmark_audit` | seconds | none |
| `flagship_donor_level` | seconds | none |
| `flagship_stability` | ~1 min | none |
| `autism_subgroup` | ~2 min | none |
| `gtex_brain` | **>25 min** | none |
| `ablation_honest` | **~77 min** | none |
| separation sweep (20 reps/step) | **~1 h 47 min** | none |
| `gate_calibration` | **~16 min** | **yes** |
| `tcga_crc` | ~10 min | **yes** |
| `cancer_r38` (BRCA) | **~59 min** | **yes** |

**Offline total ≈ 3.5 hours**, of which ~3.2 is the three long offline jobs. Add roughly
another 1.5 hours for the network packages. Run everything long detached and check the
fast packages meanwhile.

> **Do not run the cBioPortal packages concurrently.** Two large `molecular-data/fetch`
> POSTs in parallel throttle each other badly — measured here, two jobs that each finish
> alone made no progress at all for 45 minutes when run together. Run them one at a time.

> `cancer_r38` is by far the slowest and produces no output between the fetch line and
> the final result. A silent hour is expected, not a hang — measured at ~59 minutes.

---

## 3. What counts as a successful reproduction

The bar is **not the same for every result**, and applying the wrong one manufactures
disagreements. Use this table.

| Result | Package | Bar | Why this bar |
|---|---|---|---|
| Result 1 ablation | `gate_ablation` | **byte-identical** | synthetic data, fixed seed, no external input |
| Result 2 benchmark audit | `benchmark_audit` | **byte-identical** | reads a deposited CSV |
| Result 4 flagship stats | `flagship_stats` | **byte-identical** | reads a deposited partition |
| Result 4 autism arm | `results/autism_subgroup` | **byte-identical** | reads the deposited gene matrix |
| large-N calibration | `gtex_brain` | **byte-identical** | reads deposited pathway scores |
| Result 1 calibration | `gate_calibration` | **conclusion** | live cBioPortal; data may drift |
| Result 3 cancer | `cancer_r38` | **conclusion** (in practice exact, except the two bootstrap ARIs) | live cBioPortal; DL baselines vary across torch releases, but not at the pinned version |
| Gate 7 somatic | `tcga_crc` | **conclusion** | live cBioPortal |
| scoping count | `psychiatric_meta` | **byte-identical** offline | count recomputable from the deposited TSV |

**"Byte-identical"** means `sha256` of your result JSON equals the deposited one. Six
packages meet this and were verified from a clean room. If one of them differs, that is a
genuine finding — report it.

**"Conclusion"** means these must hold; the digits need not:

- `gate_calibration` — the discrete control (IDH-mut vs IDH-wt in LGG) is **certified**;
  the continuum control (LUAD immune gradient) is **not certified**. The gate's direction
  is the claim, not the p-value.
- `cancer_r38` — PSF is **competitive with, not better than**, the DL baselines (VAE-GMM
  edges PSF on enrichment, 89.1% vs 87.6%); k=5 bootstrap stability **fails** the 0.80
  bar; recovery is **metric-dependent**. Certified 2026-07-30: all four k-way ARIs and
  all PAM50 enrichments re-derived **exactly**, DEC and VAE-GMM included, so with the
  pinned torch this package is in practice much tighter than its "conclusion" bar. The
  two bootstrap-stability ARIs are the exception — they move by up to 0.008 between runs
  (0.399/0.436 vs 0.408/0.435). Quote them to two decimals.
- `tcga_crc` — the somatic association is **detected** (BRAF-V600E/KRAS/MSI strata), with
  the expected-count guard not tripped.

**Known-benign causes of small numeric differences**, to rule out before reporting:

1. **Gene column order** — ssGSEA ranks within a sample, and `log2(CPM+1)` maps every
   zero count to an exact tie (~29% of a bulk matrix). Re-parsing or re-sorting shifts
   scores. Read deposited matrices as shipped. Mean-z is invariant.
2. **torch version** — DL baselines are not deterministic across releases. Pinned at
   2.11.0 and recorded in `provenance.environment.torch`.
3. **SciPy version** — a newer SciPy returns `nan` for the LUAD dip where the deposited
   notes record 0.98. The gate decision is unchanged; it is driven by the SigClust p.
4. **Two dip p-values per calibration record** — `pc1_dip_p_standardized_pathways`
   z-scores pathway columns before PCA; `gateA.dip_pc1_p` uses the gate's centred-only
   reduction. They differ by design and neither decides the gate.

---

## 4. Order of work

1. Install and confirm `pathway_subtyping.__version__ == "0.9.0"`.
2. Run the four fast offline packages (~3 min) and diff against the deposit. If these do
   not match byte-for-byte, stop — something is wrong with the environment, and the long
   jobs will not fix it.
3. Launch `ablation_honest` and `gtex_brain` detached.
4. Run the three network packages; check `provenance.fetch.matrix_sha256` against the
   deposited value before comparing any number.
5. Run the separation sweep last, if at all — it is the longest job and supports a
   supporting figure, not a headline.
