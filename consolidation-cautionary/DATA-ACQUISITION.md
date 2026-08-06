# Data acquisition, runtime budget, and what counts as a successful reproduction

One page for a reviewer who intends to re-derive everything from source. It answers
three questions the per-package READMEs answer only in pieces: **what do I download**,
**how long will it take**, and **what result means I succeeded**.

---

## 1. What you actually need to download

**Nothing, for most of it.** Five packages are verified to reproduce byte-identically from inputs already
inside the bundle. Downloads are needed only to *re-derive* what those inputs were built from, and for the
network packages. ⚠️ Note `gate_calibration` appears on both sides: it fetches from
cBioPortal by default **but ships a deposited offline replay** (`--cache-in`), so it is one
of the six packages with an offline path. `job7` also needs a download (see §2).

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
| **NCBI GEO (raw counts)** | **GSE80655 supplementary** `GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz` (14.3 MB) | **`revision-analyses-2026-08-03/job7`** (Methods column-order test) — pass via `--expr` | none |
| cBioPortal | `brca_tcga_pan_can_atlas_2018` | **`revision-analyses-2026-08-03/job6`** (Result 3 certification statistics) — delegates to `cancer_r38`'s fetcher, so the matrix is identical by construction | none |

**The GSE80655 raw-counts row is a different acquisition route from the row above it.**
Elsewhere in this package GSE80655 is reached through *deposited pathway scores*; `job7` needs
the *gene-level counts* because it re-scores from scratch to vary gene column order. The same
warning applies with more force: **GEO revises supplementary files in place without a version
bump**, so a future download may not be byte-identical to the one used here. `job7` prints the
matrix shape, the number of matched samples and the per-set panel coverage before it reports
anything — check those three lines before comparing verdicts.

**No `mygene.info` call is required.** The symbol→Ensembl mapping for the 14-set panel is
deposited at `revision-analyses-2026-08-03/inputs/symbol_to_ensembl_scz_panel.json`
(228 of 231 symbols mapped). It was built from mygene.info and is included precisely so this
job does not depend on a third-party API remaining available. ⚠️ If you rebuild it, coverage
changes the answer: a 5-of-14-set mapping gives partition ARI 0.979 instead of 1.000.

**GSE64018 is not in this table on purpose.** It appears in the bundle's
double-dissociation material, which belongs to the companion "Stable but confounded"
paper, and supports no result reported in the manuscript.

> ⚠️ **Do not re-fetch GEO expecting a match.** GEO revises series matrices and platform
> annotation files **in place, with no version bump**. The autism arm records the SHA-256
> of both source files in `GEO_SOURCE_PROVENANCE.json` precisely so drift is detectable —
> but the reproducible route is the deposited matrix, not a re-fetch.

> ⚠️ **cBioPortal is not a versioned API.** A study can gain or lose samples between
> fetches. Every network result now records endpoint, study, UTC fetch date, post-filter
> shape, and a SHA-256 of the assembled matrix under `provenance.fetch`.
>
> ⚠️ **Correction (2026-08-07): the deposited results do NOT carry a `matrix_sha256`
> value.** The emitter that writes it (`scripts/_provenance.py`) postdates these runs, so
> the field is absent from every deposited artifact, and two `cancer_r38` results
> (`brca_pam50_validation_with_DL.json`, `cptac_brca_multiomic.json`) carry no `provenance`
> block at all. An earlier revision told reviewers to "compare `matrix_sha256` first";
> **that check cannot be performed on this deposit.** Compare the post-filter matrix shape
> and study identifier instead — a difference there means the upstream data moved, not that
> the analysis is wrong.

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
| `ablation_honest` | **~122 min** | none |
| separation sweep (20 reps/step) | **~1 h 47 min** | none |
| `gate_calibration` | **~16 min** live; offline `--cache-in` replay is longer (SigClust n_ref=150) | **optional** — offline via `--cache-in` |
| `tcga_crc` | ~10 min | **yes** |
| `cancer_r38` (BRCA) | **~59 min** | **yes** |

**Offline total ≈ 4.2 hours**, of which ~3.2 is the three long offline jobs. Add roughly
another 1.5 hours for the network packages. Run everything long detached and check the
fast packages meanwhile.

### Revision analyses (`revision-analyses-2026-08-03/`)

These were omitted from the table above until 2026-08-06. They are **the longest jobs in
the paper** — budget for them separately. Wall clocks are the `elapsed_sec` recorded in
each job's own deposited JSONL, not estimates.

| Job | Wall clock | Network | Feeds |
|---|---|---|---|
| `job1_sigclust` | **~26 min** | none | Result 5 comparator |
| `job1b_sigclust_sweep` | ~28 s | none | Result 5 separation sweep |
| `job1c_canonical_sigclust` | not recorded | none | Result 5 (needs R + CRAN `sigclust`) |
| `job1d_canonical_sweep` | not recorded | none | Result 5 (needs R + CRAN `sigclust`) |
| `job1e_sigclust_hetero` | not recorded | none | Result 5 (needs R + CRAN `sigclust`; 120 R invocations — long) |
| **`*_nrep100` variants of job1c/1d/1e** | as above | none | **the values Result 5 actually reports** — `--nrep 100` converges the k-means inside `sigclust`; the plain files are the `nrep=1` default and give TPR 27/30, FPR 1/30 |
| `job2_shrinkage` | ~1 h 57 min | none | **superseded by `job2b`** |
| `job2b_shrinkage_parity` | **~5 h 12 min** | none | Result 1 λ table — **longest job in the paper** |
| `job3b_nsweep_validated` | **~4 h 38 min** | none | Result 5 *n*-sweep |
| `job4_donor_continuous` | not recorded | none | Result 4 (needs `statsmodels`) |
| `job5_heteroscedastic_fpr` | **~4 h 31 min** | none | Result 5 headline table |
| `job6_brca_certification` | **~37 min** | **yes** | Result 3 certification table |
| `job6b_pc1_diagnostics` | seconds | none | Result 3 PC1 diagnostics |
| `job7_column_order` | ~2 s | **yes — a 14.3 MB GEO download** (`--expr` is required and the file is NOT deposited; see §1b) | Methods column-order test |

**Revision-analyses total ≈ 15.5 hours** for the recorded jobs (excluding the superseded
`job2` and the four unrecorded ones). Three jobs each exceed 4½ hours. Nothing here is
needed to *use* the screen — these are characterisation sweeps.

`job1c` / `job1d` / `job1e` require **R ≥ 4 with the CRAN `sigclust` package**; `job4`
requires `statsmodels`.

The three R jobs **probe for a working `Rscript` and a loadable `sigclust` before doing
any work**, and exit non-zero writing no output if either is missing. ⚠️ *This was added
2026-08-07 after an external review found the opposite: an earlier revision of this file
asserted they "fail loudly", but `which Rscript` was the only check, so an Rscript that
existed and errored let every call fall through to a `continue`. The job exited 0 and
wrote a complete-looking file with **zero certifications** — indistinguishable from Result
5's real finding. The claim was made without being tested; it is now true and tested.*

`job4` raises a bare `ModuleNotFoundError` if `statsmodels` is absent — it fails, and
writes nothing, but the message is a traceback rather than a written explanation.

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
| Result 1 calibration | `gate_calibration` | **offline** via `--cache-in` (byte-identity unconfirmed); **conclusion** on a live fetch | replays deposited `results/cached_inputs/`; live cBioPortal data may drift |
| Result 3 cancer | `cancer_r38` | **conclusion** (in practice exact, except the two bootstrap ARIs) | live cBioPortal; DL baselines vary across torch releases, but not at the pinned version |
| Gate 7 somatic | `tcga_crc` | **conclusion** | live cBioPortal |
| scoping count | `psychiatric_meta` | **conclusion** — NOT an offline reproduction | scripts query live NCBI E-utilities + recount3; only the headline count is re-derivable, by summing the deposited TSV, which is not regenerating the package |

**"Byte-identical"** means `sha256` of your result JSON equals the deposited one. **Five**
packages are verified to meet this from a clean room; a sixth, `gate_calibration`, ships an
offline `--cache-in` replay whose byte-identity is not independently confirmed.

> **Package count, settled 2026-08-06.** **Six** packages ship an offline path:
> `gate_ablation`, `benchmark_audit`, `flagship_stats`, `autism_subgroup`, `gtex_brain`,
> and `gate_calibration` (via `--cache-in`, replaying `results/cached_inputs/`). Five are
> verified byte-identical; `gate_calibration`'s offline path is verified but its
> byte-identity is not independently confirmed. **`psychiatric_meta` is NOT one of them** —
> its scripts query live NCBI E-utilities and the recount3 index, so only its headline
> count is re-derivable, by summing a deposited TSV. Earlier counts of "seven" and of
> "five" were both wrong.
 If one of them differs, that is a
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
  two bootstrap-stability ARIs are the exception — they move by up to 0.016 between deposited runs (0.3918/0.4338 vs 0.4076/0.4350). The cause was
  unseeded per-bootstrap GMM initialisation at the three cancer_r38 call sites, since fixed; a
  re-run will give one stable value differing slightly from both deposited figures.
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
4. Run the network packages; **note that `provenance.fetch.matrix_sha256` is absent from the
   deposited results** (see §1b), so compare the post-filter shape and study id against the
   deposited value before comparing any number.
5. Run the separation sweep last, if at all — it is the longest job and supports a
   supporting figure, not a headline.
