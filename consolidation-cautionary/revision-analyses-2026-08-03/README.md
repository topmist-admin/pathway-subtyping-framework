# Revision analyses (2026-08-03/04) — Result 5, and the new numbers in Results 1, 3 and 4

These analyses were added after an independent pre-submission review. They are **newer than
the 2026-07-30 third-party reproduction round**, so they carry the same public-data,
deterministic-seed discipline as the rest of the package but have **not** themselves been
independently reproduced. That is stated in the manuscript's Data & Code Availability section
and is repeated here so nobody infers otherwise from their presence in the deposit.

Everything here uses public data only. Seed 42 throughout — **and, since 2026-08-05, genuinely so**: five scripts previously derived per-condition seeds from `hash(<str>)`, which Python salts per process, so their outputs could not be regenerated. They now use `zlib.crc32`. Results affected by that fix were re-run and the numbers here are the re-run values. No file in this directory writes to
the rest of the deposited tree.

## Running

Scripts resolve the repository root from their own location, so they run unchanged from an
unpacked deposit:

```bash
cd <deposit-root>
python -m venv .venv && . .venv/bin/activate
pip install -r consolidation-cautionary/requirements.txt
python consolidation-cautionary/revision-analyses-2026-08-03/scripts/<job>.py --out /tmp/out.jsonl
```

Set `PSF_REPO` to override root detection. Two jobs have extra requirements — see
**Dependencies** below.

## What each job establishes, and the number to check

| job | question | key deposited result |
|---|---|---|
| **`job1_sigclust_comparator.py`** | cluster-index statistic vs bootstrap-ARI, **same null, same datasets, same seeds** | the two statistics are **indistinguishable**: TPR 30/30, FPR 0/30 for *both*. (Deposited earlier with `--skip-ari`, i.e. one arm only; re-run 2026-08-05 with both.) |
| `job1b_sigclust_separation_sweep.py` | same, across the separation sweep | cluster index certifies **0.00 at every δ** |
| **`job1c_canonical_sigclust.py`** | **canonical CRAN `sigclust`**, ablation grid | **TPR 27/30, FPR 1/30**, no abstentions |
| **`job1e_sigclust_heteroscedastic.py`** | canonical `sigclust` on **job5's heteroscedastic generators**, importing job5's `gen()` and reproducing its seeding verbatim | **0/20 at n=120, 300 and 600** — SigClust never certifies where the screen reaches 1.00. Added 2026-08-05; this column was previously quoted with no deposited source. |
| **`job1d_canonical_sigclust_sweep.py`** | canonical `sigclust`, separation sweep | **0.00 at every δ**; median p 0.615 at δ=3.0 |
| `job2_covariance_shrinkage_sweep.py` | ⛔ **SUPERSEDED — do not cite** | λ=0 gave 0.30 where the deposited sweep gives 0.55; **failed its own parity control**. Retained deliberately: it is the worked example of why job 2b patches one value inside the released code path instead of reimplementing it. |
| **`job2b_shrinkage_exact_parity.py`** | is the screen's conservatism an estimation artifact? | **parity OK at all 7 δ**; λ=1 gives 0.05 / 0.35 / 0.90 at δ=0.0 / 2.5 / 3.0 vs released 0.00 / 0.15 / 0.55 |
| `job3_nsweep_generator_diversity.py` | ⛔ **SUPERSEDED — do not cite** | its `curved_manifold` "negative" is bimodal (dip rejects unimodality in 50–70% of replicates), so its 4/5 "false certifications" were mislabelled true positives. Retained as the record of a discarded finding. |
| **`job3b_nsweep_validated_negatives.py`** | FPR and abstention vs *n*, on **dip-validated** negatives | abstention 0.40→0.13 and false certification **0.03→0.17** across n=60/120/300 |
| **`job4_donor_continuous_null.py`** | donor-level diagnosis association without the degenerate null | deposited V null has **exactly 3 atoms**; 44/48 donors split; **GEE p=0.067** |
| **`job5_heteroscedastic_fpr_confirmation.py`** | is the heteroscedastic false certification real? | **0.15 / 0.80 / 1.00** at n=120/300/600 vs matched control 0.00 / 0.00 / 0.05 |
| **`job6_brca_certification_stats.py`** | the TCGA-BRCA certification statistics | k=5: obs 0.4275, p95 0.2458, **p=0.0050 (floor)**; k=2: obs 0.8852, p95 0.6607, **p=0.0050** |
| **`job6b_pc1_diagnostics.py`** | PC1 dip + standardised gap for the two cached controls, **with the feature count recorded per cohort** | TCGA-LGG (50 features) dip 0.9963, gap **1.359**; TCGA-LUAD (9-pathway immune panel) dip 0.9739, gap 1.887 — **not on a common basis with LGG/BRCA and not comparable**. Added 2026-08-05; the LGG gap was previously quoted with no deposited source. |
| **`job7_column_order_invariance.py`** | does gene column order change Result 4? | partition **ARI 1.000** under canonical sort *and* reversal; stability 0.9540 / 0.9538 / 0.9535 |

**The two superseded jobs are kept on purpose.** Both produced clean-looking, reportable
numbers that were wrong, and both were caught by a control rather than by inspection. Deleting
them would remove the evidence that the checks did their job.

## Dependencies beyond `requirements.txt`

- **`job1c` and `job1d` require R ≥ 4 with the CRAN `sigclust` package**
  (`install.packages("sigclust")`). This is the only non-Python dependency in the paper. They
  shell out via `Rscript --vanilla` and skip cleanly with a printed error if R is absent.
- **`job3b` and `job5` require `diptest`**, which is the `discreteness` extra rather than a
  core dependency because it is GPLv2+ while the framework is MIT:
  `pip install "pathway-subtyping[discreteness]==0.9.0"`. ⚠️ **Without it, `dip_of()` returns
  NaN rather than raising**, so a missing install silently disables the validation screen that
  makes these results meaningful. Both jobs will then treat every generator as qualified.

## External inputs

- **`job6`** fetches TCGA-BRCA from the public cBioPortal API by delegating to the deposited
  `cross-domain/cancer_r38/scripts/fetch_and_validate_brca.py` — the pathway matrix is
  therefore identical to Result 3's by construction, not by reimplementation.
- **`job7`** needs the GSE80655 **raw counts** supplementary file, which is a different
  acquisition route from the cached pathway scores used elsewhere in this package:
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80655/suppl/GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz`
  (14.3 MB), passed via `--expr`. **That file is pinned**: its SHA-256 is deposited at
  `inputs/GEO_SOURCE_PROVENANCE_job7.json` and the script checks it at runtime, warning if GEO
  has revised the file in place. The script also hashes the assembled scoring matrix.
  The symbol→Ensembl mapping it needs is **deposited** at
  `inputs/symbol_to_ensembl_scz_panel.json`, so no `mygene.info` call is required; that file
  was built from mygene.info (228 of 231 panel symbols mapped) and is included precisely so
  this job does not depend on a third-party API staying up.

  ⚠️ Panel coverage matters here: an incomplete 5-of-14-set mapping gives partition ARI 0.979
  rather than 1.000. If you rebuild the mapping yourself, check the coverage line the script
  prints before reading the verdict.

## Provenance notes

- Every output is newline-delimited JSON with a leading `provenance` record carrying the
  script name and full argv, so a partial run still leaves an interpretable file.
- Local filesystem paths were removed from two `provenance` records before deposit; the
  acquisition routes are documented above instead.
- `job2b` asserts parity against the deposited separation sweep per-δ and writes `parity_ok`
  into its output. **If `parity_ok` is false, nothing else in that file is usable** — that is
  the standard job 2 failed.
