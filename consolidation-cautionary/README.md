# Reproducibility package — *"Stable but confounded: brain region, not diagnosis, drives a validation-passing molecular subtype in postmortem psychiatric brain"*

Regenerates every number in the manuscript from **public** inputs and the **public** framework release. This package is designed so an outside scientist can reproduce the paper with nothing but what is downloadable from PyPI / Codeberg / Zenodo / GEO.

Status: **draft package for the pre-post reproduction pass** (intended to become a public deposit — Codeberg subfolder + Zenodo DOI — at manuscript finalization).

---

## 1. Environment (public framework only)

```
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt        # installs pathway-subtyping==0.9.0 from PyPI
```

Do **not** import a local working tree of the framework. Verify the release:
`python -c "import pathway_subtyping as p; print(p.__version__)"` → `0.9.0`.
(v0.9.0 is a superset of v0.7.0 — the scoring, GMM, bootstrap and confound gates this
package uses are unchanged; it additionally ships the discreteness gate used by the
cautionary-framework packages. Layer A below is deterministic and version-independent.)
Provenance: PyPI https://pypi.org/project/pathway-subtyping/0.9.0/ · Codeberg (primary) https://codeberg.org/pathways/pathway-subtyping-framework (tag `v0.9.0`) · GitHub https://github.com/topmist-admin/pathway-subtyping-framework (tag `v0.9.0`) · RRID:SCR_028051.

## 2. Inputs to fetch (public)

| What | Where | Notes |
|---|---|---|
| **GSE28521** (Voineagu 2011) | GEO | `GSE28521_series_matrix.txt.gz` + `GPL6883.annot.gz` |
| **GSE64018** (Gupta/Irimia 2014) | GEO | `GSE64018_countlevel_12asd_12ctl.txt.gz`, `GSE64018_adjfpkm_12asd_12ctl.txt.gz`, series matrix |
| **GSE80655** (Ramaker 2017) | GEO | `GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz` + series matrix |
| **Pathway panels** | shipped in `panels/` | `schizophrenia_pathways.gmt` (curated 14-set, SCZ analyses) + `hallmark_200genes.gmt` (50 Hallmark sets, §3.1 autism) — also in the public repo `data/pathways/` |
| **Corrected benchmark** (§3.5) | Zenodo concept DOI 10.5281/zenodo.19323753 (resolves to the current version) | `corrected_benchmark_47datasets_v2.csv` |

Put the GEO files under a data root with `GSE28521/ GSE64018/ GSE80655/` subfolders. Record every file's SHA-256.

**⚠️ Use the 3-level `brain region` (AnCg/DLPFC/nAcc) for GSE80655** — never a 2-level "ACC" coding (see §4).

## 3. Scripts (in `scripts/`)

| Script | Manuscript coverage | Framework? |
|---|---|---|
| `reproduce_autism_subgroup_meanz.py` | §3.1 autism frontal-cortex subgroup (k=2, fails bootstrap stability) | yes (mean-z, GMM/BIC, bootstrap gate) |
| `reproduce_consolidation_paper.py` | Fig 1 (cross-diagnosis region gate) + Fig 2 (ASD/SCZ marker tables) | yes (ssGSEA, GMM/BIC, confound gate) |
| `reproduce_crossdiagnosis_region_gate.py` | Table 2 (§3.3) at k=3 **and** BIC-k + stability@BIC-k | yes |
| `double_dissociation_interaction_test.py` | §3.4 marker×diagnosis interaction (the double dissociation) | no (OLS only) |
| `scz_within_region_marker_panel.py` | §3.2 per-region SST/GABA panel (from raw GSE80655) | no (OLS only) |
| `region_crosstab_from_labels.py` | §3.2 Fig-1 headline (χ²=125.12, V=0.666) — **Layer A, exact** | no (pandas+scipy only) |

Example invocations (adjust `<DATA>` and the panel/labels paths):

```
python scripts/reproduce_autism_subgroup_meanz.py \
  --gse28521-matrix <DATA>/GSE28521/GSE28521_series_matrix.txt.gz \
  --gse28521-annot  <DATA>/GSE28521/GPL6883.annot.gz \
  --gmt hallmark_200genes.gmt --out out/

python scripts/region_crosstab_from_labels.py --labels data/partition/sample_metadata_with_subtypes.csv

python scripts/reproduce_consolidation_paper.py \
  --gse80655-expr  <DATA>/GSE80655/GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz \
  --gse80655-matrix <DATA>/GSE80655/GSE80655_series_matrix.txt.gz \
  --gse64018-counts <DATA>/GSE64018/GSE64018_countlevel_12asd_12ctl.txt.gz \
  --gse64018-adjfpkm <DATA>/GSE64018/GSE64018_adjfpkm_12asd_12ctl.txt.gz \
  --gse64018-matrix <DATA>/GSE64018/GSE64018_series_matrix.txt.gz \
  --gse28521-matrix <DATA>/GSE28521/GSE28521_series_matrix.txt.gz \
  --gse28521-annot  <DATA>/GSE28521/GPL6883.annot.gz \
  --gmt schizophrenia_pathways.gmt --out out/

python scripts/reproduce_crossdiagnosis_region_gate.py \
  --gse80655-expr <DATA>/GSE80655/...Data...txt.gz \
  --gse80655-matrix <DATA>/GSE80655/GSE80655_series_matrix.txt.gz \
  --gmt schizophrenia_pathways.gmt --out out/

python scripts/double_dissociation_interaction_test.py --data-dir <DATA>
```

(The Figure-1 scripts call the Ensembl REST API once to map panel symbols → Ensembl IDs; needs network.)

## 4. The region-coding trap (must get right)

GSE80655's derived `brain_region` field can carry a **corrupted 2-level coding** (nAcc merged into "ACC"). The manuscript headline **requires the 3-level `brain region` (AnCg/DLPFC/nAcc)** re-derived from the original GEO metadata. On the 2-level coding the same partition gives only χ²=34.6, V=0.50 (wrong). `region_crosstab_from_labels.py` warns if it sees <3 levels.

## 4b. §3.1 scoring: mean-z on the microarray, not ssGSEA (must get right)

The §3.1 autism analysis scores **GSE28521 (Illumina microarray)** against the **50 Hallmark sets** (`hallmark_200genes.gmt`, gene symbols). Use **mean-z** (`ExpressionScoringMethod.MEAN_Z`), *not* ssGSEA: on this array platform the current public release's ssGSEA variance filter culls all but ~3 Hallmark pathways (47/50 become zero-variance and are removed), leaving too few features to cluster and selecting a spurious high k. Mean-z scores all 50 sets and recovers the k=2 partition. `reproduce_autism_subgroup_meanz.py` uses mean-z and the framework's own bootstrap gate (n_bootstrap=100), so its stability ARI is directly comparable to the SCZ 0.923. Expected: k=2, disease-enriched n=11 (9/11 ASD), silhouette 0.25, bootstrap ARI **0.71 → fails (<0.80)**, random-gene-set ARI **0.237** (40 draws). **Deposited** at `results/autism_subgroup/` — `autism_subgroup_result.json`, the cluster labels, and `autism_subgroup_genes.csv.gz`, the 32x9,147 gene matrix, so the whole analysis including the random-gene-set control re-runs OFFLINE via `--genes` with no GEO access. *(All framework summaries z-normalize pathway scores; the matrix is unit-variance — §2.2.)*

## 5. Two layers of reproduction (label this in your memo)

- **Layer A — deterministic.** `region_crosstab_from_labels.py` recomputes the Fig-1 headline from the **deposited partition labels** (`data/partition/sample_metadata_with_subtypes.csv`, framework v0.3.0 / seed 42). Must reproduce **exactly**: χ²=125.12, p=4.30×10⁻²⁶, V=0.666 (0.660 Bergsma); diagnosis p=0.408; composition 2/0/37 · 45/47/0 · 1/1/8. *(Verified 2026-07-10 against a public v0.7.0 install; v0.9.0 is a superset and is the version to install today.)*
- **Layer B — from-scratch pipeline.** The framework scripts re-cluster from raw public data with the public **v0.9.0** release (the deposited run used v0.7.0). Because the deposited primary analysis used framework v0.3.0 (seed 42) and v0.7.0 differs, you will **not** reproduce the exact seed-42 partition — expected, and itself a documented result. What must hold is the *conclusion*: region V ≈ 0.67–0.72, diagnosis n.s., across seeds and k.

## ⚠️ Known reproducibility gaps (report these, don't paper over them)

1. ~~`scz_within_region_marker_panel.py` consumes processed intermediates~~ **RESOLVED 2026-07-11**: the script now builds everything from **raw public GSE80655** (Ensembl expression + covariates incl. brain pH from the series matrix) — no processed intermediates. §3.2 per-region SST/GABA numbers are reproducible from raw public data (SST AnCg *d* = −1.38, DLPFC −1.54; pan-GABA AnCg −1.20; PVALB n.s.; nAcc null).
2. **The deposited partition labels** (`data/partition/sample_metadata_with_subtypes.csv`) are shipped in this package for convenience, but the manuscript's "Data and Code Availability" implies they are publicly deposited — confirm they (and this whole package) are actually published to Codeberg/Zenodo before the paper is posted, or a reader cannot reproduce Fig 1.
3. This package is not yet a public deposit — it is the pre-post reproduction harness. Publishing it (Codeberg subfolder + Zenodo DOI, cited from the paper) is a manuscript-finalization step.

## Provenance
Every script writes intermediate CSVs + (where applicable) a `provenance.json` with input SHA-256s, seed (20260708), and package versions. Keep these with your reproduction memo.
