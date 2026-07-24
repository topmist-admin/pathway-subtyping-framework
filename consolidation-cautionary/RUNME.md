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
| **Result 1** methods validation | [`cross-domain/gate_ablation/`](cross-domain/gate_ablation/) | `gate_ablation_results.json`, `clusterer_sweep_results.md`, figure | none (synthetic) |
| **Result 1** real-data calibration (**linchpin**) | [`cross-domain/gate_calibration/`](cross-domain/gate_calibration/) | `gate_calibration.json` | cBioPortal (public, no auth) |
| **Result 2** benchmark audit | Zenodo `10.5281/zenodo.21262112` (corrected 47-dataset benchmark v2.0) | `corrected_benchmark_47datasets_v2.csv` | Zenodo |
| **Result 3** cancer worked example | [`cross-domain/cancer_r38/`](cross-domain/cancer_r38/) | `brca_pam50_validation.json`, `cptac_brca_multiomic.json` | cBioPortal |
| **Result 4** psychiatry flagship | [`README.md`](README.md) (v0.7.0 package) + [`genetic-anchoring/`](genetic-anchoring/) | see that README | GEO |
| Large-N calibration point | [`cross-domain/gtex_brain/`](cross-domain/gtex_brain/) | `gtex_brain_region_confound.json` | recount3 (R) |
| Scoping (negative result) | [`cross-domain/psychiatric_meta/`](cross-domain/psychiatric_meta/) | `track_a_recount3.tsv` | NCBI E-utilities |
| Gate-6 domain remap | [`cross-domain/`](cross-domain/) | `results/confound_remap_results.json` | none (seeded) |

---

## Headline numbers, traced to source

Every figure below is read directly from the deposited artifact named in its row.
None is restated from prose.

**Gate ablation** — `cross-domain/gate_ablation/results/gate_ablation_results.json`
- `stability_only`: TPR 1.000, FPR **0.367**, continuum certification rate 0.733
- `stability+discreteness`: TPR 0.967, FPR **0.000**, continuum rate 0.000
- Honest reading: the discreteness gate removes continuum false positives **at a
  cost of ~3% of true positives** (TPR 1.000 → 0.967). Not "TPR held".

**Gate calibration (linchpin)** — `cross-domain/gate_calibration/results/gate_calibration.json`
- Discrete positive control (pooled COAD/GBM/LUAD, n=1262): recovery of
  tumor-of-origin ARI **0.921**; Gate A **certified** (SigClust p=0.0066)
- Continuum negative control (LUAD immune infiltration, n=510): immune-score
  Hartigan dip p **0.9817** (unimodal); Gate A **not certified** (verdict
  "continuum", dip on PC1 p=0.9972)
- Both calls correct → the gate is **calibrated, not merely pessimistic**. This is
  the answer to the hostile reviewer question the whole paper hinges on.

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

**Psychiatric meta-cohort scoping** — `cross-domain/psychiatric_meta/results/track_a_recount3.tsv`
- 57 candidates → 30 bulk-RNA-seq keepers → only **6 studies / 243 samples** in
  recount3; **~595** with the GSE80655 anchor (n=352) added. Not large-N.

---

## What this bundle does NOT contain

- **CommonMind (n=986)** and any other controlled-access psychiatric cohort. Access
  is gated on controlled-access approval. When that analysis runs, controlled data cannot be
  redistributed — the package will ship the **code plus access instructions**, not
  the data, following the pattern of the packages above.
- Any claim that the framework *finds* strong validated subtypes. The paper's
  contribution is the calibrated instrument and the honest audit; see
  `cross-domain/README.md` and each package's honesty notes.

## Open items before this bundle is submission-ready

1. **Publish v0.8.0** — hard prerequisite (see blocker above). Currently untagged
   and held.
2. **Add a cautionary-bundle pin file.** `requirements.txt` covers only the
   v0.7.0 package; the cross-domain packages have no pinned environment.
3. **Zenodo deposit** of this bundle, cited from the paper.
4. **Result 2 writeup** — the corrected 47-dataset benchmark is deposited on Zenodo
   but has no package here turning it into the "reproducibility is usually low"
   audit narrative.
5. **Rewrite the abstract** — the withdrawn manuscript's abstract still carries
   retracted numbers (R²=0.889, CMS4 75.9%, ARI 0.870, bootstrap 0.923). None may
   survive into the new paper.
