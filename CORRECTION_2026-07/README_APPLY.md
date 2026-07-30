# How to apply this correction (author actions — not auto-executed)

> **Status (2026-07-29):** steps 1 and 2 are **done** — the repo carries the erratum and
> the annotated retracted model, and Zenodo **v2.0 (`10.5281/zenodo.21262112`)** was
> published on 2026-07-08 under the concept DOI `10.5281/zenodo.19323753`. The **description-only v2.1** is
> also published (2026-07-29, `10.5281/zenodo.21694795`), withdrawing the
> reproducibility claim the v2.0 description still carried — see
> [`ZENODO_v2.1_description.md`](ZENODO_v2.1_description.md). Only step 3's downstream
> updates remain. The DOI in the step-2 heading below is v1.1, the record as it stood when this
> was written.

This folder stages the correction. Nothing here has been committed, pushed, or re-deposited.
All of the following mint permanent/public changes and are the author's deliberate actions.

## 1. Repository (Codeberg / GitHub mirrors)
- Replace `research-results/benchmarks/bootstrap_threshold_calibration_47datasets_final_zenodo.csv`
  with `corrected_benchmark_47datasets_v2.csv` (or add alongside + deprecate the old).
- Replace / annotate `src/pathway_subtyping/threshold_model_real47.json` using
  `threshold_model_real47_CORRECTED.json` (keep the retraction visible).
- Add `ERRATUM_2026-07-08.md` at repo root; link it from `README.md` and `KNOWN-ISSUES.md`.
- Remove/deprecate any code path that consumes the retracted silhouette model
  (the shipped pipeline already uses a different sample-size×k table; document that too).
- Commit message suggestion: "correction: retract adaptive-threshold model (R^2=0.889 not
  reproduced); flag degenerate benchmark rows; withdraw independence claim (see ERRATUM)."
- Push deliberately, normal batch size (GitHub topmist-admin had a prior bulk-edit block).

## 2. Zenodo (DOI 10.5281/zenodo.19324360)
- Zenodo DOIs are permanent; you CANNOT delete the old version.
- Publish a NEW VERSION under the concept DOI containing: the corrected CSV, the retracted
  model file, and ERRATUM_2026-07-08.md. Use the description in `ZENODO_v2_description.md`.
- The old version remains citable but is superseded; add a note pointing to the new version.
- Optional: request a tombstone for the old version from Zenodo support only if you want it
  formally retracted rather than superseded (heavier signal, permanent tombstone).

## 3. Downstream
- Biosketch (ARP IDA v3): benchmark cite already marked "under revision" -> update to the new
  version DOI once minted.
- rs-9284565: withdrawn. rs-8913089: rebuild (see REBUILD_rs8913089 plan).
