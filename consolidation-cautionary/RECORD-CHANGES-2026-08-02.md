# Files in this record were replaced on 2026-08-02

**This record's DOI (`10.5281/zenodo.21704904`) and version (v0.9.0) are unchanged. Its files
are not.** Zenodo's edit window allows a published record's files to be replaced without
minting a new version, and no version entry is created when that happens — so this note exists
to make the replacement visible to anyone who downloads the record after the fact.

## Why replace rather than mint v0.9.1

The principal change removes **retracted text** from a deposited file. Zenodo versions are
permanent: minting v0.9.1 would have left v0.9.0 resolvable forever with the retracted sentence
still inside it. Replacing the files removes it from the citable record, which is the point of
issuing a correction at all.

## What changed — 10 of 380 files; none added, none removed

| File | Change |
|---|---|
| `consolidation-cautionary/cross-domain/gate_ablation/results/clusterer_sweep_results.md` | **The reason for this replacement.** Removes the retracted claim *"Gate A rejects the 1-D continuum for 100% of runs"* and the two-way PASS/REJECT labelling that counted abstentions as rejections. Now labels every run CERTIFY / REJECT / ABSTAIN and reports *"declines to certify … 100% of runs … 60% explicit rejection and 40% abstention."* |
| `scripts/gate_ablation_study.py` | Adds `--render-only`, which regenerates the file above **from the deposited `clusterer_sweep_raw.csv` without re-clustering**. Also corrects a false statement in the script's own output text ("the abstentions are the larger share" — the split is 9 rejections to 6 abstentions). |
| `consolidation-cautionary/scripts/build_manuscript_figures.py` | Figure 1: legend renamed to "Discreteness screen"; the error bar on the screen's false-certify rate removed (no interval is honest at n=2 testable — see the manuscript caption); panel-A rates now read from `ablation_honest.json` instead of being hardcoded. Figure 3: method labels and panel titles no longer assert a ranking. |
| `REVIEWER-CONCERN-HANDLING-MATRIX.md` | Dispositions for R1.4, R2.2, R2.6 and R3.10 updated. |
| `consolidation-cautionary/scripts/build_submission_pdf.py` | Carries commit `18f803e` (strip the internal draft banner in the builder), which post-dated the original upload. |
| `consolidation-cautionary/scripts/build_full_reviewer_bundle.py` | Carries commit `d3fd745` (withhold self-audit documents from the reviewer bundle), which post-dated the original upload. |
| `consolidation-cautionary/FULL-RERUN-REVIEWER-PLAN.md` | Dated status banner added: three of its listed blockers are verified resolved; the rest are explicitly marked as carried forward unreviewed. |
| `CHANGELOG.md`, `consolidation-cautionary/RUNME.md` | Now include the DOI line that was committed just after the original tarball was built. |
| `zenodo-metadata.json` | Regenerated with the new stamp. |

**370 of 380 files are byte-identical to the 2026-07-30 upload.** No analysis result — no
JSON, no CSV, no deposited output — changed. Every number cited in the accompanying manuscript
is unaffected, and the corrected write-up above was re-rendered from deposited data rather than
recomputed.

## Independent certification of these changes

This replacement is **not** an unreviewed edit. Two independent reviews bear on it:

1. **2026-07-30 — clean-room reproducibility review.** Verdict GO; the associated reviewer
   bundle verified at 390/390 files by SHA-256; every headline number reproduced.
2. **2026-08-01 — independent re-certification of exactly these changes.** Verdict
   **CERTIFIED — GO.** It confirmed, by recomputation rather than inspection:
   - **all 72 deposited result files (.json/.csv/.tsv) byte-identical** to the validated
     version — zero changed, so no headline number moved;
   - the corrected write-up matches the unchanged raw CSV (15/15 discrete certified; 15/15
     continuum declined; 9 explicit rejections + 6 abstentions, split 3+2 per clusterer);
   - `clusterer_sweep_raw.csv` byte-identical before and after the re-render;
   - the figures value-invariant (Figures 2 and 4 pixel-identical);
   - the retracted figures still quarantined as disclosed retractions.

   Its assessment of this specific correction: *"a correction that removes an over-claim from
   the deposit; it strengthens rather than weakens the record."*

**Two files in this archive post-date that certification** and were therefore not reviewed by
it — stated plainly so the boundary is visible:

- `consolidation-cautionary/RECORD-CHANGES-2026-08-02.md` — this note.
- `consolidation-cautionary/FULL-RERUN-REVIEWER-PLAN.md` — a dated status banner was added on
  2026-08-02 marking three of its listed blockers as verified resolved (the manuscript's version
  pin, the absence of a v0.9.0 software deposit, and a Methods description of the flagship
  scoring) and marking every other item as carried forward unreviewed. That file is a planning
  document; it is not part of any reproduction path and no reported number depends on it.

Everything else in this replacement falls inside the certified set. Anyone re-running the
verification should expect the differences listed above and no others.
