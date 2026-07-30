# Zenodo v2.1 description (paste into the new version of the record)

**PUBLISHED 2026-07-29** as version DOI `10.5281/zenodo.21694795`.
**Concept DOI:** `10.5281/zenodo.19323753` — cite this, not a version DOI. It now
resolves to v2.1 (verified).
**Supersedes:** v2.0 (`10.5281/zenodo.21262112`), v1.1 (`10.5281/zenodo.19324360`),
v1.0 (`10.5281/zenodo.19323754`). All remain accessible; Zenodo DOIs are permanent
and published records cannot be deleted, only superseded.

**The benchmark data are byte-identical to v2.0.** `corrected_benchmark_47datasets_v2.csv`
(md5 `e4c0cbae…`) and `threshold_model_real47_CORRECTED.json` (md5 `4fa35748…`) are
unchanged, so anything computed against v2.0 remains valid without re-running — including
the certified byte-identical `benchmark_audit` reproduction. `ERRATUM_2026-07-08.md` was
refreshed (md5 `0f0743c0…`) so the deposited copy names the concept DOI explicitly, as the
repo copy does; it is documentation and feeds no analysis.

---

## Paste-ready description

**Version 2.1 — description corrected (2026-07-29).** The data files are unchanged
from v2.0; this version corrects the record's own description, which overstated what
the corrected benchmark supports.

The v2.0 description closed with the claim that the benchmark "shows that
discrete-subtype reproducibility is dataset-dependent and low on most independent
cohorts — a negative methodological result." **That claim is withdrawn.** It rests on
the `bootstrap_ari_5th_percentile` column, and that column will not carry it, for two
reasons:

1. It is a 5th *percentile*, not a point estimate, so any fixed bar applied to it —
   0.5, or the 0.80 used elsewhere — is near-zero by construction. A low pass rate is
   a property of the statistic, not a finding about subtypes.
2. The column is internally inconsistent with the silhouette column. Well-separated
   clusters are trivially reproducible under resampling, so a high silhouette forces a
   high stability ARI. Yet the record contains rows such as GSE66351 (silhouette
   0.967, bootstrap ARI −0.0095) and GSE16759 (silhouette 0.958, bootstrap ARI
   −0.1394), which are not physically consistent. The column does not measure
   bootstrap stability of the reported partition.

**No reproducibility-rate claim should be drawn from this record, in either
direction.** The withdrawal is not replaced with another rate.

What the record **does** support is unchanged and remains the reason to cite it: the
adaptive bootstrap-threshold model previously reported at R²=0.889 **does not
reproduce and is retracted.** Refit on the released data it gives R²≈0.111 across all
rows, ≈0.015 on valid rows, and ≈0.001 under a stricter ground-truth screen, with the
slope reversing sign. Silhouette does not predict bootstrap-ARI reproducibility in
these data, and no silhouette-calibrated threshold is supported by them.

The corrections carried forward from v2.0 also stand: fourteen rows have degenerate
ground truth (`n_true_clusters=1`) and are flagged invalid — the empty-input
adjusted-Rand-index artifact gave 13 of them a spurious ARI=1.0, the 14th (GSE136196)
returned 0.0; one row had impossible labels; the claim that the benchmark excludes the
manuscript analysis dataset (TCGA-COAD) was incorrect and is withdrawn; and the counts
are restated as 41 unique datasets, not 47, with 40,778 samples in passing rows, not
36,551.

**Known limitation, disclosed rather than fixed.** The erratum's ground-truth rule
(`n_true_clusters > n_samples`) is too weak. It catches only GSE92332 (1,957 labels /
533 samples) and misses rows where the label count approaches the sample count —
GSE2109 (2,158/2,158), GSE5204 (79/79), GSE17537 (55/54), GSE44228 (72/69), GSE42127
(176/145), GSE66351 (190/96), GSE16759 (16/8). Seven rows still marked valid fail a
ratio screen of `n_true_clusters / n_samples >= 0.5`. The data files are left as-is so
that analyses run against v2.0 stay reproducible bit-for-bit; the stricter screen is
applied in the analysis code instead, and the retraction conclusion holds under it
(R²≈0.001). Users applying their own screen should apply the ratio test, not the
erratum rule.

See `ERRATUM_2026-07-08.md` in this record for the full correction notice.

---

## Metadata changes to make alongside the description

| Field | From | To |
|---|---|---|
| Version | `2.0` | `2.1` |
| Title | "…for Pathway Subtyping Framework v0.4.0" | "…for Pathway Subtyping Framework" — drop the stale version pin; the record outlived v0.4.0 and is now read by v0.9.0 |
| Related identifier | — | `10.5281/zenodo.18638048` (concept DOI of the software), relation **isSupplementTo** |
| Related identifier | — | `https://pypi.org/project/pathway-subtyping/0.9.0/`, relation **isSupplementedBy** |

**Do not** add a relation asserting the benchmark was *produced by* v0.9.0. It was
not — it was produced under v0.4.0 and is only *read* by the current release. The
audit that re-analyses it ships in the v0.9.0 reproduction bundle.
