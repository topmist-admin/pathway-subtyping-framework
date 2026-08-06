# Flagship donor-level inference (GSE80655) — 48 donors, not 141 samples

The Result 4 flagship claims the pathway partition of postmortem psychiatric brain
tracks **brain region**, not **diagnosis**. A reviewer flagged that the naive test
cross-tabulates over 141 samples, but those come from only **48 donors** (~3 regions
each). Region varies within donor, so the sample-level unit is fair for region; but
diagnosis is fixed per donor, so a 141-sample χ² inflates the effective n and its
p-value is invalid.

**Run:** `scripts/flagship_donor_level.py` → `results/flagship_donor_level.json`
(deterministic, no network, reads the deposited partition labels).

## Result — the finding survives correct inference

| Factor | Unit | Cramér's V (Bergsma) | Test | Verdict |
|---|---|---|---|---|
| **Brain region** | sample (n=141, valid — region varies within donor) | **0.660** | χ²=125.1 | strong, real |
| Diagnosis — naive (INVALID) | sample (n=141) | 0.000 | p=0.408 | — |
| **Diagnosis — donor-level** | donor (n=48, majority partition per donor) | 0.157 | permutation p=**0.234** | not distinguishable from chance |

The donor-level permutation test permutes diagnosis **across the 48 donors** (the
correct exchangeable unit) and compares the observed Cramér's V to the null. Observed
V=0.258 sits essentially **at the null 95th percentile (0.258)** — i.e. an effect this
size arises by chance at this donor count. The region-not-diagnosis conclusion holds
under correct inference.

## What this fixes and what it does not

- **Fixes:** the manuscript can no longer be accused of comparing a within-donor factor
  (region, effective n≈141) against a between-donor factor (diagnosis) as if both had
  141 independent observations. Both effects are now reported on their correct footing,
  and region still dominates while diagnosis is null.
- **Honest limitation it exposes:** at n=48 donors the diagnosis test cannot *exclude*
  a modest true effect — the null 95th percentile for V is ~0.26. The correct claim is
  "diagnosis is not distinguishable from chance here," not "diagnosis has no effect."
  State it that way in the manuscript.
- Both Cramér's V values are now Bergsma-corrected (the earlier draft paired a corrected
  region V with an uncorrected diagnosis V — fixed).

## Partition stability (traceability of the "passes bootstrap stability" claim)

`scripts/flagship_stability.py` → `results/flagship_stability.json`. Recomputes the
k=3 partition's mean bootstrap ARI from the cached per-sample pathway scores
(`research-results/GSE80655/pathway_scores_scz.csv`) + the deposited labels:
**0.921** (n=141, n_bootstrap=100, seed 42), passing the 0.80 bar.

⚠️ **Provenance note — three values are on record for "GSE80655 stability", and they are not
interchangeable.** An earlier revision asserted that the original 0.9234 "is this number".
**That identity is withdrawn**: `results_summary.json` records no `n_bootstrap`, so it cannot
be checked even in principle, and the two runs are under different releases.

| value | what it is | source |
|---|---|---|
| **0.9206** | bootstrap ARI, n=141, k=3, `n_bootstrap=100`, seed 42, recomputed under the **released** version | `flagship_stability.json` — **this is the value the manuscript cites** |
| 0.9234 | same quantity, **original run** under an earlier release; `n_bootstrap` not recorded | `results_summary.json` → `validation_gates[2]` |
| 0.9540 | **a different quantity** — the shipped-gene-order scoring arm of the column-order test, not a competing estimate of the above | `job7_column_order.jsonl` → `stability_shipped` |

The first two agree to ~0.003 and both clear the 0.80 bar, which is the invariant that
matters; the conclusion does not depend on which is used. The third must not be compared with
them. The withdrawal was driven by the 47-dataset benchmark and the adaptive-threshold model —
**not** by this stability computation. It is a valid GSE80655 figure and is
distinct from BOTH withdrawn figures, which are two different results: the colorectal CMS4
recovery (**75.9%**, not an ARI) and the cross-pathway-set **ARI 0.870** on GSE80655.
("Cross-platform" is a retired descriptor for the latter — it is two curated panels on the
same cohort.) Do not list it among
disavowed figures.

## Provenance

Input: `../../data/partition/sample_metadata_with_subtypes.csv` (framework v0.3.0
partition, seed 42). Donor id parsed from the `title` token (`X####_REGION_...`).
10,000 donor-level permutations, seed 42.
