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

## Provenance

Input: `../../data/partition/sample_metadata_with_subtypes.csv` (framework v0.3.0
partition, seed 42). Donor id parsed from the `title` token (`X####_REGION_...`).
10,000 donor-level permutations, seed 42.
