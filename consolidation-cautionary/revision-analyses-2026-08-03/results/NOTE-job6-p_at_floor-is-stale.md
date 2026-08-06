# Correction: `p_at_floor` is stale in `job6_brca_certification.jsonl`

**Filed 2026-08-06. The manuscript is correct; this deposited field is not.**

Both `record: "gate"` lines in `job6_brca_certification.jsonl` carry:

```json
"sg_empirical_p": 0.005, "p_at_floor": false, "p_floor": 0.004975124378109453
```

`p_at_floor: false` is **wrong**. Both p-values *are* at the floor, which is what the
manuscript states ("**0.0050** (at floor)" in the Result 3 certification table, and "both
*p*-values sit at the floor 1/(n_ref+1)").

## Why the field is wrong

The run was executed with `n_ref=200`, so the floor is `1/201 = 0.0049751…`. The gate
stores the empirical p **rounded to 4 dp** (`gate_a_discreteness_null.py`:
`sg_empirical_p=round(sg_p, 4)`), giving `0.005`. The original comparison tested the
*rounded* value against the *unrounded* floor:

```
0.005 <= 0.0049751...   ->   False
```

so a p-value sitting exactly on the floor was reported as not on it.

Because `sg_p = (#{ref >= obs} + 1) / (n_ref + 1)`, the only count that yields a stored
`0.005` is 0 exceedances: `1/201 = 0.0049751 -> 0.005`. The next attainable value is
`2/201 = 0.00995 -> 0.0100`. So `0.005` is reachable *only* at the floor, and the correct
value of the flag is **`true`**.

Compared at matched precision:

```
round(0.005, 4) <= round(1/201, 4) + 1e-12   ->   0.005 <= 0.005   ->   True
```

Cross-check against the parallel real-data record: `gate_calibration_within_study.json`
reports `sg_p 0.0066` with `sg_p_at_floor: true` at `n_ref=150` (`1/151 = 0.006623`) —
the same situation, flagged correctly, because that value does not collide with the
rounding boundary.

## Why the artifact was not regenerated

`job6_brca_certification_stats.py` was corrected on 2026-08-05 — it now compares at
matched precision, and the comment above the expression records exactly this trap. The
deposited JSONL predates that fix.

Re-running the job would require a **live cBioPortal fetch** (there is no cached BRCA
matrix in the deposit), so regenerating it offline is impossible and regenerating it
online would silently mix a 2026-08-06 data pull into a 2026-08-03 run record. The
artifact is therefore left as the true record of what that run emitted, and this note
carries the correction.

**Every other field in those two records is correct** and matches the manuscript:
`observed_bootstrap_ari` 0.4275 / 0.8852, `sg_ref_p95` 0.2458 / 0.6607,
`dip_pc1_p` 0.9925, `std_gap_pc1` 1.6529 / 0.4326.

## Note on scope

`job6` is **not** covered by the offline byte-identical reproduction suite recorded in the
submission package's `PROVENANCE` §3 — that suite covers the offline packages. This is
why the stale field survived earlier verification passes.
