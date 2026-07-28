# Time-sliced KG sensitivity — Gate K

## The question no other gate asks

Every gate in the existing battery holds the knowledge graph fixed and varies the
data: resample the cohort, permute the labels, draw from a reference
distribution. None of them varies the **knowledge graph** itself.

That leaves a whole class of artifact uninspected. Pathway definitions,
regulatory edges and drug–target links are curated products that change between
database releases. A subtype that appears under Reactome 2026 and disappears
under Reactome 2025 is a property of the knowledge base, not of the patients. It
will pass bootstrap stability, pass the discreteness gate, and pass confound
association — because all three are computed against the one KG version that
produced it.

Gate K varies the knowledge graph and asks whether the finding survives.

## Why a raw agreement number is uninterpretable

The obvious implementation is to partition the cohort under KG v1, partition it
under KG v2, and report the ARI between the two labelings. That number means
nothing on its own, for exactly the reason the pre-v0.8.0 stability gate meant
nothing on its own: there is no reference telling you what to expect.

An agreement of 0.6 between two KG versions is consistent with two opposite
readings:

- the curated change was substantive and specifically disrupted this finding; or
- this partition is fragile to **any** perturbation of comparable size, and the
  KG version is incidental.

The first is a statement about the knowledge base. The second is a statement
about the partition, and a considerably worse one — it says the result would not
survive a change of *any* kind, curated or arbitrary.

> The observed agreement is fine; without a reference distribution it cannot be
> read.

## The size-matched rewiring null

`rewire_kg` perturbs the baseline graph at random while matching the observed
diff's **per-edge-type** addition and removal counts. It changes the same number
of edges of the same types as the real curation did, but chooses which ones at
random. Node set and edge-type composition are preserved; candidate additions are
drawn from nodes already participating in that edge type, so every proposal is
schema-valid.

That gives the reading rule:

| observed vs null | verdict | meaning |
|---|---|---|
| observed ≥ `ari_min` | `robust` | survives this KG swap |
| observed ≪ null | `kg-sensitive` | the *specific* curated change matters |
| observed ≈ null | `generically-fragile` | any perturbation of this size breaks it |

The middle and bottom rows are the point of the gate. **They are not
distinguishable from the agreement number alone** — `tests/test_kg_sensitivity.py`
constructs two cases with an identical observed ARI of −0.026 that receive
opposite verdicts, separated only by where the null sits (median 1.000 versus
−0.026).

## What decides the verdict

```python
if observed_ari >= ari_min:    verdict = "robust"
elif empirical_p < alpha:      verdict = "kg-sensitive"
else:                          verdict = "generically-fragile"
```

**Both terms are live.** `ari_min` decides robustness; `empirical_p` decides how a
non-robust result is *explained*. Nothing else enters the rule — the `KGDiff` and
the optional `run_kg_regression` report are reported as context and are
explicitly **not** criteria.

This is stated plainly because Gate A shipped documenting three named criteria of
which only one decided; the gap statistic and dip test are computed, reported and
inert, and that gap between prose and `passed = obs > sg_p95 and sg_p < alpha` was
found by hostile review rather than by us. A regression test asserts that a
flagged scalar score does **not** change the partition verdict. If a future
revision wants it to, the rule above must change in the same commit as this page.

`empirical_p` is an add-one empirical p in the **low** tail: how often does a
size-matched random rewiring disrupt the partition at least as much as the real
KG change did.

## Four outcomes, not two

| verdict | testable | meaning |
|---|---|---|
| `robust` | yes | agreement at or above `ari_min` |
| `kg-sensitive` | yes | disrupted more than size-matched chance |
| `generically-fragile` | yes | disrupted no more than chance |
| `not-testable (...)` | **no** | graphs identical, partition degenerate, sample counts mismatched, or the null could not be built |

`KGSensitivityResult.testable` is the flag that keeps abstentions out of failure
rates. Gate A abstained on 28 of 30 synthetic negatives and the resulting
"FPR 0.000" was quoted against n=30 rather than the testable n=2; the `testable`
field exists so that mistake is mechanical to avoid here. **Any rate computed
over many Gate K runs must use the testable subset as its denominator.**

## Usage

```python
from pathway_subtyping import kg_timeslice_sensitivity

def partition_fn(kg, cohort):
    """Must be deterministic given its inputs — seed any internal clustering."""
    scores = score_pathways_with_kg(kg, cohort)
    return cluster(scores, k=3, seed=42)

res = kg_timeslice_sensitivity(
    kg_reactome_2025,          # v1 — baseline
    kg_reactome_2026,          # v2 — comparison
    partition_fn,
    cohort,
    n_null=50,                 # size-matched random rewirings
    ari_min=0.80,
    seed=42,
)

print(res.verdict)             # "robust" | "kg-sensitive" | "generically-fragile" | "not-testable (...)"
print(res.observed_ari, res.null_median, res.empirical_p)
print(res.diff.summary)        # what actually changed between the two graphs
```

Pass `score_fns=[...]` to get a `run_kg_regression` report alongside the
partition verdict, for the scalar view. It is reported, never decisive.

## Synthetic validation

`tests/test_kg_sensitivity.py` verifies on constructed ground truth that the gate
returns `robust` when the partition ignores the graph, `kg-sensitive` when a
single curated edge drives the partition and random rewiring rarely touches it,
`generically-fragile` when the partition responds to edge *count* so every
perturbation disrupts it equally, and abstains when the graphs are identical or a
partition is degenerate. It also checks that `rewire_kg` preserves the node set
and edge count, does not mutate its input, and reproduces exactly under a seed.

## Limits to state whenever this gate is cited

It tests the KG versions **you supply**. A `robust` verdict means the finding
survived the specific swap tested, not that it is invariant to knowledge
curation in general — two adjacent releases that barely differ will return
`robust` almost by construction, which is why `res.diff.summary` should be
reported alongside every verdict.

The null holds the *number* of changed edges fixed but not their **importance**.
Curated changes are not random: a real release concentrates edits in
actively-researched regions of the graph, so a `kg-sensitive` verdict may partly
reflect that concentration rather than the finding's fragility. Degree-preserving
or module-aware rewiring would tighten this and is the obvious next refinement.

Finally, it inherits `partition_fn`'s determinism. A partition function with
unseeded internal randomness will produce a null that measures its own noise
floor. There is no way for the gate to detect that on the caller's behalf.
