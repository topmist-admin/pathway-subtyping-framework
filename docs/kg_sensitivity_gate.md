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
random. Node set and edge-type composition are always preserved; candidate
additions are drawn from nodes already participating in that edge type, so every
proposal is schema-valid.

### Choosing the null: `rewiring=`

What *else* the null preserves sets how strict it is. This is the single most
consequential parameter in the gate — it can change the verdict on identical data.

| mode | preserves | when to use |
|---|---|---|
| `uniform` | node set, edge-type counts | comparison baseline only; the weakest null |
| `degree` **(default)** | + exact in/out-degree sequence | general use |
| `module` | + the diff's within-/cross-module split | when a release concentrated its edits in a few pathways |

**`degree`** performs Maslov–Sneppen double-edge swaps: two edges `a→b` and `c→d`
become `a→d` and `c→b`, so every node keeps its exact in- and out-degree. Hubs
stay hubs.

**Each swap changes TWO edges per side**, so a budget of *(remove, add)* affords
`min(remove, add) // 2` swaps and each completed swap consumes two units from
each side. A consequence worth knowing: only an **even** balanced budget is fully
swap-able. An odd remainder — and any imbalance, as when a release mostly grows —
falls through to a degree-*weighted* residual, which preserves the shape of the
degree distribution but not its exact values. Strict preservation is impossible
there, because adding an edge necessarily raises somebody's degree.

> This was wrong in the first two releases of this module: each swap was counted
> as one removal plus one addition, so `degree` and `module` perturbed at **2× the
> requested size**. That silently defeats the size-matching the whole argument
> rests on — an over-perturbed null pushes null agreement down, which raises
> `P(null ≤ observed)` and therefore **suppresses the `kg-sensitive` verdict**.
> Found by adversarial review, fixed 2026-07-29, pinned by
> `test_degree_mode_respects_the_requested_budget`.

**`module`** additionally constrains swaps to reproduce the observed diff's
within-module / cross-module ratio, computed by `within_module_fractions()`.
Modules come from `GENE_IN_PATHWAY` membership via `module_map_from_pathways()` —
two genes share a module when they share a pathway — so no external community
detection is needed. Use it when a release concentrated its edits in
actively-researched pathways, so that a `kg-sensitive` verdict cannot be
explained away as "the real change was concentrated and the null was not." Gate K
**abstains** if the graph has no `GENE_IN_PATHWAY` edges, rather than silently
falling back.

### Why `uniform` is not the default

`uniform` destroys the degree sequence, and that misreads changes at *both* ends
of the degree distribution:

- **Hub edges** — removing one of a hub's many edges is an ordinary consequence
  of any curation, because hubs have more edges to lose. Uniform rewiring treats
  it as rare and so overstates how special the real change was.
- **Peripheral edges** — uniform rewiring hits them constantly, making a genuinely
  targeted peripheral change look unremarkable.

`tests/test_kg_sensitivity.py::test_null_choice_can_change_the_verdict` pins the
second case: the same cohort, the same real KG change, an **identical** observed
ARI of −0.026, and opposite verdicts — `generically-fragile` under `uniform`
(p = 0.066) versus `kg-sensitive` under `degree` (p = 0.016).

**Report which null you used.** `KGSensitivityResult.rewiring` records it and
`to_dict()` serializes it, because a verdict is not interpretable without it.
Running both modes and reporting the pair is more informative than either alone.

## Reading the result

Whichever null you choose, it gives the same reading rule:

| observed vs null | verdict | meaning |
|---|---|---|
| observed ≥ `ari_min` | `robust` | survives this KG swap |
| observed ≪ null | `kg-sensitive` | the *specific* curated change matters |
| observed ≈ null | `generically-fragile` | any perturbation of this size breaks it |

The middle and bottom rows are the point of the gate, and **they are not
distinguishable from the agreement number alone**. `test_kg_sensitivity.py`
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
| `not-testable (...)` | **no** | graphs identical, partition degenerate, sample counts mismatched, the null could not be built, or **the rewiring could not actually perturb the graph** |

> **The null is measured, not assumed.** `rewire_kg` skips work it cannot do — no
> schema-valid exemplar for an edge type, or a graph too dense for free slots —
> and logs that at DEBUG. A draw can therefore come back *identical* to v1. If it
> does, every null agreement is 1.0, the observed value looks extreme by
> comparison, and the gate emits a confident and entirely spurious
> `kg-sensitive`. Gate K now measures the achieved perturbation per draw and
> **abstains** when the median is zero; `null_perturbation_requested` and
> `null_perturbation_median` are on the result and in `to_dict()`, and a null that
> achieves under half the requested change logs a warning. The commonest trigger
> is a release that *introduces* an edge type the baseline lacks.

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
    rewiring="degree",         # "uniform" | "degree" (default) | "module"
)

print(res.verdict)             # "robust" | "kg-sensitive" | "generically-fragile" | "not-testable (...)"
print(res.observed_ari, res.null_median, res.empirical_p)
print(res.diff.summary)        # what actually changed between the two graphs
print(res.rewiring)            # which null produced that verdict — always report it
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

The null holds the *number* of changed edges fixed. It does not model **which**
edges a curator would plausibly touch beyond degree and module structure — the
`degree` and `module` modes address the two largest such biases, but a release
also concentrates edits by evidence type, publication recency and organism, none
of which the null represents. A `kg-sensitive` verdict under `module` is
considerably harder to explain away than one under `uniform`, but it is still a
statement about *this* null, not about all plausible curations.

Finally, it inherits `partition_fn`'s determinism. A partition function with
unseeded internal randomness will produce a null that measures its own noise
floor. There is no way for the gate to detect that on the caller's behalf.
