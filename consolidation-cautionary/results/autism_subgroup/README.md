# §3.1 autism frontal-cortex subgroup — deposited artifacts

Supports the GSE28521 arm of manuscript **Result 4**. Before this deposit those
numbers existed only as prose; they are now reproducible offline.

| File | What it is |
|---|---|
| `autism_subgroup_result.json` | the headline numbers + full provenance |
| `autism_subgroup_labels.csv` | per-sample cluster and diagnosis |
| `autism_subgroup_genes.csv.gz` | the derived 32 × 9,147 gene matrix (with `dx`) |

## Reproduce offline — no GEO, no network

```bash
python ../../scripts/reproduce_autism_subgroup_meanz.py \
  --genes autism_subgroup_genes.csv.gz \
  --gmt ../../panels/hallmark_200genes.gmt \
  --out /tmp/autism_check
```

Expected, and verified identical to the GEO-fetched run: 32 frontal samples
(16 ASD / 16 control), 50/50 Hallmark pathways scored by mean-z, BIC **k=2**,
disease-enriched cluster **n=11 with 9/11 ASD**, silhouette **0.25**, bootstrap
ARI **0.71 — fails the 0.80 bar**, random-gene-set ARI **0.237** over 40 draws.

The gene matrix is deposited because the random-gene-set control resamples from
the full gene universe, so pathway scores alone cannot regenerate it.

## ⚠️ Two things a reviewer will otherwise trip on

**The seed here is 20260708, not 42.** That is this analysis's documented seed
(see the script's module docstring). Methods states seed 42 *with this
exception noted*. Re-running at 42 produces different numbers — it does not
check these ones.

**Do not confuse this with `research-results/GSE28521/frontal-cortex/`.** That is
a *different* analysis and it matches on four fields, so it is an easy mistake:
same n=32, same 16/16 ASD/control, same best k=2. But it scores **15 curated
autism pathways** (not 50 Hallmark) under **framework 0.3.0 at seed 42**, and
reports silhouette **0.299**, not 0.25. It belongs to the withdrawn methodology
paper, not to Result 4.

## Provenance

`autism_subgroup_result.json` records the seed, the framework version, the number
of random draws, and the SHA-256 of every input — including the two GEO files, so
anyone who re-fetches can detect upstream drift rather than silently obtaining
different numbers. GEO revises series matrices and `.annot` files in place without
a version bump, which is exactly why the derived matrix is deposited here.
