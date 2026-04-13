# How It Works

A plain-language guide to what the Pathway Subtyping Framework does, why it works, and when to use each layer. No math, no API calls — just concepts.

---

## The Problem

Many diseases (autism, schizophrenia, cancer) are **genetically heterogeneous** — hundreds of different genes can be involved, and no two patients have the same mutations. Traditional gene-by-gene analysis fails because:

- Each gene has too few patients to reach statistical significance
- Different genes can disrupt the same biological process
- Patients with different mutations can have the same clinical presentation

## The Solution: Pathway-Level Aggregation

Instead of asking "which gene is broken?", we ask "which biological process is broken?"

**Pathway scoring** groups genes by their shared function. For example, the MAPK signaling pathway includes RAS, RAF, MEK, and ERK. A mutation in ANY of these genes disrupts the same cascade. By scoring at the pathway level, we combine signal from many rare variants into one meaningful measurement per patient per pathway.

**What it computes:**

```
Patient has variants in RAF and BRAF (both in MAPK pathway)
→ Weight by severity (LoF > missense > synonymous)
→ Weight by constraint (highly constrained genes matter more)
→ Sum within MAPK pathway → Z-score across all patients
→ Result: This patient's MAPK pathway score = 2.3 (elevated)
```

For expression data, the same idea applies: instead of looking at 20,000 individual genes, we compute one score per pathway using ssGSEA, GSVA, or mean-Z methods.

## Subtype Discovery

With one score per pathway per patient, we have a matrix we can cluster. The framework uses **Gaussian Mixture Models (GMM)** by default because:

- **Soft assignments** — each patient gets a probability of belonging to each subtype, not a hard label. This reflects biological reality where boundaries between subtypes are fuzzy.
- **Automatic model selection** — BIC (Bayesian Information Criterion) selects the best number of clusters without manual tuning.
- **Handles unequal cluster sizes** — real patient subgroups are rarely equal-sized.

K-means, hierarchical, and spectral clustering are also available for comparison.

## Validation Gates: Why They Matter

Finding clusters is easy. Finding clusters that are **real** is hard. The framework includes 5 built-in validation gates that protect against overfitting:

| Gate | What It Tests | How It Works |
|------|--------------|-------------|
| **Label Shuffle** | Are clusters better than random? | Randomly permute patient labels; real clusters should score much higher (ARI < 0.15 for random) |
| **Random Gene Sets** | Do the pathways matter? | Replace real pathways with random gene sets; clusters should disappear |
| **Bootstrap Stability** | Are clusters reproducible? | Resample patients with replacement; clusters should be stable (ARI > 0.8) |
| **Ancestry Independence** | Are clusters confounded by population? | Clusters should not correlate with ancestry PCs |
| **Cross-Modal** | Do clusters replicate across data types? | Subtypes from VCF should match subtypes from expression |

If any gate fails, the subtypes may not be biologically meaningful.

## The 5-Layer Architecture

The framework is modular. Install only what you need:

### Core Layer (always installed)

**What:** GMT-based pathway scoring, GMM clustering, validation gates, simulation, visualization.

**When to use:** Every analysis. This is the foundation — it takes your data (VCF, expression, or scRNA-seq), scores pathways, finds subtypes, validates them, and generates reports.

```bash
pip install pathway-subtyping
```

### Graph Layer (`[graph]`)

**What:** Knowledge graph builder (genes, pathways, GO terms, drugs as a connected network), PPI network propagation, topology-weighted pathway scoring.

**When to use:** When you want network-aware analysis. Pathway scoring treats each gene equally; topology-weighted scoring gives more weight to hub genes (well-connected genes that affect many others). Also needed for Cytoscape network visualizations.

```bash
pip install pathway-subtyping[graph]
```

### QC Layer (`[qc]`)

**What:** 12-feature molecular quality control for manufactured cells — cascade completion, dosage gating, off-target detection, batch heterogeneity, drift monitoring, stress fingerprinting.

**When to use:** If you're engineering cells (CAR T therapy, organoids, iPSC differentiation) and need to verify they're doing what they're supposed to do at the pathway level. Not needed for observational genomics studies.

```bash
pip install pathway-subtyping[qc]
```

### GNN Layer (`[gnn]`) — Experimental

**What:** TransE/RotatE knowledge graph embeddings (learn vector representations of genes and pathways), OntologyAwareGNN (graph neural network with biological attention), gene risk classification.

**When to use:** When you want to use graph neural networks for gene-level prediction (e.g., "is this gene likely to be autism-associated?") or when you want learned gene embeddings instead of raw expression as features for clustering. Always benchmark against the simpler Core layer — GNNs don't always improve results on small biological datasets.

```bash
pip install pathway-subtyping[gnn]
```

### Autism Layer (`[autism]`) — Autism-Only

**What:** 8 curated biological rules (R1-R7) for autism variant interpretation, neuro-symbolic reasoning (combining GNN predictions with domain knowledge), and therapeutic hypothesis generation with drug safety scoring (11 pediatric safety flags).

**When to use:** Only for autism spectrum disorder datasets. The engine refuses to run on non-autism data. This layer sits on TOP of the disease-agnostic pipeline — first you find subtypes with Core, then you interpret them with the autism rules.

```bash
pip install pathway-subtyping[autism]
```

## Typical Workflow

### For a genomics researcher:

1. **Prepare data** — VCF with gene/consequence annotations, or expression matrix
2. **Choose pathways** — use provided GMT files or create your own
3. **Run pipeline** — `psf --config my_config.yaml`
4. **Check validation gates** — if they fail, your subtypes may not be real
5. **Interpret results** — pathway enrichment per subtype, gene contributions

### For a cell manufacturer:

1. **Run Core pipeline** on your batch expression data
2. **Define engineering spec** — intended pathways, excluded pathways, therapeutic windows
3. **Run QC layer** — off-target check, dosage gate, heterogeneity profiling
4. **Get release decision** — RELEASE / HOLD / REJECT with specific reasons

### For an autism researcher:

1. **Run Core pipeline** on ASD expression data
2. **Load biological context** — SFARI genes, pLI scores, CHD8 targets, SynGO
3. **Run autism rules** — which rules fire for each individual?
4. **Generate explanation** — reasoning chain with evidence citations
5. **Generate therapeutic hypotheses** — ranked drug candidates with safety flags

---

## What This Framework Does NOT Do

- **Variant calling** — it takes VCFs as input, it doesn't generate them
- **Pathway discovery** — you supply the pathways; the framework scores them
- **Clinical diagnosis** — all results are research hypotheses, not diagnoses
- **Causal inference** — it finds associations, not causation

---

For statistical details, see [METHODS.md](METHODS.md). For the full API, see [api/index.md](api/index.md).
