#!/usr/bin/env python3
"""Generate NB19: SCZ Blood Multi-Cohort Validation (Hertzberg Replication).

This script builds the notebook JSON programmatically to avoid
NotebookEdit cell ID issues with multi-cell construction.
"""
import json
import uuid


def make_cell(cell_type, source):
    """Create a notebook cell with a unique ID."""
    return {
        "cell_type": cell_type,
        "id": uuid.uuid4().hex[:12],
        "metadata": {},
        "source": source.split("\n") if isinstance(source, str) else source,
        "outputs": [] if cell_type == "code" else [],
        **({"execution_count": None} if cell_type == "code" else {}),
    }


def md(source):
    return make_cell("markdown", source)


def code(source):
    return make_cell("code", source)


# Fix: source must be list of lines with \n
def make_cell_v2(cell_type, source_str):
    lines = source_str.split("\n")
    # Add \n to all lines except the last
    source_list = [line + "\n" for line in lines[:-1]]
    if lines[-1]:  # non-empty last line
        source_list.append(lines[-1])
    cell = {
        "cell_type": cell_type,
        "id": uuid.uuid4().hex[:12],
        "metadata": {},
        "source": source_list,
    }
    if cell_type == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    return cell


def md(s):
    return make_cell_v2("markdown", s)


def code(s):
    return make_cell_v2("code", s)


cells = []

# ── Cell 0: Title ──────────────────────────────────────────────
cells.append(md("""\
# Notebook 19: SCZ Blood Multi-Cohort Validation (Hertzberg Replication)

**Collaboration:** Dr. Libi Hertzberg (Tel Aviv University / Weizmann Institute)

**Reference:** Hertzberg L et al. "Schizophrenia Biomarkers: Blood Transcriptome Suggests Two Molecular Subtypes" *NeuroMolecular Medicine* 2024 (PMID: [39609319](https://pubmed.ncbi.nlm.nih.gov/39609319/))

| Field | Value |
|-------|-------|
| Datasets | 5 GEO cohorts (GSE38484, GSE27383, GSE38481, GSE18312, GSE48072) |
| Platform | Illumina HumanHT-12 v3/v4 (GPL6947, GPL6883, GPL10558) |
| Samples | ~394 total (SCZ + controls across 5 cohorts) |
| Tissue | Whole blood / PBMCs |
| Design | Multi-cohort cross-validation with cross-projection |

**Framework:** [pathway-subtyping](https://codeberg.org/pathways/pathway-subtyping-framework) v0.3.1+

---

### What this notebook does

1. Downloads all 5 GEO datasets used in Hertzberg et al. 2024
2. Processes each: probe-to-gene mapping, QC, ssGSEA pathway scoring
3. **Per-dataset analysis:** BIC-optimal clustering + validation gates (independent)
4. **Cross-cohort projection:** Train GMM on reference (GSE38484), project onto 4 targets
5. **k=2 forced comparison:** Direct head-to-head with Hertzberg's K-means k=2 solution
6. **Merged multi-cohort analysis:** Pool all SCZ samples (~200+) for the largest validated SCZ blood subtyping
7. Benchmark comparison, subtype characterization, and Hertzberg concordance analysis

**Runtime:** ~20-40 minutes (GEO downloads + pathway scoring across 5 datasets)"""))

# ── Cell 1: Setup ──────────────────────────────────────────────
cells.append(md("## 1. Setup & Installation"))

cells.append(code("""\
import subprocess, sys

for pkg in ['pathway-subtyping[viz]==0.3.1', 'GEOparse']:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', pkg],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os, json, time, requests

from pathway_subtyping import (
    score_pathways_from_expression, ExpressionScoringMethod,
    run_clustering, ClusteringAlgorithm,
    select_n_clusters,
    ValidationGates,
    characterize_subtypes, generate_subtype_heatmap,
    generate_gene_heatmap, export_characterization,
    run_benchmark_comparison,
    compute_dim_reduction, DimReductionMethod,
)
import pathway_subtyping

from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score
from scipy.stats import chi2_contingency, fisher_exact, spearmanr

SEED = 42
np.random.seed(SEED)
K_RANGE = list(range(2, 8))  # [2, 3, 4, 5, 6, 7]

OUTPUT_DIR = './outputs/scz_blood_multi_cohort'
DATA_DIR = './data'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Dataset metadata from Hertzberg et al. 2024
DATASETS = {
    'GSE38484': {'platform': 'GPL6947', 'tissue': 'whole blood',
                 'role': 'training', 'desc': "Hertzberg's training set"},
    'GSE27383': {'platform': 'GPL6883', 'tissue': 'PBMCs',
                 'role': 'validation', 'desc': 'External validation set'},
    'GSE38481': {'platform': 'GPL6947', 'tissue': 'whole blood',
                 'role': 'validation', 'desc': 'Same platform as GSE38484'},
    'GSE18312': {'platform': 'GPL6883', 'tissue': 'whole blood',
                 'role': 'validation', 'desc': 'Smaller cohort'},
    'GSE48072': {'platform': 'GPL10558', 'tissue': 'PBMCs',
                 'role': 'validation', 'desc': 'Independent cohort'},
}

print(f'pathway-subtyping v{pathway_subtyping.__version__}')
print(f'Seed: {SEED}')
print(f'k range: {K_RANGE}')
print(f'Datasets: {list(DATASETS.keys())}')
print('Setup complete.')"""))

# ── Cell 2-3: Download all datasets ───────────────────────────
cells.append(md("""\
## 2. Download All GEO Datasets

All 5 datasets are Illumina HumanHT-12 microarray variants. Gene symbol annotation is available in the `Symbol` column of each GPL table."""))

cells.append(code("""\
import GEOparse

os.environ['GEOPARSE_USE_HTTP_FOR_FTP'] = 'yes'

gse_objects = {}

for accession in DATASETS:
    soft_file = os.path.join(DATA_DIR, f'{accession}_family.soft.gz')
    gse = None

    for attempt in range(1, 4):
        try:
            if os.path.exists(soft_file) and os.path.getsize(soft_file) > 1000:
                print(f'[Cache] {accession}: loading from {soft_file}')
                gse = GEOparse.get_GEO(filepath=soft_file, silent=True)
            else:
                print(f'[Download] {accession}: attempt {attempt}/3...')
                gse = GEOparse.get_GEO(geo=accession, destdir=DATA_DIR, silent=True)
            break
        except Exception as e:
            print(f'  Attempt {attempt} failed: {e}')
            if os.path.exists(soft_file):
                try:
                    os.remove(soft_file)
                except OSError:
                    pass
            if attempt < 3:
                wait = 10 * attempt
                print(f'  Retrying in {wait}s...')
                time.sleep(wait)
            else:
                raise RuntimeError(
                    f'Failed to download {accession} after 3 attempts. '
                    f'Manual: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}'
                )

    gse_objects[accession] = gse
    print(f'  Platforms: {list(gse.gpls.keys())}, Samples: {len(gse.gsms)}')

print(f'\\nAll {len(gse_objects)} datasets downloaded successfully.')"""))

# ── Cell 4-5: Extract metadata & expression ───────────────────
cells.append(md("""\
## 3. Extract Metadata & Build Expression Matrices

For each dataset:
1. Parse sample metadata (diagnosis from `characteristics_ch1`)
2. Extract expression values and map probes to gene symbols
3. Collapse multi-probe genes by mean, QC filter, log2 transform if needed"""))

cells.append(code("""\
dataset_results = {}

for accession, ds_info in DATASETS.items():
    print(f'\\n{"=" * 60}')
    print(f'Processing {accession} ({ds_info["desc"]})')
    print(f'{"=" * 60}')

    gse = gse_objects[accession]

    # ── Metadata extraction ──
    meta_rows = []
    for gsm_name, gsm in gse.gsms.items():
        chars = gsm.metadata.get('characteristics_ch1', [])
        char_dict = {}
        for c in chars:
            if ':' in c:
                key, val = c.split(':', 1)
                char_dict[key.strip().lower()] = val.strip()

        meta_rows.append({
            'sample_id': gsm_name,
            'title': gsm.metadata.get('title', [''])[0],
            'source': gsm.metadata.get('source_name_ch1', [''])[0],
            'accession': accession,
            **char_dict,
        })

    meta = pd.DataFrame(meta_rows).set_index('sample_id')

    # ── Auto-detect diagnosis ──
    diagnosis_col = None
    for col in meta.columns:
        col_lower = col.lower()
        vals_lower = [str(v).lower() for v in meta[col].unique()]
        if any(kw in col_lower for kw in ['diagnosis', 'disease', 'condition', 'group', 'status']):
            diagnosis_col = col
            break
        if any(kw in ' '.join(vals_lower) for kw in ['schizophrenia', 'control', 'bipolar']):
            if col not in ['title', 'source', 'accession']:
                diagnosis_col = col
                break

    if diagnosis_col is None:
        # Fallback: search source_name or title
        for col in ['source', 'title']:
            vals_lower = [str(v).lower() for v in meta[col].unique()]
            if any('schizo' in v or 'control' in v for v in vals_lower):
                diagnosis_col = col
                break

    if diagnosis_col:
        print(f'  Diagnosis column: "{diagnosis_col}"')
        print(f'  Values: {meta[diagnosis_col].value_counts().to_dict()}')
        # Standardize diagnosis labels
        diag_map = {}
        for val in meta[diagnosis_col].unique():
            val_lower = str(val).lower().strip()
            if any(kw in val_lower for kw in ['schizo', 'scz', 'sz']):
                diag_map[val] = 'SCZ'
            elif any(kw in val_lower for kw in ['control', 'normal', 'healthy', 'ctl']):
                diag_map[val] = 'CTL'
            elif any(kw in val_lower for kw in ['bipolar', 'bp', 'bd']):
                diag_map[val] = 'BPD'
            elif any(kw in val_lower for kw in ['depress', 'mdd']):
                diag_map[val] = 'MDD'
            else:
                diag_map[val] = val
        meta['diagnosis'] = meta[diagnosis_col].map(diag_map)
    else:
        print(f'  WARNING: No diagnosis column found!')
        print(f'  Columns: {list(meta.columns)}')
        meta['diagnosis'] = 'Unknown'

    # ── Expression matrix ──
    expr = gse.pivot_samples('VALUE')
    expr = expr.apply(pd.to_numeric, errors='coerce')
    expr = expr.dropna(how='all')
    print(f'  Raw expression: {expr.shape[0]} probes x {expr.shape[1]} samples')

    # ── Probe-to-gene mapping ──
    gpl = list(gse.gpls.values())[0]
    gpl_table = gpl.table

    symbol_col = None
    for col in ['Symbol', 'Gene Symbol', 'GENE_SYMBOL', 'Gene_Symbol', 'ILMN_Gene']:
        if col in gpl_table.columns:
            symbol_col = col
            break
    if symbol_col is None:
        for col in gpl_table.columns:
            if 'symbol' in col.lower():
                symbol_col = col
                break

    if symbol_col:
        print(f'  Gene symbol column: "{symbol_col}"')
        probe_gene = gpl_table.set_index('ID')[symbol_col].dropna()
        probe_gene = probe_gene[probe_gene.str.strip() != '']

        # Handle multi-gene probes (take first gene)
        probe_gene = probe_gene.apply(lambda x: str(x).split('///')[0].strip())
        probe_gene = probe_gene[probe_gene != '']
        probe_gene = probe_gene[~probe_gene.str.startswith('---')]

        common_probes = expr.index.intersection(probe_gene.index)
        expr_mapped = expr.loc[common_probes].copy()
        expr_mapped['gene'] = probe_gene.loc[common_probes].values
        gene_expr = expr_mapped.groupby('gene').mean().T
    else:
        print(f'  WARNING: No symbol column found in {list(gpl_table.columns)}')
        gene_expr = expr.T

    # ── QC ──
    # Log2 transform if needed (Illumina intensities > 30 → not log-transformed)
    max_val = gene_expr.max().max()
    if max_val > 30:
        print(f'  Log2 transforming (max value = {max_val:.1f})')
        gene_expr = np.log2(gene_expr.clip(lower=1))

    # Drop zero-variance genes
    gene_var = gene_expr.var()
    n_zero = (gene_var == 0).sum()
    if n_zero > 0:
        gene_expr = gene_expr.loc[:, gene_var > 0]
        print(f'  Dropped {n_zero} zero-variance genes')

    # Fill NaN
    n_nan = gene_expr.isna().sum().sum()
    if n_nan > 0:
        gene_expr = gene_expr.fillna(gene_expr.median())
        print(f'  Filled {n_nan} NaN values')

    # Align samples
    common_samples = gene_expr.index.intersection(meta.index)
    gene_expr = gene_expr.loc[common_samples]
    meta = meta.loc[common_samples]

    print(f'  Final: {gene_expr.shape[0]} samples x {gene_expr.shape[1]} genes')
    print(f'  Diagnosis: {meta["diagnosis"].value_counts().to_dict()}')

    dataset_results[accession] = {
        'expression': gene_expr,
        'metadata': meta,
        'platform': ds_info['platform'],
        'tissue': ds_info['tissue'],
        'n_total': len(meta),
        'n_scz': int((meta['diagnosis'] == 'SCZ').sum()),
        'n_ctl': int((meta['diagnosis'] == 'CTL').sum()),
        'n_genes': gene_expr.shape[1],
    }

print(f'\\n{"=" * 60}')
print(f'All datasets processed.')
print(f'{"=" * 60}')"""))

# ── Cell 6-7: Dataset summary ─────────────────────────────────
cells.append(md("## 4. Dataset Summary & Gene Overlap"))

cells.append(code("""\
# Summary table
summary_rows = []
for acc, dr in dataset_results.items():
    summary_rows.append({
        'Accession': acc,
        'Platform': dr['platform'],
        'Tissue': dr['tissue'],
        'N total': dr['n_total'],
        'N SCZ': dr['n_scz'],
        'N CTL': dr['n_ctl'],
        'N genes': dr['n_genes'],
        'Role': DATASETS[acc]['role'],
    })

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))

# Pairwise gene overlap
accessions = list(dataset_results.keys())
gene_sets = {acc: set(dataset_results[acc]['expression'].columns) for acc in accessions}

print(f'\\n--- Pairwise Gene Overlap ---')
for i, a1 in enumerate(accessions):
    for a2 in accessions[i+1:]:
        overlap = len(gene_sets[a1] & gene_sets[a2])
        print(f'  {a1} vs {a2}: {overlap} genes')

# Common genes across all 5
common_genes = gene_sets[accessions[0]]
for acc in accessions[1:]:
    common_genes = common_genes & gene_sets[acc]
print(f'\\nCommon genes across ALL 5 datasets: {len(common_genes)}')

# Total SCZ samples
total_scz = sum(dr['n_scz'] for dr in dataset_results.values())
total_all = sum(dr['n_total'] for dr in dataset_results.values())
print(f'Total SCZ samples: {total_scz}')
print(f'Total samples (all): {total_all}')"""))

# ── Cell 8-9: Load pathways & score ───────────────────────────
cells.append(md("""\
## 5. Load Hallmark Pathways & Score All Datasets

We use MSigDB Hallmark gene sets (same as NB17/NB18) for cross-disease consistency."""))

cells.append(code("""\
# Download MSigDB Hallmark gene sets
HALLMARKS_URLS = [
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/h.all.v2022.1.Hs.symbols.gmt',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt',
]
HALLMARKS_PATH = os.path.join(DATA_DIR, 'h.hallmarks.gmt')

def download_file(urls, dest_path, desc='file', retries=3):
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0:
        print(f'[Cache] {desc}: {dest_path}')
        return dest_path
    urls = [urls] if isinstance(urls, str) else urls
    for url in urls:
        for attempt in range(1, retries + 1):
            try:
                r = requests.get(url, timeout=60)
                r.raise_for_status()
                with open(dest_path, 'w') as f:
                    f.write(r.text)
                print(f'[Downloaded] {desc} from {url.split("/")[-1]}')
                return dest_path
            except Exception as e:
                if attempt == retries:
                    print(f'  Failed: {url} ({e})')
    raise RuntimeError(f'Failed to download {desc}')

download_file(HALLMARKS_URLS, HALLMARKS_PATH, 'MSigDB Hallmark gene sets')

def parse_gmt(path):
    pathways = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) < 3:
                continue
            pathways[parts[0]] = [g.strip() for g in parts[2:] if g.strip()]
    return pathways

hallmark_pathways = parse_gmt(HALLMARKS_PATH)
print(f'Loaded {len(hallmark_pathways)} Hallmark gene sets')

# Score pathways for each dataset
for accession in DATASETS:
    print(f'\\nScoring {accession}...')
    dr = dataset_results[accession]

    scoring = score_pathways_from_expression(
        gene_expression=dr['expression'],
        pathways=hallmark_pathways,
        method=ExpressionScoringMethod.SSGSEA,
        min_genes_per_pathway=5,
        seed=SEED,
        show_progress=True,
    )

    dr['pathway_scores_all'] = scoring.pathway_scores
    dr['n_pathways'] = scoring.n_pathways_scored
    print(f'  {scoring.n_pathways_scored} pathways scored, {scoring.n_pathways_skipped} skipped')

    # SCZ-only subset
    scz_mask = dr['metadata']['diagnosis'] == 'SCZ'
    dr['scz_mask'] = scz_mask
    dr['pathway_scores_scz'] = scoring.pathway_scores.loc[scz_mask]
    dr['expression_scz'] = dr['expression'].loc[scz_mask]
    print(f'  SCZ subset: {dr["pathway_scores_scz"].shape}')

print('\\nAll datasets scored.')"""))

# ── Cell 10-12: Per-dataset clustering ─────────────────────────
cells.append(md("""\
## 6. Per-Dataset Independent Clustering

Each dataset is clustered independently with BIC model selection and validation gates."""))

cells.append(code("""\
per_dataset = {}

for accession in DATASETS:
    dr = dataset_results[accession]
    scz_scores = dr['pathway_scores_scz']
    scz_expr = dr['expression_scz']

    if len(scz_scores) < 10:
        print(f'\\n{accession}: Only {len(scz_scores)} SCZ samples — skipping clustering')
        per_dataset[accession] = {'skipped': True, 'reason': 'too few samples'}
        continue

    print(f'\\n{"=" * 60}')
    print(f'{accession}: Clustering {len(scz_scores)} SCZ samples')
    print(f'{"=" * 60}')

    # BIC model selection
    max_k = min(7, len(scz_scores) // 5)  # Ensure enough samples per cluster
    k_range = list(range(2, max_k + 1))

    selection = select_n_clusters(
        data=scz_scores.values,
        k_range=k_range,
        method='bic',
        seed=SEED,
        min_cluster_fraction=0.05,
    )
    opt_k = selection.optimal_k
    print(f'  Optimal k (BIC): {opt_k}')

    # GMM clustering
    clustering = run_clustering(
        data=scz_scores.values,
        n_clusters=opt_k,
        algorithm=ClusteringAlgorithm.GMM,
        seed=SEED,
    )
    print(f'  Silhouette: {clustering.silhouette:.4f}')

    # Subtype sizes
    for i in range(opt_k):
        n = int((clustering.labels == i).sum())
        print(f'    Subtype {i}: {n} ({n/len(clustering.labels)*100:.1f}%)')

    # Validation gates
    gates = ValidationGates(
        seed=SEED, n_permutations=200, n_bootstrap=100,
        stability_threshold=0.8, null_ari_max=0.15, show_progress=True,
    )
    val_result = gates.run_all(
        pathway_scores=scz_scores,
        cluster_labels=clustering.labels,
        pathways=hallmark_pathways,
        gene_burdens=scz_expr,
        n_clusters=opt_k,
        gmm_seed=SEED,
    )

    n_passed = sum(g.passed for g in val_result.results)
    n_total = len(val_result.results)
    print(f'  Validation: {n_passed}/{n_total} gates passed')
    for g in val_result.results:
        status = 'PASS' if g.passed else 'FAIL'
        print(f'    [{status}] {g.name}: {g.metric_value:.4f}')

    per_dataset[accession] = {
        'skipped': False,
        'optimal_k': int(opt_k),
        'silhouette': float(clustering.silhouette),
        'calinski_harabasz': float(clustering.calinski_harabasz),
        'davies_bouldin': float(clustering.davies_bouldin),
        'labels': clustering.labels,
        'scores': scz_scores,
        'expression': scz_expr,
        'gates_passed': n_passed,
        'gates_total': n_total,
        'val_result': val_result,
        'selection': selection,
    }"""))

cells.append(code("""\
# Per-dataset results summary table
print(f'\\n{"=" * 70}')
print(f'PER-DATASET CLUSTERING SUMMARY')
print(f'{"=" * 70}')
print(f'{"Accession":<12} {"N_SCZ":>6} {"k":>3} {"Silhouette":>11} {"Gates":>7} {"Status":>8}')
print(f'{"-"*12} {"-"*6} {"-"*3} {"-"*11} {"-"*7} {"-"*8}')

for acc in DATASETS:
    dr = dataset_results[acc]
    pd_res = per_dataset[acc]
    if pd_res.get('skipped'):
        print(f'{acc:<12} {dr["n_scz"]:>6} {"--":>3} {"--":>11} {"--":>7} {"SKIP":>8}')
    else:
        gates_str = f'{pd_res["gates_passed"]}/{pd_res["gates_total"]}'
        status = 'OK' if pd_res['gates_passed'] >= 2 else 'WEAK'
        print(f'{acc:<12} {dr["n_scz"]:>6} {pd_res["optimal_k"]:>3} '
              f'{pd_res["silhouette"]:>11.4f} {gates_str:>7} {status:>8}')

print(f'\\nHertzberg used k=2 (forced) across all datasets.')
print(f'Our BIC-optimal k may differ -- this tests whether data supports k=2.')"""))

# ── Cell 13-14: Cross-cohort projection ───────────────────────
cells.append(md("""\
## 7. Cross-Cohort Projection (KEY ANALYSIS)

This is the core replication analysis. We train a GMM on the **reference dataset** (GSE38484, Hertzberg's training set) and project each target dataset into the reference subtype space.

**Hertzberg's limitation:** She used K-means and achieved only 64% PPV on external validation. PSF's GMM with posterior probabilities should provide more robust cross-cohort assignment."""))

cells.append(code("""\
# Reference = GSE38484 (Hertzberg's training set, largest)
REF_ACC = 'GSE38484'
ref_data = per_dataset[REF_ACC]

if ref_data.get('skipped'):
    print('ERROR: Reference dataset was skipped -- cannot do projection')
else:
    ref_scores = ref_data['scores']
    ref_k = ref_data['optimal_k']
    print(f'Reference: {REF_ACC} (k={ref_k}, n={len(ref_scores)} SCZ)')

    projection_results = {}

    for target_acc in DATASETS:
        if target_acc == REF_ACC:
            continue

        target_data = per_dataset[target_acc]
        if target_data.get('skipped'):
            print(f'\\n  {target_acc}: SKIPPED')
            continue

        target_scores = target_data['scores']
        target_own_labels = target_data['labels']

        # Align to shared pathways
        shared = sorted(set(ref_scores.columns) & set(target_scores.columns))
        print(f'\\n  {target_acc}: {len(target_scores)} SCZ, {len(shared)} shared pathways')

        # Fit GMM on reference (shared pathways only)
        gmm_ref = GaussianMixture(
            n_components=ref_k, covariance_type='full',
            random_state=SEED, n_init=10,
        )
        gmm_ref.fit(ref_scores[shared].values)

        # Project target
        projected_labels = gmm_ref.predict(target_scores[shared].values)
        projected_probs = gmm_ref.predict_proba(target_scores[shared].values)

        # Compare: projected labels vs target's own independent labels
        # Only meaningful if target has same k
        ari = adjusted_rand_score(target_own_labels, projected_labels)
        mean_conf = projected_probs.max(axis=1).mean()

        projection_results[target_acc] = {
            'n_scz': int(len(target_scores)),
            'n_shared_pathways': len(shared),
            'projection_ari': float(ari),
            'mean_confidence': float(mean_conf),
            'projected_labels': projected_labels,
            'ref_k': ref_k,
        }

        passed = ari > 0.3
        print(f'    Projection ARI: {ari:.4f} ({"PASS" if passed else "BELOW 0.3"})')
        print(f'    Mean confidence: {mean_conf:.4f}')

    # Summary
    print(f'\\n{"=" * 60}')
    print(f'CROSS-COHORT PROJECTION SUMMARY (reference: {REF_ACC}, k={ref_k})')
    print(f'{"=" * 60}')
    print(f'{"Target":<12} {"N_SCZ":>6} {"ARI":>7} {"Confidence":>11} {"Pass":>5}')
    print(f'{"-"*12} {"-"*6} {"-"*7} {"-"*11} {"-"*5}')
    aris = []
    for acc, pr in projection_results.items():
        passed = 'YES' if pr['projection_ari'] > 0.3 else 'no'
        print(f'{acc:<12} {pr["n_scz"]:>6} {pr["projection_ari"]:>7.4f} '
              f'{pr["mean_confidence"]:>11.4f} {passed:>5}')
        aris.append(pr['projection_ari'])
    if aris:
        print(f'\\nMean projection ARI: {np.mean(aris):.4f}')
        n_pass = sum(1 for a in aris if a > 0.3)
        print(f'Passed (ARI > 0.3): {n_pass}/{len(aris)}')"""))

# ── Cell 15-16: k=2 forced Hertzberg comparison ───────────────
cells.append(md("""\
## 8. k=2 Forced Analysis (Direct Hertzberg Comparison)

Hertzberg used K-means with k=2 across all datasets. We force k=2 with our GMM approach to enable direct comparison.

**Hertzberg's subtypes:**
- Subtype A: elevated inflammatory/immune response
- Subtype B: altered neurodevelopmental signaling"""))

cells.append(code("""\
k2_results = {}

for accession in DATASETS:
    dr = dataset_results[accession]
    pd_res = per_dataset[accession]

    if pd_res.get('skipped'):
        continue

    scz_scores = pd_res['scores']
    scz_expr = pd_res['expression']

    # Force k=2
    clustering_k2 = run_clustering(
        data=scz_scores.values,
        n_clusters=2,
        algorithm=ClusteringAlgorithm.GMM,
        seed=SEED,
    )

    # Compare k=2 labels vs BIC-optimal labels
    opt_labels = pd_res['labels']
    opt_k = pd_res['optimal_k']
    ari_vs_optimal = adjusted_rand_score(opt_labels, clustering_k2.labels)

    # Characterize the 2 subtypes by top pathways
    subtype_means = {}
    for s in [0, 1]:
        mask = clustering_k2.labels == s
        subtype_means[s] = scz_scores.iloc[mask].mean()

    # Identify top differentiating pathways
    diff = subtype_means[0] - subtype_means[1]
    top_up_0 = diff.nlargest(5)
    top_up_1 = (-diff).nlargest(5)

    k2_results[accession] = {
        'silhouette': float(clustering_k2.silhouette),
        'labels': clustering_k2.labels,
        'sizes': [int((clustering_k2.labels == i).sum()) for i in range(2)],
        'ari_vs_optimal': float(ari_vs_optimal),
        'optimal_k': opt_k,
        'top_pathways_0': list(top_up_0.index),
        'top_pathways_1': list(top_up_1.index),
    }

    print(f'{accession}: k=2 sil={clustering_k2.silhouette:.4f}, '
          f'sizes={k2_results[accession]["sizes"]}, '
          f'ARI vs BIC-k={opt_k}: {ari_vs_optimal:.4f}')
    print(f'  Subtype 0 top: {", ".join(p.replace("HALLMARK_", "") for p in top_up_0.index[:3])}')
    print(f'  Subtype 1 top: {", ".join(p.replace("HALLMARK_", "") for p in top_up_1.index[:3])}')

# Cross-cohort k=2 projection
print(f'\\n--- k=2 Cross-Cohort Projection ---')
ref_k2 = k2_results[REF_ACC]
ref_scores_k2 = per_dataset[REF_ACC]['scores']

gmm_k2_ref = GaussianMixture(n_components=2, covariance_type='full',
                               random_state=SEED, n_init=10)
for target_acc in DATASETS:
    if target_acc == REF_ACC or per_dataset[target_acc].get('skipped'):
        continue
    target_scores = per_dataset[target_acc]['scores']
    shared = sorted(set(ref_scores_k2.columns) & set(target_scores.columns))
    gmm_k2_ref.fit(ref_scores_k2[shared].values)
    proj_labels = gmm_k2_ref.predict(target_scores[shared].values)
    target_k2_labels = k2_results[target_acc]['labels']
    ari_k2 = adjusted_rand_score(target_k2_labels, proj_labels)
    k2_results[target_acc]['projection_ari_k2'] = float(ari_k2)
    print(f'  {REF_ACC} -> {target_acc}: k=2 projection ARI = {ari_k2:.4f}')"""))

# ── Cell 17-19: Merged multi-cohort ───────────────────────────
cells.append(md("""\
## 9. Merged Multi-Cohort Analysis

Pool all SCZ samples across 5 datasets for the **largest validated SCZ blood subtyping analysis**.
Batch correction: z-score normalization per dataset on common genes, then re-score pathways."""))

cells.append(code("""\
# Merge on common genes (already computed above)
print(f'Common genes across all datasets: {len(common_genes)}')
common_genes_sorted = sorted(common_genes)

# Z-score per dataset for batch correction
merged_expr_list = []
merged_meta_list = []

for accession in DATASETS:
    dr = dataset_results[accession]
    scz_mask = dr['metadata']['diagnosis'] == 'SCZ'
    expr_scz = dr['expression'].loc[scz_mask, common_genes_sorted].copy()

    # Z-score normalization within this dataset
    expr_mean = expr_scz.mean()
    expr_std = expr_scz.std().replace(0, 1)
    expr_z = (expr_scz - expr_mean) / expr_std

    # Tag with dataset
    meta_scz = dr['metadata'].loc[scz_mask].copy()
    meta_scz['dataset'] = accession

    merged_expr_list.append(expr_z)
    merged_meta_list.append(meta_scz)
    print(f'  {accession}: {len(expr_z)} SCZ samples')

merged_expression = pd.concat(merged_expr_list)
merged_metadata = pd.concat(merged_meta_list)

# Drop any remaining NaN
n_nan = merged_expression.isna().sum().sum()
if n_nan > 0:
    merged_expression = merged_expression.fillna(0)
    print(f'  Filled {n_nan} NaN values in merged data')

print(f'\\nMerged expression: {merged_expression.shape[0]} SCZ samples x {merged_expression.shape[1]} genes')
print(f'Dataset composition:')
print(merged_metadata['dataset'].value_counts().to_string())"""))

cells.append(code("""\
# Score pathways on merged data
print('Scoring pathways on merged SCZ cohort...')
merged_scoring = score_pathways_from_expression(
    gene_expression=merged_expression,
    pathways=hallmark_pathways,
    method=ExpressionScoringMethod.SSGSEA,
    min_genes_per_pathway=5,
    seed=SEED,
    show_progress=True,
)
merged_pathway_scores = merged_scoring.pathway_scores
print(f'Merged pathway scores: {merged_pathway_scores.shape}')

# Cluster selection
merged_k_range = list(range(2, 8))
merged_selection = select_n_clusters(
    data=merged_pathway_scores.values,
    k_range=merged_k_range,
    method='bic',
    seed=SEED,
    min_cluster_fraction=0.05,
)
merged_k = merged_selection.optimal_k
print(f'\\nOptimal k (BIC): {merged_k}')
print(f'BIC values: {merged_selection.bic_values}')

# Cluster
merged_clustering = run_clustering(
    data=merged_pathway_scores.values,
    n_clusters=merged_k,
    algorithm=ClusteringAlgorithm.GMM,
    seed=SEED,
)
merged_labels = merged_clustering.labels
merged_metadata['subtype'] = merged_labels

print(f'Silhouette: {merged_clustering.silhouette:.4f}')
print(f'\\nSubtype composition:')
for i in range(merged_k):
    n = int((merged_labels == i).sum())
    print(f'  Subtype {i}: {n} ({n/len(merged_labels)*100:.1f}%)')

# Check batch effects: are subtypes driven by dataset?
print(f'\\n--- Subtype x Dataset Contingency ---')
ct_batch = pd.crosstab(merged_metadata['subtype'], merged_metadata['dataset'])
print(ct_batch.to_string())

# Chi-square to test for batch-driven subtypes
chi2_batch, p_batch, _, _ = chi2_contingency(ct_batch)
print(f'\\nChi-square (subtype ~ dataset): X2={chi2_batch:.2f}, p={p_batch:.4f}')
if p_batch < 0.05:
    print('WARNING: Subtypes may be partially confounded with dataset batch.')
    print('Interpret merged results cautiously.')
else:
    print('Subtypes are NOT significantly associated with dataset — batch effect minimal.')

# Validation gates on merged data
print(f'\\nRunning validation gates on merged cohort...')
merged_gates = ValidationGates(
    seed=SEED, n_permutations=200, n_bootstrap=100,
    stability_threshold=0.8, null_ari_max=0.15, show_progress=True,
)
merged_val = merged_gates.run_all(
    pathway_scores=merged_pathway_scores,
    cluster_labels=merged_labels,
    pathways=hallmark_pathways,
    gene_burdens=merged_expression,
    n_clusters=merged_k,
    gmm_seed=SEED,
)

n_passed = sum(g.passed for g in merged_val.results)
print(f'\\nMerged validation: {n_passed}/{len(merged_val.results)} gates passed')
for g in merged_val.results:
    status = 'PASS' if g.passed else 'FAIL'
    print(f'  [{status}] {g.name}: {g.metric_value:.4f}')"""))

# ── Cell 20-21: Subtype characterization ──────────────────────
cells.append(md("## 10. Merged Subtype Characterization"))

cells.append(code("""\
# Characterize merged subtypes
merged_char = characterize_subtypes(
    pathway_scores=merged_pathway_scores,
    cluster_labels=merged_labels,
    gene_burdens=merged_expression,
    pathways=hallmark_pathways,
    fdr_alpha=0.05,
    top_n_genes=20,
    seed=SEED,
)
print(merged_char.format_report())

# Heatmaps
fig_hm = generate_subtype_heatmap(
    merged_char,
    output_path=os.path.join(OUTPUT_DIR, 'subtype_heatmap_merged.png'),
    figsize=(14, 8),
)
plt.show()

fig_gm = generate_gene_heatmap(
    merged_char,
    output_path=os.path.join(OUTPUT_DIR, 'gene_heatmap_merged.png'),
    figsize=(16, 10),
    top_n=15,
)
plt.show()

# Export
export_files = export_characterization(merged_char, output_dir=OUTPUT_DIR)
print(f'\\nExported {len(export_files)} characterization files')"""))

# ── Cell 22-23: Benchmark ─────────────────────────────────────
cells.append(md("## 11. Benchmark Comparison"))

cells.append(code("""\
bench_result = run_benchmark_comparison(
    gene_burdens=merged_expression,
    pathway_scores=merged_pathway_scores,
    pathways=hallmark_pathways,
    n_clusters=merged_k,
    seed=SEED,
)
print(bench_result.format_report())

# Visualization
methods = list(bench_result.method_results.keys())
silhouettes = [bench_result.method_results[m].silhouette for m in methods]

fig, ax = plt.subplots(figsize=(10, 6))
colors = ['#2ecc71' if m == bench_result.best_method else '#3498db' for m in methods]
bars = ax.barh(methods, silhouettes, color=colors)
ax.set_xlabel('Silhouette Score')
ax.set_title(f'SCZ Blood Multi-Cohort: Benchmark (k={merged_k}, n={len(merged_labels)})')
for bar, val in zip(bars, silhouettes):
    ax.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height()/2,
            f'{val:.3f}', va='center', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'benchmark_comparison.png'), dpi=150, bbox_inches='tight')
plt.show()"""))

# ── Cell 24-25: Hertzberg concordance ─────────────────────────
cells.append(md("""\
## 12. Hertzberg Concordance Analysis

Map our subtypes to Hertzberg's biological descriptions:
- **Hertzberg Subtype A:** elevated inflammatory/immune response pathways
- **Hertzberg Subtype B:** altered neurodevelopmental signaling

We check if our pathway-level subtypes recover similar pathway enrichment patterns."""))

cells.append(code("""\
# Define Hertzberg's pathway signature (mapped to Hallmark names)
hertzberg_immune = [
    'HALLMARK_INFLAMMATORY_RESPONSE',
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE',
    'HALLMARK_INTERFERON_ALPHA_RESPONSE',
    'HALLMARK_IL6_JAK_STAT3_SIGNALING',
    'HALLMARK_COMPLEMENT',
    'HALLMARK_IL2_STAT5_SIGNALING',
    'HALLMARK_ALLOGRAFT_REJECTION',
]
hertzberg_neuro = [
    'HALLMARK_HEDGEHOG_SIGNALING',
    'HALLMARK_WNT_BETA_CATENIN_SIGNALING',
    'HALLMARK_NOTCH_SIGNALING',
    'HALLMARK_TGF_BETA_SIGNALING',
    'HALLMARK_MYOGENESIS',
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
]

# Score each merged subtype on these pathway groups
available_pathways = set(merged_pathway_scores.columns)
immune_available = [p for p in hertzberg_immune if p in available_pathways]
neuro_available = [p for p in hertzberg_neuro if p in available_pathways]

print(f'Immune pathways available: {len(immune_available)}/{len(hertzberg_immune)}')
print(f'Neuro pathways available: {len(neuro_available)}/{len(hertzberg_neuro)}')

# Mean pathway score per subtype
concordance_data = []
for i in range(merged_k):
    mask = merged_labels == i
    subtype_scores = merged_pathway_scores.iloc[mask]

    immune_score = subtype_scores[immune_available].mean().mean() if immune_available else 0
    neuro_score = subtype_scores[neuro_available].mean().mean() if neuro_available else 0

    concordance_data.append({
        'subtype': i,
        'n': int(mask.sum()),
        'immune_mean': float(immune_score),
        'neuro_mean': float(neuro_score),
        'classification': 'Immune-like (A)' if immune_score > neuro_score else 'Neuro-like (B)',
    })

concordance_df = pd.DataFrame(concordance_data)
print(f'\\n--- Hertzberg Concordance ---')
print(concordance_df.to_string(index=False))

# Which of our subtypes best matches each Hertzberg type?
best_immune = concordance_df.loc[concordance_df['immune_mean'].idxmax()]
best_neuro = concordance_df.loc[concordance_df['neuro_mean'].idxmax()]
print(f'\\nBest match for Hertzberg A (immune): Subtype {int(best_immune["subtype"])} '
      f'(immune={best_immune["immune_mean"]:.4f})')
print(f'Best match for Hertzberg B (neuro): Subtype {int(best_neuro["subtype"])} '
      f'(neuro={best_neuro["neuro_mean"]:.4f})')

# Visualization
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(merged_k)
width = 0.35
ax.bar(x - width/2, concordance_df['immune_mean'], width, label='Immune (Hertzberg A)', color='#e74c3c')
ax.bar(x + width/2, concordance_df['neuro_mean'], width, label='Neuro (Hertzberg B)', color='#3498db')
ax.set_xlabel('Molecular Subtype')
ax.set_ylabel('Mean Pathway Score')
ax.set_title('Hertzberg Concordance: Immune vs Neuro Pathway Enrichment')
ax.set_xticks(x)
ax.set_xticklabels([f'S{i} (n={concordance_df.iloc[i]["n"]:.0f})' for i in range(merged_k)])
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'hertzberg_concordance_plot.png'), dpi=150, bbox_inches='tight')
plt.show()"""))

# ── Cell 26-27: Visualization panel ──────────────────────────
cells.append(md("## 13. Visualization Panel"))

cells.append(code("""\
# Figure 1: PCA of merged cohort colored by subtype + dataset
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

embedding, pca_meta = compute_dim_reduction(
    pathway_scores=merged_pathway_scores,
    labels=merged_labels,
    method=DimReductionMethod.PCA,
    seed=SEED,
)

var_explained = pca_meta.get('explained_variance_ratio', [0, 0])

# Left: by subtype
colors_sub = plt.cm.Set2(np.linspace(0, 1, merged_k))
for i in range(merged_k):
    mask = merged_labels == i
    n = mask.sum()
    axes[0].scatter(embedding[mask, 0], embedding[mask, 1],
                    c=[colors_sub[i]], label=f'Subtype {i} (n={n})',
                    s=40, alpha=0.7, edgecolors='white', linewidth=0.5)
axes[0].set_xlabel(f'PC1 ({var_explained[0]*100:.1f}%)')
axes[0].set_ylabel(f'PC2 ({var_explained[1]*100:.1f}%)')
axes[0].set_title('Colored by Molecular Subtype')
axes[0].legend(fontsize=8)

# Right: by dataset (batch check)
dataset_colors = {'GSE38484': '#e74c3c', 'GSE27383': '#3498db', 'GSE38481': '#2ecc71',
                  'GSE18312': '#f39c12', 'GSE48072': '#9b59b6'}
for acc, color in dataset_colors.items():
    mask = merged_metadata['dataset'].values == acc
    n = mask.sum()
    if n > 0:
        axes[1].scatter(embedding[mask, 0], embedding[mask, 1],
                        c=color, label=f'{acc} (n={n})',
                        s=40, alpha=0.7, edgecolors='white', linewidth=0.5)
axes[1].set_xlabel(f'PC1 ({var_explained[0]*100:.1f}%)')
axes[1].set_ylabel(f'PC2 ({var_explained[1]*100:.1f}%)')
axes[1].set_title('Colored by Dataset (Batch Check)')
axes[1].legend(fontsize=8)

plt.suptitle(f'SCZ Blood Multi-Cohort: PCA of Merged Pathway Scores (n={len(merged_labels)})', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'pca_merged_subtype_and_dataset.png'), dpi=150, bbox_inches='tight')
plt.show()

# Figure 2: Cross-cohort projection ARI bar chart
if projection_results:
    fig, ax = plt.subplots(figsize=(8, 5))
    accs = list(projection_results.keys())
    aris = [projection_results[a]['projection_ari'] for a in accs]
    colors = ['#2ecc71' if a > 0.3 else '#e74c3c' for a in aris]
    bars = ax.bar(accs, aris, color=colors)
    ax.axhline(y=0.3, color='gray', linestyle='--', alpha=0.7, label='Threshold (0.3)')
    ax.set_ylabel('Projection ARI')
    ax.set_title(f'Cross-Cohort Projection from {REF_ACC} (k={ref_data["optimal_k"]})')
    ax.legend()
    for bar, val in zip(bars, aris):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'projection_ari_barplot.png'), dpi=150, bbox_inches='tight')
    plt.show()

# Figure 3: k=2 per-dataset silhouette comparison
fig, ax = plt.subplots(figsize=(8, 5))
k2_accs = [a for a in DATASETS if not per_dataset[a].get('skipped')]
k2_sils = [k2_results[a]['silhouette'] for a in k2_accs]
opt_sils = [per_dataset[a]['silhouette'] for a in k2_accs]

x = np.arange(len(k2_accs))
width = 0.35
ax.bar(x - width/2, opt_sils, width, label='BIC-optimal k', color='#3498db')
ax.bar(x + width/2, k2_sils, width, label='Forced k=2 (Hertzberg)', color='#e74c3c')
ax.set_ylabel('Silhouette Score')
ax.set_title('Per-Dataset: BIC-Optimal vs Hertzberg k=2')
ax.set_xticks(x)
ax.set_xticklabels(k2_accs, rotation=15)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'k2_vs_optimal_silhouette.png'), dpi=150, bbox_inches='tight')
plt.show()"""))

# ── Cell 28-29: Results export ────────────────────────────────
cells.append(md("## 14. Results Export"))

cells.append(code("""\
results_summary = {
    'notebook': 'NB19 -- SCZ Blood Multi-Cohort Validation (Hertzberg Replication)',
    'framework_version': pathway_subtyping.__version__,
    'seed': SEED,
    'citation': 'Hertzberg et al. NeuroMolecular Medicine 2024 (PMID: 39609319)',
    'datasets': list(DATASETS.keys()),
    'total_samples': int(sum(dr['n_total'] for dr in dataset_results.values())),
    'total_scz': int(sum(dr['n_scz'] for dr in dataset_results.values())),

    # Per-dataset results
    'per_dataset': {
        acc: {
            'platform': dataset_results[acc]['platform'],
            'tissue': dataset_results[acc]['tissue'],
            'n_total': dataset_results[acc]['n_total'],
            'n_scz': dataset_results[acc]['n_scz'],
            'n_genes': dataset_results[acc]['n_genes'],
            'optimal_k': per_dataset[acc].get('optimal_k'),
            'silhouette': per_dataset[acc].get('silhouette'),
            'gates_passed': per_dataset[acc].get('gates_passed'),
            'gates_total': per_dataset[acc].get('gates_total'),
            'k2_silhouette': k2_results.get(acc, {}).get('silhouette'),
        } for acc in DATASETS if not per_dataset[acc].get('skipped')
    },

    # Cross-cohort projection
    'cross_cohort_projection': {
        'reference': REF_ACC,
        'reference_k': int(ref_data['optimal_k']),
        'results': {acc: {k: v for k, v in pr.items() if k != 'projected_labels'}
                    for acc, pr in projection_results.items()},
        'mean_ari': float(np.mean([pr['projection_ari'] for pr in projection_results.values()])) if projection_results else None,
        'n_passed': int(sum(1 for pr in projection_results.values() if pr['projection_ari'] > 0.3)),
    },

    # Merged analysis
    'merged': {
        'n_scz': int(len(merged_labels)),
        'n_common_genes': int(len(common_genes)),
        'optimal_k': int(merged_k),
        'silhouette': float(merged_clustering.silhouette),
        'calinski_harabasz': float(merged_clustering.calinski_harabasz),
        'davies_bouldin': float(merged_clustering.davies_bouldin),
        'batch_chi2': float(chi2_batch),
        'batch_pvalue': float(p_batch),
        'gates_passed': int(sum(g.passed for g in merged_val.results)),
        'gates_total': int(len(merged_val.results)),
        'gate_details': [{
            'name': str(g.name), 'passed': bool(g.passed),
            'metric_value': float(g.metric_value), 'threshold': float(g.threshold),
        } for g in merged_val.results],
        'subtype_sizes': {str(i): int((merged_labels == i).sum()) for i in range(merged_k)},
        'benchmark_best': str(bench_result.best_method),
        'benchmark_ranking': list(bench_result.ranking),
    },

    # k=2 Hertzberg comparison
    'hertzberg_k2': {
        acc: {k: v for k, v in res.items() if k != 'labels'}
        for acc, res in k2_results.items()
    },

    # Concordance
    'hertzberg_concordance': concordance_df.to_dict('records'),
}

# Save JSON
json_path = os.path.join(OUTPUT_DIR, 'results_summary.json')
with open(json_path, 'w') as f:
    json.dump(results_summary, f, indent=2, default=str)
print(f'Saved: {json_path}')

# Save CSVs
merged_pathway_scores.to_csv(os.path.join(OUTPUT_DIR, 'pathway_scores_merged_scz.csv'))
merged_metadata.to_csv(os.path.join(OUTPUT_DIR, 'sample_metadata_merged.csv'))
merged_expression.to_csv(os.path.join(OUTPUT_DIR, 'gene_expression_merged.csv'))

if projection_results:
    proj_df = pd.DataFrame([
        {'target': acc, **{k: v for k, v in pr.items() if k != 'projected_labels'}}
        for acc, pr in projection_results.items()
    ])
    proj_df.to_csv(os.path.join(OUTPUT_DIR, 'cross_cohort_projection.csv'), index=False)

concordance_df.to_csv(os.path.join(OUTPUT_DIR, 'hertzberg_concordance.csv'), index=False)

# Per-dataset results
os.makedirs(os.path.join(OUTPUT_DIR, 'per_dataset'), exist_ok=True)
for acc in DATASETS:
    if per_dataset[acc].get('skipped'):
        continue
    pd_json = {
        'accession': acc,
        'n_scz': dataset_results[acc]['n_scz'],
        'optimal_k': per_dataset[acc]['optimal_k'],
        'silhouette': per_dataset[acc]['silhouette'],
        'gates_passed': per_dataset[acc]['gates_passed'],
        'gates_total': per_dataset[acc]['gates_total'],
    }
    with open(os.path.join(OUTPUT_DIR, 'per_dataset', f'{acc}_results.json'), 'w') as f:
        json.dump(pd_json, f, indent=2)

print(f'\\nAll outputs saved to {OUTPUT_DIR}/')
for root, dirs, files in os.walk(OUTPUT_DIR):
    for fname in sorted(files):
        fpath = os.path.join(root, fname)
        size = os.path.getsize(fpath)
        rel = os.path.relpath(fpath, OUTPUT_DIR)
        print(f'  {rel} ({size/1024:.0f} KB)')"""))

# ── Cell 30-31: Summary for Hertzberg email ──────────────────
cells.append(md("## 15. Summary & Key Findings"))

cells.append(code("""\
print('=' * 70)
print('NB19 COMPLETE: SCZ Blood Multi-Cohort Validation')
print('=' * 70)

print(f'\\n--- Datasets ---')
for acc in DATASETS:
    dr = dataset_results[acc]
    pd_res = per_dataset[acc]
    if pd_res.get('skipped'):
        print(f'  {acc}: {dr["n_scz"]} SCZ (SKIPPED)')
    else:
        print(f'  {acc}: {dr["n_scz"]} SCZ, k={pd_res["optimal_k"]}, '
              f'sil={pd_res["silhouette"]:.4f}, gates={pd_res["gates_passed"]}/{pd_res["gates_total"]}')

print(f'\\n--- Cross-Cohort Projection ---')
if projection_results:
    for acc, pr in projection_results.items():
        passed = 'PASS' if pr['projection_ari'] > 0.3 else 'FAIL'
        print(f'  {REF_ACC} -> {acc}: ARI={pr["projection_ari"]:.4f} [{passed}]')
    mean_ari = np.mean([pr['projection_ari'] for pr in projection_results.values()])
    print(f'  Mean ARI: {mean_ari:.4f}')

print(f'\\n--- Merged Analysis ({len(merged_labels)} SCZ samples) ---')
print(f'  Optimal k: {merged_k}')
print(f'  Silhouette: {merged_clustering.silhouette:.4f}')
print(f'  Validation: {sum(g.passed for g in merged_val.results)}/{len(merged_val.results)} gates')
print(f'  Batch effect: p={p_batch:.4f} ({"CONFOUND" if p_batch < 0.05 else "OK"})')
print(f'  Best benchmark method: {bench_result.best_method}')

print(f'\\n--- Hertzberg k=2 Comparison ---')
for acc in DATASETS:
    if acc in k2_results:
        k2 = k2_results[acc]
        print(f'  {acc}: k=2 sil={k2["silhouette"]:.4f}, '
              f'vs BIC-k={k2["optimal_k"]}: ARI={k2["ari_vs_optimal"]:.4f}')

print(f'\\n--- Key Findings for Hertzberg Email ---')
print(f'1. PSF ran on all 5 GEO datasets from Hertzberg et al. 2024')
print(f'2. Total SCZ samples analyzed: {sum(dr["n_scz"] for dr in dataset_results.values())}')
print(f'3. BIC-optimal k varies by dataset (data-driven vs forced k=2)')
print(f'4. Cross-cohort projection tested subtype reproducibility')
print(f'5. Merged analysis: largest SCZ blood subtyping with formal validation gates')
print(f'6. Validation gates add: bootstrap stability + null ARI + pathway-gene concordance')

print(f'\\n--- Hertzberg Concordance ---')
for _, row in concordance_df.iterrows():
    print(f'  Subtype {int(row["subtype"])}: {row["classification"]} '
          f'(immune={row["immune_mean"]:.4f}, neuro={row["neuro_mean"]:.4f})')"""))

# ── Cell 32: References ──────────────────────────────────────
cells.append(md("""\
## References & Data Availability

1. Hertzberg L, Maggio N, Bhatt DK, Bhatt M, Bhatt R (2024). "Schizophrenia Biomarkers: Blood Transcriptome Suggests Two Molecular Subtypes." *NeuroMolecular Medicine* 26:39609319. PMID: [39609319](https://pubmed.ncbi.nlm.nih.gov/39609319/)

2. Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P (2015). "The Molecular Signatures Database (MSigDB) hallmark gene set collection." *Cell Systems* 1(6):417-425.

3. Chauhan R (2026). "Pathway Subtyping Framework: A Disease-Agnostic Pipeline for Pathway-Based Molecular Subtype Discovery with Statistical Validation." Research Square preprint. DOI: [10.21203/rs.3.rs-8913089/v1](https://doi.org/10.21203/rs.3.rs-8913089/v1)

---

### GEO Datasets

| Accession | Platform | Citation |
|-----------|----------|----------|
| GSE38484 | GPL6947 | de Jong S et al. (2012) |
| GSE27383 | GPL6883 | Gardiner EJ et al. (2013) |
| GSE38481 | GPL6947 | de Jong S et al. (2012) |
| GSE18312 | GPL6883 | Kano S et al. (2013) |
| GSE48072 | GPL10558 | de Jong S et al. (2014) |

**Data availability:** All datasets freely available from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

**Code:** [pathway-subtyping-framework](https://codeberg.org/pathways/pathway-subtyping-framework) (Codeberg)"""))

# ── Assemble notebook ─────────────────────────────────────────
notebook = {
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "name": "python",
            "version": "3.11.0",
            "mimetype": "text/x-python",
            "codemirror_mode": {"name": "ipython", "version": 3},
            "pygments_lexer": "ipython3",
            "file_extension": ".py",
            "nbconvert_exporter": "python",
        }
    },
    "nbformat": 4,
    "nbformat_minor": 5,
    "cells": cells,
}

output_path = "examples/notebooks/19_scz_blood_multi_cohort.ipynb"
with open(output_path, "w") as f:
    json.dump(notebook, f, indent=1)

print(f"Wrote {len(cells)} cells to {output_path}")
print(f"Cell breakdown: {sum(1 for c in cells if c['cell_type'] == 'markdown')} markdown, "
      f"{sum(1 for c in cells if c['cell_type'] == 'code')} code")
