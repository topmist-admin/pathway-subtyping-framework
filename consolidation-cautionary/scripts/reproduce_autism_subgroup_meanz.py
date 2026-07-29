#!/usr/bin/env python3
"""§3.1 autism frontal-cortex subgroup — reproducible under PUBLIC v0.7.0 (mean-z).

WHY MEAN-Z (not ssGSEA): the current public framework release scores GSE28521
(Illumina microarray) with ssGSEA through a variance filter that culls all but a
handful of Hallmark pathways on this platform (47/50 become zero-variance and are
removed), leaving too few features to cluster and selecting a spurious high k.
The mean-of-per-gene-z-scores summary ("mean_z") is stable on the array platform,
scores all 50 Hallmark sets, and recovers the k=2 disease-enriched partition the
manuscript reports. All framework summaries z-normalize pathway scores, so the
matrix is unit-variance (see §2.2). Validation uses the framework's OWN routines
(ValidationGates.stability_test_bootstrap, n_bootstrap=100; select_n_clusters BIC;
run_clustering GMM), so §3.1's bootstrap ARI is apples-to-apples with the SCZ 0.923.

Canonical result (public v0.7.0, seed 20260708):
  50/50 Hallmark scored; BIC k=2; baseline n=21 (7 ASD), disease-enriched n=11
  (9/11 ASD); silhouette 0.25; bootstrap stability ARI 0.71 -> FAILS (<0.80);
  random-gene-set ARI 0.24; label-permutation ~0.001 (uninformative).
Conclusion: the subgroup fails bootstrap stability -> not a validated subtype
(exploratory). The exact control values are scoring-summary-dependent (an older
ssGSEA run gave lower stability / higher random-set ARI) -- reported honestly in
the manuscript as method-fragile; the load-bearing failure is on stability.

Public inputs only: GSE28521 series matrix + GPL6883.annot (GEO); Hallmark gmt
(gene symbols). Framework: pip install pathway-subtyping==0.8.0 (deposited run used 0.7.0).
"""
import argparse, gzip, io, json, logging, collections
import numpy as np, pandas as pd
logging.getLogger("pathway_subtyping").setLevel(logging.ERROR)
from pathway_subtyping.expression import score_pathways_from_expression, ExpressionScoringMethod
from pathway_subtyping.clustering import ClusteringAlgorithm, run_clustering, select_n_clusters
from pathway_subtyping.validation import ValidationGates
from sklearn.metrics import adjusted_rand_score, silhouette_score


def parse_series_expr(path):
    raw = gzip.open(path, "rt").read()
    tbl = raw.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    e = pd.read_csv(io.StringIO(tbl), sep="\t", index_col=0)
    e.index = e.index.astype(str).str.replace('"', ""); e.columns = [c.replace('"', "") for c in e.columns]
    titles = geos = None
    for line in raw.splitlines():
        if line.startswith("!Sample_title"): titles = [x.strip('"') for x in line.split("\t")[1:]]
        if line.startswith("!Sample_geo_accession"): geos = [x.strip('"') for x in line.split("\t")[1:]]
    meta = pd.DataFrame({"title": titles}, index=geos)
    meta["region"] = meta["title"].str[-1].map({"C": "Cerebellum", "F": "Frontal", "T": "Temporal"})
    meta["dx"] = np.where(meta["title"].str[0] == "A", "ASD", "CTL")
    return e, meta


def probe_to_symbol(annot_path):
    ann = gzip.open(annot_path, "rt", errors="ignore").read().splitlines()
    hidx = next(i for i, l in enumerate(ann) if l.split("\t")[0] == "ID")
    scol = ann[hidx].split("\t").index("Gene symbol")
    return {f[0]: f[scol] for l in ann[hidx + 1:] for f in [l.split("\t")]
            if len(f) > scol and f[0].startswith("ILMN") and f[scol]}


def parse_gmt(path):
    return {p[0]: p[2:] for l in open(path) for p in [l.rstrip("\n").split("\t")] if len(p) > 2}


def meanz(mat, pw, seed):
    r = score_pathways_from_expression(mat, pw, method=ExpressionScoringMethod.MEAN_Z,
                                       min_genes_per_pathway=2, seed=seed, show_progress=False)
    S = r.pathway_scores
    return S if S.shape[0] == mat.shape[0] else S.T


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _framework_version():
    try:
        import pathway_subtyping
        return pathway_subtyping.__version__
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # GEO inputs are optional: supply --genes to reproduce OFFLINE from the
    # deposited gene matrix instead. GEO series matrices and .annot files are
    # revised in place without a version bump, so re-fetching is inherently less
    # reproducible than reading the deposited intermediate -- which is why the
    # other packages in this bundle ship theirs.
    ap.add_argument("--gse28521-matrix", help="GSE28521_series_matrix.txt.gz (omit if --genes)")
    ap.add_argument("--gse28521-annot", help="GPL6883.annot.gz (omit if --genes)")
    ap.add_argument("--genes", help="deposited samples x genes matrix (.csv.gz) with a dx column; "
                                    "reproduces offline, no GEO access")
    ap.add_argument("--gmt", required=True, help="Hallmark gmt (gene symbols)")
    ap.add_argument("--out", default="out_autism_subgroup", help="output dir")
    ap.add_argument("--seed", type=int, default=20260708)
    ap.add_argument("--n-random", type=int, default=40, help="random-gene-set draws")
    a = ap.parse_args()
    import os; os.makedirs(a.out, exist_ok=True)

    inputs_sha = {}
    if a.genes:
        g = pd.read_csv(a.genes, index_col=0)
        dx = g.pop("dx").values
        sym = g                                             # samples x genes
        frontal = list(sym.index)
        inputs_sha["genes"] = _sha256(a.genes)
    else:
        if not (a.gse28521_matrix and a.gse28521_annot):
            ap.error("supply either --genes, or both --gse28521-matrix and --gse28521-annot")
        expr, meta = parse_series_expr(a.gse28521_matrix)
        frontal = [g for g in meta.index[meta["region"] == "Frontal"] if g in expr.columns]
        p2s = probe_to_symbol(a.gse28521_annot)
        e = expr.loc[[p for p in expr.index if p in p2s], frontal].astype(float)
        e.index = [p2s[p] for p in e.index]
        sym = e.groupby(level=0).mean().T                   # samples x genes
        dx = meta.loc[frontal, "dx"].values
        inputs_sha["gse28521_series_matrix"] = _sha256(a.gse28521_matrix)
        inputs_sha["gpl6883_annot"] = _sha256(a.gse28521_annot)
        # Deposit the derived matrix so the analysis reproduces without GEO.
        out_genes = sym.copy(); out_genes["dx"] = dx
        out_genes.to_csv(f"{a.out}/autism_subgroup_genes.csv.gz")
    inputs_sha["hallmark_gmt"] = _sha256(a.gmt)
    pw = parse_gmt(a.gmt)

    ps = meanz(sym, pw, a.seed)
    ms = select_n_clusters(ps.values, list(range(2, 7)), method="bic", seed=a.seed)
    lab = run_clustering(ps.values, 2, ClusteringAlgorithm.GMM, seed=a.seed).labels
    gates = ValidationGates(seed=a.seed, n_bootstrap=100, show_progress=False)
    stab = gates.stability_test_bootstrap(ps, lab, n_clusters=2, gmm_seed=a.seed)

    # disease-enriched = the smaller / higher-%ASD cluster
    comp = {int(c): {"n": int((lab == c).sum()),
                     "asd": int(((lab == c) & (dx == "ASD")).sum())} for c in set(lab)}
    dis = max(comp, key=lambda c: comp[c]["asd"] / comp[c]["n"])

    allg = list(sym.columns); rng = np.random.default_rng(a.seed); szs = [len(v) for v in pw.values()]; rl = []
    for _ in range(a.n_random):
        rp = {f"R{i}": list(rng.choice(allg, min(s, len(allg)), replace=False)) for i, s in enumerate(szs)}
        Sr = meanz(sym, rp, a.seed)
        rl.append(adjusted_rand_score(lab, run_clustering(Sr.values, 2, ClusteringAlgorithm.GMM, seed=a.seed).labels))

    res = dict(n_samples=len(frontal), pathways_scored=int(ps.shape[1]), bic_k=int(ms.optimal_k),
               composition=comp, disease_enriched_cluster=dis,
               silhouette=round(float(silhouette_score(ps.values, lab)), 3),
               bootstrap_ari=round(float(stab.metric_value), 3), bootstrap_pass=bool(stab.passed),
               random_gene_set_ari=round(float(np.mean(rl)), 3), seed=a.seed, method="mean_z",
               framework_version=_framework_version(), n_random_draws=int(a.n_random),
               input_sha256=inputs_sha,
               provenance=("Seed 20260708 is this analysis's documented seed (see module docstring), "
                           "NOT the seed 42 used elsewhere. Re-running at 42 gives different numbers. "
                           "Not to be confused with research-results/GSE28521/frontal-cortex/, which is "
                           "a DIFFERENT analysis: 15 curated autism pathways under framework 0.3.0 at "
                           "seed 42, silhouette 0.299."))
    pd.DataFrame([{"sample": s, "cluster": int(l), "dx": d} for s, l, d in zip(frontal, lab, dx)]
                 ).to_csv(f"{a.out}/autism_subgroup_labels.csv", index=False)
    json.dump(res, open(f"{a.out}/autism_subgroup_result.json", "w"), indent=2)

    print(f"frontal samples     : {res['n_samples']}  (16 ASD / 16 CTL)")
    print(f"pathways scored     : {res['pathways_scored']} / {len(pw)}")
    print(f"BIC-selected k      : {res['bic_k']}")
    for c in sorted(comp): print(f"  cluster {c}: n={comp[c]['n']:2d}  ASD={comp[c]['asd']}/{comp[c]['n']}")
    print(f"disease-enriched    : cluster {dis}  (n={comp[dis]['n']}, {comp[dis]['asd']}/{comp[dis]['n']} ASD)")
    print(f"silhouette          : {res['silhouette']}")
    print(f"bootstrap ARI       : {res['bootstrap_ari']}  passed>=0.80? {res['bootstrap_pass']}")
    print(f"random-gene-set ARI : {res['random_gene_set_ari']}")
    ok = (res['bic_k'] == 2 and not res['bootstrap_pass'])
    print(f"\n§3.1 reproduces (k=2 & fails bootstrap stability): {ok}")


if __name__ == "__main__":
    main()
