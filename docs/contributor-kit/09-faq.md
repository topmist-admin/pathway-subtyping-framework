# Frequently Asked Questions

## General Questions

### What is the Pathway Subtyping Framework?

The Pathway Subtyping Framework is an open-source tool for identifying molecular subtypes in genetically heterogeneous diseases. It works by:
1. Aggregating rare variant burden into biological pathways
2. Clustering patients based on pathway profiles
3. Validating subtypes with rigorous statistical tests

### Why use pathways instead of individual genes?

Genetic diseases often involve hundreds of genes, making gene-level analysis underpowered. Pathways:
- Reduce dimensionality (hundreds of genes → ~10 pathways)
- Provide biological interpretability
- Enable discovery of shared mechanisms across different mutations

### What diseases can this framework analyze?

Any genetically heterogeneous disease where:
- Multiple genes contribute to disease risk
- Genes can be grouped into biological pathways
- Cohort data (VCF + phenotypes) is available

Current/planned applications:
- Autism Spectrum Disorder
- Schizophrenia
- Epilepsy
- Intellectual Disability
- Parkinson's Disease
- Bipolar Disorder

---

## Data Questions

### How much data do I need?

| Analysis | Minimum | Recommended |
|----------|---------|-------------|
| Pilot/exploration | 100 samples | 500+ samples |
| Subtype discovery | 500 samples | 1,000+ samples |
| Phenotype correlations | 1,000 samples | 2,000+ samples |
| Robust replication | 2,000+ samples | 5,000+ samples |

### What data formats are supported?

**VCF files:**
- Standard VCF 4.x format
- Must include INFO fields: GENE, CONSEQUENCE, CADD
- Can be gzipped (.vcf.gz)

**Phenotype files:**
- CSV format
- Must have `sample_id` column matching VCF

**Pathway files:**
- GMT (Gene Matrix Transposed) format
- Tab-separated: PATHWAY_NAME, DESCRIPTION, GENE1, GENE2, ...

### Can I use GWAS data instead of sequencing?

The framework is designed for rare variant analysis from sequencing data (WES/WGS). GWAS common variants require different approaches. However, you could:
- Use rare variant burden from imputed GWAS
- Combine with sequencing data for subset of samples

### How do I handle multi-ethnic cohorts?

Population stratification can confound clustering. The framework provides built-in ancestry correction (v0.2):

**Recommended workflow:**
1. Compute ancestry PCs from your genotype data using `compute_ancestry_pcs()`
2. Adjust pathway scores with `adjust_pathway_scores()` to regress out ancestry effects
3. Run the pipeline on adjusted scores — the ancestry independence validation gate will automatically verify clusters are not confounded
4. Use `stratified_analysis()` to confirm subtypes replicate within ancestry groups

**Configuration (YAML):**
```yaml
ancestry:
  pcs_path: data/ancestry_pcs.csv
  correction: regress_out
  n_pcs: 10
```

**Programmatic usage:**
```python
from pathway_subtyping import compute_ancestry_pcs, adjust_pathway_scores

pcs = compute_ancestry_pcs(genotype_matrix, n_components=10, seed=42)
result = adjust_pathway_scores(pathway_scores, pcs)
# result.confounded_pathways lists pathways with R² > 0.1
# result.adjusted_scores contains corrected scores
```

**Additional best practices:**
- Use population-matched gnomAD frequencies for variant filtering
- Check `result.r_squared_per_pathway` to identify heavily confounded pathways
- For known ancestry groups, run `stratified_analysis()` to verify cross-group consistency

---

## Technical Questions

### Why did my pipeline fail?

Common failure reasons and solutions:

| Error | Likely Cause | Solution |
|-------|--------------|----------|
| VCF parsing error | Missing annotations | Ensure GENE/CONSEQUENCE/CADD in INFO |
| No variants after filtering | Filters too strict | Lower CADD threshold or include more consequences |
| Memory error | Large cohort | Increase RAM or process chromosomes separately |
| Clustering failed | Too few samples | Need minimum ~100 samples with variants |

### How do I choose the number of clusters (k)?

The framework automatically selects k using BIC (Bayesian Information Criterion):
- Lower BIC = better model fit
- Tests k = 2 through 8 by default
- Selects k with minimum BIC

You can override with a fixed k in the config:
```yaml
clustering:
  n_clusters: 4  # Force 4 clusters
```

### What do the validation gates mean?

| Gate | What it tests | Pass means |
|------|---------------|------------|
| **Label Shuffle** | Compares to random labels | Clusters are non-random |
| **Random Genes** | Compares to random gene sets | Pathway biology matters |
| **Bootstrap** | Cluster stability | Subtypes are reproducible |

All three must pass for valid subtypes.

### Can I add custom validation tests?

Yes, extend the `ValidationGates` class:

```python
from pathway_subtyping.validation import ValidationGates

class CustomValidationGates(ValidationGates):
    def my_custom_test(self, ...):
        # Custom validation logic
        pass
```

### How do I parallelize for large cohorts?

For cohorts > 10,000 samples:

1. **Split by chromosome:**
```bash
for chr in {1..22}; do
    psf --config config.yaml --chr $chr &
done
```

2. **Use cloud computing:**
- AWS: r5.4xlarge or larger
- Configure multiprocessing in config

3. **Reduce iterations:**
```yaml
validation:
  bootstrap_iterations: 50  # Reduce from 100
```

---

## Pathway Questions

### How do I create pathway definitions?

See [pathway-curation-guide.md](../guides/pathway-curation-guide.md) for detailed instructions.

Quick summary:
1. Review disease literature for implicated pathways
2. Get gene lists from GO, KEGG, or manual curation
3. Create GMT file with pathway names and gene lists
4. Validate against your VCF

### How many pathways should I include?

| Pathways | Trade-off |
|----------|-----------|
| 3-5 | Interpretable, may miss biology |
| 6-10 | Good balance |
| 10-20 | More comprehensive, harder to interpret |
| 20+ | Risk of overfitting |

**Recommendation**: Start with 4-6 well-characterized pathways.

### Can pathways overlap (share genes)?

Yes, some overlap is expected and acceptable. However:
- High overlap can cause correlated scores
- Consider hierarchical pathways (parent/child)
- Document overlap in your analysis

### Where do I find pathway gene lists?

| Source | URL | Notes |
|--------|-----|-------|
| Gene Ontology | geneontology.org | Biological process terms |
| KEGG | genome.jp/kegg | Curated pathways |
| Reactome | reactome.org | Signaling pathways |
| MSigDB | gsea-msigdb.org | Comprehensive collections |
| Disease-specific | PubMed literature | Manual curation |

---

## Results Questions

### My subtypes don't make biological sense. What now?

Possible issues:
1. **Pathway definitions**: May not capture relevant biology
2. **Variant filtering**: May be too strict/loose
3. **Sample size**: May need more samples
4. **Confounders**: Check for batch effects, ancestry

Try:
- Different pathway definitions
- Adjust variant filters
- Stratify by potential confounders

### How do I interpret pathway scores?

Pathway scores are z-normalized burden:
- **Score = 0**: Average burden
- **Score > 0**: Higher than average (more variants)
- **Score < 0**: Lower than average (fewer variants)

For subtypes:
- High positive score = subtype's dominant disruption
- Interpret in context of disease biology

### Can I validate subtypes in an independent cohort?

Yes, this is highly recommended:

1. Train model on discovery cohort
2. Apply to validation cohort (same pathways)
3. Compare cluster assignments
4. Check phenotype correlations replicate

### How do I report negative results?

If validation gates fail or no clear subtypes emerge:
- This is still valuable scientific information
- Document your approach thoroughly
- Report null results in publication
- Consider as pilot for future studies

---

## Access and Collaboration

### How long does data access take?

| Repository | Typical Timeline |
|------------|------------------|
| SFARI | 4-6 weeks |
| UK Biobank | 8-12 weeks |
| dbGaP | 6-8 weeks |
| Consortium data | Variable (contact PI) |

**Start early** - apply while learning the framework.

### Can I collaborate with external researchers?

Yes, with considerations:
- Check DUA for collaboration terms
- External collaborators may need own data access
- May need to update IRB protocol
- Document in data access amendment if required

### How do I contribute to the framework?

1. **Bug reports**: GitHub Issues
2. **Feature requests**: GitHub Issues or discussions
3. **Code contributions**: Pull requests
4. **Pathway contributions**: Submit GMT files
5. **Documentation**: Improve guides

### Who do I contact for help?

| Issue | Contact |
|-------|---------|
| Technical bugs | GitHub Issues |
| Data access | Repository data access committees |
| Framework questions | Project lead |
| Collaboration | Project lead |

---

## Legal and Compliance

### Is this research covered by my institution's IRB?

Typically yes, as secondary research with de-identified data. Check with your IRB office - you'll likely need:
- IRB exemption determination
- Protocol registered (even if exempt)

### Can I use results commercially?

Check your DUA carefully. Most academic data access:
- Permits non-commercial research
- Requires separate agreement for commercial use
- May restrict IP claims

### What if I discover a reportable finding?

For incidental findings with clinical implications:
- Most research DUAs do not require re-contact
- Consult your IRB and ethics board
- Follow institutional guidelines

---

## Getting More Help

### Documentation
- [GitHub README](https://github.com/topmist-admin/pathway-subtyping-framework)
- [Getting Started Guide](01-getting-started.md)
- [Technical Setup](06-technical-setup.md)

### Community
- GitHub Discussions
- Team Slack/Teams channel
- Weekly office hours (check calendar)

### Training
- Demo notebooks in `examples/notebooks/`
- Video tutorials (coming soon)
- Internal training sessions

---

*Don't see your question? Add it to GitHub Discussions or contact the project lead.*
