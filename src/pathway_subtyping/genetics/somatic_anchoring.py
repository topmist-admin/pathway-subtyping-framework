"""
Somatic genetic anchoring — subject-level alignment core (Gate 7 somatic mode).

The feature-level anchoring gate (``gwas_enrichment.py``) asks whether the *genes
defining* a subtype are enriched for germline disease-risk genes — a gene-set
over-representation test, and the mode appropriate to psychiatry, where germline
risk is diffuse and per-subject genetic signal is low.

Cancer inverts that trade-off. Tumors carry **somatic drivers, copy-number
alterations, and mutational signatures** — high per-tumor genetic signal, fixed
in the tumor genome, causally upstream of the transcriptome, and often available
at a lower access tier than germline genotypes. So the cancer anchor is a
*subject-level* test: do the tumors in one transcriptomic subtype carry a somatic
stratum (e.g. BRAF-V600E, KRAS, MSI-high) more than the others?

This module computes that alignment: for each named somatic stratum, a chi-square
test of independence against the partition plus a Cramér's V effect size,
BH-adjusted across strata. It is the exact statistic of the confound gate
(``ValidationGates.confound_association_gate``) with **inverted polarity** — there
a strong categorical association is bad (a nuisance classifier); here it is the
positive evidence, because a somatic-driver association cannot be manufactured by
any transcriptomic/technical artifact.

Confound caveat (read before trusting a pass): a somatic driver can co-vary with
**tissue of origin** in a pan-cancer pool and with **tumor purity** within a
tumor type. Run the confound gate first so a "somatic anchor" is not merely the
tissue/purity axis re-labelled. And like feature-level anchoring this is a
*specific* test — a null is weak evidence of absence, not proof.

Research use only. Not for clinical decision-making.
"""

from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass
class StratumAlignment:
    """Association between the partition and one somatic stratum.

    Attributes:
        stratum: Name of the somatic stratum (e.g. "BRAF_V600E", "MSI").
        chi2: Chi-square statistic (None if the table was degenerate).
        p_value: Raw chi-square p-value (None if degenerate).
        cramers_v: Bergsma-corrected Cramér's V effect size in [0, 1].
        dof: Degrees of freedom of the contingency table.
        n: Number of samples with a non-missing stratum value.
    """

    stratum: str
    chi2: Any
    p_value: Any
    cramers_v: float
    dof: int
    n: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stratum": self.stratum,
            "chi2": (round(self.chi2, 4) if isinstance(self.chi2, float) else self.chi2),
            "p_value": self.p_value,
            "cramers_v": round(self.cramers_v, 4),
            "dof": self.dof,
            "n": self.n,
        }


def somatic_alignment(
    cluster_labels,
    somatic_strata: Dict[str, Any],
    cramers_v_min: float = 0.30,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """Test each somatic stratum for association with the partition.

    Args:
        cluster_labels: Array-like of cluster assignments (n_samples,).
        somatic_strata: Mapping of stratum name -> per-sample categorical value
            (array-like of length n_samples; e.g. driver carrier status,
            MSI-high/low, a mutational-signature class). Same-length as labels.
        cramers_v_min: Minimum Cramér's V for a stratum to count as an anchor
            (0.30 ~ medium effect).
        alpha: Significance level for the BH-adjusted chi-square p-values.

    Returns:
        Dict with keys: ``per_stratum`` (name -> stat dict), ``anchored_strata``
        (list), ``best_stratum``, ``max_cramers_v``, ``passed`` (True if any
        stratum is significantly associated AND V >= cramers_v_min).
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2_contingency

    from ..statistical_rigor import benjamini_hochberg
    from ..validation import cramers_v as _cramers_v

    labels = np.asarray(list(cluster_labels))
    n_total = len(labels)

    per_stratum: Dict[str, Any] = {}
    names: List[str] = []
    pvals: List[float] = []
    for name, values in somatic_strata.items():
        vals = np.asarray(list(values))
        if len(vals) != n_total:
            raise ValueError(f"Somatic stratum '{name}' length {len(vals)} != n_samples {n_total}")
        # Drop samples missing this stratum (None/NaN) so each test uses its own n.
        mask = pd.notna(pd.Series(vals))
        lab = labels[mask.to_numpy()]
        v = vals[mask.to_numpy()]
        table = pd.crosstab(pd.Series(lab, name="cluster"), pd.Series(v, name=name)).values
        if table.shape[0] < 2 or table.shape[1] < 2:
            per_stratum[name] = StratumAlignment(
                stratum=name,
                chi2=None,
                p_value=None,
                cramers_v=0.0,
                dof=0,
                n=int(mask.sum()),
            ).to_dict()
            per_stratum[name]["note"] = "degenerate table (single cluster or single stratum level)"
            continue
        chi2, p, dof, _ = chi2_contingency(table, correction=False)
        entry = StratumAlignment(
            stratum=name,
            chi2=float(chi2),
            p_value=float(p),
            cramers_v=float(_cramers_v(table)),
            dof=int(dof),
            n=int(mask.sum()),
        ).to_dict()
        per_stratum[name] = entry
        names.append(name)
        pvals.append(float(p))

    if pvals:
        qvals = benjamini_hochberg(np.asarray(pvals), alpha=alpha)
        for i, name in enumerate(names):
            per_stratum[name]["p_adjusted"] = float(qvals[i])
            per_stratum[name]["significant"] = bool(qvals[i] < alpha)

    anchored: List[str] = []
    best_v = 0.0
    best_name = None
    for name, info in per_stratum.items():
        v = info.get("cramers_v", 0.0) or 0.0
        if v > best_v:
            best_v = v
            best_name = name
        if info.get("significant", False) and v >= cramers_v_min:
            per_stratum[name]["anchored"] = True
            anchored.append(name)
        else:
            per_stratum[name]["anchored"] = False

    return {
        "per_stratum": per_stratum,
        "anchored_strata": anchored,
        "best_stratum": best_name,
        "max_cramers_v": best_v,
        "passed": len(anchored) > 0,
    }
