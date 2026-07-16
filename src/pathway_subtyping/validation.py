"""
Validation Gates for the Pathway Subtyping Framework.

Implements mandatory negative controls and stability tests to validate
clustering results before accepting them as meaningful.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from tqdm import tqdm

from .utils.seed import get_rng

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Result from a single validation test."""

    name: str
    passed: bool
    metric_name: str
    metric_value: float
    threshold: float
    comparison: str  # ">" or "<" or "=="
    details: Dict[str, Any] = field(default_factory=dict)

    @property
    def status(self) -> str:
        return "PASS" if self.passed else "FAIL"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "metric": self.metric_name,
            "value": round(self.metric_value, 4),
            "threshold": self.threshold,
            "comparison": self.comparison,
            "details": self.details,
        }


@dataclass
class ValidationGatesResult:
    """Aggregated result from all validation gates."""

    results: List[ValidationResult]
    all_passed: bool
    summary: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "all_passed": self.all_passed,
            "summary": self.summary,
            "tests": [r.to_dict() for r in self.results],
        }


class ValidationGates:
    """
    Validation gates for clustering results.

    Implements:
    - Negative Control 1: Label shuffle (expect ARI ~0)
    - Negative Control 2: Random gene sets (expect no enrichment)
    - Stability Test: Bootstrap resampling (expect ARI >= threshold)
    """

    def __init__(
        self,
        seed: Optional[int] = 42,
        n_permutations: int = 100,
        n_bootstrap: int = 50,
        stability_threshold: float = 0.8,
        null_ari_max: float = 0.15,
        show_progress: bool = True,
    ):
        """
        Initialize validation gates.

        Args:
            seed: Random seed for reproducibility
            n_permutations: Number of permutations for null tests
            n_bootstrap: Number of bootstrap iterations for stability
            stability_threshold: Minimum ARI for stability test (default: 0.8)
            null_ari_max: Maximum expected ARI under null hypothesis (default: 0.15)
            show_progress: Show tqdm progress bars for long-running loops
        """
        self.seed = seed
        self.n_permutations = n_permutations
        self.n_bootstrap = n_bootstrap
        self.stability_threshold = stability_threshold
        self.null_ari_max = null_ari_max
        self.show_progress = show_progress
        self.rng = get_rng(seed, "validation")

    def run_all(
        self,
        pathway_scores: pd.DataFrame,
        cluster_labels: np.ndarray,
        pathways: Dict[str, List[str]],
        gene_burdens: pd.DataFrame,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
        ancestry_pcs=None,
        per_modality_scores: Optional[Dict[str, pd.DataFrame]] = None,
        confounds: Optional[Dict[str, Any]] = None,
        genetic_anchoring: Optional[Dict[str, Any]] = None,
    ) -> ValidationGatesResult:
        """
        Run all validation gates.

        Args:
            pathway_scores: DataFrame of pathway scores (samples x pathways)
            cluster_labels: Array of cluster assignments
            pathways: Dict mapping pathway names to gene lists
            gene_burdens: DataFrame of gene burdens (samples x genes)
            n_clusters: Number of clusters used
            gmm_seed: Seed for GMM clustering
            ancestry_pcs: Optional AncestryPCs for ancestry independence gate
            per_modality_scores: Optional dict mapping modality label to pathway
                score DataFrame for cross-modal concordance gate
            confounds: Optional dict mapping confound name (e.g. 'region',
                'batch', 'diagnosis') to per-sample categorical values, for the
                confound association gate. Diagnosis-like keys are treated as
                biology-of-interest and never fail the gate.
            genetic_anchoring: Optional kwargs for the Gate 7 genetic-anchoring
                gate, forwarded to ``genetic_anchoring_gate`` (requires at least
                ``subtype_gene_sets``, ``risk_genes``, ``background_universe``).
                Gate 7 is positive-evidence and low-power; when supplied it is
                appended like any other gate.

        Returns:
            ValidationGatesResult with all test outcomes
        """
        logger.info("Running validation gates...")

        results = []

        # Negative Control 1: Label shuffle
        nc1_result = self.negative_control_label_shuffle(
            pathway_scores, cluster_labels, n_clusters, gmm_seed
        )
        results.append(nc1_result)
        logger.info(
            f"  - {nc1_result.name}: {nc1_result.status} (ARI = {nc1_result.metric_value:.3f})"
        )

        # Negative Control 2: Random gene sets
        nc2_result = self.negative_control_random_gene_sets(
            gene_burdens, pathways, cluster_labels, n_clusters, gmm_seed
        )
        results.append(nc2_result)
        logger.info(
            f"  - {nc2_result.name}: {nc2_result.status} (ARI = {nc2_result.metric_value:.3f})"
        )

        # Stability Test: Bootstrap
        stability_result = self.stability_test_bootstrap(
            pathway_scores, cluster_labels, n_clusters, gmm_seed
        )
        results.append(stability_result)
        logger.info(
            f"  - {stability_result.name}: {stability_result.status} "
            f"(ARI = {stability_result.metric_value:.3f})"
        )

        # Negative Control 3: Ancestry Independence (if ancestry PCs available)
        if ancestry_pcs is not None:
            nc3_result = self.negative_control_ancestry_independence(cluster_labels, ancestry_pcs)
            results.append(nc3_result)
            logger.info(
                f"  - {nc3_result.name}: {nc3_result.status} "
                f"(p_min = {nc3_result.metric_value:.4f})"
            )

        # Gate 5: Cross-Modal Concordance (if multi-omic data available)
        if per_modality_scores is not None and len(per_modality_scores) >= 2:
            cm_result = self.cross_modal_concordance_gate(
                per_modality_scores=per_modality_scores,
                cluster_labels=cluster_labels,
                fused_sample_ids=list(pathway_scores.index),
                n_clusters=n_clusters,
            )
            results.append(cm_result)
            logger.info(
                f"  - {cm_result.name}: {cm_result.status} " f"(ARI = {cm_result.metric_value:.3f})"
            )

        # Gate 6: Confound Association (if confound annotations available)
        if confounds:
            confound_result = self.confound_association_gate(cluster_labels, confounds)
            results.append(confound_result)
            logger.info(
                f"  - {confound_result.name}: {confound_result.status} "
                f"(max nuisance Cramer's V = {confound_result.metric_value:.3f})"
            )

        # Gate 7: Genetic Anchoring (if subtype gene sets + risk catalogue available)
        if genetic_anchoring:
            anchoring_result = self.genetic_anchoring_gate(**genetic_anchoring)
            results.append(anchoring_result)
            logger.info(
                f"  - {anchoring_result.name}: {anchoring_result.status} "
                f"(max enrichment fold = {anchoring_result.metric_value:.2f})"
            )

        # Aggregate results
        all_passed = all(r.passed for r in results)
        n_passed = sum(1 for r in results if r.passed)
        n_total = len(results)

        if all_passed:
            summary = f"All {n_total} validation gates PASSED"
        else:
            failed = [r.name for r in results if not r.passed]
            summary = f"{n_passed}/{n_total} validation gates passed. FAILED: {', '.join(failed)}"

        return ValidationGatesResult(
            results=results,
            all_passed=all_passed,
            summary=summary,
        )

    def negative_control_label_shuffle(
        self,
        pathway_scores: pd.DataFrame,
        original_labels: np.ndarray,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationResult:
        """
        Negative Control 1: Label shuffle test.

        Shuffles cluster labels randomly and re-runs clustering. The clustering
        should NOT recover the shuffled labels (ARI should be ~0).

        If clustering CAN recover shuffled labels, it suggests the clusters
        are an artifact of the method rather than real structure in the data.
        """
        ari_values = []

        for i in tqdm(
            range(self.n_permutations),
            desc="Label shuffle",
            disable=not self.show_progress,
        ):
            # Shuffle labels
            shuffled_labels = self.rng.permutation(original_labels)

            # Re-cluster
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=(gmm_seed + i) if gmm_seed else None,
                reg_covar=1e-6,  # Regularization for numerical stability
            )
            gmm.fit(pathway_scores.values)

            # Skip non-converged fits
            if not gmm.converged_:
                continue

            new_labels = gmm.predict(pathway_scores.values)

            # Compare to shuffled (should NOT match)
            ari = adjusted_rand_score(shuffled_labels, new_labels)
            ari_values.append(ari)

        # Handle edge case where no GMM fits converged
        if not ari_values:
            logger.warning("Label shuffle test: No GMM fits converged")
            return ValidationResult(
                name="Negative Control 1: Label Shuffle",
                passed=False,
                metric_name="mean_null_ARI",
                metric_value=1.0,  # Worst case - fail the test
                threshold=self.null_ari_max,
                comparison="<",
                details={
                    "error": "No GMM fits converged during permutation testing",
                    "n_permutations_attempted": self.n_permutations,
                    "n_successful": 0,
                },
            )

        mean_ari = np.mean(ari_values)
        std_ari = np.std(ari_values)

        # PASS if mean ARI is close to 0 (clustering doesn't recover random labels)
        passed = mean_ari < self.null_ari_max

        return ValidationResult(
            name="Negative Control 1: Label Shuffle",
            passed=passed,
            metric_name="mean_null_ARI",
            metric_value=mean_ari,
            threshold=self.null_ari_max,
            comparison="<",
            details={
                "std_ari": round(std_ari, 4),
                "n_permutations": self.n_permutations,
                "n_successful": len(ari_values),
                "interpretation": "Clustering should NOT recover shuffled labels",
            },
        )

    def negative_control_random_gene_sets(
        self,
        gene_burdens: pd.DataFrame,
        real_pathways: Dict[str, List[str]],
        original_labels: np.ndarray,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationResult:
        """
        Negative Control 2: Random gene sets test.

        Replaces biological pathway definitions with random gene sets of the
        same sizes. Clustering on random gene sets should NOT recover the
        original cluster structure (ARI should be ~0).

        If clustering CAN recover structure with random gene sets, it suggests
        the results are not driven by meaningful biological pathways.
        """
        all_genes = list(gene_burdens.columns)
        ari_values = []

        for i in tqdm(
            range(self.n_permutations),
            desc="Random gene sets",
            disable=not self.show_progress,
        ):
            # Create random pathways with same sizes
            random_pathway_scores = {}

            for pathway_name, pathway_genes in real_pathways.items():
                # How many genes from this pathway are in our burden data?
                n_genes_needed = min(len(pathway_genes), len(all_genes))
                if n_genes_needed < 2:
                    continue

                # Sample random genes
                random_genes = list(self.rng.choice(all_genes, size=n_genes_needed, replace=False))

                # Compute score for random pathway
                common_genes = [g for g in random_genes if g in gene_burdens.columns]
                if len(common_genes) >= 2:
                    random_pathway_scores[pathway_name] = gene_burdens[common_genes].mean(axis=1)

            if len(random_pathway_scores) < 2:
                continue

            # Build random pathway score matrix
            random_scores_df = pd.DataFrame(random_pathway_scores)
            random_scores_df = (random_scores_df - random_scores_df.mean()) / random_scores_df.std()
            random_scores_df = random_scores_df.fillna(0)

            # Cluster on random pathway scores
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=(gmm_seed + i + 1000) if gmm_seed else None,
                reg_covar=1e-6,  # Regularization for numerical stability
            )
            gmm.fit(random_scores_df.values)

            # Skip non-converged fits
            if not gmm.converged_:
                continue

            random_labels = gmm.predict(random_scores_df.values)

            # Compare to original clustering
            ari = adjusted_rand_score(original_labels, random_labels)
            ari_values.append(ari)

        if not ari_values:
            # Edge case: couldn't run any permutations
            return ValidationResult(
                name="Negative Control 2: Random Gene Sets",
                passed=False,
                metric_name="mean_random_ARI",
                metric_value=1.0,
                threshold=self.null_ari_max,
                comparison="<",
                details={"error": "Could not compute random gene set permutations"},
            )

        mean_ari = np.mean(ari_values)
        std_ari = np.std(ari_values)

        # PASS if mean ARI is close to 0 (random gene sets don't replicate structure)
        passed = mean_ari < self.null_ari_max

        return ValidationResult(
            name="Negative Control 2: Random Gene Sets",
            passed=passed,
            metric_name="mean_random_ARI",
            metric_value=mean_ari,
            threshold=self.null_ari_max,
            comparison="<",
            details={
                "std_ari": round(std_ari, 4),
                "n_permutations": len(ari_values),
                "interpretation": "Random gene sets should NOT replicate cluster structure",
            },
        )

    def stability_test_bootstrap(
        self,
        pathway_scores: pd.DataFrame,
        original_labels: np.ndarray,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> ValidationResult:
        """
        Stability Test: Bootstrap resampling.

        Resamples the data with replacement and re-runs clustering. Stable
        clusters should be recovered across bootstrap samples (high ARI).

        If clustering is unstable across bootstrap samples, it suggests the
        clusters are not robust features of the data.
        """
        n_samples = len(pathway_scores)
        ari_values = []

        for i in tqdm(
            range(self.n_bootstrap),
            desc="Bootstrap stability",
            disable=not self.show_progress,
        ):
            # Bootstrap sample (with replacement)
            bootstrap_idx = self.rng.choice(n_samples, size=n_samples, replace=True)

            # Get unique indices for clustering (avoid duplicate samples confusing GMM)
            unique_idx = np.unique(bootstrap_idx)
            if len(unique_idx) < n_clusters * 2:
                continue  # Not enough unique samples

            bootstrap_scores = pathway_scores.iloc[unique_idx]
            bootstrap_original_labels = original_labels[unique_idx]

            # Re-cluster on bootstrap sample
            gmm = GaussianMixture(
                n_components=n_clusters,
                covariance_type="full",
                n_init=5,
                random_state=(gmm_seed + i + 2000) if gmm_seed else None,
                reg_covar=1e-6,  # Regularization for numerical stability
            )
            try:
                gmm.fit(bootstrap_scores.values)

                # Skip non-converged fits
                if not gmm.converged_:
                    continue

                bootstrap_labels = gmm.predict(bootstrap_scores.values)

                # Compare to original labels for these samples
                ari = adjusted_rand_score(bootstrap_original_labels, bootstrap_labels)
                ari_values.append(ari)
            except Exception:
                continue  # Skip failed fits

        if not ari_values:
            logger.warning("Bootstrap stability test: No valid bootstrap iterations completed")
            return ValidationResult(
                name="Stability Test: Bootstrap",
                passed=False,
                metric_name="mean_bootstrap_ARI",
                metric_value=0.0,
                threshold=self.stability_threshold,
                comparison=">=",
                details={
                    "error": "No valid bootstrap iterations completed",
                    "n_bootstrap_attempted": self.n_bootstrap,
                    "n_successful": 0,
                },
            )

        mean_ari = np.mean(ari_values)
        std_ari = np.std(ari_values)

        # PASS if mean ARI >= threshold (clustering is stable)
        passed = mean_ari >= self.stability_threshold

        return ValidationResult(
            name="Stability Test: Bootstrap",
            passed=passed,
            metric_name="mean_bootstrap_ARI",
            metric_value=mean_ari,
            threshold=self.stability_threshold,
            comparison=">=",
            details={
                "std_ari": round(std_ari, 4),
                "n_bootstrap": len(ari_values),
                "interpretation": "Clustering should be stable across bootstrap samples",
            },
        )

    def negative_control_ancestry_independence(
        self,
        cluster_labels: np.ndarray,
        ancestry_pcs,
    ) -> ValidationResult:
        """
        Negative Control 3: Ancestry Independence Test.

        Tests whether discovered clusters are confounded with ancestry
        PCs. If clusters correlate significantly with ancestry, the
        subtypes may reflect population structure rather than biology.

        Pass criterion: No ancestry PC significantly associated with
        clusters after Bonferroni correction (p > 0.05 / n_PCs).
        """
        from .ancestry import check_ancestry_independence

        report = check_ancestry_independence(
            cluster_labels=cluster_labels,
            ancestry_pcs=ancestry_pcs,
            significance_threshold=0.05,
        )

        # Use minimum p-value as the test statistic
        min_pval = (
            min(report.independence_test_pvalues.values())
            if report.independence_test_pvalues
            else 1.0
        )
        n_pcs = len(report.independence_test_pvalues)
        bonferroni_threshold = 0.05 / n_pcs if n_pcs > 0 else 0.05

        return ValidationResult(
            name="Negative Control 3: Ancestry Independence",
            passed=report.overall_independence_passed,
            metric_name="min_ancestry_pvalue",
            metric_value=min_pval,
            threshold=bonferroni_threshold,
            comparison=">",
            details={
                "n_pcs_tested": n_pcs,
                "pvalues": {k: round(v, 6) for k, v in report.independence_test_pvalues.items()},
                "interpretation": "Clusters should NOT correlate with ancestry PCs",
            },
        )

    def cross_modal_concordance_gate(
        self,
        per_modality_scores: Dict[str, pd.DataFrame],
        cluster_labels: np.ndarray,
        fused_sample_ids: List[str],
        n_clusters: int,
    ) -> ValidationResult:
        """
        Gate 5: Cross-Modal Concordance (optional).

        Tests whether subtypes are consistent across data modalities.
        Clusters each modality independently and measures agreement.
        Only runs when multi-omic data is available (>= 2 modalities).
        """
        from .cross_modal_validation import cross_modal_concordance

        result = cross_modal_concordance(
            per_modality_scores=per_modality_scores,
            cluster_labels=cluster_labels,
            fused_sample_ids=fused_sample_ids,
            n_clusters=n_clusters,
            n_permutations=self.n_permutations,
            seed=self.seed,
            show_progress=self.show_progress,
        )

        return ValidationResult(
            name="Gate 5: Cross-Modal Concordance",
            passed=result.gate_passed,
            metric_name="mean_concordance_ARI",
            metric_value=result.mean_concordance_ari,
            threshold=result.null_ari_95th,
            comparison=">",
            details={
                "mean_concordance_nmi": round(result.mean_concordance_nmi, 4),
                "mean_transfer_ari": round(result.mean_transfer_ari, 4),
                "null_ari_95th": round(result.null_ari_95th, 4),
                "n_modality_pairs": len(result.pair_results),
                "pair_results": [p.to_dict() for p in result.pair_results],
                "interpretation": ("Subtypes should be consistent across data modalities"),
            },
        )

    def confound_association_gate(
        self,
        cluster_labels: np.ndarray,
        confounds: Dict[str, Any],
        cramers_v_max: float = 0.30,
        alpha: float = 0.05,
    ) -> ValidationResult:
        """
        Confound Association Gate (mandatory when confounds are available).

        Tests whether the partition is explained by a technical or anatomical
        confound (e.g. brain region, sequencing batch) rather than the biology
        of interest. For each named confound, computes a chi-square test of
        independence against the cluster labels and the Cramer's V effect size.

        This is the gate whose absence let an anatomy artifact pass the full
        validation battery: on GSE80655 a k=3 partition passed every stability
        control (bootstrap ARI ~0.92) yet corresponded almost exactly to brain
        region (Cramer's V ~0.67, p ~4e-26) while being independent of diagnosis
        (p ~0.41). A stability-passing partition can still be a confound
        classifier; this gate is the check that catches it.

        Pass criterion: NO confound is both statistically associated
        (BH-adjusted p < ``alpha``) AND non-trivially so (Cramer's V >=
        ``cramers_v_max``). Diagnosis is treated as the biology of interest, not
        a nuisance confound, and is reported for context but never fails the gate
        (a partition SHOULD track diagnosis).

        Args:
            cluster_labels: Array of cluster assignments (n_samples,).
            confounds: Mapping of confound name -> per-sample categorical values
                (array-like of length n_samples). Diagnosis, if present, should
                be keyed 'diagnosis' to be treated as biology-of-interest.
            cramers_v_max: Effect-size threshold above which an association is
                considered a confound (0.30 ~ medium; 0.10 small, 0.50 large).
            alpha: Significance level for the (BH-adjusted) chi-square p-values.

        Returns:
            ValidationResult. metric_value is the maximum Cramer's V across the
            *nuisance* confounds (the worst offender); passed is False if any
            nuisance confound is significant AND exceeds cramers_v_max.
        """
        from scipy.stats import chi2_contingency

        from .statistical_rigor import benjamini_hochberg

        labels = np.asarray(cluster_labels)
        n = len(labels)
        biology_keys = {"diagnosis", "dx", "disease", "condition", "phenotype"}

        per_confound: Dict[str, Any] = {}
        names: List[str] = []
        pvals: List[float] = []
        for name, values in confounds.items():
            vals = np.asarray(list(values))
            if len(vals) != n:
                logger.warning(
                    "Confound '%s' length %d != n_samples %d; skipping.",
                    name,
                    len(vals),
                    n,
                )
                continue
            # Contingency table cluster x confound-level.
            df = pd.crosstab(pd.Series(labels, name="cluster"), pd.Series(vals, name=name))
            table = df.values
            if table.shape[0] < 2 or table.shape[1] < 2:
                # Nothing to test (single cluster or single confound level).
                per_confound[name] = {
                    "chi2": None,
                    "p_value": None,
                    "cramers_v": 0.0,
                    "dof": 0,
                    "is_biology_of_interest": name.lower() in biology_keys,
                    "note": "degenerate table (single cluster or single level)",
                }
                continue
            chi2, p, dof, _ = chi2_contingency(table, correction=False)
            v = cramers_v(table)
            per_confound[name] = {
                "chi2": round(float(chi2), 4),
                "p_value": float(p),
                "cramers_v": round(float(v), 4),
                "dof": int(dof),
                "is_biology_of_interest": name.lower() in biology_keys,
            }
            names.append(name)
            pvals.append(float(p))

        # BH-adjust p-values across the confounds actually tested.
        # benjamini_hochberg returns q-values (adjusted p), not a boolean mask.
        if pvals:
            qvals = benjamini_hochberg(np.asarray(pvals), alpha=alpha)
            for i, name in enumerate(names):
                per_confound[name]["p_adjusted"] = float(qvals[i])
                per_confound[name]["significant"] = bool(qvals[i] < alpha)

        # Evaluate the gate over NUISANCE confounds only.
        worst_v = 0.0
        worst_name = None
        failing = []
        for name, info in per_confound.items():
            if info.get("is_biology_of_interest"):
                continue
            v = info.get("cramers_v", 0.0) or 0.0
            sig = info.get("significant", False)
            if v > worst_v:
                worst_v = v
                worst_name = name
            if sig and v >= cramers_v_max:
                failing.append(name)

        passed = len(failing) == 0
        interpretation = (
            "Partition is not explained by a nuisance confound."
            if passed
            else (
                f"Partition is confounded by: {', '.join(failing)} "
                f"(Cramer's V >= {cramers_v_max}). The subtype may be an "
                "artifact of this variable, not biology."
            )
        )

        return ValidationResult(
            name="Confound Association Gate",
            passed=passed,
            metric_name="max_nuisance_cramers_v",
            metric_value=worst_v,
            threshold=cramers_v_max,
            comparison="<",
            details={
                "worst_confound": worst_name,
                "failing_confounds": failing,
                "per_confound": {
                    k: {
                        kk: (round(vv, 6) if isinstance(vv, float) else vv)
                        for kk, vv in v.items()
                    }
                    for k, v in per_confound.items()
                },
                "cramers_v_max": cramers_v_max,
                "alpha": alpha,
                "interpretation": interpretation,
            },
        )

    def genetic_anchoring_gate(
        self,
        subtype_gene_sets: Dict[Any, Any],
        risk_genes: Any,
        background_universe: Any,
        reference_universe: Optional[Any] = None,
        alpha: float = 0.05,
        min_fold: float = 1.0,
    ) -> ValidationResult:
        """
        Gate 7: Genetic Anchoring Gate (feature-level).

        Tests whether the genes defining each subtype are over-represented for
        disease genetic-risk genes, by a hypergeometric over-representation test
        against a background-matched gene universe (BH-adjusted across subtypes).

        Unlike every other gate in the battery, this is a POSITIVE-EVIDENCE gate.
        The negative controls and the confound gate can only certify that a
        partition is *reproducible*; they can never show it is *real* or *causal*.
        A germline-variant enrichment is the exception: it cannot be manufactured
        by any postmortem or technical confound (PMI, RIN, dissection, batch,
        platform), because the variants are fixed at conception, upstream of all
        of them. So a subtype whose defining genes carry disease genetic risk is
        genetically implicated — necessary, not sufficient, but confound-immune.

        The null must be background-matched, not genome-wide. Brain-expressed
        genes are already enriched for brain-disease genetic risk, so a
        region-identity subtype could show spurious genome-wide enrichment. The
        gate decides on ``background_universe`` (e.g. brain-expressed genes);
        ``reference_universe`` (e.g. genome-wide protein-coding) is reported per
        subtype for contrast but never decides the gate.

        This is a SPECIFIC, LOW-POWER test: a non-anchored result is weak evidence
        of absence, never proof a subtype lacks genetic grounding. A pass reports
        which subtype axes are anchored; it does not upgrade an unstable or
        confounded partition (Gate 7 is complementary to, not a substitute for,
        the discreteness and confound gates).

        Pass criterion: at least one subtype's defining genes are BOTH
        significantly over-represented for risk genes (BH-adjusted p < ``alpha``)
        AND enriched in direction (fold > ``min_fold``) under the
        background-matched null.

        Args:
            subtype_gene_sets: Mapping of subtype label -> iterable of gene IDs
                defining that subtype. Identifiers must share a namespace with
                ``risk_genes`` and the universes (e.g. Ensembl gene IDs, version
                suffixes stripped).
            risk_genes: Disease genetic-risk reference gene IDs.
            background_universe: Background-matched universe of testable genes
                (the discriminating null, e.g. brain-expressed genes).
            reference_universe: Optional genome-wide universe reported alongside
                each subtype for contrast; does not affect the gate decision.
            alpha: Significance level for the BH-adjusted hypergeometric p-values.
            min_fold: Minimum fold enrichment for a subtype to count as anchored.

        Returns:
            ValidationResult. metric_value is the maximum fold enrichment across
            subtypes under the matched null; passed is True if any subtype is
            significantly anchored.
        """
        from .genetics.gwas_enrichment import hypergeometric_enrichment
        from .statistical_rigor import benjamini_hochberg

        risk = set(risk_genes)
        background = set(background_universe)
        reference = set(reference_universe) if reference_universe is not None else None

        labels: List[Any] = []
        matched: Dict[Any, Any] = {}
        pvals: List[float] = []
        for label, genes in subtype_gene_sets.items():
            res = hypergeometric_enrichment(
                genes, risk, background, label=str(label), null="background-matched"
            )
            matched[label] = res
            labels.append(label)
            # NaN p (undefined test) is treated as non-significant (p = 1.0).
            pvals.append(res.p_value if np.isfinite(res.p_value) else 1.0)

        qvals = benjamini_hochberg(np.asarray(pvals), alpha=alpha) if pvals else np.array([])

        per_subtype: Dict[str, Any] = {}
        anchored: List[str] = []
        best_fold = 0.0
        best_label: Optional[str] = None
        for i, label in enumerate(labels):
            res = matched[label]
            q = float(qvals[i])
            fold = res.fold if np.isfinite(res.fold) else 0.0
            is_anchored = bool(q < alpha and fold > min_fold and res.risk_hits > 0)
            entry = res.to_dict()
            entry["p_adjusted"] = q
            entry["anchored"] = is_anchored
            if reference is not None:
                entry["reference_null"] = hypergeometric_enrichment(
                    subtype_gene_sets[label], risk, reference,
                    label=str(label), null="genome-wide",
                ).to_dict()
            per_subtype[str(label)] = entry
            if fold > best_fold:
                best_fold = fold
                best_label = str(label)
            if is_anchored:
                anchored.append(str(label))

        passed = len(anchored) > 0
        interpretation = (
            f"Genetically anchored subtype axis/axes: {', '.join(anchored)} "
            "(defining genes over-represented for disease genetic risk under a "
            "background-matched null — confound-immune positive evidence)."
            if passed
            else (
                "No subtype's defining genes were significantly over-represented "
                "for disease genetic risk under the background-matched null. This "
                "is a specific, low-power test — a null is weak evidence of "
                "absence, not proof the partition lacks genetic grounding."
            )
        )

        return ValidationResult(
            name="Gate 7: Genetic Anchoring",
            passed=passed,
            metric_name="max_enrichment_fold",
            metric_value=best_fold,
            threshold=min_fold,
            comparison=">",
            details={
                "anchored_subtypes": anchored,
                "best_subtype": best_label,
                "risk_genes_n": len(risk),
                "background_universe_n": len(background),
                "reference_universe_n": (len(reference) if reference is not None else None),
                "per_subtype": per_subtype,
                "alpha": alpha,
                "min_fold": min_fold,
                "gate_polarity": "positive_evidence",
                "interpretation": interpretation,
            },
        )


def cramers_v(
    contingency: np.ndarray, bias_correction: bool = True
) -> float:
    """
    Cramér's V effect size for a contingency table (association strength between
    two categorical variables, 0 = independent, 1 = perfect association).

    Args:
        contingency: r x c contingency table of counts.
        bias_correction: Apply the Bergsma (2013) small-sample bias correction.

    Returns:
        Cramér's V in [0, 1]; 0.0 for degenerate tables (a single row/column,
        or no samples).
    """
    from scipy.stats import chi2_contingency

    table = np.asarray(contingency, dtype=float)
    n = table.sum()
    if n <= 0 or table.shape[0] < 2 or table.shape[1] < 2:
        return 0.0

    chi2 = chi2_contingency(table, correction=False)[0]
    phi2 = chi2 / n
    r, c = table.shape
    if bias_correction:
        phi2 = max(0.0, phi2 - (c - 1) * (r - 1) / (n - 1))
        r = r - (r - 1) ** 2 / (n - 1)
        c = c - (c - 1) ** 2 / (n - 1)
    denom = min(r - 1, c - 1)
    if denom <= 0:
        return 0.0
    return float(np.sqrt(phi2 / denom))


def format_validation_report(result: ValidationGatesResult) -> str:
    """Format validation results for display in report."""
    lines = [
        "## Validation Gates",
        "",
        f"**Overall Status:** {'PASS' if result.all_passed else 'FAIL'}",
        "",
        f"{result.summary}",
        "",
        "| Test | Status | Metric | Value | Threshold |",
        "|------|--------|--------|-------|-----------|",
    ]

    for test in result.results:
        status_icon = "✓" if test.passed else "✗"
        lines.append(
            f"| {test.name} | {status_icon} {test.status} | {test.metric_name} | "
            f"{test.metric_value:.3f} | {test.comparison} {test.threshold} |"
        )

    lines.extend(
        [
            "",
            "### Interpretation",
            "",
            "- **Negative Control 1 (Label Shuffle)**: Tests if clustering can recover "
            "randomly shuffled labels. PASS means clustering does NOT find spurious patterns.",
            "",
            "- **Negative Control 2 (Random Gene Sets)**: Tests if clusters are driven by "
            "biological pathways vs. random gene groupings. PASS means biological pathways matter.",
            "",
            "- **Stability Test (Bootstrap)**: Tests if clusters are robust to resampling. "
            "PASS means clusters are stable features of the data.",
            "",
            "- **Cross-Modal Concordance (Gate 5)**: Tests if subtypes discovered from "
            "fused data also appear when clustering each modality independently. "
            "PASS means subtypes reflect shared biology, not modality-specific artifacts.",
            "",
        ]
    )

    return "\n".join(lines)
