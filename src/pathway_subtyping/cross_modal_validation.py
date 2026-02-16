"""
Cross-Modal Validation Gate for the Pathway Subtyping Framework.

Tests whether molecular subtypes are consistent across different data modalities
(e.g., VCF variant burden vs bulk RNA-seq expression). If subtypes discovered
from one data type also appear when clustering a completely different data type,
they likely reflect real biology rather than modality-specific artifacts.

This implements Validation Gate 5 for multi-omic analyses.

Research use only. Not for clinical decision-making.
"""

import itertools
import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.mixture import GaussianMixture
from tqdm import tqdm

from .utils.seed import get_rng

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class CrossModalPairResult:
    """Concordance metrics for a single modality pair."""

    modality_a: str
    modality_b: str
    concordance_ari: float
    concordance_nmi: float
    transfer_ari_a_to_b: float
    transfer_ari_b_to_a: float
    n_shared_samples: int
    n_shared_pathways: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "modality_a": self.modality_a,
            "modality_b": self.modality_b,
            "concordance_ari": round(self.concordance_ari, 4),
            "concordance_nmi": round(self.concordance_nmi, 4),
            "transfer_ari_a_to_b": round(self.transfer_ari_a_to_b, 4),
            "transfer_ari_b_to_a": round(self.transfer_ari_b_to_a, 4),
            "n_shared_samples": self.n_shared_samples,
            "n_shared_pathways": self.n_shared_pathways,
        }


@dataclass
class SingleCellCompositionResult:
    """Result of testing whether bulk subtypes have distinct cell-type compositions."""

    f_statistic: float
    p_value: float
    n_subtypes: int
    n_cell_types: int
    per_cell_type_pvalues: Dict[str, float]
    significant_cell_types: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "f_statistic": round(self.f_statistic, 4),
            "p_value": round(self.p_value, 6),
            "n_subtypes": self.n_subtypes,
            "n_cell_types": self.n_cell_types,
            "per_cell_type_pvalues": {
                k: round(v, 6) for k, v in self.per_cell_type_pvalues.items()
            },
            "significant_cell_types": list(self.significant_cell_types),
        }


@dataclass
class CrossModalValidationResult:
    """Full cross-modal validation result."""

    pair_results: List[CrossModalPairResult]
    mean_concordance_ari: float
    mean_concordance_nmi: float
    mean_transfer_ari: float
    null_ari_95th: float
    gate_passed: bool
    single_cell_composition: Optional[SingleCellCompositionResult] = None

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "pair_results": [p.to_dict() for p in self.pair_results],
            "mean_concordance_ari": round(self.mean_concordance_ari, 4),
            "mean_concordance_nmi": round(self.mean_concordance_nmi, 4),
            "mean_transfer_ari": round(self.mean_transfer_ari, 4),
            "null_ari_95th": round(self.null_ari_95th, 4),
            "gate_passed": self.gate_passed,
        }
        if self.single_cell_composition is not None:
            result["single_cell_composition"] = self.single_cell_composition.to_dict()
        return result

    def format_report(self) -> str:
        """Format a human-readable report."""
        lines = [
            "# Cross-Modal Validation Report",
            "",
            f"**Gate Status:** {'PASS' if self.gate_passed else 'FAIL'}",
            f"**Mean Concordance ARI:** {self.mean_concordance_ari:.4f}",
            f"**Mean Concordance NMI:** {self.mean_concordance_nmi:.4f}",
            f"**Mean Transfer ARI:** {self.mean_transfer_ari:.4f}",
            f"**Null ARI 95th Percentile:** {self.null_ari_95th:.4f}",
            "",
            "## Pairwise Results",
            "",
            "| Modality A | Modality B | Concordance ARI | Concordance NMI "
            "| Transfer A→B | Transfer B→A | Shared Samples |",
            "|------------|------------|-----------------|-----------------|"
            "--------------|--------------|----------------|",
        ]

        for p in self.pair_results:
            lines.append(
                f"| {p.modality_a} | {p.modality_b} | {p.concordance_ari:.4f} "
                f"| {p.concordance_nmi:.4f} | {p.transfer_ari_a_to_b:.4f} "
                f"| {p.transfer_ari_b_to_a:.4f} | {p.n_shared_samples} |"
            )

        if self.single_cell_composition is not None:
            sc = self.single_cell_composition
            lines.extend(
                [
                    "",
                    "## Single-Cell Composition Test",
                    "",
                    f"**F-statistic:** {sc.f_statistic:.4f}",
                    f"**p-value:** {sc.p_value:.6f}",
                    f"**Significant cell types:** {', '.join(sc.significant_cell_types) or 'None'}",
                ]
            )

        lines.extend(
            [
                "",
                "## Interpretation",
                "",
                "Cross-modal concordance tests whether subtypes discovered from fused "
                "multi-omic data replicate when clustering each modality independently. "
                "A concordance ARI significantly above the permutation null indicates "
                "that subtypes reflect shared biology across data types, not "
                "modality-specific artifacts.",
            ]
        )

        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        return [
            "Hubert L, Arabie P. Comparing partitions. " "J Classif. 1985;2(1):193-218.",
            "Strehl A, Ghosh J. Cluster ensembles — a knowledge reuse framework "
            "for combining multiple partitions. JMLR. 2002;3:583-617.",
            "Ritchie MD, et al. Methods of integrating data to uncover "
            "genotype-phenotype interactions. Nat Rev Genet. 2015;16(2):85-97.",
        ]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _transfer_ari(
    scores_train: pd.DataFrame,
    scores_test: pd.DataFrame,
    labels_test: np.ndarray,
    n_clusters: int,
    seed: Optional[int] = None,
) -> float:
    """Train GMM on scores_train, predict on scores_test, compute ARI vs labels_test."""
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type="full",
        n_init=5,
        random_state=seed,
        reg_covar=1e-6,
    )
    gmm.fit(scores_train.values)

    if not gmm.converged_:
        logger.warning("[CrossModal] GMM did not converge during transfer validation")
        return 0.0

    predicted = gmm.predict(scores_test.values)
    return adjusted_rand_score(labels_test, predicted)


def _cluster_modality(
    scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Cluster a single modality's scores with GMM."""
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type="full",
        n_init=5,
        random_state=seed,
        reg_covar=1e-6,
    )
    gmm.fit(scores.values)

    if not gmm.converged_:
        logger.warning("[CrossModal] GMM did not converge during independent clustering")

    return gmm.predict(scores.values)


def _compute_pairwise_concordance(
    scores_a: pd.DataFrame,
    scores_b: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> Tuple[float, float, float, float]:
    """Compute concordance and transfer metrics for one modality pair.

    Returns:
        (concordance_ari, concordance_nmi, transfer_a_to_b, transfer_b_to_a)
    """
    # Independent clustering on each modality's own features
    labels_a = _cluster_modality(scores_a, n_clusters, seed=seed)
    labels_b = _cluster_modality(scores_b, n_clusters, seed=seed)

    concordance_ari = adjusted_rand_score(labels_a, labels_b)
    concordance_nmi = normalized_mutual_info_score(labels_a, labels_b)

    # Transfer validation: requires shared feature space
    shared_pathways = sorted(set(scores_a.columns) & set(scores_b.columns))

    if len(shared_pathways) >= 2:
        aligned_a = scores_a[shared_pathways]
        aligned_b = scores_b[shared_pathways]

        transfer_a_to_b = _transfer_ari(aligned_a, aligned_b, labels_b, n_clusters, seed=seed)
        transfer_b_to_a = _transfer_ari(aligned_b, aligned_a, labels_a, n_clusters, seed=seed)
    else:
        logger.warning(
            "[CrossModal] No shared pathways for transfer validation; " "using concordance only"
        )
        transfer_a_to_b = 0.0
        transfer_b_to_a = 0.0

    return concordance_ari, concordance_nmi, transfer_a_to_b, transfer_b_to_a


def _compute_null_distribution(
    per_modality_scores: Dict[str, pd.DataFrame],
    n_clusters: int,
    n_permutations: int,
    rng: np.random.RandomState,
    show_progress: bool = True,
) -> np.ndarray:
    """Compute null ARI distribution by permuting cluster labels.

    For each permutation: cluster each modality independently, then shuffle
    one set of labels and compute ARI. This establishes the baseline
    concordance expected by chance.
    """
    labels = list(per_modality_scores.keys())
    if len(labels) < 2:
        return np.array([0.0])

    # Use first two modalities for null calibration
    key_a, key_b = labels[0], labels[1]
    scores_a = per_modality_scores[key_a]
    scores_b = per_modality_scores[key_b]

    # Find shared samples
    shared_samples = sorted(set(scores_a.index) & set(scores_b.index))
    if len(shared_samples) < n_clusters * 3:
        return np.array([0.0])

    aligned_a = scores_a.loc[shared_samples]
    aligned_b = scores_b.loc[shared_samples]

    null_aris = []

    for i in tqdm(
        range(n_permutations),
        desc="Null calibration",
        disable=not show_progress,
    ):
        # Cluster each modality independently
        seed_i = int(rng.randint(0, 2**31))
        labels_a = _cluster_modality(aligned_a, n_clusters, seed=seed_i)
        labels_b = _cluster_modality(aligned_b, n_clusters, seed=seed_i + 1)

        # Permute one set of labels
        shuffled_b = rng.permutation(labels_b)
        ari = adjusted_rand_score(labels_a, shuffled_b)
        null_aris.append(ari)

    return np.array(null_aris)


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def cross_modal_concordance(
    per_modality_scores: Dict[str, pd.DataFrame],
    cluster_labels: np.ndarray,
    fused_sample_ids: List[str],
    n_clusters: int,
    n_permutations: int = 100,
    concordance_threshold: Optional[float] = None,
    seed: Optional[int] = None,
    show_progress: bool = True,
) -> CrossModalValidationResult:
    """Test whether subtypes are consistent across data modalities.

    For each pair of modalities:
    1. Find shared samples
    2. Cluster each modality independently (same k)
    3. Compute concordance ARI/NMI between independent clusterings
    4. Compute bidirectional transfer ARI (train on A, predict B and vice versa)

    The gate passes if mean concordance ARI exceeds the 95th percentile of
    a permutation null distribution.

    Args:
        per_modality_scores: Dict mapping modality label to pathway score DataFrame.
        cluster_labels: Labels from fused clustering (used for reference).
        fused_sample_ids: Sample IDs from the fused analysis.
        n_clusters: Number of clusters used.
        n_permutations: Number of permutations for null calibration.
        concordance_threshold: Override for gate threshold. If None,
            uses 95th percentile of permutation null.
        seed: Random seed for reproducibility.
        show_progress: Show tqdm progress bars.

    Returns:
        CrossModalValidationResult with concordance metrics and gate decision.

    Raises:
        ValueError: If fewer than 2 modalities provided.
    """
    if len(per_modality_scores) < 2:
        raise ValueError(
            "[CrossModal] Need at least 2 modalities for cross-modal validation, "
            f"got {len(per_modality_scores)}"
        )

    rng = get_rng(seed, "cross_modal_validation")
    modality_labels = list(per_modality_scores.keys())

    logger.info(
        f"[CrossModal] Running cross-modal concordance for "
        f"{len(modality_labels)} modalities: {modality_labels}"
    )

    # Compute pairwise concordance
    pair_results: List[CrossModalPairResult] = []

    for label_a, label_b in itertools.combinations(modality_labels, 2):
        scores_a = per_modality_scores[label_a]
        scores_b = per_modality_scores[label_b]

        shared_samples = sorted(set(scores_a.index) & set(scores_b.index))
        shared_pathways = sorted(set(scores_a.columns) & set(scores_b.columns))
        n_shared = len(shared_samples)

        if n_shared < n_clusters * 3:
            logger.warning(
                f"[CrossModal] Skipping {label_a} vs {label_b}: only "
                f"{n_shared} shared samples (need >= {n_clusters * 3})"
            )
            continue

        aligned_a = scores_a.loc[shared_samples]
        aligned_b = scores_b.loc[shared_samples]

        gmm_seed = int(rng.randint(0, 2**31))
        ari, nmi, transfer_ab, transfer_ba = _compute_pairwise_concordance(
            aligned_a, aligned_b, n_clusters, seed=gmm_seed
        )

        pair_results.append(
            CrossModalPairResult(
                modality_a=label_a,
                modality_b=label_b,
                concordance_ari=ari,
                concordance_nmi=nmi,
                transfer_ari_a_to_b=transfer_ab,
                transfer_ari_b_to_a=transfer_ba,
                n_shared_samples=n_shared,
                n_shared_pathways=len(shared_pathways),
            )
        )

        logger.info(
            f"[CrossModal] {label_a} vs {label_b}: ARI={ari:.4f}, "
            f"NMI={nmi:.4f}, transfer={transfer_ab:.4f}/{transfer_ba:.4f}"
        )

    # Compute null distribution
    null_aris = _compute_null_distribution(
        per_modality_scores, n_clusters, n_permutations, rng, show_progress
    )
    null_95th = float(np.percentile(null_aris, 95)) if len(null_aris) > 0 else 0.0

    # Compute means
    if pair_results:
        mean_ari = float(np.mean([p.concordance_ari for p in pair_results]))
        mean_nmi = float(np.mean([p.concordance_nmi for p in pair_results]))
        transfer_values = []
        for p in pair_results:
            transfer_values.append(p.transfer_ari_a_to_b)
            transfer_values.append(p.transfer_ari_b_to_a)
        mean_transfer = float(np.mean(transfer_values))
    else:
        mean_ari = 0.0
        mean_nmi = 0.0
        mean_transfer = 0.0

    # Gate decision
    threshold = concordance_threshold if concordance_threshold is not None else null_95th
    gate_passed = mean_ari > threshold

    logger.info(
        f"[CrossModal] Gate {'PASSED' if gate_passed else 'FAILED'}: "
        f"mean_ARI={mean_ari:.4f} vs threshold={threshold:.4f}"
    )

    return CrossModalValidationResult(
        pair_results=pair_results,
        mean_concordance_ari=mean_ari,
        mean_concordance_nmi=mean_nmi,
        mean_transfer_ari=mean_transfer,
        null_ari_95th=null_95th,
        gate_passed=gate_passed,
    )


def single_cell_composition_test(
    bulk_labels: np.ndarray,
    cell_type_fractions: pd.DataFrame,
    significance_threshold: float = 0.05,
) -> SingleCellCompositionResult:
    """Test whether bulk-defined subtypes have distinct cell-type compositions.

    Uses one-way ANOVA per cell type across subtypes with Bonferroni correction.

    Args:
        bulk_labels: Array of subtype labels (one per sample).
        cell_type_fractions: DataFrame (samples x cell_types) with proportions.
        significance_threshold: Significance level (default: 0.05).

    Returns:
        SingleCellCompositionResult with per-cell-type p-values.

    Raises:
        ValueError: If inputs have incompatible shapes.
    """
    if len(bulk_labels) != len(cell_type_fractions):
        raise ValueError(
            f"[CrossModal] bulk_labels length ({len(bulk_labels)}) does not match "
            f"cell_type_fractions rows ({len(cell_type_fractions)})"
        )

    unique_labels = np.unique(bulk_labels)
    n_subtypes = len(unique_labels)
    cell_types = list(cell_type_fractions.columns)
    n_cell_types = len(cell_types)
    bonferroni_threshold = (
        significance_threshold / n_cell_types if n_cell_types > 0 else significance_threshold
    )

    per_cell_type_pvalues: Dict[str, float] = {}
    significant_cell_types: List[str] = []
    f_statistics: List[float] = []

    for ct in cell_types:
        fractions = cell_type_fractions[ct].values
        groups = [fractions[bulk_labels == label] for label in unique_labels]

        # Need at least 2 groups with data
        groups = [g for g in groups if len(g) > 0]
        if len(groups) < 2:
            per_cell_type_pvalues[ct] = 1.0
            continue

        f_stat, p_val = stats.f_oneway(*groups)

        # Handle NaN from constant data
        if np.isnan(f_stat):
            f_stat = 0.0
            p_val = 1.0

        per_cell_type_pvalues[ct] = p_val
        f_statistics.append(f_stat)

        if p_val < bonferroni_threshold:
            significant_cell_types.append(ct)

    # Overall F-statistic: mean of per-cell-type F-stats
    overall_f = float(np.mean(f_statistics)) if f_statistics else 0.0
    # Overall p-value: minimum across cell types
    overall_p = float(min(per_cell_type_pvalues.values())) if per_cell_type_pvalues else 1.0

    logger.info(
        f"[CrossModal] Single-cell composition test: "
        f"{len(significant_cell_types)}/{n_cell_types} cell types significant"
    )

    return SingleCellCompositionResult(
        f_statistic=overall_f,
        p_value=overall_p,
        n_subtypes=n_subtypes,
        n_cell_types=n_cell_types,
        per_cell_type_pvalues=per_cell_type_pvalues,
        significant_cell_types=significant_cell_types,
    )


def generate_synthetic_multimodal_data(
    n_samples: int = 200,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    n_modalities: int = 2,
    shared_signal_strength: float = 1.0,
    modality_noise: float = 0.5,
    seed: Optional[int] = None,
) -> Tuple[Dict[str, pd.DataFrame], np.ndarray]:
    """Generate synthetic multi-modal data with planted shared subtype structure.

    Each modality observes the same underlying subtypes but with independent
    noise. Useful for testing cross-modal validation and threshold calibration.

    Args:
        n_samples: Number of samples.
        n_pathways: Number of pathways per modality.
        n_subtypes: Number of planted subtypes.
        n_modalities: Number of modalities to generate.
        shared_signal_strength: Strength of shared subtype signal.
        modality_noise: Standard deviation of per-modality noise.
        seed: Random seed.

    Returns:
        Tuple of (per_modality_scores dict, true_labels array).
    """
    rng = get_rng(seed, "cross_modal_synthetic")

    # Assign samples to subtypes
    true_labels = np.repeat(np.arange(n_subtypes), n_samples // n_subtypes)
    # Handle remainder
    remainder = n_samples - len(true_labels)
    if remainder > 0:
        true_labels = np.concatenate([true_labels, rng.randint(0, n_subtypes, remainder)])

    sample_ids = [f"SAMPLE_{i:04d}" for i in range(n_samples)]

    # Generate shared subtype centroids
    centroids = rng.normal(0, shared_signal_strength, (n_subtypes, n_pathways))

    per_modality_scores: Dict[str, pd.DataFrame] = {}

    for m in range(n_modalities):
        label = f"modality_{m + 1}"
        pathway_names = [f"PATHWAY_{p + 1}" for p in range(n_pathways)]

        # Build data: centroid[subtype] + per-modality noise
        data = np.zeros((n_samples, n_pathways))
        for i in range(n_samples):
            data[i] = centroids[true_labels[i]] + rng.normal(0, modality_noise, n_pathways)

        df = pd.DataFrame(data, index=sample_ids, columns=pathway_names)
        per_modality_scores[label] = df

    logger.info(
        f"[CrossModal] Generated synthetic data: {n_modalities} modalities, "
        f"{n_samples} samples, {n_pathways} pathways, {n_subtypes} subtypes"
    )

    return per_modality_scores, true_labels
