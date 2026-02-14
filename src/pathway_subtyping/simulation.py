"""
Simulation Framework for the Pathway Subtyping Framework.

Generates synthetic data with known ground truth for validation:
- Planted subtype structure with configurable effect sizes
- Variable noise levels
- Type I error rate estimation
- Power analysis

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


@dataclass
class SimulationConfig:
    """
    Configuration for synthetic data generation.

    Attributes:
        n_samples: Number of samples to generate
        n_pathways: Number of pathways
        n_genes_per_pathway: Average genes per pathway
        n_subtypes: Number of planted subtypes
        effect_size: Cohen's d effect size for subtype differences
        noise_level: Standard deviation of background noise
        subtype_proportions: Relative sizes of subtypes (default: equal)
        seed: Random seed for reproducibility
    """

    n_samples: int = 200
    n_pathways: int = 15
    n_genes_per_pathway: int = 20
    n_subtypes: int = 3
    effect_size: float = 1.0  # Cohen's d
    noise_level: float = 1.0
    subtype_proportions: Optional[List[float]] = None
    seed: Optional[int] = 42

    # Ancestry/population stratification (optional)
    n_ancestry_groups: int = 0  # 0 = no ancestry simulation
    ancestry_effect_size: float = 0.5  # Effect of ancestry on pathway scores
    ancestry_confounding: float = 0.0  # Correlation between subtypes and ancestry (0-1)

    def __post_init__(self):
        if self.subtype_proportions is None:
            # Equal proportions by default
            self.subtype_proportions = [1.0 / self.n_subtypes] * self.n_subtypes

        # Normalize proportions
        total = sum(self.subtype_proportions)
        self.subtype_proportions = [p / total for p in self.subtype_proportions]


@dataclass
class SimulatedData:
    """
    Container for simulated data with ground truth.

    Attributes:
        gene_burdens: DataFrame of gene burdens (samples x genes)
        pathway_scores: DataFrame of pathway scores (samples x pathways)
        pathways: Dictionary mapping pathway names to gene lists
        true_labels: Array of true subtype labels
        config: SimulationConfig used to generate data
        subtype_pathway_effects: Which pathways differ between subtypes
    """

    gene_burdens: pd.DataFrame
    pathway_scores: pd.DataFrame
    pathways: Dict[str, List[str]]
    true_labels: np.ndarray
    config: SimulationConfig
    subtype_pathway_effects: Dict[int, List[str]] = field(default_factory=dict)
    ancestry_labels: Optional[np.ndarray] = None
    ancestry_pcs: Optional[pd.DataFrame] = None

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "n_samples": len(self.true_labels),
            "n_pathways": len(self.pathways),
            "n_subtypes": len(np.unique(self.true_labels)),
            "config": {
                "effect_size": self.config.effect_size,
                "noise_level": self.config.noise_level,
                "seed": self.config.seed,
            },
            "subtype_sizes": {
                int(k): int(v) for k, v in zip(*np.unique(self.true_labels, return_counts=True))
            },
        }
        if self.ancestry_labels is not None:
            result["ancestry"] = {
                "n_groups": len(np.unique(self.ancestry_labels)),
                "ancestry_effect_size": self.config.ancestry_effect_size,
                "ancestry_confounding": self.config.ancestry_confounding,
            }
        return result


def generate_synthetic_data(config: SimulationConfig) -> SimulatedData:
    """
    Generate synthetic data with planted subtype structure.

    Creates pathway scores where:
    - A subset of pathways have different means between subtypes
    - Effect size controls the magnitude of differences
    - Noise is added to all pathways

    Args:
        config: SimulationConfig specifying data generation parameters

    Returns:
        SimulatedData with gene burdens, pathway scores, and ground truth
    """
    rng = np.random.RandomState(config.seed)

    # Assign samples to subtypes
    n_per_subtype = []
    remaining = config.n_samples

    for i, prop in enumerate(config.subtype_proportions[:-1]):
        n_i = int(round(config.n_samples * prop))
        n_per_subtype.append(n_i)
        remaining -= n_i
    n_per_subtype.append(remaining)

    true_labels = np.concatenate([np.full(n, i) for i, n in enumerate(n_per_subtype)])

    # Shuffle to avoid ordering effects
    shuffle_idx = rng.permutation(config.n_samples)
    true_labels = true_labels[shuffle_idx]

    # Generate pathways
    pathways = {}
    all_genes = []

    for i in range(config.n_pathways):
        n_genes = max(5, int(rng.normal(config.n_genes_per_pathway, 5)))
        pathway_genes = [f"GENE_{i}_{j}" for j in range(n_genes)]
        pathways[f"PATHWAY_{i}"] = pathway_genes
        all_genes.extend(pathway_genes)

    all_genes = list(set(all_genes))

    # Generate gene burdens (baseline)
    gene_burdens = pd.DataFrame(
        rng.exponential(scale=0.5, size=(config.n_samples, len(all_genes))),
        columns=all_genes,
        index=[f"SAMPLE_{i}" for i in range(config.n_samples)],
    )

    # Select pathways that will differ between subtypes
    # Each subtype has ~1/3 of pathways elevated
    n_effect_pathways = max(3, config.n_pathways // 3)
    pathway_names = list(pathways.keys())
    subtype_pathway_effects = {}

    for subtype in range(config.n_subtypes):
        # Select random pathways for this subtype
        start_idx = (subtype * n_effect_pathways) % config.n_pathways
        effect_pathways = []

        for i in range(n_effect_pathways):
            idx = (start_idx + i) % config.n_pathways
            effect_pathways.append(pathway_names[idx])

        subtype_pathway_effects[subtype] = effect_pathways

        # Add effect to genes in these pathways for this subtype's samples
        subtype_mask = true_labels == subtype

        for pathway in effect_pathways:
            for gene in pathways[pathway]:
                if gene in gene_burdens.columns:
                    # Add effect (scaled by noise_level to maintain signal-to-noise)
                    effect_magnitude = config.effect_size * config.noise_level
                    gene_burdens.loc[subtype_mask, gene] += rng.normal(
                        effect_magnitude, 0.1, size=np.sum(subtype_mask)
                    )

    # Add noise
    noise = rng.normal(0, config.noise_level, size=gene_burdens.shape)
    gene_burdens = gene_burdens + noise

    # Ensure non-negative
    gene_burdens = gene_burdens.clip(lower=0)

    # Simulate ancestry/population structure (optional)
    ancestry_labels = None
    ancestry_pcs_df = None

    if config.n_ancestry_groups > 0:
        ancestry_labels, ancestry_pcs_df = _simulate_ancestry(
            rng, config, gene_burdens, true_labels
        )

    # Compute pathway scores
    pathway_scores = pd.DataFrame(index=gene_burdens.index)

    for pathway_name, pathway_genes in pathways.items():
        common_genes = [g for g in pathway_genes if g in gene_burdens.columns]
        if common_genes:
            pathway_scores[pathway_name] = gene_burdens[common_genes].mean(axis=1)

    # Z-score normalize
    pathway_scores = (pathway_scores - pathway_scores.mean()) / pathway_scores.std()
    pathway_scores = pathway_scores.fillna(0)

    return SimulatedData(
        gene_burdens=gene_burdens,
        pathway_scores=pathway_scores,
        pathways=pathways,
        true_labels=true_labels,
        config=config,
        subtype_pathway_effects=subtype_pathway_effects,
        ancestry_labels=ancestry_labels,
        ancestry_pcs=ancestry_pcs_df,
    )


@dataclass
class RecoveryResult:
    """
    Result of ground truth recovery analysis.

    Attributes:
        ari: Adjusted Rand Index between predicted and true labels
        nmi: Normalized Mutual Information
        n_clusters_predicted: Number of clusters found
        n_clusters_true: True number of subtypes
        correct_k: Whether correct number of clusters was found
    """

    ari: float
    nmi: float
    n_clusters_predicted: int
    n_clusters_true: int
    correct_k: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ari": round(self.ari, 4),
            "nmi": round(self.nmi, 4),
            "n_clusters_predicted": self.n_clusters_predicted,
            "n_clusters_true": self.n_clusters_true,
            "correct_k": self.correct_k,
        }


def evaluate_recovery(
    predicted_labels: np.ndarray,
    true_labels: np.ndarray,
) -> RecoveryResult:
    """
    Evaluate how well clustering recovered the true subtype structure.

    Args:
        predicted_labels: Predicted cluster assignments
        true_labels: True subtype labels

    Returns:
        RecoveryResult with ARI, NMI, and other metrics
    """
    ari = adjusted_rand_score(true_labels, predicted_labels)
    nmi = normalized_mutual_info_score(true_labels, predicted_labels)

    n_pred = len(np.unique(predicted_labels))
    n_true = len(np.unique(true_labels))

    return RecoveryResult(
        ari=ari,
        nmi=nmi,
        n_clusters_predicted=n_pred,
        n_clusters_true=n_true,
        correct_k=(n_pred == n_true),
    )


def run_gmm_clustering(
    pathway_scores: pd.DataFrame,
    n_clusters: int,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Run GMM clustering on pathway scores.

    Args:
        pathway_scores: DataFrame of pathway scores
        n_clusters: Number of clusters
        seed: Random seed

    Returns:
        Array of cluster assignments
    """
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type="full",
        n_init=10,
        random_state=seed,
        reg_covar=1e-6,
    )

    gmm.fit(pathway_scores.values)
    return gmm.predict(pathway_scores.values)


@dataclass
class TypeIErrorResult:
    """
    Result of Type I error rate estimation.

    Attributes:
        null_ari_mean: Mean ARI under null hypothesis (random data)
        null_ari_std: Standard deviation of null ARI
        type_i_rate: Proportion of null simulations with ARI > threshold
        threshold: ARI threshold used
        n_simulations: Number of null simulations
    """

    null_ari_mean: float
    null_ari_std: float
    type_i_rate: float
    threshold: float
    n_simulations: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "null_ari_mean": round(self.null_ari_mean, 4),
            "null_ari_std": round(self.null_ari_std, 4),
            "type_i_error_rate": round(self.type_i_rate, 4),
            "threshold": self.threshold,
            "n_simulations": self.n_simulations,
        }


def estimate_type_i_error(
    n_samples: int = 100,
    n_pathways: int = 15,
    n_clusters: int = 3,
    ari_threshold: float = 0.15,
    n_simulations: int = 100,
    seed: Optional[int] = None,
) -> TypeIErrorResult:
    """
    Estimate Type I error rate (false positive clustering).

    Generates random data with no true structure and measures how often
    clustering finds spurious structure (ARI > threshold).

    Args:
        n_samples: Number of samples per simulation
        n_pathways: Number of pathways
        n_clusters: Number of clusters to fit
        ari_threshold: ARI threshold for "finding structure"
        n_simulations: Number of null simulations
        seed: Random seed

    Returns:
        TypeIErrorResult with error rate estimate
    """
    rng = np.random.RandomState(seed)
    null_aris = []

    for i in range(n_simulations):
        # Generate random data with no structure
        null_config = SimulationConfig(
            n_samples=n_samples,
            n_pathways=n_pathways,
            n_subtypes=1,  # No true subtypes
            effect_size=0.0,
            noise_level=1.0,
            seed=seed + i if seed else rng.randint(0, 10000),
        )

        null_data = generate_synthetic_data(null_config)

        # Run clustering
        predicted = run_gmm_clustering(
            null_data.pathway_scores,
            n_clusters=n_clusters,
            seed=seed + i if seed else None,
        )

        # Generate random "true" labels for comparison
        random_labels = rng.randint(0, n_clusters, size=n_samples)

        # Compute ARI (should be ~0 for random vs random)
        ari = adjusted_rand_score(random_labels, predicted)
        null_aris.append(ari)

    null_aris = np.array(null_aris)

    # Type I error: proportion with ARI > threshold
    type_i_rate = np.mean(null_aris > ari_threshold)

    return TypeIErrorResult(
        null_ari_mean=np.mean(null_aris),
        null_ari_std=np.std(null_aris),
        type_i_rate=type_i_rate,
        threshold=ari_threshold,
        n_simulations=n_simulations,
    )


@dataclass
class PowerAnalysisResult:
    """
    Result of power analysis.

    Attributes:
        effect_sizes: Effect sizes tested
        recovery_rates: ARI values at each effect size
        power_at_threshold: Power (proportion ARI > threshold) at each effect size
        threshold: ARI threshold for "successful" recovery
        recommended_effect_size: Minimum effect size for 80% power
    """

    effect_sizes: List[float]
    recovery_rates: Dict[float, List[float]]  # effect_size -> list of ARIs
    power_at_threshold: Dict[float, float]
    threshold: float
    recommended_effect_size: Optional[float]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "effect_sizes": self.effect_sizes,
            "power_at_threshold": {str(k): round(v, 3) for k, v in self.power_at_threshold.items()},
            "threshold": self.threshold,
            "recommended_effect_size": self.recommended_effect_size,
        }


def run_power_analysis(
    n_samples: int = 100,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    effect_sizes: Optional[List[float]] = None,
    ari_threshold: float = 0.8,
    n_simulations_per_effect: int = 50,
    seed: Optional[int] = None,
) -> PowerAnalysisResult:
    """
    Run power analysis across effect sizes.

    Determines what effect size is needed for reliable recovery.

    Args:
        n_samples: Number of samples
        n_pathways: Number of pathways
        n_subtypes: Number of subtypes
        effect_sizes: Effect sizes to test (default: 0.1 to 2.0)
        ari_threshold: ARI threshold for "successful" recovery
        n_simulations_per_effect: Simulations per effect size
        seed: Random seed

    Returns:
        PowerAnalysisResult with power curves
    """
    if effect_sizes is None:
        effect_sizes = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]

    recovery_rates = {}
    power_at_threshold = {}

    logger.info(f"Running power analysis for {len(effect_sizes)} effect sizes...")

    for effect in effect_sizes:
        logger.info(f"  Effect size = {effect}")
        aris = []

        for i in range(n_simulations_per_effect):
            config = SimulationConfig(
                n_samples=n_samples,
                n_pathways=n_pathways,
                n_subtypes=n_subtypes,
                effect_size=effect,
                noise_level=1.0,
                seed=seed + i if seed else None,
            )

            sim_data = generate_synthetic_data(config)

            # Cluster
            predicted = run_gmm_clustering(
                sim_data.pathway_scores,
                n_clusters=n_subtypes,
                seed=seed + i if seed else None,
            )

            # Evaluate
            result = evaluate_recovery(predicted, sim_data.true_labels)
            aris.append(result.ari)

        recovery_rates[effect] = aris
        power_at_threshold[effect] = np.mean(np.array(aris) >= ari_threshold)

    # Find recommended effect size (first to achieve 80% power)
    recommended = None
    for effect in effect_sizes:
        if power_at_threshold[effect] >= 0.8:
            recommended = effect
            break

    return PowerAnalysisResult(
        effect_sizes=effect_sizes,
        recovery_rates=recovery_rates,
        power_at_threshold=power_at_threshold,
        threshold=ari_threshold,
        recommended_effect_size=recommended,
    )


@dataclass
class SampleSizeAnalysisResult:
    """
    Result of sample size analysis.

    Attributes:
        sample_sizes: Sample sizes tested
        power_by_size: Power at each sample size
        recommended_n: Minimum n for 80% power
    """

    sample_sizes: List[int]
    power_by_size: Dict[int, float]
    effect_size: float
    threshold: float
    recommended_n: Optional[int]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_sizes": self.sample_sizes,
            "power_by_size": {str(k): round(v, 3) for k, v in self.power_by_size.items()},
            "effect_size": self.effect_size,
            "threshold": self.threshold,
            "recommended_n": self.recommended_n,
        }


def run_sample_size_analysis(
    effect_size: float = 1.0,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    sample_sizes: Optional[List[int]] = None,
    ari_threshold: float = 0.8,
    n_simulations: int = 50,
    seed: Optional[int] = None,
) -> SampleSizeAnalysisResult:
    """
    Analyze power as a function of sample size.

    Args:
        effect_size: Effect size to use
        n_pathways: Number of pathways
        n_subtypes: Number of subtypes
        sample_sizes: Sample sizes to test
        ari_threshold: ARI threshold for success
        n_simulations: Simulations per sample size
        seed: Random seed

    Returns:
        SampleSizeAnalysisResult with power by sample size
    """
    if sample_sizes is None:
        sample_sizes = [30, 50, 75, 100, 150, 200, 300, 500]

    power_by_size = {}

    logger.info(f"Running sample size analysis for {len(sample_sizes)} sizes...")

    for n in sample_sizes:
        logger.info(f"  n = {n}")
        aris = []

        for i in range(n_simulations):
            config = SimulationConfig(
                n_samples=n,
                n_pathways=n_pathways,
                n_subtypes=n_subtypes,
                effect_size=effect_size,
                noise_level=1.0,
                seed=seed + i if seed else None,
            )

            sim_data = generate_synthetic_data(config)

            predicted = run_gmm_clustering(
                sim_data.pathway_scores,
                n_clusters=n_subtypes,
                seed=seed + i if seed else None,
            )

            result = evaluate_recovery(predicted, sim_data.true_labels)
            aris.append(result.ari)

        power_by_size[n] = np.mean(np.array(aris) >= ari_threshold)

    # Find recommended sample size
    recommended = None
    for n in sample_sizes:
        if power_by_size[n] >= 0.8:
            recommended = n
            break

    return SampleSizeAnalysisResult(
        sample_sizes=sample_sizes,
        power_by_size=power_by_size,
        effect_size=effect_size,
        threshold=ari_threshold,
        recommended_n=recommended,
    )


def validate_framework(
    n_samples: int = 200,
    n_pathways: int = 15,
    n_subtypes: int = 3,
    effect_size: float = 1.0,
    n_runs: int = 10,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Run comprehensive framework validation on simulated data.

    Args:
        n_samples: Number of samples
        n_pathways: Number of pathways
        n_subtypes: Number of subtypes
        effect_size: Effect size for planted structure
        n_runs: Number of validation runs
        seed: Random seed

    Returns:
        Dictionary with validation metrics
    """
    logger.info("Running comprehensive framework validation...")

    # Run multiple simulations
    aris = []
    nmis = []
    correct_k = 0

    for i in range(n_runs):
        config = SimulationConfig(
            n_samples=n_samples,
            n_pathways=n_pathways,
            n_subtypes=n_subtypes,
            effect_size=effect_size,
            seed=seed + i if seed else None,
        )

        sim_data = generate_synthetic_data(config)

        predicted = run_gmm_clustering(
            sim_data.pathway_scores,
            n_clusters=n_subtypes,
            seed=seed + i if seed else None,
        )

        result = evaluate_recovery(predicted, sim_data.true_labels)
        aris.append(result.ari)
        nmis.append(result.nmi)
        if result.correct_k:
            correct_k += 1

    # Estimate Type I error
    type_i = estimate_type_i_error(
        n_samples=n_samples,
        n_pathways=n_pathways,
        n_clusters=n_subtypes,
        n_simulations=50,
        seed=seed,
    )

    return {
        "n_runs": n_runs,
        "config": {
            "n_samples": n_samples,
            "n_pathways": n_pathways,
            "n_subtypes": n_subtypes,
            "effect_size": effect_size,
        },
        "recovery": {
            "mean_ari": round(np.mean(aris), 4),
            "std_ari": round(np.std(aris), 4),
            "mean_nmi": round(np.mean(nmis), 4),
            "correct_k_rate": round(correct_k / n_runs, 4),
        },
        "type_i_error": type_i.to_dict(),
    }


def _simulate_ancestry(
    rng: np.random.RandomState,
    config: SimulationConfig,
    gene_burdens: pd.DataFrame,
    true_labels: np.ndarray,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Add simulated population structure to gene burden data.

    Creates ancestry groups and adds ancestry-correlated effects to
    gene burdens. If ancestry_confounding > 0, ancestry groups are
    correlated with true subtypes (creating the confounding problem
    that ancestry correction is designed to address).

    Args:
        rng: Random state for reproducibility.
        config: SimulationConfig with ancestry parameters.
        gene_burdens: DataFrame of gene burdens to modify in-place.
        true_labels: True subtype assignments.

    Returns:
        Tuple of (ancestry_labels, ancestry_pcs DataFrame).
    """
    n_samples = config.n_samples

    # Assign ancestry groups
    if config.ancestry_confounding > 0 and config.n_ancestry_groups >= config.n_subtypes:
        # Confounded: ancestry partially correlated with subtypes
        ancestry_labels = np.zeros(n_samples, dtype=int)
        for i in range(n_samples):
            if rng.random() < config.ancestry_confounding:
                ancestry_labels[i] = int(true_labels[i]) % config.n_ancestry_groups
            else:
                ancestry_labels[i] = rng.randint(0, config.n_ancestry_groups)
    else:
        # Random (non-confounded) ancestry assignment
        ancestry_labels = rng.randint(0, config.n_ancestry_groups, size=n_samples)

    # Add ancestry-specific effects to gene burdens
    n_ancestry_genes = max(5, len(gene_burdens.columns) // 4)

    for group in range(config.n_ancestry_groups):
        mask = ancestry_labels == group
        if not np.any(mask):
            continue

        ancestry_genes = rng.choice(
            gene_burdens.columns,
            size=min(n_ancestry_genes, len(gene_burdens.columns)),
            replace=False,
        )
        for gene in ancestry_genes:
            gene_burdens.loc[mask, gene] += rng.normal(
                config.ancestry_effect_size * (group + 1),
                0.1,
                size=int(np.sum(mask)),
            )

    # Generate simulated ancestry PCs
    n_pcs = min(10, config.n_ancestry_groups * 2)
    pc_values = np.zeros((n_samples, n_pcs))

    for i in range(n_pcs):
        for group in range(config.n_ancestry_groups):
            mask = ancestry_labels == group
            if not np.any(mask):
                continue
            pc_values[mask, i] = rng.normal(
                (group - config.n_ancestry_groups / 2) * (1.0 / (i + 1)),
                0.3,
                size=int(np.sum(mask)),
            )

    ancestry_pcs_df = pd.DataFrame(
        pc_values,
        index=gene_burdens.index,
        columns=[f"PC{j+1}" for j in range(n_pcs)],
    )

    return ancestry_labels, ancestry_pcs_df


# ---------------------------------------------------------------------------
# Expression Data Simulation
# ---------------------------------------------------------------------------


@dataclass
class ExpressionSimulationConfig:
    """
    Configuration for synthetic expression data generation.

    Attributes:
        n_samples: Number of samples to generate.
        n_genes: Total number of genes in expression matrix.
        n_pathways: Number of pathways.
        n_genes_per_pathway: Genes per pathway.
        n_subtypes: Number of planted subtypes.
        effect_size: Log2 fold-change for subtype-specific pathways.
        noise_level: Standard deviation of background noise.
        base_expression_mean: Mean expression in log2-scale.
        base_expression_std: Std of baseline expression.
        dropout_rate: Fraction of zeros (RNA-seq zero-inflation).
        subtype_proportions: Relative sizes of subtypes.
        seed: Random seed for reproducibility.
    """

    n_samples: int = 200
    n_genes: int = 500
    n_pathways: int = 15
    n_genes_per_pathway: int = 20
    n_subtypes: int = 3
    effect_size: float = 1.5
    noise_level: float = 1.0
    base_expression_mean: float = 6.0
    base_expression_std: float = 2.0
    dropout_rate: float = 0.1
    subtype_proportions: Optional[List[float]] = None
    seed: Optional[int] = 42

    def __post_init__(self):
        if self.subtype_proportions is None:
            self.subtype_proportions = [1.0 / self.n_subtypes] * self.n_subtypes
        total = sum(self.subtype_proportions)
        self.subtype_proportions = [p / total for p in self.subtype_proportions]


@dataclass
class SimulatedExpressionData:
    """Container for simulated expression data with ground truth."""

    gene_expression: pd.DataFrame
    pathways: Dict[str, List[str]]
    true_labels: np.ndarray
    config: ExpressionSimulationConfig
    subtype_pathway_effects: Dict[int, List[str]] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_samples": len(self.gene_expression),
            "n_genes": len(self.gene_expression.columns),
            "n_pathways": len(self.pathways),
            "n_subtypes": self.config.n_subtypes,
            "effect_size": self.config.effect_size,
            "subtype_pathway_effects": {
                str(k): v for k, v in self.subtype_pathway_effects.items()
            },
        }


def generate_synthetic_expression_data(
    config: ExpressionSimulationConfig,
) -> SimulatedExpressionData:
    """
    Generate synthetic expression data with planted subtype structure.

    Creates a gene expression matrix where:
    - Baseline expression ~ Normal(base_mean, base_std) in log2 space
    - Subtype-specific pathways have differential expression
    - Effect size controls log2 fold-change magnitude
    - Dropout simulates zero-inflation in RNA-seq

    Args:
        config: Simulation configuration.

    Returns:
        SimulatedExpressionData with expression matrix, pathways,
        true labels, and subtype-pathway effect mapping.
    """
    rng = np.random.RandomState(config.seed)

    # 1. Assign samples to subtypes
    n_per_subtype = []
    remaining = config.n_samples
    for i, prop in enumerate(config.subtype_proportions[:-1]):
        n_i = int(round(config.n_samples * prop))
        n_per_subtype.append(n_i)
        remaining -= n_i
    n_per_subtype.append(remaining)

    true_labels = np.concatenate(
        [np.full(n, i) for i, n in enumerate(n_per_subtype)]
    )
    shuffle_idx = rng.permutation(config.n_samples)
    true_labels = true_labels[shuffle_idx]

    # 2. Generate pathway definitions
    all_genes = [f"GENE_{j}" for j in range(config.n_genes)]
    pathways = {}
    gene_pool = list(range(config.n_genes))
    rng.shuffle(gene_pool)

    for i in range(config.n_pathways):
        n_genes = min(
            config.n_genes_per_pathway,
            len(gene_pool) - i * config.n_genes_per_pathway,
        )
        start = i * config.n_genes_per_pathway
        end = start + n_genes
        if end > len(gene_pool):
            # Wrap around for extra pathways
            idxs = gene_pool[start:] + gene_pool[: end - len(gene_pool)]
        else:
            idxs = gene_pool[start:end]
        pathways[f"PATHWAY_{i}"] = [all_genes[j] for j in idxs]

    # 3. Generate baseline expression (log2-scale)
    sample_ids = [f"SAMPLE_{i:04d}" for i in range(config.n_samples)]
    expression = rng.normal(
        config.base_expression_mean,
        config.base_expression_std,
        (config.n_samples, config.n_genes),
    )

    # 4. Add subtype-specific pathway effects
    n_effect_pathways = max(2, config.n_pathways // config.n_subtypes)
    subtype_effects: Dict[int, List[str]] = {}

    pathway_names = list(pathways.keys())
    for subtype in range(config.n_subtypes):
        start = subtype * n_effect_pathways
        end = min(start + n_effect_pathways, len(pathway_names))
        effect_pw = pathway_names[start:end]
        subtype_effects[subtype] = effect_pw

        subtype_mask = true_labels == subtype
        for pw_name in effect_pw:
            for gene_name in pathways[pw_name]:
                gene_idx = all_genes.index(gene_name)
                expression[subtype_mask, gene_idx] += rng.normal(
                    config.effect_size, 0.2, size=int(np.sum(subtype_mask))
                )

    # 5. Add noise
    expression += rng.normal(0, config.noise_level * 0.3, expression.shape)

    # 6. Apply dropout (zero-inflation)
    if config.dropout_rate > 0:
        dropout_mask = rng.random(expression.shape) < config.dropout_rate
        expression[dropout_mask] = 0.0

    # 7. Clip negatives (log2 values shouldn't be very negative)
    expression = np.clip(expression, 0, None)

    gene_expression = pd.DataFrame(
        expression, index=sample_ids, columns=all_genes
    )

    logger.info(
        f"[Simulation] Generated synthetic expression data: "
        f"{config.n_samples} samples, {config.n_genes} genes, "
        f"{config.n_subtypes} subtypes"
    )

    return SimulatedExpressionData(
        gene_expression=gene_expression,
        pathways=pathways,
        true_labels=true_labels,
        config=config,
        subtype_pathway_effects=subtype_effects,
    )
