"""
Demo Pipeline Orchestrator

Provides the main pipeline for running the pathway subtyping framework
from data loading through subtype clustering and output generation.
"""

import hashlib
import json
import logging
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import yaml

from .data_quality import (
    DataQualityReport,
    VCFDataQualityError,
    load_vcf_with_quality_check,
)
from .utils.seed import get_rng, set_global_seed
from .characterization import (
    CharacterizationResult,
    characterize_subtypes,
    export_characterization,
    generate_subtype_heatmap,
)
from .validation import ValidationGates, ValidationGatesResult

# Configure logging
logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Configuration for the demo pipeline."""

    name: str = "demo_run"
    output_dir: str = "outputs/demo_run"
    seed: Optional[int] = 42
    verbose: bool = True

    # Input paths
    vcf_path: str = ""
    phenotype_path: str = ""
    pathway_db: str = ""

    # Clustering
    n_clusters: Optional[int] = None
    n_clusters_range: List[int] = field(default_factory=lambda: [2, 8])

    # Input type
    input_type: str = "vcf"  # "vcf" or "expression"

    # Expression-specific config
    expression_path: str = ""
    expression_input_type: str = "tpm"  # "counts", "tpm", "fpkm", "log2"
    expression_scoring_method: str = "ssgsea"  # "mean_z", "ssgsea", "gsva"
    ssgsea_alpha: float = 0.25

    # Ancestry correction (optional)
    ancestry_pcs_path: Optional[str] = None
    ancestry_correction: Optional[str] = None  # "regress_out", "covariate_aware"
    ancestry_n_pcs: int = 10

    # Validation gates
    validation_run_gates: bool = True
    validation_n_permutations: int = 100
    validation_n_bootstrap: int = 50
    validation_stability_threshold: Optional[float] = None  # None = auto-calibrate
    validation_null_ari_max: Optional[float] = None  # None = auto-calibrate
    validation_calibrate: bool = True
    validation_alpha: float = 0.05

    # Output settings
    disclaimer: str = "Research use only. Not medical advice."

    @classmethod
    def from_yaml(cls, yaml_path: str) -> "PipelineConfig":
        """Load configuration from YAML file."""
        with open(yaml_path, "r") as f:
            config_dict = yaml.safe_load(f)

        # Flatten nested config
        pipeline = config_dict.get("pipeline", {})
        data = config_dict.get("data", {})
        clustering = config_dict.get("clustering", {})
        output = config_dict.get("output", {})
        ancestry = config_dict.get("ancestry", {})
        validation = config_dict.get("validation", {})

        return cls(
            name=pipeline.get("name", "demo_run"),
            output_dir=pipeline.get("output_dir", "outputs/demo_run"),
            seed=pipeline.get("seed", 42),
            verbose=pipeline.get("verbose", True),
            vcf_path=data.get("vcf_path", ""),
            phenotype_path=data.get("phenotype_path", ""),
            pathway_db=data.get("pathway_db", ""),
            n_clusters=clustering.get("n_clusters"),
            n_clusters_range=clustering.get("n_clusters_range", [2, 8]),
            input_type=data.get("input_type", "vcf"),
            expression_path=data.get("expression_path", ""),
            expression_input_type=data.get("expression_input_type", "tpm"),
            expression_scoring_method=data.get(
                "expression_scoring_method", "ssgsea"
            ),
            ssgsea_alpha=float(data.get("ssgsea_alpha", 0.25)),
            ancestry_pcs_path=ancestry.get("pcs_path"),
            ancestry_correction=ancestry.get("correction"),
            ancestry_n_pcs=ancestry.get("n_pcs", 10),
            validation_run_gates=validation.get("run_gates", True),
            validation_n_permutations=validation.get("n_permutations", 100),
            validation_n_bootstrap=validation.get("n_bootstrap", 50),
            validation_stability_threshold=validation.get("stability_threshold"),
            validation_null_ari_max=validation.get("null_ari_max"),
            validation_calibrate=validation.get("calibrate", True),
            validation_alpha=float(validation.get("alpha", 0.05)),
            disclaimer=output.get("disclaimer", "Research use only."),
        )


class DemoPipeline:
    """
    Main pipeline orchestrator for the pathway subtyping framework.

    This provides a simplified end-to-end pipeline that:
    1. Loads VCF, phenotype, and pathway data
    2. Computes gene burden scores
    3. Aggregates to pathway scores
    4. Clusters samples into subtypes
    5. Generates outputs (tables, figures, reports)
    """

    def __init__(self, config: PipelineConfig):
        """Initialize the pipeline with configuration."""
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.rng = None

        # Data holders
        self.variants_df: Optional[pd.DataFrame] = None
        self.phenotypes_df: Optional[pd.DataFrame] = None
        self.pathways: Dict[str, List[str]] = {}
        self.gene_burdens: Optional[pd.DataFrame] = None
        self.pathway_scores: Optional[pd.DataFrame] = None
        self.cluster_assignments: Optional[pd.DataFrame] = None
        self.n_clusters: int = 0

        # Validation results
        self.validation_result: Optional[ValidationGatesResult] = None

        # Data quality report
        self.data_quality_report: Optional[DataQualityReport] = None

        # Ancestry correction
        self.ancestry_pcs = None
        self.ancestry_adjustment = None
        self.ancestry_report = None

        # Expression data (for expression input mode)
        self.gene_expression: Optional[pd.DataFrame] = None
        self.expression_scoring_result = None

        # Threshold calibration
        self.calibrated_thresholds = None

        # Characterization
        self.characterization_result: Optional[CharacterizationResult] = None

        # Timing
        self.start_time: Optional[datetime] = None
        self.end_time: Optional[datetime] = None

    def setup(self) -> None:
        """Set up output directories and logging."""
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "figures").mkdir(exist_ok=True)

        # Set up logging
        log_file = self.output_dir / "pipeline.log"
        logging.basicConfig(
            level=logging.INFO if self.config.verbose else logging.WARNING,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout),
            ],
        )

        # Set seed for reproducibility
        if self.config.seed is not None:
            set_global_seed(self.config.seed)
            self.rng = get_rng(self.config.seed, "pipeline")
            logger.info(f"Set random seed: {self.config.seed}")

        logger.info(f"Output directory: {self.output_dir}")

    def load_data(self) -> None:
        """Load input data based on input type."""
        logger.info("Loading input data...")

        if self.config.input_type == "expression":
            self._load_expression()
        else:
            self._load_vcf()

        # Load phenotypes
        self._load_phenotypes()

        # Load pathways
        self._load_pathways()

        if self.config.input_type == "expression":
            logger.info(
                f"Loaded: {len(self.gene_expression.columns)} genes, "
                f"{len(self.phenotypes_df)} samples, "
                f"{len(self.pathways)} pathways"
            )
        else:
            logger.info(
                f"Loaded: {len(self.variants_df)} variants, "
                f"{len(self.phenotypes_df)} samples, "
                f"{len(self.pathways)} pathways"
            )

    def _load_vcf(self) -> None:
        """Load and parse VCF file with data quality checking.

        Features (v0.2):
        - Multi-allelic variant support
        - Graceful handling of missing annotations
        - Data quality reporting with warnings
        - User-friendly error messages with fix guidance
        """
        vcf_path = Path(self.config.vcf_path)
        if not vcf_path.exists():
            raise FileNotFoundError(
                f"VCF file not found: {vcf_path}\n\n"
                f"Please check:\n"
                f"  1. The file path is correct in your config\n"
                f"  2. You have read permissions for the file\n"
                f"  3. The file hasn't been moved or deleted\n\n"
                f"Config path: {self.config.vcf_path}"
            )

        try:
            # Use the new data quality module for robust VCF parsing
            self.variants_df, self.genotypes_df, self.samples, self.data_quality_report = (
                load_vcf_with_quality_check(
                    str(vcf_path),
                    strict=False,  # Allow pipeline to continue with warnings
                    min_gene_coverage=50.0,
                )
            )

            # Log data quality summary
            if self.data_quality_report.multi_allelic_variants > 0:
                logger.info(
                    f"Expanded {self.data_quality_report.multi_allelic_variants} "
                    f"multi-allelic variants to "
                    f"{self.data_quality_report.multi_allelic_expanded} bi-allelic records"
                )

            if self.data_quality_report.gene_coverage < 80:
                logger.warning(
                    f"Gene annotation coverage is {self.data_quality_report.gene_coverage:.1f}%. "
                    f"Consider running VEP or ANNOVAR for better results."
                )

            if self.data_quality_report.warnings:
                for warning in self.data_quality_report.warnings:
                    logger.warning(f"Data quality: {warning}")

            # Check if data is usable
            if not self.data_quality_report.is_usable:
                raise VCFDataQualityError(
                    "VCF data quality is insufficient for analysis",
                    self.data_quality_report,
                    [
                        "Annotate variants with gene symbols using VEP or ANNOVAR",
                        "Use the annotation helper: python scripts/annotate_vcf.py",
                        "See troubleshooting: docs/troubleshooting.md#data-issues",
                    ],
                )

        except VCFDataQualityError:
            # Re-raise data quality errors with full context
            raise
        except Exception as e:
            # Wrap other errors with helpful context
            raise RuntimeError(
                f"Failed to parse VCF file: {vcf_path}\n\n"
                f"Error: {e}\n\n"
                f"Common causes:\n"
                f"  1. File is corrupted or truncated\n"
                f"  2. File is not in valid VCF format\n"
                f"  3. Encoding issues (try UTF-8)\n\n"
                f"To validate your VCF:\n"
                f"  bcftools view -h {vcf_path} | head -20\n"
                f'  python -c "from pathway_subtyping.data_quality '
                f"import validate_vcf_for_pipeline; "
                f"validate_vcf_for_pipeline('{vcf_path}')\""
            ) from e

    def _load_expression(self) -> None:
        """Load expression matrix for expression input mode."""
        from .expression import ExpressionInputType, load_expression_matrix

        expr_path = Path(self.config.expression_path)
        if not expr_path.exists():
            raise FileNotFoundError(
                f"Expression file not found: {expr_path}\n\n"
                f"Expected: CSV/TSV with samples as rows and genes as columns\n"
                f"(or genes as rows — auto-detected)\n\n"
                f"Config path: {self.config.expression_path}"
            )

        input_type = ExpressionInputType(self.config.expression_input_type)
        self.gene_expression, self.data_quality_report_expr = (
            load_expression_matrix(str(expr_path), input_type=input_type)
        )
        self.samples = list(self.gene_expression.index)

    def _load_phenotypes(self) -> None:
        """Load phenotype CSV file."""
        pheno_path = Path(self.config.phenotype_path)
        if not pheno_path.exists():
            raise FileNotFoundError(f"Phenotype file not found: {pheno_path}")

        self.phenotypes_df = pd.read_csv(pheno_path)
        self.phenotypes_df.set_index("sample_id", inplace=True)

    def _load_pathways(self) -> None:
        """Load pathway GMT file."""
        gmt_path = Path(self.config.pathway_db)
        if not gmt_path.exists():
            raise FileNotFoundError(f"Pathway file not found: {gmt_path}")

        self.pathways = {}
        with open(gmt_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    pathway_name = parts[0]
                    # parts[1] is description, parts[2:] are genes
                    genes = parts[2:]
                    self.pathways[pathway_name] = genes

    def _load_ancestry_pcs(self) -> None:
        """Load pre-computed ancestry principal components."""
        from .ancestry import AncestryPCs

        pcs_path = Path(self.config.ancestry_pcs_path)
        if not pcs_path.exists():
            raise FileNotFoundError(
                f"Ancestry PCs file not found: {pcs_path}\n\n"
                f"Expected format: CSV with sample_id as index, PC1..PCn as columns.\n"
                f"Generate PCs from genotype data using PLINK:\n"
                f"  plink2 --bfile your_data --pca 10 --out ancestry_pcs"
            )

        pcs_df = pd.read_csv(pcs_path, index_col=0)
        n_components = len(pcs_df.columns)

        self.ancestry_pcs = AncestryPCs(
            components=pcs_df,
            explained_variance_ratio=np.zeros(n_components),
            n_components=n_components,
            sample_ids=list(pcs_df.index),
        )
        logger.info(f"[Ancestry] Loaded {n_components} ancestry PCs for {len(pcs_df)} samples")

    def _adjust_for_ancestry(self) -> None:
        """Apply ancestry correction to pathway scores."""
        from .ancestry import AncestryMethod, adjust_pathway_scores

        method = AncestryMethod(self.config.ancestry_correction)
        self.ancestry_adjustment = adjust_pathway_scores(
            pathway_scores=self.pathway_scores,
            ancestry_pcs=self.ancestry_pcs,
            method=method,
            n_pcs=self.config.ancestry_n_pcs,
        )
        self.pathway_scores = self.ancestry_adjustment.adjusted_scores
        logger.info(
            f"[Ancestry] Applied correction: {method.value} "
            f"({self.ancestry_adjustment.n_pcs_used} PCs, "
            f"{len(self.ancestry_adjustment.highly_confounded_pathways)} confounded pathways)"
        )

    def compute_gene_burdens(self) -> None:
        """Compute gene-level burden scores for each sample."""
        logger.info("Computing gene burden scores...")

        # Get unique genes
        genes = self.variants_df["gene"].unique()

        # Initialize burden matrix
        burden_data = {}

        for gene in genes:
            if not gene:
                continue

            # Get variants in this gene
            gene_variants = self.variants_df[self.variants_df["gene"] == gene]

            for sample in self.samples:
                # Sum weighted burden across variants
                burden = 0.0
                for idx, var in gene_variants.iterrows():
                    gt = self.genotypes_df.loc[idx, sample]
                    if gt > 0:
                        # Weight by consequence and CADD
                        weight = 1.0
                        if "frameshift" in var["consequence"] or "stop" in var["consequence"]:
                            weight = 1.0  # LoF
                        elif "missense" in var["consequence"]:
                            weight = 0.5 if var["cadd"] > 25 else 0.1
                        else:
                            weight = 0.1

                        # Normalize CADD score (cap at 40, handle missing)
                        cadd_score = var["cadd"]
                        if cadd_score <= 0:
                            # Missing CADD: use consequence-based default instead of 0
                            # This prevents silent data loss when CADD is unavailable
                            if "frameshift" in var["consequence"] or "stop" in var["consequence"]:
                                cadd_score = 35.0  # High impact default
                            elif "missense" in var["consequence"]:
                                cadd_score = 20.0  # Moderate impact default
                            else:
                                cadd_score = 10.0  # Low impact default
                            logger.debug(
                                f"Missing CADD for {var['id']}, using default {cadd_score}"
                            )
                        # Cap CADD at 40 to prevent normalization issues
                        cadd_normalized = min(cadd_score, 40.0) / 40.0
                        burden += gt * weight * cadd_normalized

                if gene not in burden_data:
                    burden_data[gene] = {}
                burden_data[gene][sample] = burden

        self.gene_burdens = pd.DataFrame(burden_data).fillna(0)
        logger.info(f"Computed burdens for {len(burden_data)} genes")

    def compute_pathway_scores(self) -> None:
        """Aggregate gene burdens to pathway scores."""
        logger.info("Computing pathway scores...")

        pathway_scores = {}

        for pathway_name, pathway_genes in self.pathways.items():
            # Find genes in this pathway that we have burden data for
            common_genes = [g for g in pathway_genes if g in self.gene_burdens.columns]

            if len(common_genes) < 2:
                continue

            # Aggregate: mean burden across pathway genes
            pathway_scores[pathway_name] = self.gene_burdens[common_genes].mean(axis=1)

        self.pathway_scores = pd.DataFrame(pathway_scores)

        # Check for and handle zero-variance pathways before normalization
        pathway_stds = self.pathway_scores.std()
        zero_variance_pathways = pathway_stds[pathway_stds == 0].index.tolist()

        if zero_variance_pathways:
            logger.warning(
                f"Removing {len(zero_variance_pathways)} zero-variance pathway(s): "
                f"{zero_variance_pathways[:5]}{'...' if len(zero_variance_pathways) > 5 else ''}"
            )
            # Remove zero-variance pathways (they provide no discriminative power)
            self.pathway_scores = self.pathway_scores.drop(columns=zero_variance_pathways)

        if self.pathway_scores.empty or len(self.pathway_scores.columns) < 2:
            raise ValueError(
                f"Insufficient pathways after filtering: "
                f"{len(self.pathway_scores.columns)} remaining. "
                f"Need at least 2 pathways with non-zero variance "
                f"for clustering. Check that your VCF contains "
                f"variants in pathway genes."
            )

        # Z-score normalize with numerical stability (epsilon prevents div-by-zero)
        means = self.pathway_scores.mean()
        stds = self.pathway_scores.std()
        # Add small epsilon only where std is very small (near zero but not exactly)
        stds = stds.replace(0, 1e-10)  # Safety net (should not occur after filtering)
        self.pathway_scores = (self.pathway_scores - means) / stds

        logger.info(
            f"Computed scores for {len(self.pathway_scores.columns)} pathways "
            f"(removed {len(zero_variance_pathways)} zero-variance)"
        )

    def compute_expression_pathway_scores(self) -> None:
        """Compute pathway scores from expression data."""
        from .expression import ExpressionScoringMethod, score_pathways_from_expression

        method = ExpressionScoringMethod(self.config.expression_scoring_method)
        result = score_pathways_from_expression(
            self.gene_expression,
            self.pathways,
            method=method,
            alpha=self.config.ssgsea_alpha,
            seed=self.config.seed,
        )

        self.pathway_scores = result.pathway_scores
        self.expression_scoring_result = result

        # For compatibility with validation gates and characterization,
        # set gene_burdens to the expression matrix
        self.gene_burdens = self.gene_expression

        logger.info(
            f"[Expression] Scored {result.n_pathways_scored} pathways "
            f"via {method.value} (skipped {result.n_pathways_skipped})"
        )

    def cluster_samples(self) -> None:
        """Cluster samples into subtypes using GMM."""
        logger.info("Clustering samples into subtypes...")

        from sklearn.metrics import silhouette_score
        from sklearn.mixture import GaussianMixture

        X = self.pathway_scores.values

        # Determine optimal number of clusters using BIC
        if self.config.n_clusters is None:
            best_bic = np.inf
            best_k = 2
            for k in range(self.config.n_clusters_range[0], self.config.n_clusters_range[1] + 1):
                gmm = GaussianMixture(
                    n_components=k,
                    covariance_type="full",
                    n_init=10,
                    random_state=self.config.seed,
                    reg_covar=1e-6,  # Regularization for numerical stability
                )
                gmm.fit(X)
                if not gmm.converged_:
                    logger.warning(f"GMM with k={k} did not converge during BIC search")
                    continue  # Skip non-converged models
                bic = gmm.bic(X)
                if bic < best_bic:
                    best_bic = bic
                    best_k = k
            n_clusters = best_k
            logger.info(f"Selected {n_clusters} clusters via BIC")
        else:
            n_clusters = self.config.n_clusters

        self.n_clusters = n_clusters

        # Fit final model
        gmm = GaussianMixture(
            n_components=n_clusters,
            covariance_type="full",
            n_init=10,
            random_state=self.config.seed,
            reg_covar=1e-6,  # Regularization for numerical stability
        )
        gmm.fit(X)

        # Verify convergence - critical for reliable results
        if not gmm.converged_:
            logger.warning(
                f"GMM did not converge after {gmm.n_iter_} iterations. "
                f"Results may be unreliable. Consider: (1) reducing n_clusters, "
                f"(2) increasing max_iter, or (3) checking for data issues."
            )

        labels = gmm.predict(X)
        probs = gmm.predict_proba(X)

        # Create cluster labels based on top pathways
        cluster_labels = self._label_clusters(labels)

        # Build assignments dataframe
        self.cluster_assignments = pd.DataFrame(
            {
                "sample_id": self.pathway_scores.index,
                "cluster_id": labels,
                "cluster_label": [cluster_labels[label] for label in labels],
                "confidence": probs.max(axis=1),
            }
        )

        # Add planted subtype for validation if available
        if "planted_subtype" in self.phenotypes_df.columns:
            self.cluster_assignments["planted_subtype"] = self.cluster_assignments["sample_id"].map(
                self.phenotypes_df["planted_subtype"]
            )

        # Compute silhouette score
        if n_clusters > 1:
            sil_score = silhouette_score(X, labels)
            logger.info(f"Silhouette score: {sil_score:.3f}")

        logger.info(f"Assigned {len(labels)} samples to {n_clusters} clusters")

    def _label_clusters(self, labels: np.ndarray) -> Dict[int, str]:
        """Assign biological labels to clusters based on top pathways."""
        cluster_labels = {}

        for cluster_id in np.unique(labels):
            # Get samples in this cluster
            cluster_mask = labels == cluster_id
            cluster_scores = self.pathway_scores.iloc[cluster_mask].mean()

            # Get top pathway
            top_pathway = cluster_scores.idxmax()

            # Map to biological label (generic mapping)
            top_pathway_upper = top_pathway.upper()
            if "SYNAPTIC" in top_pathway_upper or "GLUTAMAT" in top_pathway_upper:
                label = "synaptic"
            elif "CHROMATIN" in top_pathway_upper or "HISTONE" in top_pathway_upper:
                label = "chromatin"
            elif (
                "ION_CHANNEL" in top_pathway_upper
                or "SODIUM" in top_pathway_upper
                or "POTASSIUM" in top_pathway_upper
            ):
                label = "ion_channel"
            elif "DOPAMINE" in top_pathway_upper:
                label = "dopamine"
            elif "GABA" in top_pathway_upper:
                label = "gaba"
            elif "IMMUNE" in top_pathway_upper or "COMPLEMENT" in top_pathway_upper:
                label = "immune"
            elif "MTOR" in top_pathway_upper or "PI3K" in top_pathway_upper:
                label = "mtor"
            else:
                label = f"subtype_{cluster_id}"

            cluster_labels[cluster_id] = label

        return cluster_labels

    def run_validation_gates(self) -> None:
        """Run validation gates to verify clustering quality."""
        logger.info("Running validation gates...")

        # Determine thresholds: explicit config > auto-calibration > defaults
        stability_threshold = self.config.validation_stability_threshold
        null_ari_max = self.config.validation_null_ari_max

        if (stability_threshold is None or null_ari_max is None) and self.config.validation_calibrate:
            from .threshold_calibration import calibrate_thresholds

            n_samples = len(self.pathway_scores)
            ct = calibrate_thresholds(
                n_samples=n_samples,
                n_clusters=self.n_clusters,
                n_pathways=len(self.pathways),
                alpha=self.config.validation_alpha,
            )
            self.calibrated_thresholds = ct
            logger.info(
                f"[Calibration] Auto-calibrated thresholds for n={n_samples}, "
                f"k={self.n_clusters}: null_ari={ct.null_ari_threshold:.4f}, "
                f"stability={ct.stability_threshold:.4f}"
            )

            if stability_threshold is None:
                stability_threshold = ct.stability_threshold
            if null_ari_max is None:
                null_ari_max = ct.null_ari_threshold
        else:
            # Use explicit values or legacy defaults
            if stability_threshold is None:
                stability_threshold = 0.8
            if null_ari_max is None:
                null_ari_max = 0.15

        validator = ValidationGates(
            seed=self.config.seed,
            n_permutations=self.config.validation_n_permutations,
            n_bootstrap=self.config.validation_n_bootstrap,
            stability_threshold=stability_threshold,
            null_ari_max=null_ari_max,
        )

        cluster_labels = self.cluster_assignments["cluster_id"].values

        self.validation_result = validator.run_all(
            pathway_scores=self.pathway_scores,
            cluster_labels=cluster_labels,
            pathways=self.pathways,
            gene_burdens=self.gene_burdens,
            n_clusters=self.n_clusters,
            gmm_seed=self.config.seed,
            ancestry_pcs=self.ancestry_pcs,
        )

        # Store ancestry report if available
        if self.ancestry_pcs is not None:
            from .ancestry import check_ancestry_independence

            self.ancestry_report = check_ancestry_independence(cluster_labels, self.ancestry_pcs)

        # Log summary
        status = "PASS" if self.validation_result.all_passed else "FAIL"
        logger.info(f"Validation Gates: {status}")
        logger.info(f"  {self.validation_result.summary}")

    def characterize(self) -> None:
        """Run subtype characterization analysis."""
        logger.info("Running subtype characterization...")

        cluster_labels = self.cluster_assignments["cluster_id"].values

        # Build cluster names from label mapping
        cluster_names = {}
        for _, row in self.cluster_assignments.drop_duplicates("cluster_id").iterrows():
            cluster_names[int(row["cluster_id"])] = row["cluster_label"]

        # Confidence scores (if available)
        confidence = None
        if "confidence" in self.cluster_assignments.columns:
            confidence = self.cluster_assignments["confidence"].values

        self.characterization_result = characterize_subtypes(
            pathway_scores=self.pathway_scores,
            cluster_labels=cluster_labels,
            gene_burdens=self.gene_burdens,
            pathways=self.pathways if self.pathways else None,
            cluster_names=cluster_names,
            confidence_scores=confidence,
            seed=self.config.seed,
        )

        # Generate characterization heatmap
        heatmap_path = self.output_dir / "figures" / "subtype_heatmap.png"
        generate_subtype_heatmap(
            self.characterization_result, output_path=str(heatmap_path)
        )

        # Export characterization CSV files
        char_dir = self.output_dir / "characterization"
        export_characterization(self.characterization_result, str(char_dir))

        logger.info("Subtype characterization complete")

    def generate_outputs(self) -> None:
        """Generate all output artifacts."""
        logger.info("Generating outputs...")

        # 1. Pathway scores table
        self._save_pathway_scores()

        # 2. Cluster assignments
        self._save_cluster_assignments()

        # 3. Summary figure
        self._generate_summary_figure()

        # 4. Reports (JSON and Markdown)
        self._generate_reports()

        # 5. Run metadata
        self._save_metadata()

        logger.info(f"All outputs saved to: {self.output_dir}")

    def _save_pathway_scores(self) -> None:
        """Save pathway scores to CSV."""
        output_path = self.output_dir / "pathway_scores.csv"
        self.pathway_scores.to_csv(output_path)
        logger.info(f"Saved pathway scores: {output_path}")

    def _save_cluster_assignments(self) -> None:
        """Save cluster assignments to CSV."""
        output_path = self.output_dir / "subtype_assignments.csv"
        self.cluster_assignments.to_csv(output_path, index=False)
        logger.info(f"Saved cluster assignments: {output_path}")

    def _generate_summary_figure(self) -> None:
        """Generate summary visualization."""
        try:
            # Set non-interactive backend before importing pyplot
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import seaborn as sns
            from sklearn.decomposition import PCA

            # Ensure we have a clean figure state
            plt.close("all")
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))

            # 1. PCA of samples colored by cluster
            pca = PCA(n_components=2, random_state=self.config.seed)
            X_pca = pca.fit_transform(self.pathway_scores.values)

            scatter = axes[0].scatter(
                X_pca[:, 0],
                X_pca[:, 1],
                c=self.cluster_assignments["cluster_id"],
                cmap="Set2",
                alpha=0.7,
                s=50,
            )
            axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
            axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
            axes[0].set_title("Sample Clustering (PCA)")
            plt.colorbar(scatter, ax=axes[0], label="Cluster")

            # 2. Pathway scores heatmap (top pathways)
            top_pathways = self.pathway_scores.var().nlargest(10).index
            cluster_order = self.cluster_assignments.sort_values("cluster_id")["sample_id"]
            heatmap_data = self.pathway_scores.loc[cluster_order, top_pathways]

            sns.heatmap(
                heatmap_data.T,
                cmap="RdBu_r",
                center=0,
                ax=axes[1],
                xticklabels=False,
                yticklabels=True,
            )
            axes[1].set_title("Pathway Scores by Sample")
            axes[1].set_xlabel("Samples (ordered by cluster)")

            # 3. Cluster distribution
            cluster_counts = self.cluster_assignments["cluster_label"].value_counts()
            axes[2].bar(range(len(cluster_counts)), cluster_counts.values, color="steelblue")
            axes[2].set_xticks(range(len(cluster_counts)))
            axes[2].set_xticklabels(cluster_counts.index, rotation=45, ha="right")
            axes[2].set_ylabel("Number of Samples")
            axes[2].set_title("Cluster Distribution")

            plt.tight_layout()
            output_path = self.output_dir / "figures" / "summary.png"
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close()

            logger.info(f"Saved summary figure: {output_path}")

        except ImportError as e:
            logger.warning(f"Could not generate figure (missing dependency): {e}")
        except Exception as e:
            logger.warning(f"Could not generate figure: {e}")

    def _generate_reports(self) -> None:
        """Generate JSON and Markdown reports."""
        # Compute ground truth validation metrics (if planted subtypes available)
        ground_truth_validation = {}
        if "planted_subtype" in self.cluster_assignments.columns:
            from sklearn.metrics import adjusted_rand_score

            ari = adjusted_rand_score(
                self.cluster_assignments["planted_subtype"],
                self.cluster_assignments["cluster_label"],
            )
            ground_truth_validation["adjusted_rand_index"] = round(ari, 4)
            ground_truth_validation["ari_threshold_met"] = ari > 0.7

        # Validation gates results
        validation_gates = {}
        if self.validation_result:
            validation_gates = self.validation_result.to_dict()

        # Data quality results
        data_quality = {}
        if self.data_quality_report:
            data_quality = self.data_quality_report.to_dict()

        # Ancestry results
        ancestry_info = {}
        if self.ancestry_adjustment:
            ancestry_info["adjustment"] = self.ancestry_adjustment.to_dict()
        if self.ancestry_report:
            ancestry_info["independence_test"] = self.ancestry_report.to_dict()

        # Build input-specific report sections
        if self.config.input_type == "expression":
            input_files = {
                "expression": self.config.expression_path,
                "phenotypes": self.config.phenotype_path,
                "pathways": self.config.pathway_db,
            }
            summary = {
                "input_type": "expression",
                "scoring_method": self.config.expression_scoring_method,
                "n_samples": len(self.samples),
                "n_genes": (
                    len(self.gene_expression.columns)
                    if self.gene_expression is not None
                    else 0
                ),
                "n_pathways": len(self.pathways),
                "n_clusters": self.cluster_assignments["cluster_id"].nunique(),
            }
        else:
            input_files = {
                "vcf": self.config.vcf_path,
                "phenotypes": self.config.phenotype_path,
                "pathways": self.config.pathway_db,
            }
            summary = {
                "input_type": "vcf",
                "n_variants": len(self.variants_df),
                "n_samples": len(self.samples),
                "n_genes": (
                    len(self.gene_burdens.columns)
                    if self.gene_burdens is not None
                    else 0
                ),
                "n_pathways": len(self.pathways),
                "n_clusters": self.cluster_assignments["cluster_id"].nunique(),
            }

        # Threshold calibration results
        calibration_info = {}
        if self.calibrated_thresholds:
            calibration_info = self.calibrated_thresholds.to_dict()

        # JSON report
        report = {
            "pipeline_name": self.config.name,
            "timestamp": datetime.now().isoformat(),
            "seed": self.config.seed,
            "input_files": input_files,
            "summary": summary,
            "clusters": self.cluster_assignments["cluster_label"].value_counts().to_dict(),
            "data_quality": data_quality,
            "ancestry_correction": ancestry_info,
            "ground_truth_validation": ground_truth_validation,
            "validation_gates": validation_gates,
            "threshold_calibration": calibration_info,
            "characterization": (
                self.characterization_result.to_dict()
                if self.characterization_result
                else {}
            ),
            "disclaimer": self.config.disclaimer,
        }

        json_path = self.output_dir / "report.json"
        with open(json_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"Saved JSON report: {json_path}")

        # Markdown report
        md_content = self._generate_markdown_report(report)
        md_path = self.output_dir / "report.md"
        with open(md_path, "w") as f:
            f.write(md_content)
        logger.info(f"Saved Markdown report: {md_path}")

    def _generate_markdown_report(self, report: Dict[str, Any]) -> str:
        """Generate human-readable Markdown report."""
        lines = [
            "# Pathway Subtyping Framework - Analysis Report",
            "",
            f"**Pipeline:** {report['pipeline_name']}",
            f"**Date:** {report['timestamp']}",
            f"**Seed:** {report['seed']}",
            "",
            "---",
            "",
            "## Input Summary",
            "",
            "| Metric | Value |",
            "|--------|-------|",
        ]

        if report["summary"].get("input_type") == "expression":
            lines.extend([
                f"| Input Type | expression ({report['summary'].get('scoring_method', 'ssgsea')}) |",
                f"| Samples | {report['summary']['n_samples']} |",
                f"| Genes | {report['summary']['n_genes']} |",
                f"| Pathways | {report['summary']['n_pathways']} |",
            ])
        else:
            lines.extend([
                f"| Variants | {report['summary'].get('n_variants', 'N/A')} |",
                f"| Samples | {report['summary']['n_samples']} |",
                f"| Genes | {report['summary']['n_genes']} |",
                f"| Pathways | {report['summary']['n_pathways']} |",
            ])

        lines.extend([
            "",
            "## Clustering Results",
            "",
            f"**Number of clusters:** {report['summary']['n_clusters']}",
            "",
            "| Cluster | Samples |",
            "|---------|---------|",
        ])

        for cluster, count in report["clusters"].items():
            lines.append(f"| {cluster} | {count} |")

        # Data Quality section (v0.2)
        data_quality = report.get("data_quality", {})
        if data_quality:
            annotation_cov = data_quality.get("annotation_coverage", {})
            multi_allelic = data_quality.get("multi_allelic", {})
            is_usable = data_quality.get("is_usable", True)
            warnings = data_quality.get("warnings", [])

            status_str = "PASS" if is_usable else "FAIL"
            lines.extend(
                [
                    "",
                    "## Data Quality",
                    "",
                    f"**Status:** {status_str}",
                    "",
                    "### Annotation Coverage",
                    "",
                    "| Field | Coverage |",
                    "|-------|----------|",
                    f"| GENE | {annotation_cov.get('gene', 'N/A')} |",
                    f"| CONSEQUENCE | {annotation_cov.get('consequence', 'N/A')} |",
                    f"| CADD | {annotation_cov.get('cadd', 'N/A')} |",
                ]
            )

            if multi_allelic.get("original", 0) > 0:
                lines.extend(
                    [
                        "",
                        f"**Multi-allelic variants:** {multi_allelic.get('original', 0)} "
                        f"(expanded to {multi_allelic.get('expanded', 0)} bi-allelic records)",
                    ]
                )

            if warnings:
                lines.extend(["", "### Warnings", ""])
                for warning in warnings:
                    lines.append(f"- {warning}")

        # Ancestry Correction section
        ancestry = report.get("ancestry_correction", {})
        if ancestry:
            lines.append("")
            if self.ancestry_adjustment:
                lines.append(self.ancestry_adjustment.format_report())
            if self.ancestry_report:
                lines.append("")
                lines.append(self.ancestry_report.format_report())

        # Validation Gates section
        validation_gates = report.get("validation_gates", {})
        if validation_gates:
            all_passed = validation_gates.get("all_passed", False)
            summary = validation_gates.get("summary", "")
            tests = validation_gates.get("tests", [])

            status_str = "PASS" if all_passed else "FAIL"
            lines.extend(
                [
                    "",
                    "## Validation Gates",
                    "",
                    f"**Overall Status:** {status_str}",
                    "",
                    f"{summary}",
                    "",
                    "| Test | Status | Metric | Value | Threshold |",
                    "|------|--------|--------|-------|-----------|",
                ]
            )

            for test in tests:
                status_icon = "PASS" if test["status"] == "PASS" else "FAIL"
                lines.append(
                    f"| {test['name']} | {status_icon} | {test['metric']} | "
                    f"{test['value']:.3f} | {test['comparison']} {test['threshold']} |"
                )

            lines.extend(
                [
                    "",
                    "### Interpretation",
                    "",
                    "- **Negative Control 1 (Label Shuffle)**: Clustering should NOT recover "
                    "randomly shuffled labels. PASS means no spurious patterns.",
                    "- **Negative Control 2 (Random Gene Sets)**: Clusters should be driven by "
                    "biological pathways, not random genes. PASS means pathways matter.",
                    "- **Stability Test (Bootstrap)**: Clusters should be robust to resampling. "
                    "PASS means stable features.",
                ]
            )

        # Threshold calibration section
        calibration = report.get("threshold_calibration", {})
        if calibration:
            cal_method = calibration.get("calibration_method", "unknown")
            cal_interp = calibration.get("interpolated", False)
            method_str = cal_method
            if cal_interp:
                method_str += " (interpolated)"
            lines.extend(
                [
                    "",
                    "### Threshold Calibration",
                    "",
                    f"Thresholds were calibrated using the **{method_str}** method "
                    f"for n={calibration.get('n_samples', '?')}, "
                    f"k={calibration.get('n_clusters', '?')}.",
                    "",
                    f"- Null ARI threshold: {calibration.get('null_ari_threshold', '?')}",
                    f"- Stability threshold: {calibration.get('stability_threshold', '?')}",
                    f"- Significance level: {calibration.get('alpha', 0.05)}",
                ]
            )

        # Subtype characterization
        if self.characterization_result:
            lines.append("")
            lines.append(self.characterization_result.format_report())

        # Ground truth validation (planted subtypes)
        lines.extend(["", "## Ground Truth Validation", ""])

        ground_truth = report.get("ground_truth_validation", {})
        if ground_truth:
            ari = ground_truth.get("adjusted_rand_index", "N/A")
            threshold_met = ground_truth.get("ari_threshold_met", False)
            status = "PASS" if threshold_met else "FAIL"
            lines.extend(
                [
                    "| Metric | Value | Status |",
                    "|--------|-------|--------|",
                    f"| Adjusted Rand Index | {ari} | {status} |",
                    "",
                    "*Threshold: ARI > 0.7 for planted subtype recovery*",
                ]
            )
        else:
            lines.append("*No planted ground truth available for validation*")

        lines.extend(
            [
                "",
                "---",
                "",
                "## Disclaimer",
                "",
                report["disclaimer"],
                "",
                "---",
                "",
                "*Generated by Pathway Subtyping Framework v0.2*",
            ]
        )

        return "\n".join(lines)

    def _save_metadata(self) -> None:
        """Save run metadata for reproducibility."""

        # Compute file hashes
        def file_hash(path: str) -> str:
            if not Path(path).exists():
                return "N/A"
            with open(path, "rb") as f:
                return hashlib.sha256(f.read()).hexdigest()[:16]

        metadata = {
            "run_id": f"{self.config.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            "config_file": "configs/demo.yaml",
            "timestamp": datetime.now().isoformat(),
            "git_commit": self._get_git_commit(),
            "random_seed": self.config.seed,
            "input_type": self.config.input_type,
            "input_files": (
                {
                    "expression": self.config.expression_path,
                    "expression_hash": file_hash(self.config.expression_path),
                    "phenotypes": self.config.phenotype_path,
                    "phenotypes_hash": file_hash(self.config.phenotype_path),
                    "pathways": self.config.pathway_db,
                    "pathways_hash": file_hash(self.config.pathway_db),
                }
                if self.config.input_type == "expression"
                else {
                    "vcf": self.config.vcf_path,
                    "vcf_hash": file_hash(self.config.vcf_path),
                    "phenotypes": self.config.phenotype_path,
                    "phenotypes_hash": file_hash(self.config.phenotype_path),
                    "pathways": self.config.pathway_db,
                    "pathways_hash": file_hash(self.config.pathway_db),
                }
            ),
            "runtime_seconds": (
                (self.end_time - self.start_time).total_seconds()
                if self.end_time and self.start_time
                else None
            ),
        }

        yaml_path = self.output_dir / "run_metadata.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(metadata, f, default_flow_style=False)
        logger.info(f"Saved metadata: {yaml_path}")

    def _get_git_commit(self) -> str:
        """Get current git commit hash."""
        try:
            import subprocess

            result = subprocess.run(
                ["git", "rev-parse", "--short", "HEAD"],
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent.parent,
            )
            return result.stdout.strip() if result.returncode == 0 else "unknown"
        except Exception:
            return "unknown"

    def run(self) -> None:
        """Execute the full pipeline."""
        self.start_time = datetime.now()

        logger.info("=" * 60)
        logger.info("Pathway Subtyping Framework - Demo Pipeline")
        logger.info("=" * 60)

        try:
            self.setup()
            self.load_data()

            if self.config.input_type == "expression":
                self.compute_expression_pathway_scores()
            else:
                self.compute_gene_burdens()
                self.compute_pathway_scores()

            # Ancestry correction (optional)
            if self.config.ancestry_pcs_path:
                self._load_ancestry_pcs()
            if self.config.ancestry_correction and self.ancestry_pcs:
                self._adjust_for_ancestry()

            self.cluster_samples()
            self.run_validation_gates()

            # Characterization is best-effort — never block the core pipeline
            try:
                self.characterize()
            except Exception as e:
                logger.warning(
                    f"[Characterization] Skipped due to error: {e}"
                )

            self.generate_outputs()

            self.end_time = datetime.now()
            runtime = (self.end_time - self.start_time).total_seconds()

            logger.info("=" * 60)
            logger.info(f"Pipeline completed successfully in {runtime:.1f} seconds")
            logger.info(f"Outputs: {self.output_dir}")
            logger.info("=" * 60)

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise
