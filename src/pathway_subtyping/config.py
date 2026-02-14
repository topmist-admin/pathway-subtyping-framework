"""
Configuration loading and validation.

Provides utility functions for loading pipeline configuration from YAML files.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

logger = logging.getLogger(__name__)


class ConfigValidationError(ValueError):
    """Exception raised for configuration validation errors."""

    def __init__(
        self, message: str, field: Optional[str] = None, suggestions: Optional[List[str]] = None
    ):
        self.field = field
        self.suggestions = suggestions or []
        super().__init__(self._format_message(message))

    def _format_message(self, message: str) -> str:
        lines = [message]
        if self.field:
            lines.append(f"Field: {self.field}")
        if self.suggestions:
            lines.append("\nSuggested fixes:")
            for i, s in enumerate(self.suggestions, 1):
                lines.append(f"  {i}. {s}")
        return "\n".join(lines)


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load pipeline configuration from YAML file.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If YAML is invalid
    """
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(path) as f:
        config = yaml.safe_load(f)

    if config is None:
        raise ConfigValidationError(
            "Config file is empty or invalid",
            suggestions=["Check the file contains valid YAML", "Ensure proper indentation"],
        )

    return config


def validate_config(config: Dict[str, Any], check_files: bool = True) -> bool:
    """
    Validate pipeline configuration.

    Args:
        config: Configuration dictionary
        check_files: If True, verify input files exist

    Returns:
        True if valid

    Raises:
        ConfigValidationError: If configuration is invalid
    """
    # Check required sections
    required_sections = ["pipeline", "data"]
    for section in required_sections:
        if section not in config:
            raise ConfigValidationError(
                f"Missing required config section: {section}",
                field=section,
                suggestions=[f"Add a '{section}:' section to your config file"],
            )

    # Validate data section
    data = config.get("data", {})
    _validate_data_section(data, check_files)

    # Validate pipeline section
    pipeline = config.get("pipeline", {})
    _validate_pipeline_section(pipeline)

    # Validate clustering section (optional)
    clustering = config.get("clustering", {})
    if clustering:
        _validate_clustering_section(clustering)

    # Validate ancestry section (optional)
    ancestry = config.get("ancestry", {})
    if ancestry:
        _validate_ancestry_section(ancestry, check_files)

    # Validate validation section (optional)
    validation = config.get("validation", {})
    if validation:
        _validate_validation_section(validation)

    # Validate variant_qc section (optional)
    variant_qc = config.get("variant_qc", {})
    if variant_qc:
        _validate_variant_qc_section(variant_qc)

    return True


def _validate_data_section(data: Dict[str, Any], check_files: bool) -> None:
    """Validate the data section of config."""
    input_type = data.get("input_type", "vcf")

    if input_type == "expression":
        required_fields = ["expression_path", "phenotype_path", "pathway_db"]
    else:
        required_fields = ["vcf_path", "phenotype_path", "pathway_db"]

    for field in required_fields:
        if not data.get(field):
            raise ConfigValidationError(
                f"Missing required field: data.{field}",
                field=f"data.{field}",
                suggestions=[f"Add '{field}: /path/to/file' under the 'data:' section"],
            )

        if check_files:
            path = Path(data[field])
            if not path.exists():
                raise ConfigValidationError(
                    f"File not found: {data[field]}",
                    field=f"data.{field}",
                    suggestions=[
                        "Check the file path is correct",
                        "Use absolute paths if relative paths don't work",
                        f"Verify the file exists: ls -la {data[field]}",
                    ],
                )

    # Validate expression-specific fields
    if input_type == "expression":
        valid_input_types = ["counts", "tpm", "fpkm", "log2"]
        expr_input_type = data.get("expression_input_type", "tpm")
        if expr_input_type not in valid_input_types:
            raise ConfigValidationError(
                f"Invalid expression_input_type: {expr_input_type}",
                field="data.expression_input_type",
                suggestions=[f"Use one of: {', '.join(valid_input_types)}"],
            )

        valid_methods = ["mean_z", "ssgsea", "gsva"]
        scoring_method = data.get("expression_scoring_method", "ssgsea")
        if scoring_method not in valid_methods:
            raise ConfigValidationError(
                f"Invalid expression_scoring_method: {scoring_method}",
                field="data.expression_scoring_method",
                suggestions=[f"Use one of: {', '.join(valid_methods)}"],
            )


def _validate_pipeline_section(pipeline: Dict[str, Any]) -> None:
    """Validate the pipeline section of config."""
    # Validate seed if provided
    seed = pipeline.get("seed")
    if seed is not None and not isinstance(seed, int):
        raise ConfigValidationError(
            f"Invalid seed value: {seed} (must be an integer)",
            field="pipeline.seed",
            suggestions=["Use an integer value like: seed: 42"],
        )

    # Validate output_dir if provided
    output_dir = pipeline.get("output_dir")
    if output_dir:
        output_path = Path(output_dir)
        parent = output_path.parent
        if not parent.exists():
            logger.warning(f"Output directory parent does not exist: {parent}")


def _validate_clustering_section(clustering: Dict[str, Any]) -> None:
    """Validate the clustering section of config."""
    # Validate n_clusters if provided
    n_clusters = clustering.get("n_clusters")
    if n_clusters is not None:
        if not isinstance(n_clusters, int) or n_clusters < 2:
            raise ConfigValidationError(
                f"Invalid n_clusters value: {n_clusters} (must be integer >= 2)",
                field="clustering.n_clusters",
                suggestions=["Use an integer >= 2, e.g., n_clusters: 4"],
            )

    # Validate n_clusters_range if provided
    n_clusters_range = clustering.get("n_clusters_range")
    if n_clusters_range is not None:
        if not isinstance(n_clusters_range, list) or len(n_clusters_range) != 2:
            raise ConfigValidationError(
                f"Invalid n_clusters_range: {n_clusters_range} (must be [min, max])",
                field="clustering.n_clusters_range",
                suggestions=["Use format: n_clusters_range: [2, 8]"],
            )
        min_k, max_k = n_clusters_range
        if not isinstance(min_k, int) or not isinstance(max_k, int):
            raise ConfigValidationError(
                "n_clusters_range values must be integers",
                field="clustering.n_clusters_range",
            )
        if min_k < 2:
            raise ConfigValidationError(
                f"Minimum clusters must be >= 2, got {min_k}",
                field="clustering.n_clusters_range",
            )
        if max_k < min_k:
            raise ConfigValidationError(
                f"Maximum clusters ({max_k}) must be >= minimum ({min_k})",
                field="clustering.n_clusters_range",
            )


def _validate_ancestry_section(data: Dict[str, Any], check_files: bool) -> None:
    """Validate the ancestry section of config."""
    correction = data.get("correction")
    valid_methods = ["regress_out", "covariate_aware", "stratified"]
    if correction is not None and correction not in valid_methods:
        raise ConfigValidationError(
            f"Invalid ancestry correction method: {correction}",
            field="ancestry.correction",
            suggestions=[f"Use one of: {', '.join(valid_methods)}"],
        )

    pcs_path = data.get("pcs_path")
    if pcs_path and check_files:
        path = Path(pcs_path)
        if not path.exists():
            raise ConfigValidationError(
                f"Ancestry PCs file not found: {pcs_path}",
                field="ancestry.pcs_path",
                suggestions=[
                    "Provide a CSV with sample IDs as index and PC columns",
                    "Generate PCs from genotype data: plink2 --bfile data --pca 10",
                ],
            )

    n_pcs = data.get("n_pcs", 10)
    if not isinstance(n_pcs, int) or n_pcs < 1:
        raise ConfigValidationError(
            f"Invalid n_pcs value: {n_pcs} (must be positive integer)",
            field="ancestry.n_pcs",
            suggestions=["Use a positive integer, typically 10"],
        )


def _validate_validation_section(data: Dict[str, Any]) -> None:
    """Validate the validation section of config."""
    # Validate threshold values
    for field_name in ["stability_threshold", "null_ari_max"]:
        value = data.get(field_name)
        if value is not None:
            if not isinstance(value, (int, float)) or value < 0 or value > 1:
                raise ConfigValidationError(
                    f"Invalid {field_name}: {value} (must be a number in [0, 1])",
                    field=f"validation.{field_name}",
                    suggestions=[f"Use a value between 0 and 1, e.g., {field_name}: 0.8"],
                )

    # Validate alpha
    alpha = data.get("alpha")
    if alpha is not None:
        if not isinstance(alpha, (int, float)) or alpha <= 0 or alpha >= 1:
            raise ConfigValidationError(
                f"Invalid alpha: {alpha} (must be in (0, 1))",
                field="validation.alpha",
                suggestions=["Use a value between 0 and 1 (exclusive), e.g., alpha: 0.05"],
            )

    # Validate n_permutations
    n_perm = data.get("n_permutations")
    if n_perm is not None:
        if not isinstance(n_perm, int) or n_perm < 1:
            raise ConfigValidationError(
                f"Invalid n_permutations: {n_perm} (must be positive integer)",
                field="validation.n_permutations",
                suggestions=["Use a positive integer, e.g., n_permutations: 100"],
            )

    # Validate n_bootstrap
    n_boot = data.get("n_bootstrap")
    if n_boot is not None:
        if not isinstance(n_boot, int) or n_boot < 1:
            raise ConfigValidationError(
                f"Invalid n_bootstrap: {n_boot} (must be positive integer)",
                field="validation.n_bootstrap",
                suggestions=["Use a positive integer, e.g., n_bootstrap: 50"],
            )


def _validate_variant_qc_section(data: Dict[str, Any]) -> None:
    """Validate the variant_qc section of config."""
    # Validate min_qual
    min_qual = data.get("min_qual")
    if min_qual is not None:
        if not isinstance(min_qual, (int, float)) or min_qual < 0:
            raise ConfigValidationError(
                f"Invalid min_qual: {min_qual} (must be a non-negative number)",
                field="variant_qc.min_qual",
                suggestions=["Use a non-negative number, e.g., min_qual: 30.0"],
            )

    # Validate min_call_rate
    min_call_rate = data.get("min_call_rate")
    if min_call_rate is not None:
        if not isinstance(min_call_rate, (int, float)) or min_call_rate < 0 or min_call_rate > 1:
            raise ConfigValidationError(
                f"Invalid min_call_rate: {min_call_rate} (must be in [0, 1])",
                field="variant_qc.min_call_rate",
                suggestions=["Use a value between 0 and 1, e.g., min_call_rate: 0.9"],
            )

    # Validate hwe_p_threshold
    hwe_p = data.get("hwe_p_threshold")
    if hwe_p is not None:
        if not isinstance(hwe_p, (int, float)) or hwe_p < 0 or hwe_p > 1:
            raise ConfigValidationError(
                f"Invalid hwe_p_threshold: {hwe_p} (must be in [0, 1])",
                field="variant_qc.hwe_p_threshold",
                suggestions=["Use a value between 0 and 1, e.g., hwe_p_threshold: 1e-6"],
            )

    # Validate max_maf
    max_maf = data.get("max_maf")
    if max_maf is not None:
        if not isinstance(max_maf, (int, float)) or max_maf < 0 or max_maf > 1:
            raise ConfigValidationError(
                f"Invalid max_maf: {max_maf} (must be in [0, 1])",
                field="variant_qc.max_maf",
                suggestions=["Use a value between 0 and 1, e.g., max_maf: 0.01"],
            )

    # Validate min_gq
    min_gq = data.get("min_gq")
    if min_gq is not None:
        if not isinstance(min_gq, int) or min_gq < 0:
            raise ConfigValidationError(
                f"Invalid min_gq: {min_gq} (must be a non-negative integer)",
                field="variant_qc.min_gq",
                suggestions=["Use a non-negative integer, e.g., min_gq: 20"],
            )

    # Validate min_dp
    min_dp = data.get("min_dp")
    if min_dp is not None:
        if not isinstance(min_dp, int) or min_dp < 0:
            raise ConfigValidationError(
                f"Invalid min_dp: {min_dp} (must be a non-negative integer)",
                field="variant_qc.min_dp",
                suggestions=["Use a non-negative integer, e.g., min_dp: 10"],
            )


def validate_gmt_file(gmt_path: str) -> Dict[str, List[str]]:
    """
    Validate and load a GMT (Gene Matrix Transposed) pathway file.

    Args:
        gmt_path: Path to GMT file

    Returns:
        Dictionary of pathway_name -> gene_list

    Raises:
        FileNotFoundError: If file doesn't exist
        ConfigValidationError: If GMT format is invalid
    """
    path = Path(gmt_path)
    if not path.exists():
        raise FileNotFoundError(f"GMT file not found: {gmt_path}")

    pathways = {}
    errors = []

    with open(path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                errors.append(
                    f"Line {line_num}: Expected at least 3 tab-separated fields, got {len(parts)}"
                )
                continue

            pathway_name = parts[0]
            # parts[1] is description (ignored)
            genes = [g.strip() for g in parts[2:] if g.strip()]

            if not pathway_name:
                errors.append(f"Line {line_num}: Empty pathway name")
                continue

            if len(genes) < 2:
                errors.append(f"Line {line_num}: Pathway '{pathway_name}' has fewer than 2 genes")
                continue

            if pathway_name in pathways:
                errors.append(f"Line {line_num}: Duplicate pathway name '{pathway_name}'")
                continue

            pathways[pathway_name] = genes

    if errors:
        raise ConfigValidationError(
            f"GMT file has {len(errors)} error(s):\n" + "\n".join(errors[:5]),
            field="pathway_db",
            suggestions=[
                "GMT format: PATHWAY_NAME<TAB>DESCRIPTION<TAB>GENE1<TAB>GENE2<TAB>...",
                "Each line must have at least 3 tab-separated fields",
                "Each pathway must have at least 2 genes",
            ],
        )

    if not pathways:
        raise ConfigValidationError(
            "GMT file contains no valid pathways",
            field="pathway_db",
            suggestions=[
                "Check the file is not empty",
                "Ensure proper GMT format with tab separators",
            ],
        )

    logger.info(f"Loaded {len(pathways)} pathways from GMT file")
    return pathways
