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

    return True


def _validate_data_section(data: Dict[str, Any], check_files: bool) -> None:
    """Validate the data section of config."""
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
