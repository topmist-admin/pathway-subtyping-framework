"""
Public Dataset Validation Module for the Pathway Subtyping Framework.

Downloads and integrates open-access genomic datasets to validate
pathway definitions and framework behavior against real biological data.

Data sources:
- ClinVar gene-specific summary (NCBI FTP, no auth required)
- Reactome pathway definitions (reactome.org, no auth required)

References:
- Landrum MJ et al. ClinVar: improvements to accessing data.
  Nucleic Acids Res. 2020;48(D1):D835-D844.
- Gillespie M et al. The Reactome Pathway Knowledgebase 2022.
  Nucleic Acids Res. 2022;50(D1):D649-D653.

Research use only. Not for clinical decision-making.
"""

import io
import logging
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .characterization import characterize_subtypes
from .clustering import run_clustering
from .simulation import SimulatedData, SimulationConfig

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================


@dataclass
class DatasetInfo:
    """Metadata about a downloadable public dataset."""

    name: str
    url: str
    description: str
    file_format: str  # "tsv", "gmt.zip"
    license_note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "url": self.url,
            "description": self.description,
            "file_format": self.file_format,
            "license_note": self.license_note,
        }


@dataclass
class ClinVarGeneSummary:
    """
    Parsed ClinVar gene-level data for one gene.

    Attributes:
        symbol: HGNC gene symbol
        gene_id: NCBI Gene ID
        n_pathogenic: Count of pathogenic submissions
        n_likely_pathogenic: Count of likely pathogenic submissions
        n_uncertain: Count of VUS submissions
        n_benign: Count of benign submissions
        n_likely_benign: Count of likely benign submissions
    """

    symbol: str
    gene_id: int
    n_pathogenic: int
    n_likely_pathogenic: int
    n_uncertain: int
    n_benign: int
    n_likely_benign: int

    @property
    def total_pathogenic(self) -> int:
        """Total pathogenic + likely pathogenic submissions."""
        return self.n_pathogenic + self.n_likely_pathogenic

    def to_dict(self) -> Dict[str, Any]:
        return {
            "symbol": str(self.symbol),
            "gene_id": int(self.gene_id),
            "n_pathogenic": int(self.n_pathogenic),
            "n_likely_pathogenic": int(self.n_likely_pathogenic),
            "n_uncertain": int(self.n_uncertain),
            "n_benign": int(self.n_benign),
            "n_likely_benign": int(self.n_likely_benign),
            "total_pathogenic": int(self.total_pathogenic),
        }


@dataclass
class PathwayCoverageResult:
    """
    Result from validating pathway gene coverage against ClinVar.

    Attributes:
        pathway_name: Name of the curated pathway
        total_genes: Total genes in the pathway
        genes_in_clinvar: How many pathway genes appear in ClinVar
        genes_with_pathogenic: How many have pathogenic variant submissions
        coverage_fraction: genes_in_clinvar / total_genes
        pathogenic_fraction: genes_with_pathogenic / total_genes
        gene_details: Per-gene detail list
    """

    pathway_name: str
    total_genes: int
    genes_in_clinvar: int
    genes_with_pathogenic: int
    coverage_fraction: float
    pathogenic_fraction: float
    gene_details: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pathway_name": str(self.pathway_name),
            "total_genes": int(self.total_genes),
            "genes_in_clinvar": int(self.genes_in_clinvar),
            "genes_with_pathogenic": int(self.genes_with_pathogenic),
            "coverage_fraction": float(self.coverage_fraction),
            "pathogenic_fraction": float(self.pathogenic_fraction),
            "gene_details": self.gene_details,
        }


@dataclass
class BiologicalPlausibilityResult:
    """
    Result from biological plausibility checks on clustering output.

    Attributes:
        disease_name: Disease being validated
        n_subtypes: Number of discovered subtypes
        n_enriched_pathways: Number of significantly enriched pathways
        pathway_gene_clinvar_overlap: Fraction of enriched pathway genes in ClinVar
        subtypes_biologically_distinct: Whether subtypes have distinct pathway profiles
        details: Additional detail dict
    """

    disease_name: str
    n_subtypes: int
    n_enriched_pathways: int
    pathway_gene_clinvar_overlap: float
    subtypes_biologically_distinct: bool
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "disease_name": str(self.disease_name),
            "n_subtypes": int(self.n_subtypes),
            "n_enriched_pathways": int(self.n_enriched_pathways),
            "pathway_gene_clinvar_overlap": float(self.pathway_gene_clinvar_overlap),
            "subtypes_biologically_distinct": bool(self.subtypes_biologically_distinct),
            "details": self.details,
        }


@dataclass
class ValidationReport:
    """
    Complete validation report combining all checks.

    Attributes:
        timestamp: When validation was run
        data_sources: List of datasets used
        pathway_coverage: Per-pathway coverage results
        reactome_cross_ref: Cross-reference results against Reactome
        synthetic_validation: Synthetic data pipeline results
        biological_plausibility: Plausibility check result
        overall_pass: Whether all critical checks passed
        warnings: Any warnings generated
    """

    timestamp: str
    data_sources: List[DatasetInfo] = field(default_factory=list)
    pathway_coverage: List[PathwayCoverageResult] = field(default_factory=list)
    reactome_cross_ref: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    synthetic_validation: Dict[str, Any] = field(default_factory=dict)
    biological_plausibility: Optional[BiologicalPlausibilityResult] = None
    overall_pass: bool = False
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "timestamp": str(self.timestamp),
            "data_sources": [ds.to_dict() for ds in self.data_sources],
            "pathway_coverage": [pc.to_dict() for pc in self.pathway_coverage],
            "reactome_cross_ref": self.reactome_cross_ref,
            "synthetic_validation": self.synthetic_validation,
            "biological_plausibility": (
                self.biological_plausibility.to_dict()
                if self.biological_plausibility
                else None
            ),
            "overall_pass": bool(self.overall_pass),
            "warnings": list(self.warnings),
        }
        return result

    def format_report(self) -> str:
        """Generate a publication-ready markdown validation report."""
        lines = [
            "# Pathway Subtyping Framework — Validation Report",
            "",
            f"**Generated:** {self.timestamp}",
            f"**Overall status:** {'PASS' if self.overall_pass else 'NEEDS REVIEW'}",
            "",
        ]

        # Data sources
        if self.data_sources:
            lines.append("## Data Sources")
            lines.append("")
            for ds in self.data_sources:
                lines.append(f"- **{ds.name}**: {ds.description}")
                lines.append(f"  License: {ds.license_note}")
            lines.append("")

        # Pathway coverage
        if self.pathway_coverage:
            lines.append("## Pathway Coverage (ClinVar)")
            lines.append("")
            lines.append(
                "| Pathway | Genes | In ClinVar | With Pathogenic | Coverage | Pathogenic % |"
            )
            lines.append("|---------|-------|------------|-----------------|----------|--------------|")
            for pc in self.pathway_coverage:
                lines.append(
                    f"| {pc.pathway_name} | {pc.total_genes} "
                    f"| {pc.genes_in_clinvar} | {pc.genes_with_pathogenic} "
                    f"| {pc.coverage_fraction:.0%} | {pc.pathogenic_fraction:.0%} |"
                )
            lines.append("")

        # Reactome cross-reference
        if self.reactome_cross_ref:
            lines.append("## Reactome Cross-Reference")
            lines.append("")
            lines.append("| Curated Pathway | Best Reactome Match | Jaccard |")
            lines.append("|-----------------|---------------------|---------|")
            for pathway, match in self.reactome_cross_ref.items():
                best = match.get("best_match", "No match")
                jaccard = match.get("jaccard", 0.0)
                lines.append(f"| {pathway} | {best} | {jaccard:.3f} |")
            lines.append("")

        # Synthetic validation
        if self.synthetic_validation:
            lines.append("## Synthetic Data Validation")
            lines.append("")
            sv = self.synthetic_validation
            lines.append(f"- Samples: {sv.get('n_samples', 'N/A')}")
            lines.append(f"- Planted subtypes: {sv.get('n_subtypes', 'N/A')}")
            lines.append(f"- Recovered clusters: {sv.get('n_clusters_found', 'N/A')}")
            lines.append(f"- ARI vs true labels: {sv.get('ari', 'N/A')}")
            lines.append("")

        # Biological plausibility
        if self.biological_plausibility:
            bp = self.biological_plausibility
            lines.append("## Biological Plausibility")
            lines.append("")
            lines.append(f"- Disease: {bp.disease_name}")
            lines.append(f"- Subtypes found: {bp.n_subtypes}")
            lines.append(f"- Enriched pathways: {bp.n_enriched_pathways}")
            lines.append(
                f"- ClinVar overlap of enriched pathway genes: "
                f"{bp.pathway_gene_clinvar_overlap:.0%}"
            )
            lines.append(
                f"- Subtypes biologically distinct: "
                f"{'Yes' if bp.subtypes_biologically_distinct else 'No'}"
            )
            lines.append("")

        # Warnings
        if self.warnings:
            lines.append("## Warnings")
            lines.append("")
            for w in self.warnings:
                lines.append(f"- {w}")
            lines.append("")

        lines.append("---")
        lines.append("*Research use only. Not for clinical decision-making.*")

        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return citations for data sources used in validation."""
        return [
            (
                "Landrum MJ et al. ClinVar: improvements to accessing data. "
                "Nucleic Acids Res. 2020;48(D1):D835-D844."
            ),
            (
                "Gillespie M et al. The Reactome Pathway Knowledgebase 2022. "
                "Nucleic Acids Res. 2022;50(D1):D649-D653."
            ),
        ]


# =============================================================================
# DATASET REGISTRY
# =============================================================================


DATASETS = {
    "clinvar_gene_summary": DatasetInfo(
        name="ClinVar Gene-Specific Summary",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
            "gene_specific_summary.txt"
        ),
        description="Gene-level summary of ClinVar variant submissions per gene",
        file_format="tsv",
        license_note="Public domain (NCBI)",
    ),
    "reactome_pathways": DatasetInfo(
        name="Reactome Pathways GMT",
        url="https://reactome.org/download/current/ReactomePathways.gmt.zip",
        description="Reactome pathway gene sets in GMT format",
        file_format="gmt.zip",
        license_note="CC BY 4.0 (Reactome)",
    ),
}

DEFAULT_CACHE_DIR = Path(__file__).parent.parent.parent / "data" / "validation_cache"


# =============================================================================
# DOWNLOAD & CACHE
# =============================================================================


def download_dataset(
    dataset_key: str,
    cache_dir: Optional[Path] = None,
    force: bool = False,
    timeout: int = 60,
) -> Path:
    """
    Download a public dataset with caching.

    Uses urllib.request (stdlib) to avoid adding a requests dependency.
    Files are cached in cache_dir and not re-downloaded unless force=True.

    Args:
        dataset_key: Key from DATASETS dict (e.g., "clinvar_gene_summary")
        cache_dir: Directory for cached downloads (default: data/validation_cache/)
        force: If True, re-download even if cached
        timeout: Download timeout in seconds

    Returns:
        Path to the downloaded (and possibly decompressed) file

    Raises:
        KeyError: If dataset_key is not recognized
        urllib.error.URLError: If download fails
    """
    if dataset_key not in DATASETS:
        raise KeyError(
            f"Unknown dataset key: '{dataset_key}'. "
            f"Available: {list(DATASETS.keys())}"
        )

    dataset = DATASETS[dataset_key]
    if cache_dir is None:
        cache_dir = DEFAULT_CACHE_DIR

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Determine local filename
    url_filename = dataset.url.split("/")[-1]
    local_path = cache_dir / url_filename

    # Check cache
    if local_path.exists() and local_path.stat().st_size > 0 and not force:
        logger.info("[ValidationDatasets] Using cached %s: %s", dataset_key, local_path)
        return _decompress_if_needed(local_path, dataset.file_format)

    # Download
    logger.info("[ValidationDatasets] Downloading %s from %s", dataset.name, dataset.url)
    try:
        urllib.request.urlretrieve(dataset.url, str(local_path))
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
        raise urllib.error.URLError(
            f"Failed to download {dataset.name}: {e}"
        ) from e

    logger.info("[ValidationDatasets] Downloaded to %s (%d bytes)", local_path, local_path.stat().st_size)
    return _decompress_if_needed(local_path, dataset.file_format)


def _decompress_if_needed(path: Path, file_format: str) -> Path:
    """
    Decompress .zip files if needed. Returns path to usable file.

    Args:
        path: Path to the downloaded file
        file_format: Expected format (e.g., "gmt.zip", "tsv")

    Returns:
        Path to the extracted or original file
    """
    if file_format.endswith(".zip"):
        extract_dir = path.parent / path.stem
        # Check if already extracted
        if extract_dir.exists():
            gmt_files = list(extract_dir.glob("*.gmt"))
            if gmt_files:
                return gmt_files[0]

        extract_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "r") as zf:
            zf.extractall(extract_dir)

        # Find the GMT file
        gmt_files = list(extract_dir.glob("*.gmt"))
        if gmt_files:
            return gmt_files[0]

        # Fallback: return first extracted file
        extracted = list(extract_dir.iterdir())
        if extracted:
            return extracted[0]

        raise FileNotFoundError(f"No files found after extracting {path}")

    return path


# =============================================================================
# PARSING
# =============================================================================


def load_clinvar_gene_summary(
    file_path: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
) -> Dict[str, ClinVarGeneSummary]:
    """
    Load ClinVar gene-specific summary into a dict keyed by gene symbol.

    Downloads if not cached. Parses the NCBI tab-delimited format.

    Args:
        file_path: Explicit path to pre-downloaded file (skips download)
        cache_dir: Cache directory for automatic download

    Returns:
        Dict mapping gene symbol (str) to ClinVarGeneSummary
    """
    if file_path is None:
        file_path = download_dataset("clinvar_gene_summary", cache_dir=cache_dir)
    else:
        file_path = Path(file_path)

    logger.info("[ValidationDatasets] Loading ClinVar gene summary from %s", file_path)

    genes = {}
    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Parse header
            if line.startswith("#"):
                # The header line starts with #Symbol
                if "Symbol" in line:
                    header = line.lstrip("#").split("\t")
                continue

            if header is None:
                # Try to detect header without #
                if "Symbol" in line or "GeneID" in line:
                    header = line.split("\t")
                    continue
                continue

            parts = line.split("\t")
            if len(parts) < len(header):
                continue

            row = dict(zip(header, parts))

            symbol = row.get("Symbol", "").strip()
            if not symbol:
                continue

            def _safe_int(val, default=0):
                try:
                    return int(val)
                except (ValueError, TypeError):
                    return default

            genes[symbol.upper()] = ClinVarGeneSummary(
                symbol=symbol,
                gene_id=_safe_int(row.get("GeneID", "0")),
                n_pathogenic=_safe_int(row.get("Number_Pathogenic", "0")),
                n_likely_pathogenic=_safe_int(row.get("Number_Likely_Pathogenic", "0")),
                n_uncertain=_safe_int(row.get("Number_Uncertain_Significance", "0")),
                n_benign=_safe_int(row.get("Number_Benign", "0")),
                n_likely_benign=_safe_int(row.get("Number_Likely_Benign", "0")),
            )

    logger.info("[ValidationDatasets] Loaded %d genes from ClinVar", len(genes))
    return genes


def load_reactome_pathways(
    file_path: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
    species: str = "Homo sapiens",
) -> Dict[str, List[str]]:
    """
    Load Reactome pathways from GMT file, filtered to specified species.

    Downloads if not cached. Reactome GMT lines have the format:
    PathwayID<TAB>Description<TAB>Gene1<TAB>Gene2<TAB>...

    Species filtering uses multiple heuristics:
    - Description field containing species name (e.g., "Homo sapiens")
    - Reactome stable ID prefix (R-HSA- for human)

    Args:
        file_path: Explicit path to pre-downloaded GMT
        cache_dir: Cache directory for automatic download
        species: Species filter (default: "Homo sapiens")

    Returns:
        Dict mapping pathway name to list of gene symbols
    """
    if file_path is None:
        file_path = download_dataset("reactome_pathways", cache_dir=cache_dir)
    else:
        file_path = Path(file_path)

    logger.info("[ValidationDatasets] Loading Reactome pathways from %s", file_path)

    # Reactome stable ID species prefixes
    species_prefixes = {
        "Homo sapiens": "R-HSA-",
        "Mus musculus": "R-MMU-",
        "Rattus norvegicus": "R-RNO-",
        "Danio rerio": "R-DRE-",
        "Drosophila melanogaster": "R-DME-",
    }
    species_prefix = species_prefixes.get(species, "")

    pathways = {}
    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            pathway_name = parts[0]
            description = parts[1]
            genes = [g.strip() for g in parts[2:] if g.strip()]

            # Filter by species using description text or ID prefix
            if species:
                desc_match = species.lower() in description.lower() if description else False
                prefix_match = pathway_name.startswith(species_prefix) if species_prefix else False

                # If neither matches and species filter is active, skip
                if not desc_match and not prefix_match:
                    continue

            if genes:
                pathways[pathway_name] = genes

    logger.info(
        "[ValidationDatasets] Loaded %d %s pathways from Reactome",
        len(pathways),
        species,
    )
    return pathways


# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================


def validate_pathway_coverage(
    pathways: Dict[str, List[str]],
    clinvar_genes: Dict[str, ClinVarGeneSummary],
) -> List[PathwayCoverageResult]:
    """
    Validate pathway gene lists against ClinVar annotations.

    Checks what fraction of pathway genes appear in ClinVar and
    how many have pathogenic variant submissions.

    Args:
        pathways: Dict mapping pathway name to gene list
        clinvar_genes: ClinVar gene summary data

    Returns:
        List of PathwayCoverageResult, one per pathway, sorted by pathogenic fraction
    """
    logger.info(
        "[ValidationDatasets] Validating %d pathways against %d ClinVar genes",
        len(pathways),
        len(clinvar_genes),
    )

    # Normalize ClinVar keys to uppercase for matching
    clinvar_upper = {k.upper(): v for k, v in clinvar_genes.items()}

    results = []
    for pathway_name, genes in pathways.items():
        gene_details = []
        in_clinvar = 0
        with_pathogenic = 0

        for gene in genes:
            gene_upper = gene.upper()
            if gene_upper in clinvar_upper:
                cv = clinvar_upper[gene_upper]
                in_clinvar += 1
                has_pathogenic = cv.total_pathogenic > 0
                if has_pathogenic:
                    with_pathogenic += 1
                gene_details.append({
                    "gene": gene,
                    "in_clinvar": True,
                    "total_pathogenic": int(cv.total_pathogenic),
                })
            else:
                gene_details.append({
                    "gene": gene,
                    "in_clinvar": False,
                    "total_pathogenic": 0,
                })

        total = len(genes) if genes else 1  # Avoid division by zero

        results.append(PathwayCoverageResult(
            pathway_name=pathway_name,
            total_genes=len(genes),
            genes_in_clinvar=in_clinvar,
            genes_with_pathogenic=with_pathogenic,
            coverage_fraction=in_clinvar / total,
            pathogenic_fraction=with_pathogenic / total,
            gene_details=gene_details,
        ))

    # Sort by pathogenic fraction (most covered first)
    results.sort(key=lambda r: r.pathogenic_fraction, reverse=True)
    return results


def validate_pathway_against_reactome(
    pathways: Dict[str, List[str]],
    reactome_pathways: Dict[str, List[str]],
    min_overlap: float = 0.1,
) -> Dict[str, Dict[str, Any]]:
    """
    Cross-reference curated pathways against Reactome definitions.

    For each curated pathway, finds the best-matching Reactome pathway
    by Jaccard similarity. This validates that our pathway definitions
    are consistent with an independent curated database.

    Args:
        pathways: Curated pathway definitions
        reactome_pathways: Reactome pathway definitions
        min_overlap: Minimum Jaccard index to report a match

    Returns:
        Dict mapping curated pathway name to match details
    """
    logger.info(
        "[ValidationDatasets] Cross-referencing %d pathways against %d Reactome pathways",
        len(pathways),
        len(reactome_pathways),
    )

    # Pre-compute Reactome gene sets (uppercase)
    reactome_sets = {
        name: set(g.upper() for g in genes)
        for name, genes in reactome_pathways.items()
    }

    results = {}
    for pathway_name, genes in pathways.items():
        curated_set = set(g.upper() for g in genes)

        best_match = None
        best_jaccard = 0.0
        best_overlap_genes = []

        for rname, rset in reactome_sets.items():
            intersection = curated_set & rset
            union = curated_set | rset
            if not union:
                continue
            jaccard = len(intersection) / len(union)

            if jaccard > best_jaccard:
                best_jaccard = jaccard
                best_match = rname
                best_overlap_genes = sorted(intersection)

        results[pathway_name] = {
            "best_match": best_match if best_jaccard >= min_overlap else None,
            "jaccard": float(best_jaccard),
            "overlap_genes": best_overlap_genes if best_jaccard >= min_overlap else [],
            "n_overlap": len(best_overlap_genes) if best_jaccard >= min_overlap else 0,
        }

    return results


def generate_disease_realistic_synthetic(
    pathways: Dict[str, List[str]],
    clinvar_genes: Dict[str, ClinVarGeneSummary],
    n_samples: int = 100,
    n_subtypes: int = 3,
    effect_size: float = 1.0,
    noise_level: float = 1.0,
    seed: Optional[int] = 42,
) -> SimulatedData:
    """
    Generate synthetic data using real gene names from curated pathways.

    Unlike generate_synthetic_data() which uses GENE_0_0 placeholders,
    this creates data with actual gene symbols from the provided pathways,
    weighted by ClinVar pathogenicity counts to be more disease-realistic.

    Args:
        pathways: Real pathway definitions (e.g., from autism_pathways.gmt)
        clinvar_genes: ClinVar gene summary for realistic weighting
        n_samples: Number of synthetic samples
        n_subtypes: Number of planted subtypes
        effect_size: Cohen's d for subtype differences
        noise_level: Standard deviation of background noise
        seed: Random seed for reproducibility

    Returns:
        SimulatedData with real gene names and pathways
    """
    rng = np.random.RandomState(seed)
    logger.info(
        "[ValidationDatasets] Generating disease-realistic synthetic data: "
        "%d samples, %d subtypes, effect_size=%.1f",
        n_samples,
        n_subtypes,
        effect_size,
    )

    # Normalize ClinVar keys
    clinvar_upper = {k.upper(): v for k, v in clinvar_genes.items()}

    # Collect all genes from pathways
    all_genes = sorted(set(g for genes in pathways.values() for g in genes))
    if not all_genes:
        raise ValueError("No genes found in provided pathways")

    pathway_names = sorted(pathways.keys())
    n_pathways = len(pathway_names)

    # Assign samples to subtypes (equal proportions)
    n_per_subtype = [n_samples // n_subtypes] * n_subtypes
    remainder = n_samples - sum(n_per_subtype)
    for i in range(remainder):
        n_per_subtype[i] += 1

    true_labels = np.concatenate([np.full(n, i) for i, n in enumerate(n_per_subtype)])
    shuffle_idx = rng.permutation(n_samples)
    true_labels = true_labels[shuffle_idx]

    # Generate gene burdens with ClinVar-informed baseline
    gene_baselines = np.zeros(len(all_genes))
    for i, gene in enumerate(all_genes):
        cv = clinvar_upper.get(gene.upper())
        if cv is not None and cv.total_pathogenic > 0:
            # Genes with more pathogenic variants get higher baseline
            gene_baselines[i] = 0.3 + min(cv.total_pathogenic, 100) / 100.0 * 0.7
        else:
            gene_baselines[i] = 0.3

    gene_burdens = np.zeros((n_samples, len(all_genes)))
    for i, baseline in enumerate(gene_baselines):
        gene_burdens[:, i] = rng.exponential(scale=baseline, size=n_samples)

    # Select effect pathways for each subtype
    n_effect_pathways = max(2, n_pathways // n_subtypes)
    subtype_pathway_effects = {}

    for subtype in range(n_subtypes):
        start = (subtype * n_effect_pathways) % n_pathways
        effect_pws = []
        for j in range(n_effect_pathways):
            effect_pws.append(pathway_names[(start + j) % n_pathways])
        subtype_pathway_effects[subtype] = effect_pws

        # Add effect to genes in effect pathways for this subtype
        subtype_mask = true_labels == subtype
        effect_magnitude = effect_size * noise_level
        for pw in effect_pws:
            for gene in pathways[pw]:
                if gene in all_genes:
                    gene_idx = all_genes.index(gene)
                    gene_burdens[subtype_mask, gene_idx] += rng.normal(
                        effect_magnitude, 0.1, size=subtype_mask.sum()
                    )

    # Add noise
    gene_burdens += rng.normal(0, noise_level * 0.3, size=gene_burdens.shape)
    gene_burdens = np.clip(gene_burdens, 0, None)

    gene_burdens_df = pd.DataFrame(
        gene_burdens,
        columns=all_genes,
        index=[f"SAMPLE_{i}" for i in range(n_samples)],
    )

    # Compute pathway scores (mean of gene burdens per pathway, then Z-score)
    pathway_scores = np.zeros((n_samples, n_pathways))
    for j, pw_name in enumerate(pathway_names):
        pw_genes = [g for g in pathways[pw_name] if g in all_genes]
        if pw_genes:
            gene_indices = [all_genes.index(g) for g in pw_genes]
            pathway_scores[:, j] = gene_burdens[:, gene_indices].mean(axis=1)

    # Z-score normalize
    for j in range(n_pathways):
        col = pathway_scores[:, j]
        std = col.std()
        if std > 0:
            pathway_scores[:, j] = (col - col.mean()) / std
        else:
            pathway_scores[:, j] = 0.0

    pathway_scores_df = pd.DataFrame(
        pathway_scores,
        columns=pathway_names,
        index=[f"SAMPLE_{i}" for i in range(n_samples)],
    )

    config = SimulationConfig(
        n_samples=n_samples,
        n_pathways=n_pathways,
        n_genes_per_pathway=len(all_genes) // max(n_pathways, 1),
        n_subtypes=n_subtypes,
        effect_size=effect_size,
        noise_level=noise_level,
        seed=seed,
    )

    return SimulatedData(
        gene_burdens=gene_burdens_df,
        pathway_scores=pathway_scores_df,
        pathways=pathways,
        true_labels=true_labels,
        config=config,
        subtype_pathway_effects=subtype_pathway_effects,
    )


def run_biological_plausibility_check(
    pathway_scores: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathways: Dict[str, List[str]],
    clinvar_genes: Dict[str, ClinVarGeneSummary],
    disease_name: str = "autism",
    seed: Optional[int] = 42,
) -> BiologicalPlausibilityResult:
    """
    Check whether clustering results are biologically plausible.

    Runs characterization on the clustering result and validates that:
    1. Enriched pathways contain genes with ClinVar pathogenic variants
    2. Subtypes are driven by distinct biological processes

    Args:
        pathway_scores: Pathway score matrix (samples x pathways)
        cluster_labels: Cluster assignments
        pathways: Pathway definitions
        clinvar_genes: ClinVar data
        disease_name: Name for reporting
        seed: Random seed

    Returns:
        BiologicalPlausibilityResult
    """
    logger.info("[ValidationDatasets] Running biological plausibility check for %s", disease_name)

    unique_labels = np.unique(cluster_labels)
    n_subtypes = len(unique_labels)

    # Run characterization
    cluster_names = {int(k): f"Subtype_{k}" for k in unique_labels}
    char_result = characterize_subtypes(
        pathway_scores=pathway_scores,
        cluster_labels=cluster_labels,
        cluster_names=cluster_names,
        seed=seed,
    )

    # Count enriched pathways (significant from any subtype)
    enriched_pathways = set()
    for profile in char_result.subtype_profiles:
        for ep in profile.enriched_pathways:
            if ep.significant:
                enriched_pathways.add(ep.pathway)

    # Check ClinVar overlap for enriched pathway genes
    clinvar_upper = {k.upper(): v for k, v in clinvar_genes.items()}
    enriched_genes = set()
    for pw_name in enriched_pathways:
        if pw_name in pathways:
            enriched_genes.update(g.upper() for g in pathways[pw_name])

    if enriched_genes:
        genes_in_clinvar = sum(1 for g in enriched_genes if g in clinvar_upper)
        clinvar_overlap = genes_in_clinvar / len(enriched_genes)
    else:
        clinvar_overlap = 0.0

    # Check if subtypes have distinct top pathways
    subtype_top_pathways = []
    for profile in char_result.subtype_profiles:
        sig_pws = [ep.pathway for ep in profile.enriched_pathways if ep.significant]
        top_pw = sig_pws[0] if sig_pws else None
        subtype_top_pathways.append(top_pw)

    # Subtypes are distinct if they have different top enriched pathways
    non_none = [p for p in subtype_top_pathways if p is not None]
    subtypes_distinct = len(set(non_none)) == len(non_none) if non_none else False

    return BiologicalPlausibilityResult(
        disease_name=disease_name,
        n_subtypes=n_subtypes,
        n_enriched_pathways=len(enriched_pathways),
        pathway_gene_clinvar_overlap=float(clinvar_overlap),
        subtypes_biologically_distinct=subtypes_distinct,
        details={
            "enriched_pathways": sorted(enriched_pathways),
            "subtype_top_pathways": subtype_top_pathways,
            "n_enriched_genes": len(enriched_genes),
            "n_genes_in_clinvar": (
                sum(1 for g in enriched_genes if g in clinvar_upper)
                if enriched_genes
                else 0
            ),
        },
    )


# =============================================================================
# FULL VALIDATION ORCHESTRATOR
# =============================================================================


def run_full_validation(
    gmt_path: Optional[str] = None,
    disease_name: str = "autism",
    n_samples: int = 100,
    n_subtypes: int = 3,
    effect_size: float = 1.0,
    cache_dir: Optional[Path] = None,
    seed: Optional[int] = 42,
    skip_download: bool = False,
) -> ValidationReport:
    """
    Run the complete validation pipeline with public datasets.

    Workflow:
    1. Download ClinVar and Reactome data (with caching)
    2. Load curated pathway GMT
    3. Validate pathway coverage against ClinVar
    4. Cross-reference against Reactome
    5. Generate disease-realistic synthetic data
    6. Run clustering pipeline on synthetic data
    7. Check biological plausibility of results
    8. Produce validation report

    Args:
        gmt_path: Path to curated pathway GMT (default: autism_pathways.gmt)
        disease_name: Disease name for labeling
        n_samples: Synthetic cohort size
        n_subtypes: Expected subtypes
        effect_size: Signal strength for synthetic data
        cache_dir: Download cache directory
        seed: Random seed
        skip_download: If True, skip ClinVar/Reactome download (offline mode)

    Returns:
        ValidationReport with all results
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    warnings = []
    data_sources = []

    # Resolve GMT path
    if gmt_path is None:
        default_gmt = Path(__file__).parent.parent.parent / "data" / "pathways" / f"{disease_name}_pathways.gmt"
        if default_gmt.exists():
            gmt_path = str(default_gmt)
        else:
            raise FileNotFoundError(
                f"No GMT file found for disease '{disease_name}'. "
                f"Expected: {default_gmt}. Pass --gmt-path explicitly."
            )

    # Load curated pathways
    logger.info("[ValidationDatasets] Loading curated pathways from %s", gmt_path)
    curated_pathways = {}
    with open(gmt_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                curated_pathways[parts[0]] = parts[2:]

    logger.info("[ValidationDatasets] Loaded %d curated pathways", len(curated_pathways))

    # --- Step 1: ClinVar ---
    clinvar_genes = {}
    if not skip_download:
        try:
            clinvar_genes = load_clinvar_gene_summary(cache_dir=cache_dir)
            data_sources.append(DATASETS["clinvar_gene_summary"])
        except Exception as e:
            warnings.append(f"ClinVar download failed: {e}")
            logger.warning("[ValidationDatasets] ClinVar download failed: %s", e)
    else:
        # Try to load from cache
        if cache_dir is None:
            cache_dir = DEFAULT_CACHE_DIR
        cached = Path(cache_dir) / "gene_specific_summary.txt"
        if cached.exists():
            clinvar_genes = load_clinvar_gene_summary(file_path=cached)
            data_sources.append(DATASETS["clinvar_gene_summary"])
        else:
            warnings.append("ClinVar data not available (offline mode, no cache)")

    # --- Step 2: Reactome ---
    reactome_pathways = {}
    if not skip_download:
        try:
            reactome_pathways = load_reactome_pathways(cache_dir=cache_dir)
            data_sources.append(DATASETS["reactome_pathways"])
        except Exception as e:
            warnings.append(f"Reactome download failed: {e}")
            logger.warning("[ValidationDatasets] Reactome download failed: %s", e)
    else:
        if cache_dir is None:
            cache_dir = DEFAULT_CACHE_DIR
        cached_dir = Path(cache_dir) / "ReactomePathways"
        gmt_files = list(cached_dir.glob("*.gmt")) if cached_dir.exists() else []
        if gmt_files:
            reactome_pathways = load_reactome_pathways(file_path=gmt_files[0])
            data_sources.append(DATASETS["reactome_pathways"])
        else:
            warnings.append("Reactome data not available (offline mode, no cache)")

    # --- Step 3: Pathway coverage ---
    pathway_coverage = []
    if clinvar_genes:
        pathway_coverage = validate_pathway_coverage(curated_pathways, clinvar_genes)

    # --- Step 4: Reactome cross-reference ---
    reactome_cross_ref = {}
    if reactome_pathways:
        reactome_cross_ref = validate_pathway_against_reactome(
            curated_pathways, reactome_pathways
        )

    # --- Step 5-6: Synthetic data + clustering ---
    synthetic_validation = {}
    bio_plausibility = None

    try:
        sim_data = generate_disease_realistic_synthetic(
            pathways=curated_pathways,
            clinvar_genes=clinvar_genes if clinvar_genes else {},
            n_samples=n_samples,
            n_subtypes=n_subtypes,
            effect_size=effect_size,
            seed=seed,
        )

        # Run clustering
        from sklearn.metrics import adjusted_rand_score

        cluster_result = run_clustering(
            sim_data.pathway_scores,
            n_clusters=n_subtypes,
            seed=seed,
        )

        ari = adjusted_rand_score(sim_data.true_labels, cluster_result.labels)

        synthetic_validation = {
            "n_samples": n_samples,
            "n_subtypes": n_subtypes,
            "effect_size": effect_size,
            "n_clusters_found": int(cluster_result.n_clusters),
            "ari": float(round(ari, 4)),
            "silhouette": float(round(cluster_result.silhouette, 4)),
        }

        # --- Step 7: Biological plausibility ---
        if clinvar_genes:
            bio_plausibility = run_biological_plausibility_check(
                pathway_scores=sim_data.pathway_scores,
                cluster_labels=cluster_result.labels,
                pathways=curated_pathways,
                clinvar_genes=clinvar_genes,
                disease_name=disease_name,
                seed=seed,
            )

    except Exception as e:
        warnings.append(f"Synthetic validation failed: {e}")
        logger.warning("[ValidationDatasets] Synthetic validation error: %s", e)

    # --- Determine overall pass ---
    overall_pass = True

    # Check coverage: at least 50% of pathways should have >50% ClinVar coverage
    if pathway_coverage:
        good_coverage = sum(1 for pc in pathway_coverage if pc.coverage_fraction >= 0.5)
        if good_coverage < len(pathway_coverage) * 0.5:
            overall_pass = False
            warnings.append(
                f"Low ClinVar coverage: only {good_coverage}/{len(pathway_coverage)} "
                "pathways have >=50% gene coverage"
            )

    # Check synthetic ARI
    if synthetic_validation and synthetic_validation.get("ari", 0) < 0.5:
        overall_pass = False
        warnings.append(
            f"Low ARI ({synthetic_validation.get('ari', 0):.3f}) on synthetic data — "
            "framework may not recover planted subtypes"
        )

    # No data at all
    if not clinvar_genes and not reactome_pathways:
        overall_pass = False

    return ValidationReport(
        timestamp=timestamp,
        data_sources=data_sources,
        pathway_coverage=pathway_coverage,
        reactome_cross_ref=reactome_cross_ref,
        synthetic_validation=synthetic_validation,
        biological_plausibility=bio_plausibility,
        overall_pass=overall_pass,
        warnings=warnings,
    )
