"""
Signaling Pathway Database Integration

Loads cell-cell signaling interaction databases (CellPhoneDB, CellChatDB)
and converts them into pathway gene sets (Dict[str, List[str]]) compatible
with the existing pathway scoring and clustering infrastructure.

Data sources:
- CellPhoneDB: Ligand-receptor interactions (auto-download from GitHub)
- CellChatDB: Ligand-receptor interactions (user-exported CSV from R)

References:
- Efremova M, et al. CellPhoneDB: inferring cell-cell communication
  from combined expression of multi-subunit ligand-receptor complexes.
  Nat Protoc. 2020;15(4):1484-1506.
- Jin S, et al. Inference and analysis of cell-cell communication using
  CellChat. Nat Commun. 2021;12(1):1088.

Research use only. Not for clinical decision-making.
"""

import csv
import logging
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CELLPHONEDB_INTERACTION_URL = (
    "https://raw.githubusercontent.com/ventolab/cellphonedb-data/"
    "master/data/interaction_input.csv"
)
CELLPHONEDB_GENE_URL = (
    "https://raw.githubusercontent.com/ventolab/cellphonedb-data/" "master/data/gene_input.csv"
)
CELLPHONEDB_COMPLEX_URL = (
    "https://raw.githubusercontent.com/ventolab/cellphonedb-data/" "master/data/complex_input.csv"
)

DEFAULT_CACHE_DIR = Path(__file__).parent.parent.parent / "data" / "validation_cache"

# ---------------------------------------------------------------------------
# Enum
# ---------------------------------------------------------------------------


class SignalingDatabase(Enum):
    """Supported signaling interaction databases."""

    CELLPHONEDB = "cellphonedb"
    CELLCHATDB = "cellchatdb"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class SignalingInteraction:
    """A single ligand-receptor interaction.

    Attributes:
        interaction_id: Unique identifier for the interaction.
        partner_a: First partner (gene symbol, UniProt ID, or complex name).
        partner_b: Second partner.
        genes_a: Resolved gene symbols for partner A.
        genes_b: Resolved gene symbols for partner B.
        pathway_name: Signaling pathway classification.
        source_database: Which database this came from.
        annotation: Additional annotation (e.g. "Secreted Signaling").
    """

    interaction_id: str
    partner_a: str
    partner_b: str
    genes_a: List[str]
    genes_b: List[str]
    pathway_name: str
    source_database: SignalingDatabase
    annotation: str = ""

    @property
    def all_genes(self) -> List[str]:
        """All unique gene symbols involved in this interaction, sorted."""
        return sorted(set(self.genes_a + self.genes_b))

    def to_dict(self) -> Dict[str, Any]:
        return {
            "interaction_id": self.interaction_id,
            "partner_a": self.partner_a,
            "partner_b": self.partner_b,
            "genes_a": list(self.genes_a),
            "genes_b": list(self.genes_b),
            "pathway_name": self.pathway_name,
            "source_database": self.source_database.value,
            "annotation": self.annotation,
        }


@dataclass
class SignalingDatabaseResult:
    """Result of loading a signaling database.

    Attributes:
        database: Which database was loaded.
        pathway_gene_sets: Pathway name -> sorted list of gene symbols.
        n_interactions: Total interactions parsed.
        n_pathways: Number of unique pathway gene sets.
        n_unique_genes: Total unique genes across all pathways.
        interactions: Full list of parsed interactions.
        warnings: Any warnings during parsing.
    """

    database: SignalingDatabase
    pathway_gene_sets: Dict[str, List[str]]
    n_interactions: int = 0
    n_pathways: int = 0
    n_unique_genes: int = 0
    interactions: List[SignalingInteraction] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "database": self.database.value,
            "n_interactions": int(self.n_interactions),
            "n_pathways": int(self.n_pathways),
            "n_unique_genes": int(self.n_unique_genes),
            "n_pathway_gene_sets": len(self.pathway_gene_sets),
            "pathway_names": sorted(self.pathway_gene_sets.keys()),
            "warnings": list(self.warnings),
        }

    def format_report(self) -> str:
        """Human-readable summary report."""
        lines = [
            "# Signaling Database Report",
            "",
            f"**Database:** {self.database.value}",
            f"**Interactions loaded:** {self.n_interactions}",
            f"**Pathway gene sets:** {self.n_pathways}",
            f"**Unique genes:** {self.n_unique_genes}",
            "",
            "## Pathway Gene Sets (top 20 by size)",
            "",
            "| Pathway | Genes |",
            "|---------|-------|",
        ]
        sorted_pws = sorted(self.pathway_gene_sets.items(), key=lambda x: len(x[1]), reverse=True)
        for name, genes in sorted_pws[:20]:
            lines.append(f"| {name} | {len(genes)} |")
        if len(sorted_pws) > 20:
            lines.append(f"| ... and {len(sorted_pws) - 20} more | |")

        if self.warnings:
            lines.extend(["", "## Warnings"])
            for w in self.warnings:
                lines.append(f"- {w}")

        lines.extend(
            [
                "",
                "---",
                "*Research use only. Not for clinical decision-making.*",
            ]
        )
        return "\n".join(lines)

    def get_citations(self) -> List[str]:
        """Return citations for the loaded database."""
        if self.database == SignalingDatabase.CELLPHONEDB:
            return [
                "Efremova M, et al. CellPhoneDB: inferring cell-cell communication "
                "from combined expression of multi-subunit ligand-receptor complexes. "
                "Nat Protoc. 2020;15(4):1484-1506."
            ]
        elif self.database == SignalingDatabase.CELLCHATDB:
            return [
                "Jin S, et al. Inference and analysis of cell-cell communication "
                "using CellChat. Nat Commun. 2021;12(1):1088."
            ]
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _download_file(
    url: str,
    cache_dir: Path,
    filename: str,
    force: bool = False,
    timeout: int = 60,
) -> Path:
    """Download a file with caching.

    Parameters
    ----------
    url : str
        URL to download.
    cache_dir : Path
        Directory for cached downloads.
    filename : str
        Local filename to save as.
    force : bool
        If True, re-download even if cached.
    timeout : int
        Download timeout in seconds.

    Returns
    -------
    Path
        Path to the downloaded (or cached) file.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    local_path = cache_dir / filename

    if local_path.exists() and local_path.stat().st_size > 0 and not force:
        logger.info("[SignalingDatabases] Using cached %s: %s", filename, local_path)
        return local_path

    logger.info("[SignalingDatabases] Downloading %s from %s", filename, url)
    try:
        urllib.request.urlretrieve(url, str(local_path))
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
        raise urllib.error.URLError(f"Failed to download {filename}: {e}") from e

    logger.info(
        "[SignalingDatabases] Downloaded %s (%d bytes)",
        local_path,
        local_path.stat().st_size,
    )
    return local_path


def _parse_cellphonedb_genes(gene_csv_path: Path) -> Dict[str, str]:
    """Parse CellPhoneDB gene_input.csv into UniProt -> gene symbol mapping.

    Parameters
    ----------
    gene_csv_path : Path
        Path to gene_input.csv.

    Returns
    -------
    Dict[str, str]
        Mapping from UniProt accession to HGNC gene symbol.
    """
    uniprot_to_gene: Dict[str, str] = {}
    with open(gene_csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            uniprot = row.get("uniprot", "").strip()
            gene_name = row.get("gene_name", "").strip()
            hgnc = row.get("hgnc_symbol", "").strip()
            symbol = hgnc if hgnc else gene_name
            if uniprot and symbol:
                uniprot_to_gene[uniprot] = symbol
    logger.info(
        "[SignalingDatabases] Parsed %d UniProt->gene mappings",
        len(uniprot_to_gene),
    )
    return uniprot_to_gene


def _parse_cellphonedb_complexes(
    complex_csv_path: Path,
    uniprot_to_gene: Dict[str, str],
) -> Dict[str, List[str]]:
    """Parse CellPhoneDB complex_input.csv into complex name -> gene list.

    Parameters
    ----------
    complex_csv_path : Path
        Path to complex_input.csv.
    uniprot_to_gene : Dict[str, str]
        UniProt -> gene symbol mapping.

    Returns
    -------
    Dict[str, List[str]]
        Mapping from complex name to list of gene symbols.
    """
    complexes: Dict[str, List[str]] = {}
    with open(complex_csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            complex_name = row.get("complex_name", "").strip()
            if not complex_name:
                continue
            genes: List[str] = []
            for col in [
                "uniprot_1",
                "uniprot_2",
                "uniprot_3",
                "uniprot_4",
                "uniprot_5",
            ]:
                uniprot = row.get(col, "").strip()
                if uniprot and uniprot in uniprot_to_gene:
                    genes.append(uniprot_to_gene[uniprot])
            if genes:
                complexes[complex_name] = genes
    logger.info("[SignalingDatabases] Parsed %d complexes", len(complexes))
    return complexes


def _resolve_partner_genes(
    partner: str,
    uniprot_to_gene: Dict[str, str],
    complexes: Dict[str, List[str]],
) -> List[str]:
    """Resolve a CellPhoneDB partner identifier to gene symbols.

    A partner can be a UniProt accession, a complex name, or a gene name.

    Parameters
    ----------
    partner : str
        Partner identifier string.
    uniprot_to_gene : Dict[str, str]
        UniProt -> gene mapping.
    complexes : Dict[str, List[str]]
        Complex name -> gene list mapping.

    Returns
    -------
    List[str]
        Resolved gene symbols (may be empty if unresolvable).
    """
    partner = partner.strip()
    if not partner:
        return []

    # Check complex lookup first
    if partner in complexes:
        return list(complexes[partner])

    # Check UniProt lookup
    if partner in uniprot_to_gene:
        return [uniprot_to_gene[partner]]

    # Check for "simple:GENENAME" or similar colon-prefixed patterns
    if ":" in partner:
        gene = partner.split(":")[-1].strip()
        if gene:
            return [gene]

    # Fallback: treat as a gene name directly
    return [partner]


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def load_cellphonedb(
    interactions_path: Optional[Path] = None,
    gene_path: Optional[Path] = None,
    complex_path: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
    force_download: bool = False,
    timeout: int = 60,
    min_genes_per_pathway: int = 2,
) -> SignalingDatabaseResult:
    """Load CellPhoneDB interactions and convert to pathway gene sets.

    Downloads interaction_input.csv, gene_input.csv, and complex_input.csv
    from the CellPhoneDB GitHub repository (or uses provided file paths).
    Resolves multi-subunit complexes to their constituent gene symbols and
    groups interactions by their ``classification`` field to produce
    pathway-like gene sets.

    Parameters
    ----------
    interactions_path : Path, optional
        Explicit path to interaction_input.csv (skips download).
    gene_path : Path, optional
        Explicit path to gene_input.csv (skips download).
    complex_path : Path, optional
        Explicit path to complex_input.csv (skips download).
    cache_dir : Path, optional
        Cache directory for downloads. Defaults to ``data/validation_cache/``.
    force_download : bool
        If True, re-download even if cached.
    timeout : int
        Download timeout in seconds.
    min_genes_per_pathway : int
        Minimum genes per pathway gene set to include (default: 2).

    Returns
    -------
    SignalingDatabaseResult
        Result with ``pathway_gene_sets`` as ``Dict[str, List[str]]``.

    Raises
    ------
    FileNotFoundError
        If explicit paths do not exist.
    ValueError
        If CSV files are empty or have missing required columns.
    urllib.error.URLError
        If download fails.
    """
    cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    warnings: List[str] = []

    # Resolve file paths (download if not provided)
    if interactions_path is not None:
        interactions_path = Path(interactions_path)
        if not interactions_path.exists():
            raise FileNotFoundError(f"Interactions file not found: {interactions_path}")
    else:
        interactions_path = _download_file(
            CELLPHONEDB_INTERACTION_URL,
            cache_dir,
            "cellphonedb_interaction_input.csv",
            force=force_download,
            timeout=timeout,
        )

    if gene_path is not None:
        gene_path = Path(gene_path)
        if not gene_path.exists():
            raise FileNotFoundError(f"Gene file not found: {gene_path}")
    else:
        gene_path = _download_file(
            CELLPHONEDB_GENE_URL,
            cache_dir,
            "cellphonedb_gene_input.csv",
            force=force_download,
            timeout=timeout,
        )

    if complex_path is not None:
        complex_path = Path(complex_path)
        if not complex_path.exists():
            raise FileNotFoundError(f"Complex file not found: {complex_path}")
    else:
        complex_path = _download_file(
            CELLPHONEDB_COMPLEX_URL,
            cache_dir,
            "cellphonedb_complex_input.csv",
            force=force_download,
            timeout=timeout,
        )

    # Parse gene mappings and complexes
    uniprot_to_gene = _parse_cellphonedb_genes(gene_path)
    complexes = _parse_cellphonedb_complexes(complex_path, uniprot_to_gene)

    # Parse interactions
    interactions: List[SignalingInteraction] = []
    n_skipped_no_classification = 0

    with open(interactions_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid CSV: {interactions_path}")
        required = {"partner_a", "partner_b"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing required columns in interactions CSV: {missing}")

        for row in reader:
            partner_a = row.get("partner_a", "").strip()
            partner_b = row.get("partner_b", "").strip()
            classification = row.get("classification", "").strip()
            interaction_id = row.get("id_cp_interaction", "").strip()
            annotation = row.get("annotation_strategy", "").strip()

            if not classification:
                n_skipped_no_classification += 1
                continue

            genes_a = _resolve_partner_genes(partner_a, uniprot_to_gene, complexes)
            genes_b = _resolve_partner_genes(partner_b, uniprot_to_gene, complexes)

            interactions.append(
                SignalingInteraction(
                    interaction_id=interaction_id,
                    partner_a=partner_a,
                    partner_b=partner_b,
                    genes_a=genes_a,
                    genes_b=genes_b,
                    pathway_name=classification,
                    source_database=SignalingDatabase.CELLPHONEDB,
                    annotation=annotation,
                )
            )

    if n_skipped_no_classification > 0:
        warnings.append(
            f"Skipped {n_skipped_no_classification} interactions " f"with empty classification"
        )

    if not interactions:
        raise ValueError(
            "No valid interactions found in CellPhoneDB data. "
            "Check that the CSV contains rows with non-empty classification."
        )

    # Convert to pathway gene sets
    pathway_gene_sets = convert_interactions_to_pathways(
        interactions,
        group_by="pathway_name",
        min_genes=min_genes_per_pathway,
        prefix="Signaling",
    )

    # Compute stats
    all_genes: Set[str] = set()
    for genes in pathway_gene_sets.values():
        all_genes.update(genes)

    result = SignalingDatabaseResult(
        database=SignalingDatabase.CELLPHONEDB,
        pathway_gene_sets=pathway_gene_sets,
        n_interactions=len(interactions),
        n_pathways=len(pathway_gene_sets),
        n_unique_genes=len(all_genes),
        interactions=interactions,
        warnings=warnings,
    )

    logger.info(
        "[SignalingDatabases] Loaded %d CellPhoneDB interactions, "
        "%d signaling pathway gene sets (%d unique genes)",
        result.n_interactions,
        result.n_pathways,
        result.n_unique_genes,
    )
    return result


def load_cellchatdb(
    file_path: Path,
    species: str = "human",
    min_genes_per_pathway: int = 2,
) -> SignalingDatabaseResult:
    """Load CellChatDB interactions from a user-exported CSV file.

    CellChatDB data is distributed as R binary ``.rda`` files, which cannot be
    read by Python without ``rpy2``. This function accepts a CSV file exported
    from R using the following commands::

        library(CellChat)
        write.csv(CellChatDB.human$interaction, "cellchatdb_human.csv",
                  row.names = FALSE)

    The CSV should contain columns: ``pathway_name``, ``ligand``, ``receptor``,
    and optionally ``interaction_name``, ``annotation``, ``evidence``.

    Parameters
    ----------
    file_path : Path
        Path to the exported CellChatDB CSV file.
    species : str
        Species label for reporting ("human", "mouse", or "zebrafish").
    min_genes_per_pathway : int
        Minimum genes per pathway gene set to include (default: 2).

    Returns
    -------
    SignalingDatabaseResult
        Result with ``pathway_gene_sets`` as ``Dict[str, List[str]]``.

    Raises
    ------
    FileNotFoundError
        If file_path does not exist.
    ValueError
        If required columns are missing or CSV is empty.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"CellChatDB file not found: {file_path}")

    warnings: List[str] = []
    interactions: List[SignalingInteraction] = []

    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid CSV: {file_path}")

        required = {"pathway_name", "ligand", "receptor"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns in CellChatDB CSV: {missing}. "
                f"Expected: pathway_name, ligand, receptor. "
                f"Export from R with: "
                f"write.csv(CellChatDB.{species}$interaction, ...)"
            )

        for row in reader:
            pathway_name = row.get("pathway_name", "").strip()
            ligand = row.get("ligand", "").strip()
            receptor = row.get("receptor", "").strip()
            interaction_name = row.get("interaction_name", "").strip()
            annotation = row.get("annotation", "").strip()

            if not pathway_name:
                continue

            genes_a = [ligand] if ligand else []
            genes_b = [receptor] if receptor else []

            interactions.append(
                SignalingInteraction(
                    interaction_id=interaction_name or f"{ligand}_{receptor}",
                    partner_a=ligand,
                    partner_b=receptor,
                    genes_a=genes_a,
                    genes_b=genes_b,
                    pathway_name=pathway_name,
                    source_database=SignalingDatabase.CELLCHATDB,
                    annotation=annotation,
                )
            )

    if not interactions:
        raise ValueError(f"No valid interactions found in CellChatDB CSV: {file_path}")

    # Convert to pathway gene sets
    pathway_gene_sets = convert_interactions_to_pathways(
        interactions,
        group_by="pathway_name",
        min_genes=min_genes_per_pathway,
        prefix="Signaling",
    )

    # Compute stats
    all_genes: Set[str] = set()
    for genes in pathway_gene_sets.values():
        all_genes.update(genes)

    result = SignalingDatabaseResult(
        database=SignalingDatabase.CELLCHATDB,
        pathway_gene_sets=pathway_gene_sets,
        n_interactions=len(interactions),
        n_pathways=len(pathway_gene_sets),
        n_unique_genes=len(all_genes),
        interactions=interactions,
        warnings=warnings,
    )

    logger.info(
        "[SignalingDatabases] Loaded %d CellChatDB interactions (%s), "
        "%d signaling pathway gene sets (%d unique genes)",
        result.n_interactions,
        species,
        result.n_pathways,
        result.n_unique_genes,
    )
    return result


def convert_interactions_to_pathways(
    interactions: List[SignalingInteraction],
    group_by: str = "pathway_name",
    min_genes: int = 2,
    prefix: str = "Signaling",
) -> Dict[str, List[str]]:
    """Convert a list of signaling interactions to pathway gene sets.

    Groups interactions by the specified field and collects the union
    of all gene symbols per group.

    Parameters
    ----------
    interactions : List[SignalingInteraction]
        List of interactions to convert.
    group_by : str
        Field to group by: ``"pathway_name"`` (default) or ``"annotation"``.
    min_genes : int
        Minimum genes per group to include (default: 2).
    prefix : str
        Prefix to add to pathway names (e.g. ``"Signaling"``).
        Set to ``""`` to disable prefixing.

    Returns
    -------
    Dict[str, List[str]]
        Mapping from pathway name to sorted list of unique gene symbols.
    """
    gene_sets: Dict[str, Set[str]] = defaultdict(set)

    for interaction in interactions:
        key = getattr(interaction, group_by, interaction.pathway_name)
        if not key:
            continue
        gene_sets[key].update(interaction.genes_a)
        gene_sets[key].update(interaction.genes_b)

    # Build final dict with prefix, filtering, and sorting
    result: Dict[str, List[str]] = {}
    for name, genes in sorted(gene_sets.items()):
        # Remove empty strings from gene sets
        genes.discard("")
        if len(genes) < min_genes:
            continue
        prefixed_name = f"{prefix}: {name}" if prefix else name
        result[prefixed_name] = sorted(genes)

    return result


def merge_signaling_databases(
    *results: SignalingDatabaseResult,
) -> Dict[str, List[str]]:
    """Merge pathway gene sets from multiple signaling database results.

    If both CellPhoneDB and CellChatDB define the same pathway name,
    their gene sets are unioned.

    Parameters
    ----------
    *results : SignalingDatabaseResult
        One or more signaling database results.

    Returns
    -------
    Dict[str, List[str]]
        Merged pathway name -> sorted list of unique gene symbols.
    """
    merged: Dict[str, Set[str]] = defaultdict(set)

    for result in results:
        for pathway_name, genes in result.pathway_gene_sets.items():
            merged[pathway_name].update(genes)

    return {name: sorted(genes) for name, genes in sorted(merged.items())}
