"""
Reference Atlas Comparison (F12).

Compares a manufacturing batch against an external reference standard
(published atlas, gold-standard batch, or healthy tissue reference)
to answer: "are we correct?" vs just "are we consistent?"

Metrics:
    - Per-cell and per-batch distance from reference in pathway-activation space
    - Nearest reference type mapping via atlas annotations
    - Deviation report for cells that fall outside atlas coverage

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class CellMapping:
    """Mapping of a single cell to its nearest reference type."""

    cell_id: str
    nearest_type: str
    distance: float
    on_target: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cell_id": self.cell_id,
            "nearest_type": self.nearest_type,
            "distance": round(self.distance, 4),
            "on_target": self.on_target,
        }


@dataclass
class AtlasComparisonResult:
    """Result of comparing a batch against a reference atlas."""

    median_distance: float
    p90_distance: float
    mean_distance: float
    n_on_target: int
    n_off_target: int
    n_total: int
    per_cell_distances: np.ndarray
    cell_mappings: List[CellMapping] = field(default_factory=list)
    type_distribution: Dict[str, int] = field(default_factory=dict)
    distance_threshold: float = 0.0
    summary: str = ""

    @property
    def passed(self) -> bool:
        return self.n_total > 0 and (self.n_on_target / self.n_total) >= 0.8

    @property
    def fraction_on_target(self) -> float:
        if self.n_total == 0:
            return 0.0
        return self.n_on_target / self.n_total

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "median_distance": round(self.median_distance, 4),
            "p90_distance": round(self.p90_distance, 4),
            "mean_distance": round(self.mean_distance, 4),
            "fraction_on_target": round(self.fraction_on_target, 4),
            "n_on_target": self.n_on_target,
            "n_off_target": self.n_off_target,
            "n_total": self.n_total,
            "type_distribution": self.type_distribution,
            "summary": self.summary,
        }


@dataclass
class AtlasEntry:
    """A reference atlas entry: profile + annotation."""

    name: str
    profile: Dict[str, float]
    annotation: str = ""
    source: str = ""


class AtlasComparator:
    """Compares batch pathway profiles against a reference atlas.

    Supports three atlas formats:
        1. Dict-based: dict of AtlasEntry objects
        2. DataFrame: pathway scores with row annotations
        3. CSV: file path to a pre-computed atlas

    Usage::

        atlas = AtlasComparator()
        atlas.add_reference("healthy_T_cell", {"HALLMARK_INFLAMMATORY_RESPONSE": 0.6, ...})
        atlas.add_reference("exhausted_T_cell", {"HALLMARK_INFLAMMATORY_RESPONSE": 0.2, ...})
        result = atlas.compare(batch_pathway_scores, expected_type="healthy_T_cell")
    """

    def __init__(
        self,
        distance_threshold: Optional[float] = None,
        seed: Optional[int] = None,
    ):
        """Initialize the atlas comparator.

        Args:
            distance_threshold: Max distance to be considered "on target."
                If None, auto-calibrated from atlas spread.
            seed: Random seed (unused, kept for API consistency).
        """
        self._entries: Dict[str, AtlasEntry] = {}
        self._distance_threshold = distance_threshold
        self._pathway_order: Optional[List[str]] = None

    @property
    def n_references(self) -> int:
        return len(self._entries)

    def add_reference(
        self,
        name: str,
        profile: Dict[str, float],
        annotation: str = "",
        source: str = "",
    ) -> None:
        """Add a reference profile to the atlas.

        Args:
            name: Unique name for this reference type (e.g., "healthy_T_cell").
            profile: Dict of pathway -> activation score.
            annotation: Description of this reference type.
            source: Source of the reference (e.g., "Tabula Sapiens, 2022").
        """
        self._entries[name] = AtlasEntry(
            name=name, profile=profile, annotation=annotation, source=source
        )
        logger.info("[QC Atlas] Added reference: %s (%d pathways)", name, len(profile))

    def load_atlas_csv(
        self,
        path: Union[str, Path],
        name_column: str = "cell_type",
        source: str = "",
    ) -> int:
        """Load reference profiles from a CSV file.

        CSV format: rows = reference types, columns = pathways,
        plus a name_column for the type label.

        Args:
            path: Path to the CSV file.
            name_column: Column containing reference type names.
            source: Source attribution for the atlas.

        Returns:
            Number of references loaded.
        """
        df = pd.read_csv(path)
        if name_column not in df.columns:
            raise ValueError(f"Name column '{name_column}' not found in CSV")

        pathway_cols = [c for c in df.columns if c != name_column]
        count = 0

        for _, row in df.iterrows():
            name = str(row[name_column])
            profile = {col: float(row[col]) for col in pathway_cols}
            self.add_reference(name, profile, source=source)
            count += 1

        logger.info("[QC Atlas] Loaded %d references from %s", count, path)
        return count

    def load_from_dataframe(
        self,
        df: pd.DataFrame,
        annotations: Optional[Dict[str, str]] = None,
        source: str = "",
    ) -> int:
        """Load reference profiles from a DataFrame.

        Args:
            df: DataFrame where index = reference type names, columns = pathways.
            annotations: Optional dict of name -> annotation.
            source: Source attribution.

        Returns:
            Number of references loaded.
        """
        annotations = annotations or {}
        count = 0

        for name in df.index:
            profile = df.loc[name].to_dict()
            annotation = annotations.get(str(name), "")
            self.add_reference(str(name), profile, annotation=annotation, source=source)
            count += 1

        return count

    def compare(
        self,
        pathway_scores: pd.DataFrame,
        expected_type: Optional[str] = None,
    ) -> AtlasComparisonResult:
        """Compare batch against the atlas.

        Args:
            pathway_scores: DataFrame (cells x pathways) of batch scores.
            expected_type: If specified, cells are "on target" only if their
                nearest reference is this type. Otherwise, on-target means
                within distance threshold of any reference.

        Returns:
            AtlasComparisonResult with per-cell distances and mappings.
        """
        if not self._entries:
            raise ValueError("Atlas is empty. Add references before comparing.")

        # Build reference matrix aligned to batch pathways
        batch_pathways = pathway_scores.columns.tolist()
        ref_names = list(self._entries.keys())
        ref_matrix = np.zeros((len(ref_names), len(batch_pathways)))

        for i, name in enumerate(ref_names):
            profile = self._entries[name].profile
            for j, pw in enumerate(batch_pathways):
                ref_matrix[i, j] = profile.get(pw, 0.0)

        batch_matrix = pathway_scores.values
        n_cells = len(batch_matrix)

        # Compute distances: each cell to each reference
        # distances shape: (n_cells, n_references)
        distances = np.zeros((n_cells, len(ref_names)))
        for i in range(len(ref_names)):
            diff = batch_matrix - ref_matrix[i]
            distances[:, i] = np.sqrt(np.sum(diff ** 2, axis=1))

        # Nearest reference per cell
        nearest_indices = np.argmin(distances, axis=1)
        nearest_distances = distances[np.arange(n_cells), nearest_indices]

        # Auto-calibrate threshold if not set
        threshold = self._distance_threshold
        if threshold is None:
            # Use median inter-reference distance as threshold
            if len(ref_names) > 1:
                inter_ref = []
                for i in range(len(ref_names)):
                    for j in range(i + 1, len(ref_names)):
                        inter_ref.append(np.linalg.norm(ref_matrix[i] - ref_matrix[j]))
                threshold = float(np.median(inter_ref))
            else:
                threshold = float(np.median(nearest_distances)) * 2.0

        # Classify on-target
        cell_mappings: List[CellMapping] = []
        type_counts: Dict[str, int] = {name: 0 for name in ref_names}
        n_on_target = 0

        for i in range(n_cells):
            nearest_name = ref_names[nearest_indices[i]]
            dist = float(nearest_distances[i])
            type_counts[nearest_name] += 1

            if expected_type is not None:
                on_target = nearest_name == expected_type and dist <= threshold
            else:
                on_target = dist <= threshold

            if on_target:
                n_on_target += 1

            cell_mappings.append(
                CellMapping(
                    cell_id=str(pathway_scores.index[i]),
                    nearest_type=nearest_name,
                    distance=dist,
                    on_target=on_target,
                )
            )

        median_dist = float(np.median(nearest_distances))
        p90_dist = float(np.percentile(nearest_distances, 90))
        mean_dist = float(np.mean(nearest_distances))
        n_off_target = n_cells - n_on_target

        summary = (
            f"Atlas comparison: {n_on_target}/{n_cells} on-target "
            f"({100 * n_on_target / max(n_cells, 1):.1f}%), "
            f"median distance={median_dist:.3f}, "
            f"P90={p90_dist:.3f}, threshold={threshold:.3f}. "
            f"{'PASS' if n_on_target / max(n_cells, 1) >= 0.8 else 'FAIL'}."
        )

        logger.info("[QC Atlas] %s", summary)

        return AtlasComparisonResult(
            median_distance=median_dist,
            p90_distance=p90_dist,
            mean_distance=mean_dist,
            n_on_target=n_on_target,
            n_off_target=n_off_target,
            n_total=n_cells,
            per_cell_distances=nearest_distances,
            cell_mappings=cell_mappings,
            type_distribution=type_counts,
            distance_threshold=threshold,
            summary=summary,
        )

    def atlas_registry(self) -> List[Dict[str, Any]]:
        """List all reference entries with metadata."""
        return [
            {
                "name": entry.name,
                "n_pathways": len(entry.profile),
                "annotation": entry.annotation,
                "source": entry.source,
            }
            for entry in self._entries.values()
        ]
