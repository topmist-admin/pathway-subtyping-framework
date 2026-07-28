"""
Versioned knowledge-graph source manifest for v0.6 F3.

Pins PSF's upstream KG sources (OmniPath, SIGNOR, Reactome) to specific
releases so a KG built today can be rebuilt bit-exactly 18 months from
now. Each ``KGSource`` entry carries:

    - A stable source ID (``omnipath``, ``signor``, ``reactome``)
    - The release label we target (e.g. ``"2025"``)
    - A canonical download URL for that release's bundle
    - A SHA256 hash of the archive, so integrity is verifiable offline
    - Release date + citation

The ``KG_SOURCES`` registry is the single source of truth; loaders that
wire these archives into ``KnowledgeGraphBuilder`` read from it. Live
network-dependent fetching is gated on an optional ``[kg_sources]`` extra.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class KGSource:
    """A pinned upstream data source for the knowledge graph.

    Attributes:
        source_id: Stable short name (``omnipath``, ``signor``, ``reactome``).
        release: Release label we target (e.g. ``"2025"``, ``"3.0"``).
        url: Canonical bundle URL for the pinned release.
        sha256: Expected SHA-256 of the archive (64-char lowercase hex).
        release_date: ISO-8601 date string for the release.
        citation: Short literature / DOI reference for the release.
        notes: Free-text notes — e.g. edge-direction fixes vs prior release.
    """

    source_id: str
    release: str
    url: str
    sha256: str
    release_date: str
    citation: str
    notes: str = ""

    def __post_init__(self) -> None:
        if len(self.sha256) != 64 or not all(c in "0123456789abcdef" for c in self.sha256):
            raise ValueError(f"{self.source_id} sha256 must be 64 lowercase hex characters")

    def verify_archive(self, path: Path) -> bool:
        """Return True if the archive at ``path`` matches the pinned hash."""
        if not path.exists():
            return False
        h = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(8192), b""):
                h.update(chunk)
        observed = h.hexdigest()
        if observed != self.sha256:
            logger.warning(
                "[KGSource:%s/%s] hash mismatch: expected %s observed %s",
                self.source_id,
                self.release,
                self.sha256,
                observed,
            )
            return False
        return True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "source_id": self.source_id,
            "release": self.release,
            "url": self.url,
            "sha256": self.sha256,
            "release_date": self.release_date,
            "citation": self.citation,
            "notes": self.notes,
        }


# --------------------------------------------------------------------------- #
# Pinned v0.6 KG sources
# --------------------------------------------------------------------------- #

# The sha256 values below are PLACEHOLDERS pinned at the snapshot the v0.6
# release-engineering step will populate. Until the real bundles are
# archived into the PSF mirror, these hashes will fail verify_archive()
# — by design, so that an unpinned build is never silently accepted.
# The release-engineering workflow is documented in
# docs/migration/v05-to-v06-kg.md.

_PLACEHOLDER_HASH = "0" * 64

KG_SOURCES_V06: Dict[str, KGSource] = {
    "omnipath": KGSource(
        source_id="omnipath",
        release="2025",
        url="https://archive.omnipathdb.org/omnipath_2025_bundle.tar.gz",
        sha256=_PLACEHOLDER_HASH,
        release_date="2025-04-01",
        citation="Turei et al. 2016; OmniPath 2025 release",
        notes=(
            "+15% curated interactions vs 2024; adds edge-direction "
            "corrections in FGF, Notch, and WNT subnetworks."
        ),
    ),
    "signor": KGSource(
        source_id="signor",
        release="3.0",
        url="https://signor.uniroma2.it/archive/signor_3_0_bundle.tar.gz",
        sha256=_PLACEHOLDER_HASH,
        release_date="2025-11-15",
        citation="Licata et al. 2020; SIGNOR 3.0 release",
        notes="Adds neuronal-signalling subnetwork (~2,000 new edges).",
    ),
    "reactome": KGSource(
        source_id="reactome",
        release="2026",
        url="https://reactome.org/download/current/ReactomeBundle_2026.tar.gz",
        sha256=_PLACEHOLDER_HASH,
        release_date="2026-03-20",
        citation="Jassal et al. 2020; Reactome 2026 release",
        notes="Refactored cell-cycle pathway hierarchy; expect minor re-IDs.",
    ),
}


KG_SOURCES_V05: Dict[str, KGSource] = {
    "omnipath": KGSource(
        source_id="omnipath",
        release="2024",
        url="https://archive.omnipathdb.org/omnipath_2024_bundle.tar.gz",
        sha256=_PLACEHOLDER_HASH,
        release_date="2024-03-15",
        citation="Turei et al. 2016; OmniPath 2024 release",
    ),
    "signor": KGSource(
        source_id="signor",
        release="2.4",
        url="https://signor.uniroma2.it/archive/signor_2_4_bundle.tar.gz",
        sha256=_PLACEHOLDER_HASH,
        release_date="2023-11-10",
        citation="Licata et al. 2020; SIGNOR 2.4 release",
    ),
    "reactome": KGSource(
        source_id="reactome",
        release="2025",
        url="https://reactome.org/download/current/ReactomeBundle_2025.tar.gz",
        sha256=_PLACEHOLDER_HASH,
        release_date="2025-03-18",
        citation="Jassal et al. 2020; Reactome 2025 release",
    ),
}


# Current-default alias. Changing this is the mechanism by which PSF
# semver-bumps its KG.
KG_SOURCES: Dict[str, KGSource] = KG_SOURCES_V06


# --------------------------------------------------------------------------- #
# Access helpers
# --------------------------------------------------------------------------- #


def list_sources(registry: Optional[Dict[str, KGSource]] = None) -> List[KGSource]:
    """Return all pinned sources in the registry (defaults to v0.6)."""
    reg = registry if registry is not None else KG_SOURCES
    return list(reg.values())


def get_source(
    source_id: str,
    registry: Optional[Dict[str, KGSource]] = None,
) -> KGSource:
    """Look up a pinned source by ID. Raises KeyError if absent."""
    reg = registry if registry is not None else KG_SOURCES
    if source_id not in reg:
        raise KeyError(
            f"no pinned source '{source_id}' in registry; " f"available: {sorted(reg.keys())}"
        )
    return reg[source_id]


def manifest_digest(registry: Optional[Dict[str, KGSource]] = None) -> str:
    """Return a single SHA-256 over the whole pinned manifest.

    Two KGs are bitwise-reproducible only if their manifest digests match.
    The digest covers source_id + release + sha256 for every source, in
    deterministic order.
    """
    reg = registry if registry is not None else KG_SOURCES
    h = hashlib.sha256()
    for source_id in sorted(reg.keys()):
        src = reg[source_id]
        line = f"{src.source_id}|{src.release}|{src.sha256}".encode("utf-8")
        h.update(line)
        h.update(b"\n")
    return h.hexdigest()
