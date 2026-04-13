"""
Drug-Pathway Mapping for Therapeutic Hypothesis Generation.

Maps disrupted pathways to candidate drugs via gene-level targeting,
with mechanism classification and ASD relevance scoring.

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger(__name__)


class DrugMechanism(Enum):
    """Drug mechanism of action."""

    ANTAGONIST = "antagonist"
    AGONIST = "agonist"
    INHIBITOR = "inhibitor"
    ACTIVATOR = "activator"
    MODULATOR = "modulator"
    BLOCKER = "blocker"
    ENHANCER = "enhancer"
    STABILIZER = "stabilizer"
    UNKNOWN = "unknown"


class DrugStatus(Enum):
    """Drug development status."""

    APPROVED = "approved"
    INVESTIGATIONAL = "investigational"
    EXPERIMENTAL = "experimental"
    WITHDRAWN = "withdrawn"
    UNKNOWN = "unknown"


@dataclass
class DrugCandidate:
    """A candidate drug for therapeutic repurposing."""

    drug_id: str
    drug_name: str
    target_genes: List[str] = field(default_factory=list)
    mechanism: str = ""
    mechanism_type: DrugMechanism = DrugMechanism.UNKNOWN
    indications: List[str] = field(default_factory=list)
    contraindications: List[str] = field(default_factory=list)
    asd_relevance_score: float = 0.0
    status: DrugStatus = DrugStatus.UNKNOWN
    pathways: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "drug_id": self.drug_id,
            "drug_name": self.drug_name,
            "target_genes": self.target_genes,
            "mechanism": self.mechanism,
            "mechanism_type": self.mechanism_type.value,
            "status": self.status.value,
            "asd_relevance_score": round(self.asd_relevance_score, 4),
            "pathways": self.pathways,
        }


@dataclass
class DrugTargetDatabase:
    """In-memory database of drug-gene-pathway associations."""

    drugs: Dict[str, DrugCandidate] = field(default_factory=dict)
    gene_to_drugs: Dict[str, Set[str]] = field(default_factory=dict)

    def add_drug(self, drug: DrugCandidate) -> None:
        self.drugs[drug.drug_id] = drug
        for gene in drug.target_genes:
            if gene not in self.gene_to_drugs:
                self.gene_to_drugs[gene] = set()
            self.gene_to_drugs[gene].add(drug.drug_id)

    def get_drug(self, drug_id: str) -> Optional[DrugCandidate]:
        return self.drugs.get(drug_id)

    def get_drugs_for_gene(self, gene_id: str) -> List[DrugCandidate]:
        drug_ids = self.gene_to_drugs.get(gene_id, set())
        return [self.drugs[did] for did in drug_ids if did in self.drugs]

    def search_drugs(
        self,
        mechanism_type: Optional[DrugMechanism] = None,
        status: Optional[DrugStatus] = None,
        min_asd_relevance: float = 0.0,
    ) -> List[DrugCandidate]:
        results = list(self.drugs.values())
        if mechanism_type:
            results = [d for d in results if d.mechanism_type == mechanism_type]
        if status:
            results = [d for d in results if d.status == status]
        results = [d for d in results if d.asd_relevance_score >= min_asd_relevance]
        return results

    def get_statistics(self) -> Dict[str, Any]:
        return {
            "n_drugs": len(self.drugs),
            "n_gene_targets": len(self.gene_to_drugs),
            "by_status": {
                s.value: sum(1 for d in self.drugs.values() if d.status == s)
                for s in DrugStatus
            },
        }


class PathwayDrugMapper:
    """Maps disrupted pathways to candidate drugs.

    Usage::

        mapper = PathwayDrugMapper(drug_db)
        candidates = mapper.map_pathway("MAPK", pathway_genes=["RAS", "RAF"])
    """

    def __init__(
        self,
        drug_db: DrugTargetDatabase,
        include_experimental: bool = False,
        min_asd_relevance: float = 0.0,
    ):
        self.drug_db = drug_db
        self.include_experimental = include_experimental
        self.min_asd_relevance = min_asd_relevance

    def map_pathway(
        self,
        pathway_id: str,
        pathway_genes: Optional[List[str]] = None,
        disrupted_genes: Optional[List[str]] = None,
    ) -> List[DrugCandidate]:
        """Find drugs targeting genes in a pathway.

        Args:
            pathway_id: Pathway identifier.
            pathway_genes: Genes in the pathway.
            disrupted_genes: Genes with damaging variants.

        Returns:
            List of DrugCandidate sorted by relevance.
        """
        genes = pathway_genes or []
        candidates: Dict[str, DrugCandidate] = {}

        for gene in genes:
            drugs = self.drug_db.get_drugs_for_gene(gene)
            for drug in drugs:
                if not self._filter_drug(drug):
                    continue
                if drug.drug_id not in candidates:
                    candidates[drug.drug_id] = drug

        result = list(candidates.values())
        result.sort(key=lambda d: -d.asd_relevance_score)
        return result

    def map_multiple_pathways(
        self,
        pathway_scores: Dict[str, float],
        pathway_genes: Optional[Dict[str, List[str]]] = None,
        min_zscore: float = 1.5,
    ) -> Dict[str, List[DrugCandidate]]:
        """Map drugs across multiple disrupted pathways.

        Args:
            pathway_scores: Pathway -> z-score.
            pathway_genes: Pathway -> gene list.
            min_zscore: Minimum absolute z-score for a pathway to be considered.

        Returns:
            Dict mapping pathway_id to list of DrugCandidate.
        """
        pw_genes = pathway_genes or {}
        results: Dict[str, List[DrugCandidate]] = {}

        for pw, score in pathway_scores.items():
            if abs(score) < min_zscore:
                continue
            genes = pw_genes.get(pw, [])
            candidates = self.map_pathway(pw, genes)
            if candidates:
                results[pw] = candidates

        return results

    def _filter_drug(self, drug: DrugCandidate) -> bool:
        if drug.asd_relevance_score < self.min_asd_relevance:
            return False
        if drug.status == DrugStatus.WITHDRAWN:
            return False
        if drug.status == DrugStatus.EXPERIMENTAL and not self.include_experimental:
            return False
        return True
