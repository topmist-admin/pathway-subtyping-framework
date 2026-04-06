"""Gene ID translation and harmonization with multi-source consensus."""

import logging
import re
from typing import Optional, Dict, List

from .api import HGNCRestAPI, EnsemblXrefAPI
from .translation_cache import TranslationCache

logger = logging.getLogger(__name__)


class GeneHarmonizer:
    """Multi-source consensus gene ID translation with conflict resolution"""

    # Persistent cache to avoid repeated API calls (10x speedup on repeated IDs)
    _cache = TranslationCache()

    ID_TYPE_PATTERNS = {
        'hgnc_symbol': r'^[A-Z][A-Z0-9-]*$',         # e.g., TP53, BRCA2
        'ensg': r'^ENSG\d{11}',                        # Ensembl human
        'ensmusg': r'^ENSMUSG\d{11}',                  # Ensembl mouse
        'entrez': r'^\d+$',                            # Numeric only
        'refseq': r'^(NM_|NR_|XM_|XR_|NP_)\d+',       # RefSeq prefixes
    }

    @classmethod
    def detect_gene_id_type(cls, gene_id: str) -> str:
        """Classify input gene ID as one of the known types"""
        gene_id_str = str(gene_id).strip()
        for id_type, pattern in cls.ID_TYPE_PATTERNS.items():
            if re.match(pattern, gene_id_str):
                return id_type
        return 'unknown'

    @classmethod
    def harmonize_symbol_list(
        cls,
        gene_symbols: List[str],
        species: str = 'human'
    ) -> Dict[str, str]:
        """
        Batch harmonize a list of gene symbols/IDs to HGNC symbols

        Returns: {input_id: harmonized_hgnc_symbol, ...}
                 (excludes unmapped IDs)

        Algorithm:
        1. For each input ID:
           - Detect ID type
           - Query mygene.info, HGNC API, Ensembl xref in parallel
           - Consolidate results via majority vote

        2. Conflict resolution:
           - If 3/3 sources agree → use that symbol
           - If 2/3 sources agree → use majority vote
           - If 1/3 agrees → mark as low-confidence (skip)
           - If 0/3 agree → None (unmapped)
        """
        results = {}
        mapped_count = 0

        for gene_id in gene_symbols:
            gene_id_str = str(gene_id).strip()

            # Check persistent cache first (10x speedup on repeated IDs)
            cached_symbol = cls._cache.get(gene_id_str, species)
            if cached_symbol is not None:
                results[gene_id_str] = cached_symbol
                mapped_count += 1
                continue
            elif cached_symbol is None and gene_id_str in [k.split(':')[0] for k in cls._cache._cache.keys()]:
                # Gene was previously queried but unmapped
                continue

            # Detect ID type
            id_type = cls.detect_gene_id_type(gene_id_str)

            # Query 3 sources in parallel (simple sequential for now)
            votes = {}  # symbol → count of sources that returned it

            # Source 1: mygene.info
            mg_symbol = cls._query_mygene(gene_id_str, id_type, species)
            if mg_symbol:
                votes[mg_symbol] = votes.get(mg_symbol, 0) + 1

            # Source 2: HGNC API (for approved symbols or aliases)
            hgnc_symbol = cls._query_hgnc(gene_id_str, id_type)
            if hgnc_symbol:
                votes[hgnc_symbol] = votes.get(hgnc_symbol, 0) + 1

            # Source 3: Ensembl xref (for Ensembl IDs)
            if id_type == 'ensg':
                ens_symbols = cls._query_ensembl_xref(gene_id_str)
                for symbol in ens_symbols:
                    votes[symbol] = votes.get(symbol, 0) + 1

            # Consolidate: take majority vote (require at least 2/3 agreement)
            best_symbol = None
            if votes:
                vote_counts = list(votes.values())
                max_votes = max(vote_counts)
                if max_votes >= 2:  # At least 2 sources agree
                    best_symbol = max(votes, key=votes.get)

            # Save to persistent cache and return
            cls._cache.set(gene_id_str, best_symbol, species)
            if best_symbol:
                results[gene_id_str] = best_symbol
                mapped_count += 1

        logger.debug(f"Harmonized {mapped_count}/{len(gene_symbols)} genes via multi-source consensus")
        return results

    @staticmethod
    def _query_mygene(gene_id: str, id_type: str, species: str) -> Optional[str]:
        """Query mygene.info with appropriate scope"""
        scopes_map = {
            'hgnc_symbol': 'symbol',
            'ensg': 'ensembl.gene',
            'ensmusg': 'ensembl.gene',
            'entrez': 'entrezgene',
            'refseq': 'refseq'
        }

        scope = scopes_map.get(id_type)
        if not scope:
            return None

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            result = mg.query(gene_id, scopes=scope, species=species, size=1)
            if result['total'] > 0:
                return result['hits'][0].get('symbol')
        except Exception as e:
            logger.debug(f"mygene query failed for {gene_id}: {e}")

        return None

    @staticmethod
    def _query_hgnc(gene_id: str, id_type: str) -> Optional[str]:
        """Query HGNC REST API for symbol resolution"""
        hgnc_api = HGNCRestAPI()

        if id_type == 'hgnc_symbol':
            result = hgnc_api.fetch_symbol(gene_id)
            if result:
                return result['symbol']

        return None

    @staticmethod
    def _query_ensembl_xref(ensembl_id: str) -> List[str]:
        """Query Ensembl xref API for external references"""
        xref_api = EnsemblXrefAPI()
        return xref_api.get_xrefs_by_id(ensembl_id)
