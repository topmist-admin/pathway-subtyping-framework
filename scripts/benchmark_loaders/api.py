"""API clients for gene nomenclature and cross-reference services."""

import logging
from typing import Optional, Dict, List
import requests

logger = logging.getLogger(__name__)


class HGNCRestAPI:
    """Query official HGNC (Hugo Gene Nomenclature Committee) database"""

    BASE_URL = "https://rest.genenames.org"

    @classmethod
    def fetch_symbol(cls, symbol: str) -> Optional[Dict]:
        """
        Query HGNC for a symbol by checking three endpoints in order:
        1. /fetch?hgnc_symbol=SYMBOL (approved symbol)
        2. /fetch?alias_name=SYMBOL (is it an alias?)
        3. /fetch?prev_name=SYMBOL (is it a withdrawn symbol?)

        Returns: {'symbol': 'APPROVED', 'status': 'approved'|'alias'|'withdrawn'}
                or None if not found
        """
        endpoints = [
            f"{cls.BASE_URL}/fetch?hgnc_symbol={symbol}",
            f"{cls.BASE_URL}/fetch?alias_name={symbol}",
            f"{cls.BASE_URL}/fetch?prev_name={symbol}"
        ]

        for endpoint in endpoints:
            try:
                resp = requests.get(endpoint, timeout=5)
                if resp.status_code == 200:
                    data = resp.json()
                    if data.get('response', {}).get('numFound', 0) > 0:
                        doc = data['response']['docs'][0]
                        return {
                            'symbol': doc.get('hgnc_symbol', doc.get('symbol')),
                            'status': 'approved' if 'hgnc_symbol' in endpoint else 'alias_or_withdrawn'
                        }
            except Exception as e:
                logger.debug(f"HGNC lookup failed for {symbol}: {e}")

        return None


class EnsemblXrefAPI:
    """Query Ensembl external references for ID resolution"""

    BASE_URL = "https://rest.ensembl.org"

    @classmethod
    def get_xrefs_by_id(cls, ensembl_id: str) -> List[str]:
        """
        For an Ensembl gene ID (ENSG...), return all external gene references
        including: gene symbol, Entrez, RefSeq, HGNC, CCDS, etc.

        Endpoint: /xrefs/id/{id}?external_db=GO,HGNC,Entrez_Gene,RefSeq_peptide

        Return: [list of external gene symbols/IDs]
        """
        url = f"{cls.BASE_URL}/xrefs/id/{ensembl_id}?external_db=HGNC,Entrez_Gene"
        try:
            resp = requests.get(url, timeout=5)
            if resp.status_code == 200:
                xrefs = resp.json()
                # Extract all HGNC symbols (should be primary hit)
                symbols = [x['display_id'] for x in xrefs if x.get('db_display_name') == 'HGNC']
                return symbols
        except Exception as e:
            logger.debug(f"Ensembl xref lookup failed for {ensembl_id}: {e}")

        return []
