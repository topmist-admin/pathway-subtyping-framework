"""Gene ID translators for Phase 1 improvements.

Each translator handles a specific gene ID type and provides:
- detection: Can this translator handle this gene ID?
- translate: Convert to HGNC symbol
- translate_batch: Bulk conversion with caching
"""

import logging
import re
import time
from typing import Optional, List, Dict, Tuple

logger = logging.getLogger(__name__)

# Simple in-memory cache for API results during a single run
_TRANSLATION_CACHE = {}


class TranslatorBase:
    """Abstract base for gene ID translators."""

    @classmethod
    def can_handle(cls, gene_id: str) -> bool:
        """Check if this translator can handle the gene ID."""
        raise NotImplementedError

    @classmethod
    def translate_one(cls, gene_id: str) -> Optional[str]:
        """Translate a single gene ID to HGNC symbol."""
        raise NotImplementedError

    @classmethod
    def translate_batch(cls, gene_ids: List[str]) -> Dict[str, Optional[str]]:
        """Translate multiple gene IDs with caching."""
        results = {}
        for gene_id in gene_ids:
            result = cls.translate_one(gene_id)
            if result:
                results[gene_id] = result
        return results


class EntrezGeneTranslator(TranslatorBase):
    """
    Translate Entrez gene IDs (numeric) to HGNC symbols.

    Example: 7157 → TP53, 672 → BRCA1

    Methods:
    1. mygene.info batch API (primary) — querymany() for 1000+ IDs per request
    2. NCBI Entrez API (fallback)
    3. In-memory cache
    """

    ENTREZ_PATTERN = re.compile(r'^\d+$')  # Numeric only

    @classmethod
    def can_handle(cls, gene_id: str) -> bool:
        """Check if gene_id is numeric (Entrez format)."""
        return cls.ENTREZ_PATTERN.match(str(gene_id).strip()) is not None

    @classmethod
    def translate_one(cls, gene_id: str) -> Optional[str]:
        """
        Translate Entrez ID to HGNC symbol.

        Tries mygene.info first, falls back to NCBI if needed.
        """
        gene_id_str = str(gene_id).strip()

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            result = mg.query(gene_id_str, scopes='entrezgene', species='human', size=1)
            if result.get('total', 0) > 0:
                symbol = result['hits'][0].get('symbol')
                if symbol:
                    logger.debug(f"Entrez {gene_id_str} → {symbol} (via mygene)")
                    return symbol
        except Exception as e:
            logger.debug(f"mygene.info failed for Entrez {gene_id_str}: {e}")

        # Fallback: Try NCBI Entrez API
        return cls._query_ncbi_entrez(gene_id_str)

    @classmethod
    def translate_batch(cls, gene_ids: List[str]) -> Dict[str, Optional[str]]:
        """
        Batch translate Entrez IDs using mygene.info querymany().

        Much faster than translate_one() in a loop:
        - Single API call for up to 1000 IDs
        - ~50 IDs/second vs 1 ID every 0.5s (single query)
        - Expected 10x+ speedup
        """
        gene_ids_unique = list(set(str(gid).strip() for gid in gene_ids))

        # Prepopulate results with cached items
        results = {gid: _TRANSLATION_CACHE[gid] for gid in gene_ids_unique if gid in _TRANSLATION_CACHE}

        # Only query for uncached items
        uncached = [gid for gid in gene_ids_unique if gid not in _TRANSLATION_CACHE]
        if not uncached:
            return results

        logger.info(f"Translating {len(uncached)} Entrez IDs via mygene.info batch API...")

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            # Note: Don't call set_caching() - it requires extra libraries even to disable
            # Just use default behavior which is generally fast enough

            # Reduce batch size to 500 for faster/more reliable queries
            # Large batches (1000+) can timeout on slow networks
            batch_size = 500
            for i in range(0, len(uncached), batch_size):
                batch = uncached[i:i+batch_size]
                batch_num = i // batch_size + 1
                logger.debug(f"  Batch {batch_num}: querying {len(batch)} IDs (timeout 30s)...")

                try:
                    # Set request timeout to avoid infinite hangs
                    batch_results = mg.querymany(
                        batch,
                        scopes='entrezgene',
                        species='human',
                        verbose=False,
                        timeout=30  # 30 second timeout per batch
                    )
                    batch_success = 0
                    for hit in batch_results:
                        entrez_id = hit.get('query')
                        symbol = hit.get('symbol')
                        if entrez_id and symbol:
                            _TRANSLATION_CACHE[entrez_id] = symbol
                            results[entrez_id] = symbol
                            batch_success += 1
                    logger.debug(f"  Batch {batch_num}: {batch_success}/{len(batch)} translated")
                    # Rate limit: be respectful to the API
                    if i + batch_size < len(uncached):
                        time.sleep(0.5)  # Increased delay to avoid rate limiting
                except Exception as e:
                    logger.warning(f"  Batch {batch_num} query failed ({str(e)[:80]}). Trying medium batches (100 IDs)...")
                    # Fall back to medium-sized batches instead of one-at-a-time (much faster)
                    medium_batch_size = 100
                    medium_success = 0
                    for j in range(0, len(batch), medium_batch_size):
                        medium_batch = batch[j:j+medium_batch_size]
                        try:
                            batch_results = mg.querymany(
                                medium_batch,
                                scopes='entrezgene',
                                species='human',
                                verbose=False,
                                timeout=15  # Shorter timeout for fallback
                            )
                            for hit in batch_results:
                                entrez_id = hit.get('query')
                                symbol = hit.get('symbol')
                                if entrez_id and symbol:
                                    _TRANSLATION_CACHE[entrez_id] = symbol
                                    results[entrez_id] = symbol
                                    medium_success += 1
                            time.sleep(0.2)  # Rate limit between medium batches
                        except Exception as e2:
                            logger.debug(f"    Medium batch also failed: {str(e2)[:60]}. Skipping remaining {len(medium_batch)} IDs")
                            # Skip remaining IDs rather than one-at-a-time (prevents hour-long hangs)
                    if medium_success > 0:
                        logger.info(f"  Batch {batch_num}: Medium-batch fallback recovered {medium_success} IDs")

        except ImportError:
            logger.warning("mygene library not available. Skipping batch translation (mygene import failed)")
            # Don't fall back to one-at-a-time - just skip
            # The translators module requires mygene for batch operations

        logger.info(f"✅ Translated {len(results)}/{len(gene_ids_unique)} Entrez IDs")
        return results

    @staticmethod
    def _query_ncbi_entrez(gene_id: str) -> Optional[str]:
        """Query NCBI Entrez API for gene symbol (fallback)."""
        try:
            from Bio import Entrez
            Entrez.email = "research@example.com"

            # Fetch gene record
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
            from xml.etree import ElementTree as ET
            root = ET.parse(handle).getroot()

            # Extract gene symbol
            for elem in root.iter():
                if elem.tag == 'Gene-ref':
                    for child in elem:
                        if child.tag == 'Gene-ref_locus':
                            symbol = child.text
                            if symbol:
                                logger.debug(f"Entrez {gene_id} → {symbol} (via NCBI)")
                                return symbol

        except Exception as e:
            logger.debug(f"NCBI Entrez API failed for {gene_id}: {e}")

        return None


class GencodTranslator(TranslatorBase):
    """
    Translate Gencode versioned gene IDs to HGNC symbols.

    Example: ENSG00000141510.11 → TP53, ENSG00000012048.20 → BRCA1

    Gencode format: ENSG00000XXXXXX.VERSION (human)
    This translator extracts the Ensembl ID and translates via Ensembl API.
    """

    GENCODE_PATTERN = re.compile(r'^(ENSG\d{11})\.\d+$')  # ENSG + 11 digits + .version

    @classmethod
    def can_handle(cls, gene_id: str) -> bool:
        """Check if gene_id is Gencode versioned format."""
        return cls.GENCODE_PATTERN.match(str(gene_id).strip()) is not None

    @classmethod
    def translate_one(cls, gene_id: str) -> Optional[str]:
        """
        Translate Gencode versioned ID to HGNC symbol.

        Strips version and uses Ensembl xref API.
        """
        gene_id_str = str(gene_id).strip()

        # Extract unversioned Ensembl ID
        match = cls.GENCODE_PATTERN.match(gene_id_str)
        if not match:
            return None

        ensembl_id = match.group(1)

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            result = mg.query(ensembl_id, scopes='ensembl.gene', species='human', size=1)
            if result.get('total', 0) > 0:
                symbol = result['hits'][0].get('symbol')
                if symbol:
                    logger.debug(f"Gencode {gene_id_str} → {symbol} (via mygene)")
                    return symbol
        except Exception as e:
            logger.debug(f"Gencode translation failed for {gene_id_str}: {e}")

        return None

    @classmethod
    def translate_batch(cls, gene_ids: List[str]) -> Dict[str, Optional[str]]:
        """
        Batch translate Gencode IDs using mygene.info querymany().

        Single API call for up to 1000 IDs vs one call per ID.
        Expected 10x+ speedup.
        """
        gene_ids_unique = list(set(str(gid).strip() for gid in gene_ids))

        # Prepopulate results with cached items
        results = {gid: _TRANSLATION_CACHE[gid] for gid in gene_ids_unique if gid in _TRANSLATION_CACHE}

        # Only query for uncached items
        uncached = [gid for gid in gene_ids_unique if gid not in _TRANSLATION_CACHE]
        if not uncached:
            return results

        logger.info(f"Translating {len(uncached)} Gencode IDs via mygene.info batch API...")

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            # Note: Don't call set_caching() - it requires extra libraries even to disable
            # Just use default behavior which is generally fast enough

            # Extract unversioned Ensembl IDs
            ensembl_mapping = {}  # gencode_id -> ensembl_id
            reverse_mapping = {}  # ensembl_id -> [gencode_ids]
            for gid in uncached:
                match = cls.GENCODE_PATTERN.match(gid)
                if match:
                    ensembl_id = match.group(1)
                    ensembl_mapping[gid] = ensembl_id
                    if ensembl_id not in reverse_mapping:
                        reverse_mapping[ensembl_id] = []
                    reverse_mapping[ensembl_id].append(gid)

            if not ensembl_mapping:
                return results

            ensembl_ids = list(reverse_mapping.keys())

            # Query in batches (reduced to 500 for faster/more reliable queries)
            batch_size = 500
            for i in range(0, len(ensembl_ids), batch_size):
                batch = ensembl_ids[i:i+batch_size]
                batch_num = i // batch_size + 1
                logger.debug(f"  Batch {batch_num}: querying {len(batch)} Ensembl IDs (timeout 30s)...")

                try:
                    batch_results = mg.querymany(
                        batch,
                        scopes='ensembl.gene',
                        species='human',
                        verbose=False,
                        timeout=30  # 30 second timeout per batch
                    )
                    batch_success = 0
                    for hit in batch_results:
                        ensembl_id = hit.get('query')
                        symbol = hit.get('symbol')
                        if ensembl_id and symbol:
                            # Map back to original Gencode IDs (one to many)
                            for gencode_id in reverse_mapping.get(ensembl_id, []):
                                _TRANSLATION_CACHE[gencode_id] = symbol
                                results[gencode_id] = symbol
                                batch_success += 1
                    logger.debug(f"  Batch {batch_num}: {batch_success} Gencode IDs translated")
                    if i + batch_size < len(ensembl_ids):
                        time.sleep(0.5)  # Increased delay to avoid rate limiting
                except Exception as e:
                    logger.warning(f"  Batch query failed ({str(e)[:60]}). Trying medium batches (100 IDs)...")
                    # Fall back to medium-sized batches instead of one-at-a-time
                    medium_batch_size = 100
                    medium_success = 0
                    for j in range(0, len(batch), medium_batch_size):
                        medium_batch = batch[j:j+medium_batch_size]
                        try:
                            batch_results = mg.querymany(
                                medium_batch,
                                scopes='ensembl.gene',
                                species='human',
                                verbose=False,
                                timeout=15
                            )
                            for hit in batch_results:
                                ensembl_id = hit.get('query')
                                symbol = hit.get('symbol')
                                if ensembl_id and symbol:
                                    for gencode_id in reverse_mapping.get(ensembl_id, []):
                                        _TRANSLATION_CACHE[gencode_id] = symbol
                                        results[gencode_id] = symbol
                                        medium_success += 1
                            time.sleep(0.2)
                        except Exception as e2:
                            logger.debug(f"    Medium batch failed: {str(e2)[:60]}. Skipping {len(medium_batch)} IDs")
                    if medium_success > 0:
                        logger.info(f"  Batch {batch_num}: Medium-batch fallback recovered {medium_success} Gencode IDs")

        except ImportError:
            logger.warning("mygene not available. Skipping Gencode translation...")

        logger.info(f"✅ Translated {len(results)}/{len(gene_ids_unique)} Gencode IDs")
        return results


class RefSeqTranslator(TranslatorBase):
    """
    Translate RefSeq accession numbers to HGNC symbols.

    Example: NM_000546 → TP53, NM_007294 → BRCA1

    RefSeq formats:
    - NM_: mRNA
    - NR_: non-coding RNA
    - NP_: protein
    - XM_: predicted mRNA
    - XR_: predicted non-coding RNA
    - XP_: predicted protein
    - NG_: genomic DNA
    """

    REFSEQ_PATTERN = re.compile(r'^(NM_|NR_|XM_|XR_|NP_|XP_|NG_)\d+(\.\d+)?$')

    @classmethod
    def can_handle(cls, gene_id: str) -> bool:
        """Check if gene_id is RefSeq format."""
        return cls.REFSEQ_PATTERN.match(str(gene_id).strip()) is not None

    @classmethod
    def translate_one(cls, gene_id: str) -> Optional[str]:
        """
        Translate RefSeq accession to HGNC symbol.

        Uses mygene.info with refseq scope.
        """
        gene_id_str = str(gene_id).strip()

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            result = mg.query(gene_id_str, scopes='refseq', species='human', size=1)
            if result.get('total', 0) > 0:
                symbol = result['hits'][0].get('symbol')
                if symbol:
                    logger.debug(f"RefSeq {gene_id_str} → {symbol} (via mygene)")
                    return symbol
        except Exception as e:
            logger.debug(f"RefSeq translation failed for {gene_id_str}: {e}")

        return None

    @classmethod
    def translate_batch(cls, gene_ids: List[str]) -> Dict[str, Optional[str]]:
        """
        Batch translate RefSeq IDs using mygene.info querymany().

        Single API call for up to 1000 IDs vs one call per ID.
        Expected 10x+ speedup.
        """
        gene_ids_unique = list(set(str(gid).strip() for gid in gene_ids))

        # Prepopulate results with cached items
        results = {gid: _TRANSLATION_CACHE[gid] for gid in gene_ids_unique if gid in _TRANSLATION_CACHE}

        # Only query for uncached items
        uncached = [gid for gid in gene_ids_unique if gid not in _TRANSLATION_CACHE]
        if not uncached:
            return results

        logger.info(f"Translating {len(uncached)} RefSeq IDs via mygene.info batch API...")

        try:
            import mygene
            mg = mygene.MyGeneInfo()
            # Note: Don't call set_caching() - it requires extra libraries even to disable
            # Just use default behavior which is generally fast enough

            batch_size = 500  # Reduced from 1000 for faster/more reliable queries
            for i in range(0, len(uncached), batch_size):
                batch = uncached[i:i+batch_size]
                batch_num = i // batch_size + 1
                logger.debug(f"  Batch {batch_num}: querying {len(batch)} IDs (timeout 30s)...")

                try:
                    batch_results = mg.querymany(
                        batch,
                        scopes='refseq',
                        species='human',
                        verbose=False,
                        timeout=30  # 30 second timeout per batch
                    )
                    batch_success = 0
                    for hit in batch_results:
                        refseq_id = hit.get('query')
                        symbol = hit.get('symbol')
                        if refseq_id and symbol:
                            _TRANSLATION_CACHE[refseq_id] = symbol
                            results[refseq_id] = symbol
                            batch_success += 1
                    logger.debug(f"  Batch {batch_num}: {batch_success}/{len(batch)} translated")
                    if i + batch_size < len(uncached):
                        time.sleep(0.5)  # Increased delay to avoid rate limiting
                except Exception as e:
                    logger.warning(f"  Batch query failed ({str(e)[:60]}). Trying medium batches (100 IDs)...")
                    # Fall back to medium-sized batches instead of one-at-a-time
                    medium_batch_size = 100
                    medium_success = 0
                    for j in range(0, len(batch), medium_batch_size):
                        medium_batch = batch[j:j+medium_batch_size]
                        try:
                            batch_results = mg.querymany(
                                medium_batch,
                                scopes='refseq',
                                species='human',
                                verbose=False,
                                timeout=15
                            )
                            for hit in batch_results:
                                refseq_id = hit.get('query')
                                symbol = hit.get('symbol')
                                if refseq_id and symbol:
                                    _TRANSLATION_CACHE[refseq_id] = symbol
                                    results[refseq_id] = symbol
                                    medium_success += 1
                            time.sleep(0.2)
                        except Exception as e2:
                            logger.debug(f"    Medium batch failed: {str(e2)[:60]}. Skipping {len(medium_batch)} IDs")
                    if medium_success > 0:
                        logger.info(f"  Batch {batch_num}: Medium-batch fallback recovered {medium_success} RefSeq IDs")

        except ImportError:
            logger.warning("mygene not available. Skipping RefSeq translation...")

        logger.info(f"✅ Translated {len(results)}/{len(gene_ids_unique)} RefSeq IDs")
        return results


# Translator registry (order matters: test in this order)
# RefSeqTranslator first (more specific patterns)
# Then EST/Gencode/Entrez (broader patterns)
TRANSLATORS = [
    RefSeqTranslator,
    GencodTranslator,
    EntrezGeneTranslator,
]


def detect_and_translate(gene_id: str) -> Optional[Tuple[str, Optional[str]]]:
    """
    Detect gene ID type and translate to HGNC symbol.

    Returns:
        (translator_name, hgnc_symbol) or (None, None) if no translator matched
    """
    gene_id_str = str(gene_id).strip()

    for translator_class in TRANSLATORS:
        if translator_class.can_handle(gene_id_str):
            symbol = translator_class.translate_one(gene_id_str)
            return (translator_class.__name__, symbol)

    return (None, None)


def batch_translate(gene_ids: List[str]) -> Dict[str, Optional[str]]:
    """
    Translate a batch of mixed gene ID types.

    Groups by translator class for efficient batch processing.
    """
    results = {}

    # Group gene IDs by translator
    groups = {translator: [] for translator in TRANSLATORS}
    unhandled = []

    for gene_id in gene_ids:
        gene_id_str = str(gene_id).strip()
        handled = False
        for translator in TRANSLATORS:
            if translator.can_handle(gene_id_str):
                groups[translator].append(gene_id_str)
                handled = True
                break
        if not handled:
            unhandled.append(gene_id_str)

    # Process each group with its translator
    for translator, ids in groups.items():
        if ids:
            batch_results = translator.translate_batch(ids)
            results.update(batch_results)

    if unhandled:
        logger.debug(f"No translator found for {len(unhandled)} gene IDs: {unhandled[:3]}...")

    return results
