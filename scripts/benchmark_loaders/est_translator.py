"""EST (Expressed Sequence Tag) probe translation to HGNC symbols.

EST probes are short DNA sequences without direct HGNC annotation.
They require multi-step mapping:
  EST sequence → UniGene ID → Gene symbol

Supports platforms:
- GPL3427 (TIGR 40K EST array)
- GPL9603 (other EST arrays)
- Custom EST collections
"""

import logging
import re
from typing import Optional, Dict, List, Tuple
import json
import time

logger = logging.getLogger(__name__)


class UniGeneMapper:
    """
    Map EST sequences and UniGene IDs to HGNC symbols.

    Uses UniGene database which maps:
    - UniGene cluster ID (Hs.123456) → Gene symbol
    - EST accession → UniGene cluster

    Local cache: ~/.cache/pathway_subtyping/unigene_mapping.json
    """

    CACHE_FILE = None  # Will be set dynamically

    # UniGene cluster pattern: Hs.123456, Mm.45678, etc. (needs at least 5 digits)
    UNIGENE_PATTERN = re.compile(r'^[A-Za-z]{2}\.(\d{5,})$', re.IGNORECASE)

    @classmethod
    def _get_cache_file(cls):
        """Get or create cache file path."""
        if cls.CACHE_FILE is None:
            from pathlib import Path
            cache_dir = Path.home() / ".cache" / "pathway_subtyping"
            cache_dir.mkdir(parents=True, exist_ok=True)
            cls.CACHE_FILE = cache_dir / "unigene_mapping.json"
        return cls.CACHE_FILE

    @classmethod
    def _load_cache(cls) -> Dict:
        """Load UniGene mapping cache from disk."""
        cache_file = cls._get_cache_file()
        if not cache_file.exists():
            return {}
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            logger.debug(f"Failed to load UniGene cache: {e}")
            return {}

    @classmethod
    def _save_cache(cls, cache: Dict):
        """Save UniGene mapping cache to disk."""
        cache_file = cls._get_cache_file()
        try:
            with open(cache_file, 'w') as f:
                json.dump(cache, f, indent=2)
        except IOError as e:
            logger.debug(f"Failed to save UniGene cache: {e}")

    @classmethod
    def unigene_to_symbol(cls, unigene_id: str) -> Optional[str]:
        """
        Map UniGene cluster ID to HGNC symbol.

        Args:
            unigene_id: UniGene cluster ID (e.g., "Hs.123456")

        Returns:
            HGNC symbol or None if not found
        """
        unigene_id_str = str(unigene_id).strip()

        # Check cache first
        cache = cls._load_cache()
        if unigene_id_str in cache:
            return cache[unigene_id_str]

        # Query mygene.info with unigene scope
        try:
            import mygene
            mg = mygene.MyGeneInfo()
            # Try querying by UniGene ID
            result = mg.query(unigene_id_str, scopes='unigene', species='human', size=1)
            if result.get('total', 0) > 0:
                symbol = result['hits'][0].get('symbol')
                if symbol:
                    cache[unigene_id_str] = symbol
                    cls._save_cache(cache)
                    logger.debug(f"UniGene {unigene_id_str} → {symbol}")
                    return symbol
        except Exception as e:
            logger.debug(f"UniGene query failed for {unigene_id_str}: {e}")

        # Not found in cache or API
        cache[unigene_id_str] = None
        cls._save_cache(cache)
        return None

    @classmethod
    def is_unigene_id(cls, gene_id: str) -> bool:
        """Check if this looks like a UniGene ID."""
        return cls.UNIGENE_PATTERN.match(str(gene_id).strip()) is not None


class ESTTranslator:
    """
    Translate EST (Expressed Sequence Tag) probes to HGNC symbols.

    EST probes appear on microarrays (GPL3427, GPL9603) as short DNA sequences
    or EST IDs. They require multi-step mapping via UniGene database.

    Key challenge: Each probe ID needs to be mapped to UniGene, then to symbol.
    This requires platform-specific annotation files.

    Methods:
    1. Direct HGNC annotation (if available)
    2. UniGene mapping (Hs.123456 → symbol)
    3. mygene.info with fuzzy matching
    4. Local platform annotation cache
    """

    # EST array platform patterns
    # TIGR format: Hs_20001 (human), Mm_12345 (mouse), Rn_99999 (rat), etc.
    # Only match known organism codes to avoid matching RefSeq (NM_, NP_, etc.)
    TIGR_40K_PATTERN = re.compile(r'^(Hs|Mm|Rn|Ss|Bt|Cf|Gg|Xl|Dr|Ce|Dm|At)[_\.]\d+$', re.IGNORECASE)
    # EST_ACCESSION_PATTERN not used (too broad)

    @classmethod
    def can_handle(cls, gene_id: str) -> bool:
        """Check if this looks like an EST probe ID."""
        gene_id_str = str(gene_id).strip()
        # Check if it matches TIGR format or is already a UniGene ID
        # NOTE: We don't check EST_ACCESSION_PATTERN as it's too broad
        # and can incorrectly match RefSeq accessions
        return (
            cls.TIGR_40K_PATTERN.match(gene_id_str) is not None
            or UniGeneMapper.is_unigene_id(gene_id_str)
        )

    @classmethod
    def translate_one(
        cls,
        probe_id: str,
        platform: Optional[str] = None,
        annotation_file: Optional[str] = None
    ) -> Optional[str]:
        """
        Translate EST probe to HGNC symbol.

        Args:
            probe_id: EST probe ID (e.g., "Hs_20001", "TIGR:EST_ID", "Hs.123456")
            platform: Platform name for context (GPL3427, GPL9603, etc.)
            annotation_file: Path to platform-specific annotation file

        Returns:
            HGNC symbol or None if unmapped
        """
        probe_id_str = str(probe_id).strip()

        # Step 1: Check if it's already a UniGene ID
        if UniGeneMapper.is_unigene_id(probe_id_str):
            symbol = UniGeneMapper.unigene_to_symbol(probe_id_str)
            if symbol:
                return symbol

        # Step 2: Try to map EST ID → UniGene using mygene.info
        try:
            import mygene
            mg = mygene.MyGeneInfo()

            # Try querying EST ID directly
            result = mg.query(probe_id_str, scopes='refseq,symbol', species='human', size=1)
            if result.get('total', 0) > 0:
                symbol = result['hits'][0].get('symbol')
                if symbol:
                    logger.debug(f"EST {probe_id_str} → {symbol} (via mygene direct)")
                    return symbol
        except Exception as e:
            logger.debug(f"mygene EST query failed for {probe_id_str}: {e}")

        # Step 3: Try platform-specific annotation if provided
        if annotation_file:
            symbol = cls._query_annotation_file(probe_id_str, annotation_file)
            if symbol:
                logger.debug(f"EST {probe_id_str} → {symbol} (via annotation file)")
                return symbol

        # Step 4: Try fuzzy matching via mygene
        try:
            import mygene
            mg = mygene.MyGeneInfo()

            # Fuzzy matching on the EST ID as a query
            result = mg.querymany([probe_id_str], scopes='symbol', species='human', size=1)
            for hit in result:
                if 'symbol' in hit:
                    symbol = hit['symbol']
                    if symbol:
                        logger.debug(f"EST {probe_id_str} → {symbol} (via mygene fuzzy)")
                        return symbol
        except Exception as e:
            logger.debug(f"mygene fuzzy match failed for {probe_id_str}: {e}")

        return None

    @staticmethod
    def _query_annotation_file(probe_id: str, annotation_file: str) -> Optional[str]:
        """
        Query platform-specific annotation file for EST probe mapping.

        Expected format (TSV/CSV):
        probe_id    gene_symbol
        Hs_20001    TP53
        Hs_20002    BRCA1
        ...

        Returns:
            HGNC symbol or None if not found
        """
        try:
            import pandas as pd

            # Try loading the annotation file
            df = pd.read_csv(annotation_file, sep='\t', index_col=0)

            if probe_id in df.index:
                symbol = df.loc[probe_id, 'gene_symbol']
                if pd.notna(symbol):
                    return str(symbol).strip()

        except Exception as e:
            logger.debug(f"Failed to query annotation file {annotation_file}: {e}")

        return None

    @classmethod
    def translate_batch(
        cls,
        probe_ids: List[str],
        platform: Optional[str] = None,
        annotation_file: Optional[str] = None
    ) -> Dict[str, Optional[str]]:
        """
        Translate multiple EST probes efficiently.

        Batches UniGene lookups together for efficiency.
        """
        results = {}

        # Separate by likely UniGene IDs vs raw EST IDs
        unigene_ids = []
        raw_ids = []

        for probe_id in probe_ids:
            probe_id_str = str(probe_id).strip()
            if UniGeneMapper.is_unigene_id(probe_id_str):
                unigene_ids.append(probe_id_str)
            else:
                raw_ids.append(probe_id_str)

        # Batch translate UniGene IDs
        for unigene_id in unigene_ids:
            symbol = UniGeneMapper.unigene_to_symbol(unigene_id)
            if symbol:
                results[unigene_id] = symbol

        # Translate raw EST IDs using batched mygene queries (not one-at-a-time)
        if raw_ids:
            # Try batch mygene query first for all raw EST IDs
            try:
                import mygene
                mg = mygene.MyGeneInfo()

                # Batch query all raw EST IDs together
                batch_results = mg.querymany(raw_ids, scopes='refseq,symbol', species='human', verbose=False)
                batch_mapped = 0
                for hit in batch_results:
                    est_id = hit.get('query')
                    symbol = hit.get('symbol')
                    if est_id and symbol and est_id not in results:
                        results[est_id] = symbol
                        batch_mapped += 1
                logger.debug(f"Batch mygene query mapped {batch_mapped}/{len(raw_ids)} EST IDs")

                # For unmapped ones, try annotation file
                unmapped_ids = [eid for eid in raw_ids if eid not in results]
                if unmapped_ids and annotation_file:
                    for est_id in unmapped_ids:
                        symbol = cls._query_annotation_file(est_id, annotation_file)
                        if symbol:
                            results[est_id] = symbol

            except Exception as e:
                logger.warning(f"Batch mygene EST query failed: {e}. Falling back to per-ID lookup...")
                # Fallback: use translate_one only for unmapped IDs
                for est_id in raw_ids:
                    if est_id not in results:
                        symbol = cls.translate_one(est_id, platform, annotation_file)
                        if symbol:
                            results[est_id] = symbol

        return results


class PlatformAnnotationLoader:
    """Load and cache platform-specific annotation files for EST arrays."""

    KNOWN_PLATFORMS = {
        'GPL3427': {
            'name': 'TIGR Arabidopsis Genome Array',
            'type': 'EST',
            'description': 'EST-based microarray',
            'probe_type': 'EST',
        },
        'GPL9603': {
            'name': 'EST array',
            'type': 'EST',
            'description': 'Custom EST array',
            'probe_type': 'EST',
        },
    }

    @classmethod
    def is_est_platform(cls, gpl_id: str) -> bool:
        """Check if platform is EST-based."""
        return gpl_id in cls.KNOWN_PLATFORMS and cls.KNOWN_PLATFORMS[gpl_id]['type'] == 'EST'

    @classmethod
    def get_platform_info(cls, gpl_id: str) -> Optional[Dict]:
        """Get metadata about a platform."""
        return cls.KNOWN_PLATFORMS.get(gpl_id)
