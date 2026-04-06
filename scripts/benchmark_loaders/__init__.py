"""Modular benchmark data loading framework."""

from .errors import DataLoaderError
from .api import HGNCRestAPI, EnsemblXrefAPI
from .translation_cache import TranslationCache
from .gene_translation import GeneHarmonizer
from .translators import (
    TranslatorBase,
    EntrezGeneTranslator,
    GencodTranslator,
    RefSeqTranslator,
    detect_and_translate,
    batch_translate,
)
from .est_translator import (
    ESTTranslator,
    UniGeneMapper,
    PlatformAnnotationLoader,
)
from .loaders import (
    GDCDataLoader,
    GEODataLoader,
    SingleCellDataLoader,
    CustomDataLoader,
    RobustGEOLoader,
    BenchmarkDataLoaderFactory,
)

__all__ = [
    # Errors
    'DataLoaderError',
    # APIs
    'HGNCRestAPI',
    'EnsemblXrefAPI',
    # Caching
    'TranslationCache',
    # Gene Translation
    'GeneHarmonizer',
    # Translators (Phase 1)
    'TranslatorBase',
    'EntrezGeneTranslator',
    'GencodTranslator',
    'RefSeqTranslator',
    'ESTTranslator',
    'UniGeneMapper',
    'PlatformAnnotationLoader',
    'detect_and_translate',
    'batch_translate',
    # Data Loaders
    'GDCDataLoader',
    'GEODataLoader',
    'SingleCellDataLoader',
    'CustomDataLoader',
    'RobustGEOLoader',
    'BenchmarkDataLoaderFactory',
]
