"""
Data loaders for 47 benchmark datasets.

Handles downloading and parsing datasets from:
- GDC API (TCGA)
- GEO (Gene Expression Omnibus)
- CellxGene (single-cell)
- ArrayExpress
- BioConductor
- Custom sources (microbial, plant, etc.)

**DEPRECATED:** This module is maintained for backward compatibility.
New code should import directly from scripts.benchmark_loaders.
"""

import logging

logger = logging.getLogger(__name__)

# Retry settings (kept for backward compatibility)
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds

# Import all public classes from modular components (Phase 3.1 refactoring)
# Support both relative imports (when used as package) and direct imports (when added to sys.path)
try:
    from .benchmark_loaders import (
        DataLoaderError,
        HGNCRestAPI,
        EnsemblXrefAPI,
        GeneHarmonizer,
    )
    from .benchmark_loaders.loaders import (
        GDCDataLoader,
        GEODataLoader,
        SingleCellDataLoader,
        CustomDataLoader,
        RobustGEOLoader,
        BenchmarkDataLoaderFactory,
    )
except ImportError:
    # Fallback for direct imports (when module is imported via sys.path)
    from benchmark_loaders import (
        DataLoaderError,
        HGNCRestAPI,
        EnsemblXrefAPI,
        GeneHarmonizer,
    )
    from benchmark_loaders.loaders import (
        GDCDataLoader,
        GEODataLoader,
        SingleCellDataLoader,
        CustomDataLoader,
        RobustGEOLoader,
        BenchmarkDataLoaderFactory,
    )

# Re-export for backward compatibility
__all__ = [
    'DataLoaderError',
    'HGNCRestAPI',
    'EnsemblXrefAPI',
    'GeneHarmonizer',
    'GDCDataLoader',
    'GEODataLoader',
    'SingleCellDataLoader',
    'CustomDataLoader',
    'RobustGEOLoader',
    'BenchmarkDataLoaderFactory',
    'MAX_RETRIES',
    'RETRY_DELAY',
]
