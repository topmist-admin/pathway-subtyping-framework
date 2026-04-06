"""
Local GDC TCGA loader for manually downloaded datasets.

Use this instead of GDCDataLoader when files are downloaded manually to avoid API timeouts.

Directory structure expected:
  tcga_data/
    TCGA-BRCA/
      *.tsv files
    TCGA-LUAD/
      *.tsv files
    ... etc
"""

import logging
import os
import glob
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple

# Import with fallback for different import contexts
try:
    from .errors import DataLoaderError
except ImportError:
    from errors import DataLoaderError

logger = logging.getLogger(__name__)


class LocalGDCDataLoader:
    """Load TCGA expression data from locally downloaded files."""

    def __init__(self, tcga_data_dir: str = "tcga_data"):
        """
        Initialize loader with path to locally downloaded TCGA data.

        Args:
            tcga_data_dir: Path to directory containing TCGA-* subdirectories
        """
        self.base_dir = Path(tcga_data_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"TCGA data directory not found: {tcga_data_dir}")

    @classmethod
    def load(cls, project_id: str, data_dir: str = "tcga_data") -> Tuple[np.ndarray, np.ndarray]:
        """
        Load TCGA expression matrix from locally downloaded files.

        Args:
            project_id: TCGA project ID (e.g., 'TCGA-BRCA')
            data_dir: Path to directory containing downloaded TCGA data

        Returns:
            (expression_matrix, phenotype_labels)
        """
        loader = cls(data_dir)
        return loader._load_project(project_id)

    def _load_project(self, project_id: str) -> Tuple[np.ndarray, np.ndarray]:
        """Load all expression files for a project."""
        project_dir = self.base_dir / project_id

        if not project_dir.exists():
            raise DataLoaderError(f"Project directory not found: {project_dir}")

        # Find all TSV files
        tsv_files = sorted(glob.glob(str(project_dir / "*.tsv")))
        if not tsv_files:
            raise DataLoaderError(f"No TSV files found in {project_dir}")

        logger.info(f"Loading {project_id} from local files...")
        logger.info(f"  Found {len(tsv_files)} TSV files")

        # Parse all files
        expressions = []
        failed_count = 0

        for i, tsv_path in enumerate(tsv_files):
            if (i + 1) % 50 == 0:
                logger.debug(f"  Processed {i+1}/{len(tsv_files)} files")

            try:
                df = self._parse_tsv_file(tsv_path)
                if df is not None:
                    expressions.append(df)
                else:
                    failed_count += 1
            except Exception as e:
                logger.warning(f"  Failed to parse {Path(tsv_path).name}: {str(e)[:80]}")
                failed_count += 1
                continue

        if not expressions:
            raise DataLoaderError(f"Could not parse any expression files from {project_id}")

        logger.info(f"  Successfully parsed {len(expressions)} files ({failed_count} failed)")

        # Aggregate expression matrix
        logger.info(f"  Aggregating expression matrix...")
        expression_matrix = pd.concat(expressions, axis=1).fillna(0).groupby(level=0).sum()

        # Log2 transform
        expression_matrix = np.log2(expression_matrix + 1)

        # Transpose to (samples, genes) format
        expression_matrix = expression_matrix.T

        # Create phenotype labels (simple: based on project)
        # In production, would use actual clinical metadata
        n_samples = expression_matrix.shape[0]
        phenotype_labels = np.zeros(n_samples, dtype=int)

        logger.info(f"  Expression matrix: {expression_matrix.shape}")
        logger.info(f"  Samples: {n_samples}, Genes: {expression_matrix.shape[1]}")

        return expression_matrix.values, phenotype_labels

    def _parse_tsv_file(self, tsv_path: str) -> pd.DataFrame:
        """
        Parse a TCGA TSV expression file.

        TCGA augmented_star files have format:
        # gene-model: GENCODE v36
        gene_id | gene_name | gene_type | unstranded | stranded_first | stranded_second | ...

        We extract: gene_id (column 0) and unstranded counts (column 3)
        Skip metadata rows (N_unmapped, N_multimapping, etc)
        """
        try:
            # Read file, skip comment lines, skip metadata rows
            df = pd.read_csv(tsv_path, sep='\t', index_col=0, header=0, comment='#', skip_blank_lines=True)

            # Filter out metadata rows (start with 'N_')
            df = df[~df.index.str.startswith('N_')]

            # Extract unstranded column (index 3 = column 4 in 1-indexed, but column 3 after index_col=0)
            if 'unstranded' in df.columns:
                expr_df = df[['unstranded']].copy()
            else:
                # Fallback: take first numeric column
                numeric_cols = df.select_dtypes(include=['number'])
                if numeric_cols.empty:
                    return None
                expr_df = numeric_cols.iloc[:, [0]].copy()

            file_id = Path(tsv_path).stem
            expr_df.columns = [file_id]
            expr_df.index.name = 'gene_id'

            return expr_df

        except Exception as e:
            logger.debug(f"  Could not parse {Path(tsv_path).name}: {e}")
            return None


# Convenience function for use in benchmark loaders
def load_tcga_project(project_id: str, data_dir: str = "tcga_data") -> Tuple[np.ndarray, np.ndarray]:
    """Load TCGA project from local files."""
    return LocalGDCDataLoader.load(project_id, data_dir)
