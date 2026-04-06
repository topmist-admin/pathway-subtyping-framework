"""
Data loaders for 47 benchmark datasets.

Handles downloading and parsing datasets from:
- GDC API (TCGA)
- GEO (Gene Expression Omnibus)
- CellxGene (single-cell)
- ArrayExpress
- BioConductor
- Custom sources (microbial, plant, etc.)
"""

import gzip
import json
import logging
import os
import re
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import URLError
from urllib.request import urlopen, urlretrieve

import numpy as np
import pandas as pd
import requests

from .errors import DataLoaderError
from .gene_translation import GeneHarmonizer
from .translators import detect_and_translate, batch_translate
from .est_translator import ESTTranslator

logger = logging.getLogger(__name__)

# Retry settings
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds


class GDCDataLoader:
    """Load TCGA datasets from NCI GDC API."""

    GDC_API_URL = "https://api.gdc.cancer.gov"
    PAGE_SIZE = 200

    @classmethod
    def load(cls, project_id: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load TCGA expression matrix and clinical data from GDC.

        Args:
            project_id: GDC project ID (e.g., 'TCGA-BRCA')

        Returns:
            (gene_expression, phenotype_labels)
        """
        logger.info(f"Loading {project_id} from GDC API...")

        try:
            # Get sample metadata
            logger.info(f"  Fetching sample metadata...")
            samples = cls._fetch_samples(project_id)
            if not samples:
                raise DataLoaderError(f"No samples found for {project_id}")

            sample_ids = [s['id'] for s in samples]
            logger.info(f"  Found {len(sample_ids)} samples")

            # Get expression files (RNA-seq)
            logger.info(f"  Fetching RNA-seq files...")
            files = cls._fetch_expression_files(sample_ids)
            if not files:
                raise DataLoaderError(f"No expression files found for {project_id}")

            logger.info(f"  Found {len(files)} expression files")

            # Download and parse expression matrix
            logger.info(f"  Downloading and parsing expression data...")
            expression_df = cls._download_expression_matrix(files)

            # Get phenotype labels (e.g., PAM50 subtype, diagnosis)
            phenotype_labels = cls._extract_phenotypes(samples, project_id)

            # Align samples
            common_samples = list(set(expression_df.index) & set(phenotype_labels.index))
            if not common_samples:
                raise DataLoaderError("No common samples between expression and phenotype data")

            expression_matrix = expression_df.loc[common_samples].values
            phenotypes = phenotype_labels.loc[common_samples].values

            logger.info(f"  Expression matrix: {expression_matrix.shape} (samples × genes)")
            logger.info(f"  Phenotypes: {len(phenotypes)} labels")

            return expression_matrix, phenotypes

        except Exception as e:
            raise DataLoaderError(f"Failed to load {project_id}: {e}")

    @classmethod
    def _fetch_samples(cls, project_id: str) -> list:
        """Fetch sample metadata from GDC."""
        params = {
            'expand': 'diagnoses,samples',
            'format': 'JSON',
            'size': cls.PAGE_SIZE,
            'filters': json.dumps({
                'op': '=',
                'content': {'field': 'project.project_id', 'value': project_id}
            })
        }

        url = f"{cls.GDC_API_URL}/cases?{cls._build_query_string(params)}"

        for attempt in range(MAX_RETRIES):
            try:
                with urlopen(url, timeout=60) as response:  # GDC API can be slow with SSL handshake
                    data = json.loads(response.read().decode('utf-8'))
                    cases = data.get('data', {}).get('hits', [])
                    samples = []
                    for case in cases:
                        for sample in case.get('samples', []):
                            samples.append({'id': sample['sample_id']})
                    return samples
            except URLError as e:
                if attempt < MAX_RETRIES - 1:
                    # Exponential backoff: 5s, 10s, 20s
                    backoff = RETRY_DELAY * (2 ** attempt)
                    logger.warning(f"  Retry {attempt+1}/{MAX_RETRIES} in {backoff}s: {e}")
                    time.sleep(backoff)
                else:
                    raise

        return []

    @classmethod
    def _fetch_expression_files(cls, sample_ids: list) -> list:
        """
        Fetch RNA-seq expression files for samples.

        Batches queries to avoid HTTP 414 (URI Too Long) errors.
        With 600+ samples, passing all sample IDs in a single filter
        creates a URL that exceeds browser/server limits (~8KB).

        Solution: Query in smaller batches (50 samples) with longer timeout.
        The GDC API can be slow with large queries, so we use:
        - Smaller batch size (50 instead of 100) for faster responses
        - Longer timeout (60s instead of 30s) for SSL handshake
        - Exponential backoff for retries
        """
        all_files = []
        batch_size = 50  # Smaller batches for more reliable responses

        # Batch sample IDs to avoid URI-too-long errors
        for batch_idx in range(0, len(sample_ids), batch_size):
            batch = sample_ids[batch_idx:batch_idx + batch_size]
            logger.debug(f"  Fetching files for batch {batch_idx//batch_size + 1}: {len(batch)} samples")

            params = {
                'format': 'JSON',
                'size': cls.PAGE_SIZE,
                'filters': json.dumps({
                    'op': 'and',
                    'content': [
                        {'op': 'in', 'content': {'field': 'cases.samples.sample_id', 'value': batch}},
                        {'op': '=', 'content': {'field': 'data_category', 'value': 'Transcriptome Profiling'}},
                        {'op': '=', 'content': {'field': 'data_type', 'value': 'Gene Expression Quantification'}},
                        # Note: workflow_type filter removed (TCGA uses "augmented_star_gene_counts" files)
                        # Different projects use different metadata structures, so we filter by filename instead
                    ]
                })
            }

            url = f"{cls.GDC_API_URL}/files?{cls._build_query_string(params)}"

            for attempt in range(MAX_RETRIES):
                try:
                    with urlopen(url, timeout=60) as response:  # GDC API can be slow with SSL handshake
                        data = json.loads(response.read().decode('utf-8'))
                        batch_files = data.get('data', {}).get('hits', [])
                        all_files.extend(batch_files)
                        logger.debug(f"    Found {len(batch_files)} files in this batch")
                        break  # Success, move to next batch
                except URLError as e:
                    if attempt < MAX_RETRIES - 1:
                        # Exponential backoff: 5s, 10s, 20s (helps with SSL handshake timeouts)
                        backoff = RETRY_DELAY * (2 ** attempt)
                        logger.warning(f"    Batch retry {attempt+1}/{MAX_RETRIES} in {backoff}s: {str(e)[:60]}")
                        time.sleep(backoff)
                    else:
                        logger.error(f"    Failed to fetch batch after {MAX_RETRIES} retries")
                        raise

        logger.info(f"  Total expression files found: {len(all_files)}")
        return all_files

    @classmethod
    def _download_expression_matrix(cls, files: list, max_files: int = 500) -> pd.DataFrame:
        """
        Download and aggregate expression files into matrix.

        Filters files by name to ensure valid gene expression data:
        - Includes: gene_counts, fpkm, rsem, tpm, isoforms
        - Excludes: isoform (single form, not aggregated), junctions, fusions, mirnas
        """
        expressions = []

        # Filter files to valid expression formats
        valid_files = []
        invalid_count = 0
        for f in files:
            file_name = f['file_name'].lower()
            # Accept gene-level counts/quantifications
            if any(x in file_name for x in ['gene_counts', 'fpkm', 'rsem', 'tpm']):
                # But skip isoform-level data
                if 'isoform' not in file_name:
                    valid_files.append(f)
                    continue
            # Also accept augmented_star files (TCGA standard)
            if 'augmented_star' in file_name and 'gene' in file_name:
                valid_files.append(f)
                continue
            invalid_count += 1

        if invalid_count > 0:
            logger.debug(f"  Filtered out {invalid_count} non-gene-level files")

        for i, file_info in enumerate(valid_files[:max_files]):
            file_id = file_info['id']
            file_name = file_info['file_name']

            try:
                logger.debug(f"    Downloading file {i+1}/{min(len(valid_files), max_files)}: {file_name}")

                # Download file
                download_url = f"{cls.GDC_API_URL}/data/{file_id}"
                with tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp:
                    urlretrieve(download_url, tmp.name)
                    tmp_path = tmp.name

                # Parse TSV
                # TCGA augmented_star files have: gene_id | gene_name | unstranded | stranded_first | stranded_second | ...
                # We only need gene_id (column 0) and unstranded counts (column 2)
                try:
                    df = pd.read_csv(tmp_path, sep='\t', index_col=0, usecols=[0, 2], header=0)
                    df.index.name = 'gene_id'
                    # Name the sample column after the file ID for uniqueness
                    df.columns = [file_id]
                    expressions.append(df)
                except (ValueError, IndexError):
                    # Fallback: try reading just first and last column (for different file formats)
                    df = pd.read_csv(tmp_path, sep='\t', index_col=0, header=0)
                    # Take only the numeric columns (skip gene_name and other metadata)
                    numeric_cols = df.select_dtypes(include=['number']).iloc[:, 0]
                    expr_df = pd.DataFrame(numeric_cols)
                    expr_df.index.name = 'gene_id'
                    expr_df.columns = [file_id]
                    expressions.append(expr_df)

                os.unlink(tmp_path)

            except Exception as e:
                logger.warning(f"    Failed to parse {file_name}: {e}")
                continue

        if not expressions:
            raise DataLoaderError("No expression files could be parsed")

        # Aggregate: sum counts across files (handles duplicates)
        expression_matrix = pd.concat(expressions, axis=1).fillna(0).groupby(level=0).sum()

        # Log2 transform
        expression_matrix = np.log2(expression_matrix + 1)

        return expression_matrix

    @classmethod
    def _extract_phenotypes(cls, samples: list, project_id: str) -> pd.Series:
        """Extract phenotype labels from samples."""
        # For now, use simple mapping based on project_id
        # In production, would fetch actual phenotype data from GDC clinical module

        phenotype_map = {
            'TCGA-BRCA': {'Luminal A': 0, 'Luminal B': 1, 'HER2+': 2, 'Basal': 3},
            'TCGA-LUAD': {'Adenocarcinoma': 0, 'Squamous': 1},
            'TCGA-OV': {'Immunoreactive': 0, 'Differentiated': 1, 'Mesenchymal': 2, 'Proliferative': 3},
            'TCGA-GBM': {'Classical': 0, 'Mesenchymal': 1, 'Proneural': 2},
        }

        if project_id not in phenotype_map:
            # Default: random assignment for testing
            logger.warning(f"  No phenotype mapping for {project_id}; using random labels")
            n_samples = len(samples)
            return pd.Series(np.random.randint(0, 2, n_samples))

        # Assign labels based on mapping
        labels_dict = phenotype_map[project_id]
        n_samples = len(samples)
        n_classes = len(labels_dict)

        # For demonstration: assign proportionally
        labels = np.random.choice(list(labels_dict.values()), n_samples, p=[1/n_classes]*n_classes)

        sample_ids = [s['id'] for s in samples]
        return pd.Series(labels, index=sample_ids)

    @staticmethod
    def _build_query_string(params: dict) -> str:
        """Build URL query string from parameter dict."""
        from urllib.parse import urlencode
        return urlencode(params)


class GEODataLoader:
    """Load GEO datasets (microarray and RNA-seq)."""

    @classmethod
    def load(cls, accession_id: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load GEO dataset using GEOparse.

        Args:
            accession_id: GEO accession (e.g., 'GSE2109')

        Returns:
            (gene_expression, phenotype_labels)
        """
        logger.info(f"Loading {accession_id} from GEO...")

        try:
            import GEOparse
        except ImportError:
            raise DataLoaderError("GEOparse not installed. Install with: pip install GEOparse")

        try:
            logger.info(f"  Downloading from GEO (this may take a minute)...")
            os.environ['GEOPARSE_USE_HTTP_FOR_FTP'] = 'yes'  # Use HTTP instead of FTP for GEO

            gse = GEOparse.get_GEO(geo=accession_id, destdir=tempfile.gettempdir(), silent=True)

            # Extract expression matrix from GSM samples
            logger.info(f"  Extracting expression matrix from {len(gse.gsms)} samples...")

            # Get list of samples and their expression data
            expression_list = []
            sample_names = []
            all_genes = set()

            # First pass: collect all gene/probe IDs
            for sample_name, gsm in gse.gsms.items():
                sample_names.append(sample_name)
                gene_ids = gsm.table.index.tolist()
                all_genes.update(gene_ids)

            all_genes = sorted(list(all_genes))
            logger.info(f"  Found {len(all_genes)} unique genes across all samples")

            # Second pass: build aligned expression matrix
            expression_arrays = []
            for sample_name, gsm in gse.gsms.items():
                # Create aligned expression vector for this sample
                sample_expr = pd.Series(index=all_genes, dtype=float)

                # Robust value column detection
                value_col = cls._detect_value_column(gsm.table)
                if value_col is None:
                    logger.warning(f"    Could not find value column in {sample_name}; skipping sample")
                    # Fill with NaN for this sample, will be dropped later
                    sample_expr[:] = np.nan
                else:
                    sample_expr[gsm.table.index] = gsm.table[value_col].astype(float)

                sample_expr = sample_expr.fillna(0)
                expression_arrays.append(sample_expr.values)

            # Convert to dataframe (genes × samples)
            expression_df = pd.DataFrame(
                np.array(expression_arrays).T,
                index=all_genes,
                columns=sample_names
            )

            if expression_df is None or expression_df.empty:
                raise DataLoaderError(f"Could not extract expression matrix from {accession_id}")

            # Remove rows with all zeros and all-NaN rows
            expression_df = expression_df.loc[(expression_df != 0).any(axis=1)]
            expression_df = expression_df.dropna(how='all')

            # Impute remaining NaN values (from samples that failed to load)
            # Use column-wise median imputation (preserves sample structure)
            logger.debug(f"  Imputing NaN values...")
            for col in expression_df.columns:
                if expression_df[col].isna().any():
                    median_val = expression_df[col].median()
                    if pd.isna(median_val):
                        # If entire column is NaN, use global median
                        valid_vals = expression_df.values[~np.isnan(expression_df.values)]
                        median_val = np.nanmedian(valid_vals) if len(valid_vals) > 0 else 0
                    expression_df[col].fillna(median_val, inplace=True)

            # Log2 transform if needed (check if already log-transformed)
            max_val = expression_df.max().max()
            if max_val > 50:
                logger.info(f"  Log2 transforming (max value: {max_val:.1f})...")
                expression_df = np.log2(expression_df + 1)

            # Extract phenotypes from GSM metadata
            logger.info(f"  Extracting phenotype labels...")
            phenotypes = cls._extract_phenotypes_from_gse(gse, accession_id)

            expression_matrix = expression_df.values
            logger.info(f"  Expression matrix: {expression_matrix.shape}")

            return expression_matrix, phenotypes

        except Exception as e:
            raise DataLoaderError(f"Failed to load {accession_id}: {e}")

    @staticmethod
    def _detect_value_column(table_df: pd.DataFrame) -> Optional[str]:
        """
        Robustly detect expression value column.

        Tries in priority order: 'VALUE', 'value', 'VALUE (log2)',
        'expression', then first numeric column, then mean of numeric columns.

        Args:
            table_df: GEOparse sample table

        Returns:
            Column name if found, None if no suitable column exists
        """
        # Priority list of column names
        priority_names = ['VALUE', 'value', 'VALUE (log2)', 'expression', 'Expression']

        # Check priority names first
        for col in priority_names:
            if col in table_df.columns:
                return col

        # Find first numeric column
        numeric_cols = table_df.select_dtypes(include=[np.number]).columns.tolist()
        if numeric_cols:
            return numeric_cols[0]

        # Last resort: try to convert first column to numeric
        if len(table_df.columns) > 0:
            first_col = table_df.columns[0]
            try:
                pd.to_numeric(table_df[first_col], errors='coerce')
                non_nan_count = table_df[first_col].notna().sum()
                if non_nan_count > len(table_df) * 0.5:  # At least 50% valid
                    return first_col
            except (ValueError, TypeError):
                pass

        return None

    @staticmethod
    def _extract_phenotypes_from_gse(gse, accession_id: str) -> np.ndarray:
        """Extract phenotype labels from GSE metadata."""
        # Common phenotype column names (in priority order)
        phenotype_columns = [
            'characteristics_ch1', 'characteristics_ch2',  # Most common
            'disease state', 'diagnosis', 'disease', 'phenotype',
            'group', 'treatment', 'cell type', 'cell_type',
            'subtype', 'sample type', 'source name'
        ]

        # Try to extract from characteristics
        all_characteristics = {}
        for sample_name, gsm in gse.gsms.items():
            if 'characteristics_ch1' in gsm.metadata:
                chars = gsm.metadata['characteristics_ch1']
                if isinstance(chars, list):
                    for char_str in chars:
                        # Parse "key: value" format
                        if ':' in char_str:
                            key, value = char_str.split(':', 1)
                            key = key.strip().lower()
                            value = value.strip()
                            if key not in all_characteristics:
                                all_characteristics[key] = []
                            all_characteristics[key].append((sample_name, value))

        # Find the characteristic with most variance (likely phenotype)
        best_char = None
        best_variance = 0
        for key, values in all_characteristics.items():
            if len(values) == len(gse.gsms):  # Must have value for all samples
                unique_vals = len(set(v[1] for v in values))
                if unique_vals > best_variance and unique_vals < len(gse.gsms):  # Not all unique, not all same
                    best_char = key
                    best_variance = unique_vals

        if best_char:
            logger.debug(f"  Found phenotype characteristic: '{best_char}' ({best_variance} classes)")
            phenotype_map = {}
            phenotype_labels = []
            label_counter = 0

            for sample_name in gse.gsms.keys():
                for char_sample, value in all_characteristics[best_char]:
                    if char_sample == sample_name:
                        if value not in phenotype_map:
                            phenotype_map[value] = label_counter
                            label_counter += 1
                        phenotype_labels.append(phenotype_map[value])
                        break
                else:
                    phenotype_labels.append(0)  # Default

            if len(phenotype_labels) == len(gse.gsms):
                return np.array(phenotype_labels)

        # Fallback 1: Try to extract from title with pattern matching
        logger.warning(f"  Could not extract phenotype; trying sample titles with pattern matching...")
        phenotype_labels = GEODataLoader._extract_phenotypes_from_titles(gse)
        if phenotype_labels is not None and len(set(phenotype_labels)) > 1:
            return phenotype_labels

        # Fallback 2: Simple title-based assignment (unique titles → unique labels)
        logger.warning(f"  Pattern matching failed; using simple title-based assignment...")
        titles = []
        for gsm in gse.gsms.values():
            if 'title' in gsm.metadata:
                title = gsm.metadata['title']
                if isinstance(title, list):
                    titles.append(title[0])
                else:
                    titles.append(str(title))
            else:
                titles.append('unknown')

        unique_titles = list(set(titles))
        title_map = {title: i for i, title in enumerate(unique_titles)}
        return np.array([title_map[t] for t in titles])

    @staticmethod
    def _extract_phenotypes_from_titles(gse) -> Optional[np.ndarray]:
        """
        Extract phenotypes from GSM titles using pattern matching.

        Looks for common patterns like:
        - 'disease' vs 'control', 'case' vs 'control'
        - 'treated' vs 'untreated'
        - disease abbreviations (SCZ, ASD, CTL, PD, AD)
        - numbered cell types (CD4+, CD8+, type 1, type 2)

        Returns:
            Phenotype array if patterns found, None otherwise
        """
        phenotype_patterns = {
            'disease': ['disease', 'patient', 'affected', 'case', 'scz', 'asd', 'pd', 'ad', 'ms'],
            'control': ['control', 'ctrl', 'healthy', 'normal', 'ctl', 'unaffected'],
            'treated': ['treated', 'stimulated', 'activated'],
            'untreated': ['untreated', 'unstimulated', 'control'],
            'type1': ['type.1', 'type_1', 't1', 'type1', 'cd8', 'effector'],
            'type2': ['type.2', 'type_2', 't2', 'type2', 'cd4', 'regulatory', 'treg'],
        }

        titles = []
        for gsm in gse.gsms.values():
            if 'title' in gsm.metadata:
                title = gsm.metadata['title']
                if isinstance(title, list):
                    titles.append(title[0].lower())
                else:
                    titles.append(str(title).lower())
            else:
                titles.append('')

        # Try to classify titles into categories
        classifications = []
        for title in titles:
            category = None
            for cat, patterns in phenotype_patterns.items():
                if any(pattern in title for pattern in patterns):
                    category = cat
                    break

            if category:
                classifications.append(category)
            else:
                classifications.append(None)

        # Check if we got meaningful classifications
        unique_classes = set(c for c in classifications if c is not None)
        if len(unique_classes) > 1:
            # Map to numeric labels
            class_list = sorted(list(unique_classes))
            class_map = {c: i for i, c in enumerate(class_list)}
            class_map[None] = 0  # Default to 0 if unclassified

            return np.array([class_map[c] for c in classifications])

        return None


class SingleCellDataLoader:
    """Load single-cell RNA-seq datasets."""

    @classmethod
    def load_cellxgene(cls, dataset_id: str) -> Tuple[np.ndarray, np.ndarray]:
        """Load from CellxGene portal."""
        logger.info(f"Loading single-cell dataset {dataset_id} from CellxGene...")

        try:
            import cellxgene_census
        except ImportError:
            raise DataLoaderError("cellxgene-census not installed. Install with: pip install cellxgene-census")

        try:
            # Open CellxGene census
            with cellxgene_census.open_soma() as census:
                # This is a simplified version; actual implementation would be more complex
                logger.warning("  CellxGene loading: placeholder implementation")
                # Would need to query by dataset ID, aggregate by cell type, etc.
                raise DataLoaderError("CellxGene loader not yet fully implemented")

        except Exception as e:
            raise DataLoaderError(f"Failed to load from CellxGene: {e}")

    @classmethod
    def load_h5ad(cls, file_path: str) -> Tuple[np.ndarray, np.ndarray]:
        """Load from h5ad (AnnData) file."""
        logger.info(f"Loading single-cell dataset from {file_path}...")

        try:
            import scanpy as sc
        except ImportError:
            raise DataLoaderError("scanpy not installed. Install with: pip install scanpy")

        try:
            adata = sc.read_h5ad(file_path)

            # Get expression matrix (X) and cell type labels
            expression_matrix = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

            # Try to find cell type annotation
            cell_type_key = None
            for key in ['cell_type', 'celltype', 'cell_type_annotation', 'obs_names']:
                if key in adata.obs.columns:
                    cell_type_key = key
                    break

            if cell_type_key:
                cell_types = adata.obs[cell_type_key].values
                # Convert to numeric
                unique_types = list(set(cell_types))
                type_map = {t: i for i, t in enumerate(unique_types)}
                phenotypes = np.array([type_map[t] for t in cell_types])
            else:
                logger.warning("  No cell type annotation found; using random labels")
                phenotypes = np.random.randint(0, 5, len(adata))

            return expression_matrix, phenotypes

        except Exception as e:
            raise DataLoaderError(f"Failed to load h5ad file: {e}")


class CustomDataLoader:
    """Load custom/specialized datasets."""

    @classmethod
    def load_microbial_16s(cls, otu_table_path: str, taxonomy_path: str) -> Tuple[np.ndarray, np.ndarray]:
        """Load 16S rRNA OTU table (microbial community data)."""
        logger.info(f"Loading 16S OTU table from {otu_table_path}...")

        try:
            otu_df = pd.read_csv(otu_table_path, sep='\t', index_col=0)
            taxonomy_df = pd.read_csv(taxonomy_path, sep='\t', index_col=0)

            # OTU matrix: samples × OTUs
            otu_matrix = otu_df.values

            # Extract phyla from taxonomy
            phyla = taxonomy_df['Phylum'].values
            unique_phyla = list(set(phyla))
            phyla_map = {p: i for i, p in enumerate(unique_phyla)}
            phenotypes = np.array([phyla_map.get(p, 0) for p in phyla])

            logger.info(f"  OTU matrix: {otu_matrix.shape}")
            logger.info(f"  Unique phyla: {len(unique_phyla)}")

            return otu_matrix, phenotypes

        except Exception as e:
            raise DataLoaderError(f"Failed to load 16S data: {e}")

    @classmethod
    def load_generic_csv(cls, csv_path: str, phenotype_column: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Load generic CSV expression matrix."""
        logger.info(f"Loading expression matrix from {csv_path}...")

        try:
            df = pd.read_csv(csv_path, index_col=0)

            expression_matrix = df.values

            if phenotype_column and phenotype_column in df.columns:
                phenotypes = df[phenotype_column].values
                # Convert to numeric if needed
                if not np.issubdtype(phenotypes.dtype, np.number):
                    unique_vals = list(set(phenotypes))
                    val_map = {v: i for i, v in enumerate(unique_vals)}
                    phenotypes = np.array([val_map[v] for v in phenotypes])
            else:
                logger.warning("  No phenotype column found; using random labels")
                phenotypes = np.random.randint(0, 3, len(df))

            logger.info(f"  Expression matrix: {expression_matrix.shape}")

            return expression_matrix, phenotypes

        except Exception as e:
            raise DataLoaderError(f"Failed to load CSV: {e}")


class RobustGEOLoader:
    """
    Robust GEO dataset loader with:
    - GPL probe-to-symbol mapping (critical for gene overlap)
    - Flexible value column detection (RNA-seq, microarray, single-cell)
    - Supplementary file fallback (for single-cell with minimal SOFT tables)
    - Proper NaN imputation (fixes numpy ndarray.median() bug)
    """

    PRIORITY_VALUE_COLUMNS = [
        'VALUE', 'value', 'VALUE (log2)', 'VALUE_LOG2',
        'TPM', 'tpm', 'FPKM', 'fpkm', 'RPKM', 'rpkm',
        'count', 'Count', 'counts', 'Counts',
        'CPM', 'cpm', 'expression', 'Expression',
        'log2_expr', 'log2_count', 'normalized_count',
        'signal', 'Signal', 'intensity', 'Intensity',
    ]

    SYMBOL_COLUMNS = [
        # Tier 1: Standard and most common
        'Gene Symbol', 'Symbol', 'GENE_SYMBOL', 'Gene_Symbol', 'gene_symbol',
        'ILMN_Gene', 'GeneName', 'gene_name', 'GENE', 'gene', 'Name',
        # Tier 2: Affymetrix formats
        'gene_assignment', 'mrna_assignment', 'category',
        # Tier 3: Non-standard / legacy platforms (FIX C)
        'Gene title', 'Gene ID', 'Entrez_Gene_ID', 'ENTREZ_GENE_ID',
        'ORF', 'GB_LIST', 'GB_ACC', 'GENOME_ACC', 'SPOT_ID',
        'Description', 'transcript_id', 'RefSeq_ID',
    ]

    @classmethod
    def load(cls, accession_id: str) -> Tuple:
        """
        Load GEO dataset robustly with probe-to-symbol mapping.

        Args:
            accession_id: GEO accession (e.g., 'GSE2109')

        Returns:
            (gene_expression as DataFrame, phenotype_labels as ndarray)
            DataFrame has samples as rows, gene symbols as columns
        """
        logger.info(f"Loading {accession_id} from GEO (RobustGEOLoader)...")

        try:
            import GEOparse
        except ImportError:
            raise DataLoaderError("GEOparse not installed. Install with: pip install GEOparse")

        try:
            logger.info(f"  Downloading from GEO (this may take a minute)...")
            os.environ['GEOPARSE_USE_HTTP_FOR_FTP'] = 'yes'

            gse = GEOparse.get_GEO(geo=accession_id, destdir=tempfile.gettempdir(), silent=True)

            if not gse.gsms:
                raise DataLoaderError(f"No samples found in {accession_id}")

            logger.info(f"  Extracting expression matrix from {len(gse.gsms)} samples...")

            # Build probe → gene symbol mapping from GPL platform
            probe_to_symbol = cls._get_probe_to_symbol(gse)
            if probe_to_symbol:
                logger.info(f"  Probe-to-symbol mapping: {len(probe_to_symbol)} probes mapped")
            else:
                logger.warning(f"  Probe-to-symbol mapping: EMPTY (no GPL platform tables found or parsed)")

            # Build expression matrix (SOFT table with fallback to supplementary files)
            expression_df = cls._build_expression_matrix(gse, accession_id, probe_to_symbol)

            if expression_df is None or expression_df.empty:
                raise DataLoaderError(f"Could not extract expression matrix from {accession_id}")

            logger.info(f"  Expression matrix: {expression_df.shape[0]} samples × {expression_df.shape[1]} genes")

            # Extract phenotype labels
            phenotypes = RobustGEOLoader._extract_phenotypes_from_gse(gse, accession_id)

            # Return DataFrame (with gene symbols as columns) and phenotype labels
            return expression_df, phenotypes

        except Exception as e:
            raise DataLoaderError(f"Failed to load {accession_id}: {e}")

    @staticmethod
    def _extract_phenotypes_from_gse(gse, accession_id: str) -> np.ndarray:
        """Extract phenotype labels from GSE metadata (copied from GEODataLoader to avoid cls issues)."""
        # Common phenotype column names (in priority order)
        phenotype_columns = [
            'characteristics_ch1', 'characteristics_ch2',  # Most common
            'disease state', 'diagnosis', 'disease', 'phenotype',
            'group', 'treatment', 'cell type', 'cell_type',
            'subtype', 'sample type', 'source name'
        ]

        # Try to extract from characteristics
        all_characteristics = {}
        for sample_name, gsm in gse.gsms.items():
            if 'characteristics_ch1' in gsm.metadata:
                chars = gsm.metadata['characteristics_ch1']
                if isinstance(chars, list):
                    for char_str in chars:
                        # Parse "key: value" format
                        if ':' in char_str:
                            key, value = char_str.split(':', 1)
                            key = key.strip().lower()
                            value = value.strip()
                            if key not in all_characteristics:
                                all_characteristics[key] = []
                            all_characteristics[key].append((sample_name, value))

        # Find the characteristic with most variance (likely phenotype)
        best_char = None
        best_variance = 0
        for key, values in all_characteristics.items():
            if len(values) == len(gse.gsms):  # Must have value for all samples
                unique_vals = len(set(v[1] for v in values))
                if unique_vals > best_variance and unique_vals < len(gse.gsms):  # Not all unique, not all same
                    best_char = key
                    best_variance = unique_vals

        if best_char:
            logger.debug(f"  Found phenotype characteristic: '{best_char}' ({best_variance} classes)")
            phenotype_map = {}
            phenotype_labels = []
            label_counter = 0

            for sample_name in gse.gsms.keys():
                for char_sample, value in all_characteristics[best_char]:
                    if char_sample == sample_name:
                        if value not in phenotype_map:
                            phenotype_map[value] = label_counter
                            label_counter += 1
                        phenotype_labels.append(phenotype_map[value])
                        break
                else:
                    phenotype_labels.append(0)  # Default

            if len(phenotype_labels) == len(gse.gsms):
                return np.array(phenotype_labels)

        # Fallback 1: Try to extract from title with pattern matching
        logger.warning(f"  Could not extract phenotype; trying sample titles with pattern matching...")
        phenotype_labels = RobustGEOLoader._extract_phenotypes_from_titles(gse)
        if phenotype_labels is not None and len(set(phenotype_labels)) > 1:
            return phenotype_labels

        # Fallback 2: Simple title-based assignment (unique titles → unique labels)
        logger.warning(f"  Pattern matching failed; using simple title-based assignment...")
        titles = []
        for gsm in gse.gsms.values():
            if 'title' in gsm.metadata:
                title = gsm.metadata['title']
                if isinstance(title, list):
                    titles.append(title[0])
                else:
                    titles.append(str(title))
            else:
                titles.append('unknown')

        unique_titles = list(set(titles))
        title_map = {title: i for i, title in enumerate(unique_titles)}
        return np.array([title_map[t] for t in titles])

    @staticmethod
    def _extract_phenotypes_from_titles(gse) -> Optional[np.ndarray]:
        """
        Extract phenotypes from GSM titles using pattern matching.

        Looks for common patterns like:
        - 'disease' vs 'control', 'case' vs 'control'
        - 'treated' vs 'untreated'
        - disease abbreviations (SCZ, ASD, CTL, PD, AD)
        - numbered cell types (CD4+, CD8+, type 1, type 2)

        Returns:
            Phenotype array if patterns found, None otherwise
        """
        phenotype_patterns = {
            'disease': ['disease', 'patient', 'affected', 'case', 'scz', 'asd', 'pd', 'ad', 'ms'],
            'control': ['control', 'ctrl', 'healthy', 'normal', 'ctl', 'unaffected'],
            'treated': ['treated', 'stimulated', 'activated'],
            'untreated': ['untreated', 'unstimulated', 'control'],
            'type1': ['type.1', 'type_1', 't1', 'type1', 'cd8', 'effector'],
            'type2': ['type.2', 'type_2', 't2', 'type2', 'cd4', 'regulatory', 'treg'],
        }

        titles = []
        for gsm in gse.gsms.values():
            if 'title' in gsm.metadata:
                title = gsm.metadata['title']
                if isinstance(title, list):
                    titles.append(title[0].lower())
                else:
                    titles.append(str(title).lower())
            else:
                titles.append('')

        # Try to classify titles into categories
        classifications = []
        for title in titles:
            category = None
            for cat, patterns in phenotype_patterns.items():
                if any(pattern in title for pattern in patterns):
                    category = cat
                    break

            if category:
                classifications.append(category)
            else:
                classifications.append(None)

        # Check if we got meaningful classifications
        unique_classes = set(c for c in classifications if c is not None)
        if len(unique_classes) > 1:
            # Map to numeric labels
            class_list = sorted(list(unique_classes))
            class_map = {c: i for i, c in enumerate(class_list)}
            class_map[None] = 0  # Default to 0 if unclassified

            return np.array([class_map[c] for c in classifications])

        return None

    @classmethod
    def _get_probe_to_symbol(cls, gse) -> dict:
        """
        Build probe ID → gene symbol mapping from GPL platform tables.

        Handles multiple GPL formats:
        - Standard 'Gene Symbol' column
        - Affymetrix 'gene_assignment' with ' // ' delimiters
        - Multi-gene probes with ' /// ' delimiter (take first gene only)

        Returns:
            {probe_id: gene_symbol} dict
        """
        probe_to_symbol = {}

        if not gse.gpls:
            logger.warning("  No GPL platform tables found; probe-to-symbol mapping unavailable")
            return probe_to_symbol

        for gpl_id, gpl in gse.gpls.items():
            if not hasattr(gpl, 'table') or gpl.table is None or gpl.table.empty:
                logger.debug(f"    GPL {gpl_id} has no table; skipping")
                continue

            logger.debug(f"    Processing GPL {gpl_id}...")
            table = gpl.table

            # FIX 3: Find symbol column with case-insensitive search
            symbol_col = None

            # Exact case match first (preferred)
            for col_name in cls.SYMBOL_COLUMNS:
                if col_name in table.columns:
                    symbol_col = col_name
                    logger.debug(f"      Found symbol column (exact match): {col_name}")
                    break

            # Fall back to case-insensitive search
            if symbol_col is None:
                col_lower = [c.lower() for c in table.columns]
                for target_col in cls.SYMBOL_COLUMNS:
                    target_lower = target_col.lower()
                    if target_lower in col_lower:
                        symbol_col = table.columns[col_lower.index(target_lower)]
                        logger.debug(f"      Found symbol column (case-insensitive): {symbol_col}")
                        break

            if symbol_col is None:
                logger.debug(f"    No symbol column found in GPL {gpl_id} (tried: {cls.SYMBOL_COLUMNS[:3]}...)")
                continue

            # Build mapping
            for probe_id in table.index:
                try:
                    symbol_raw = table.loc[probe_id, symbol_col]

                    if pd.isna(symbol_raw) or symbol_raw == '' or symbol_raw is None:
                        continue

                    symbol_raw = str(symbol_raw).strip()

                    # Handle Affymetrix gene_assignment format: "A1BG // 23 // ..." → take index 1
                    if ' // ' in symbol_raw:
                        parts = symbol_raw.split(' // ')
                        if len(parts) > 1:
                            symbol_raw = parts[1].strip()
                        else:
                            symbol_raw = parts[0].strip()

                    # Handle multi-gene probes (` /// ` separator) — take first gene only
                    if ' /// ' in symbol_raw:
                        symbol_raw = symbol_raw.split(' /// ')[0].strip()

                    # Skip if empty, numeric-only, or NaN
                    if not symbol_raw or symbol_raw.lower() == 'nan' or symbol_raw.isdigit():
                        continue

                    probe_to_symbol[str(probe_id)] = symbol_raw

                except (ValueError, KeyError, AttributeError):
                    continue

        # FIX 3 (IMPROVED): Always attempt translator enhancement as final validation step
        # Even if we have many mappings, translators may improve quality/coverage

        before_count = len(probe_to_symbol)
        before_genes = set(probe_to_symbol.values()) if probe_to_symbol else set()

        if len(probe_to_symbol) == 0:
            # No GPL mappings: try translators on probe IDs directly
            logger.warning(f"  WARNING: No probes mapped from GPL tables; using translators on probe IDs")
            probe_to_symbol = cls._enhance_probe_mapping_with_translators(gse, {})
        else:
            # Have some mappings: try to enhance/validate using translators
            logger.debug(f"  Built probe-to-symbol mapping: {len(probe_to_symbol)} probes mapped")

            # Always invoke translators to fill gaps and validate coverage
            enhanced = cls._enhance_probe_mapping_with_translators(gse, probe_to_symbol)

            # Only use enhanced result if it actually improved coverage
            after_count = len(enhanced)
            after_genes = set(enhanced.values()) if enhanced else set()

            if after_count > before_count or len(after_genes) > len(before_genes):
                coverage_improved = after_count - before_count
                gene_improved = len(after_genes) - len(before_genes)
                logger.info(f"  Translators improved mapping: +{coverage_improved} probes, +{gene_improved} unique genes")
                probe_to_symbol = enhanced
            else:
                logger.debug(f"  Translators didn't improve mapping ({before_count} → {after_count} probes)")
                # Keep original mapping

        return probe_to_symbol

    @classmethod
    def _enhance_probe_mapping_with_translators(cls, gse, probe_to_symbol: dict) -> dict:
        """
        Enhance probe-to-symbol mapping using Phase 1 gene translators.

        When GPL table doesn't provide complete annotations, try translating
        probe IDs using Entrez, RefSeq, EST, and Gencode translators.

        Args:
            gse: GEO Series object
            probe_to_symbol: Partial probe → symbol mapping from GPL tables

        Returns:
            Enhanced probe_to_symbol mapping
        """
        if not gse.gpls:
            return probe_to_symbol

        unmapped_probes = []

        # Collect unmapped probe IDs from GPL tables
        for gpl_id, gpl in gse.gpls.items():
            if not hasattr(gpl, 'table') or gpl.table is None or gpl.table.empty:
                continue

            for probe_id in gpl.table.index:
                probe_id_str = str(probe_id).strip()
                if probe_id_str not in probe_to_symbol:
                    unmapped_probes.append(probe_id_str)

        if not unmapped_probes:
            logger.debug(f"  No unmapped probes; skipping translator enhancement")
            return probe_to_symbol

        logger.debug(f"  Enhancing {len(unmapped_probes)} unmapped probes using translators...")

        # Try translators in priority order
        # EST first (solves GPL3427, GPL9603)
        est_probes = [p for p in unmapped_probes if ESTTranslator.can_handle(p)]
        if est_probes:
            logger.debug(f"    EST translator: {len(est_probes)} probes")
            est_results = ESTTranslator.translate_batch(est_probes)
            probe_to_symbol.update(est_results)

        # Generic translators for remaining unmapped
        remaining = [p for p in unmapped_probes if p not in probe_to_symbol]
        if remaining:
            logger.debug(f"    Generic translators: {len(remaining)} probes")
            translated = batch_translate(remaining)
            probe_to_symbol.update(translated)

        enhanced_count = len(probe_to_symbol) - len(set(probe_to_symbol) & set(['']))
        logger.debug(f"  Translator enhancement: {len(unmapped_probes)} unmapped → {enhanced_count} total mapped")

        return probe_to_symbol

    @classmethod
    def _detect_value_column(cls, table_df: pd.DataFrame) -> Optional[str]:
        """
        Robustly detect expression value column.

        Tries in priority order: standard microarray columns, then RNA-seq TPM/FPKM,
        then first numeric column.

        Args:
            table_df: GEOparse sample table

        Returns:
            Column name if found, None if no suitable column exists
        """
        # Check priority names first
        for col in cls.PRIORITY_VALUE_COLUMNS:
            if col in table_df.columns:
                return col

        # Find first numeric column
        numeric_cols = table_df.select_dtypes(include=[np.number]).columns.tolist()
        if numeric_cols:
            return numeric_cols[0]

        # Last resort: try to convert first column to numeric
        if len(table_df.columns) > 0:
            first_col = table_df.columns[0]
            try:
                pd.to_numeric(table_df[first_col], errors='coerce')
                non_nan_count = table_df[first_col].notna().sum()
                if non_nan_count > len(table_df) * 0.5:  # At least 50% valid
                    return first_col
            except (ValueError, TypeError):
                pass

        return None

    @classmethod
    def _try_soft_table(cls, gse, probe_to_symbol: dict) -> Optional[pd.DataFrame]:
        """
        Build expression DataFrame from SOFT table using probe-to-symbol mapping.

        Returns:
            DataFrame with samples as rows, gene symbols as columns (or None if failure)
        """
        try:
            sample_names = []
            expression_arrays = []

            for sample_name, gsm in gse.gsms.items():
                sample_names.append(sample_name)

                # Detect value column
                value_col = cls._detect_value_column(gsm.table)
                if value_col is None:
                    logger.debug(f"    {sample_name}: no value column found; skipping")
                    expression_arrays.append({})
                    continue

                # Build sample expression dict: gene_symbol → value
                sample_expr = {}
                for probe_id, value in gsm.table[value_col].items():
                    try:
                        # Map probe to symbol (fallback to probe ID if not found)
                        gene_symbol = probe_to_symbol.get(str(probe_id), str(probe_id))
                        sample_expr[gene_symbol] = float(value)
                    except (ValueError, TypeError):
                        pass

                expression_arrays.append(sample_expr)

            # Align all genes across samples
            all_genes = set()
            for expr_dict in expression_arrays:
                all_genes.update(expr_dict.keys())

            all_genes_before_filter = len(all_genes)

            # Filter out unmapped probes (numeric-only or empty)
            all_genes = sorted([g for g in all_genes if g and not g.isdigit() and not g.replace('.', '').isdigit()])

            logger.debug(f"    Gene symbols: {all_genes_before_filter} → {len(all_genes)} after filtering pure-numeric probes")

            if len(all_genes) < 10:
                logger.debug(f"    Too few gene symbols ({len(all_genes)}); likely mapping failure")
                return None

            # Build matrix: samples × genes
            matrix_data = []
            for expr_dict in expression_arrays:
                row = [expr_dict.get(gene, 0.0) for gene in all_genes]
                matrix_data.append(row)

            expression_df = pd.DataFrame(matrix_data, columns=all_genes, index=sample_names)

            return expression_df

        except Exception as e:
            logger.debug(f"    _try_soft_table failed: {e}")
            return None

    @classmethod
    def _download_supplementary_files_direct(cls, gse, accession_id: str, dest_dir: str) -> int:
        """
        Download supplementary files directly from GEO metadata URLs.

        GEOparse's download_supplementary_files() often fails silently.
        This method uses the supplementary_files list from metadata directly.

        Returns:
            Number of files downloaded successfully
        """
        count = 0

        # Extract supplementary file URLs from metadata
        if not hasattr(gse, 'metadata') or 'supplementary_file' not in gse.metadata:
            return 0

        supp_urls = gse.metadata.get('supplementary_file', [])
        if isinstance(supp_urls, str):
            supp_urls = [supp_urls]

        for url in supp_urls:
            try:
                # Convert FTP URLs to HTTPS (keeping ftp.ncbi subdomain which works over HTTPS)
                http_url = url.replace('ftp://', 'https://')
                filename = url.split('/')[-1]
                dest_path = Path(dest_dir) / filename

                # Skip if already exists
                if dest_path.exists():
                    logger.debug(f"      Supplementary file already exists: {filename}")
                    count += 1
                    continue

                # Try to download (with timeout)
                logger.debug(f"      Downloading supplementary file: {filename}...")
                response = requests.get(http_url, timeout=30, allow_redirects=True)
                if response.status_code == 200:
                    with open(dest_path, 'wb') as f:
                        f.write(response.content)
                    logger.debug(f"        Downloaded: {filename}")
                    count += 1
                else:
                    logger.debug(f"        HTTP {response.status_code}: {filename}")
            except Exception as e:
                logger.debug(f"      Failed to download {filename}: {str(e)[:100]}")

        return count

    @classmethod
    def _try_supplementary(cls, gse, accession_id: str) -> Optional[pd.DataFrame]:
        """
        Fallback to GSE supplementary files for single-cell datasets.

        Tries CSV, TSV, and MTX (10x format) files.
        Searches multiple possible download locations.

        Returns:
            DataFrame with samples as rows, genes as columns (or None if failure)
        """
        try:
            logger.debug(f"    Attempting supplementary file fallback...")

            # Create temporary directory for supplementary files
            gse_supp_dir = Path(tempfile.gettempdir()) / f"Supp_{accession_id}"
            gse_supp_dir.mkdir(exist_ok=True)

            # First, try direct download from metadata URLs (more reliable than GEOparse's method)
            direct_count = cls._download_supplementary_files_direct(gse, accession_id, str(gse_supp_dir))
            if direct_count > 0:
                logger.debug(f"    Downloaded {direct_count} supplementary files directly from GEO")

            # Also attempt GEOparse download (as fallback)
            for attempt in range(3):
                try:
                    gse.download_supplementary_files(
                        directory=tempfile.gettempdir(),
                        download_sra=False
                    )
                    break
                except Exception as e:
                    if attempt < 2:
                        logger.debug(f"    GEOparse download attempt {attempt+1}/3 failed (non-critical)")
                        time.sleep(2)
                    else:
                        logger.debug(f"    GEOparse download failed after 3 attempts (continuing with direct downloads)")

            import glob
            import scipy.io
            import tarfile
            import zipfile

            # FIX 1: Try GSE-level supplementary directory FIRST (where direct downloads go)
            # This is more reliable than per-sample directories for bulk expression matrices
            if gse_supp_dir.exists() and any(gse_supp_dir.glob('*')):
                logger.debug(f"    Trying GSE-level supplementary directory: {gse_supp_dir.name}")

                # Extract archives first
                archive_files = list(gse_supp_dir.glob('*.tar')) + list(gse_supp_dir.glob('*.tar.gz')) + \
                               list(gse_supp_dir.glob('*.tgz')) + list(gse_supp_dir.glob('*.zip'))
                for archive_file in archive_files:
                    try:
                        logger.debug(f"      Extracting archive: {archive_file.name}...")
                        if str(archive_file).endswith(('.tar', '.tar.gz', '.tgz')):
                            with tarfile.open(archive_file, 'r:*') as tar:
                                tar.extractall(path=gse_supp_dir)
                        elif str(archive_file).endswith('.zip'):
                            with zipfile.ZipFile(archive_file, 'r') as zf:
                                zf.extractall(path=gse_supp_dir)
                    except Exception as e:
                        logger.debug(f"      Failed to extract {archive_file.name}: {e}")

                # Try CSV/TSV files (most common for bulk expression matrices)
                for pattern in ['*.csv.gz', '*.csv', '*.tsv.gz', '*.tsv', '*.txt.gz', '*.txt']:
                    files = [f for f in gse_supp_dir.glob(pattern)
                            if 'barcodes' not in f.name and 'features' not in f.name and 'genes' not in f.name
                            and 'readme' not in f.name.lower()]
                    if files:
                        csv_file = files[0]
                        logger.debug(f"      Loading {csv_file.name} from GSE level...")
                        try:
                            # Try standard CSV/TSV loading first
                            if str(csv_file).endswith('.gz'):
                                df = pd.read_csv(csv_file, compression='gzip', index_col=0, sep='\t')
                            else:
                                df = pd.read_csv(csv_file, index_col=0, sep='\t')

                            # Extract only numeric columns (handles metadata columns)
                            numeric_cols = df.select_dtypes(include=['number']).columns
                            if len(numeric_cols) > 0:
                                df = df[numeric_cols]

                            # Transpose if genes are rows
                            if df.shape[0] > df.shape[1]:
                                df = df.T

                            # Validate data
                            if df.shape[0] >= 2 and df.shape[1] >= 5:
                                logger.debug(f"      Successfully loaded GSE-level file: {df.shape[0]} samples × {df.shape[1]} genes")
                                return df
                        except Exception as e1:
                            # Fallback: Try loading with error_bad_lines handling
                            try:
                                logger.debug(f"      Standard parsing failed, trying with on_bad_lines='skip'...")
                                if str(csv_file).endswith('.gz'):
                                    df = pd.read_csv(csv_file, compression='gzip', index_col=0, sep='\t',
                                                   on_bad_lines='skip', engine='python')
                                else:
                                    df = pd.read_csv(csv_file, index_col=0, sep='\t',
                                                   on_bad_lines='skip', engine='python')

                                # Extract only numeric columns (handles metadata columns)
                                numeric_cols = df.select_dtypes(include=['number']).columns
                                if len(numeric_cols) > 0:
                                    df = df[numeric_cols]

                                # Transpose if genes are rows
                                if df.shape[0] > df.shape[1]:
                                    df = df.T

                                # Validate data
                                if df.shape[0] >= 2 and df.shape[1] >= 5:
                                    logger.debug(f"      Successfully loaded {csv_file.name} (with error recovery): {df.shape[0]} samples × {df.shape[1]} genes")
                                    return df
                            except Exception as e2:
                                logger.debug(f"      Failed to load {csv_file.name}: {str(e1)[:80]} / {str(e2)[:80]}")
                            continue

            # If GSE-level files didn't work, fall back to per-sample directories
            # FIX 2: Search for supplementary sample directories created by GEOparse
            # GEOparse downloads per-sample files into directories like Supp_GSM2836573_Description/
            temp_path = Path(tempfile.gettempdir())

            # FIX 0: Scope GSM directories to ONLY those belonging to this study (prevent cross-contamination)
            gsm_ids = set(gse.gsms.keys())
            all_supp_dirs = list(temp_path.glob("Supp_GSM*"))
            supp_sample_dirs = [
                d for d in all_supp_dirs
                if any(d.name.startswith(f"Supp_{gsm_id}") for gsm_id in gsm_ids)
            ]

            if not supp_sample_dirs:
                logger.debug(f"    No supplementary sample directories found (Supp_GSM* pattern)")
                return None

            logger.debug(f"    Found {len(supp_sample_dirs)} supplementary sample directories")

            # Try to load MTX (10x format) from any sample directory
            for supp_dir in supp_sample_dirs:
                # Extract any archives (.tar, .tar.gz, .zip, .rar, .7z) in the directory first
                archive_files = list(supp_dir.glob('*.tar')) + list(supp_dir.glob('*.tar.gz')) + \
                               list(supp_dir.glob('*.tgz')) + list(supp_dir.glob('*.zip'))
                for archive_file in archive_files:
                    try:
                        logger.debug(f"    Extracting archive: {archive_file.name}...")
                        if str(archive_file).endswith(('.tar', '.tar.gz', '.tgz')):
                            with tarfile.open(archive_file, 'r:*') as tar:
                                tar.extractall(path=supp_dir)
                        elif str(archive_file).endswith('.zip'):
                            with zipfile.ZipFile(archive_file, 'r') as zf:
                                zf.extractall(path=supp_dir)
                    except Exception as e:
                        logger.debug(f"    Failed to extract {archive_file.name}: {e}")
                        continue
                if not supp_dir.is_dir():
                    continue

                # Try MTX (10x format) - most common for single-cell data
                mtx_files = list(supp_dir.glob('*.mtx.gz')) + list(supp_dir.glob('*.mtx'))
                if mtx_files:
                    mtx_file = mtx_files[0]
                    logger.debug(f"    Loading {mtx_file.name} (10x format)...")
                    try:
                        mtx = scipy.io.mmread(mtx_file).T.tocsr()  # Transpose: genes → samples

                        # Find barcodes and features files (may have sample-specific prefixes)
                        barcodes_files = list(supp_dir.glob('*barcodes*')) + list(supp_dir.glob('*barcode*'))
                        features_files = list(supp_dir.glob('*genes*')) + list(supp_dir.glob('*features*'))

                        barcodes_file = barcodes_files[0] if barcodes_files else None
                        features_file = features_files[0] if features_files else None

                        if barcodes_file and features_file:
                            barcodes = pd.read_csv(barcodes_file, compression='gzip' if str(barcodes_file).endswith('.gz') else None, header=None, sep='\t')
                            features = pd.read_csv(features_file, compression='gzip' if str(features_file).endswith('.gz') else None, header=None, sep='\t')

                            df = pd.DataFrame.sparse.from_spmatrix(mtx, index=barcodes[0].values, columns=features[1].values)
                            logger.debug(f"    Successfully loaded 10x MTX: {df.shape[0]} samples × {df.shape[1]} genes")
                            return df.astype(float)
                        else:
                            logger.debug(f"    MTX found but missing barcodes/features files")
                    except Exception as e:
                        logger.debug(f"    Failed to load MTX format: {e}")
                        continue

                # Try CSV/TSV files if MTX not available (but exclude barcodes/features files)
                for pattern in ['*.csv.gz', '*.csv', '*.tsv.gz', '*.tsv']:
                    files = [f for f in supp_dir.glob(pattern)
                            if 'barcodes' not in f.name and 'features' not in f.name and 'genes' not in f.name]
                    if files:
                        csv_file = files[0]
                        logger.debug(f"    Loading {csv_file.name}...")
                        try:
                            if str(csv_file).endswith('.gz'):
                                df = pd.read_csv(csv_file, compression='gzip', index_col=0)
                            else:
                                df = pd.read_csv(csv_file, index_col=0)

                            # Transpose if genes are rows
                            if df.shape[0] > df.shape[1]:
                                df = df.T

                            if df.shape[0] >= 2:  # Minimum 2 samples
                                logger.debug(f"    Successfully loaded CSV/TSV: {df.shape[0]} samples × {df.shape[1]} genes")
                                return df
                        except Exception as e:
                            logger.debug(f"    Failed to load {csv_file.name}: {e}")
                            continue

                # Try generic .gz compressed files (auto-detect TSV/CSV after decompression)
                # Skip non-expression file types (.gtf, .gff, .bam, .vcf, .CEL, .bedgraph, .sam, .fastq)
                for pattern in ['*.gz']:
                    files = [f for f in supp_dir.glob(pattern)
                            if 'barcodes' not in f.name and 'features' not in f.name and 'genes' not in f.name
                            and not any(f.name.lower().endswith(ext) for ext in ['.gtf.gz', '.gff.gz', '.bam.gz', '.vcf.gz', '.cel.gz', '.sam.gz', '.fastq.gz', '.bedgraph.gz'])]
                    if files:
                        gz_file = files[0]
                        logger.debug(f"    Loading {gz_file.name} (generic .gz format)...")
                        try:
                            import gzip
                            with gzip.open(gz_file, 'rt') as f:
                                # Try reading as TSV first (most common for expression matrices)
                                try:
                                    df = pd.read_csv(f, sep='\t', index_col=0)
                                except (pd.errors.ParserError, ValueError):
                                    # Fall back to CSV if TSV parsing fails
                                    f.seek(0)
                                    df = pd.read_csv(f, sep=',', index_col=0)

                            # Validate: ensure we have numeric data (convert to float, drop non-numeric)
                            if df.empty or df.shape[0] < 2 or df.shape[1] < 5:
                                logger.debug(f"    Loaded .gz has insufficient data: {df.shape}")
                                continue

                            # Try to convert to numeric; drop columns that can't be converted
                            df = df.apply(pd.to_numeric, errors='coerce')

                            # Drop rows with all NaN (likely header rows or metadata)
                            df = df.dropna(how='all')

                            # Drop columns with all NaN (failed conversions)
                            df = df.dropna(axis=1, how='all')

                            if df.empty or df.shape[0] < 2 or df.shape[1] < 5:
                                logger.debug(f"    After type conversion, .gz has insufficient numeric data: {df.shape}")
                                continue

                            # Transpose if genes are rows
                            if df.shape[0] > df.shape[1]:
                                df = df.T

                            logger.debug(f"    Successfully loaded .gz: {df.shape[0]} samples × {df.shape[1]} genes")
                            return df
                        except Exception as e:
                            logger.debug(f"    Failed to load {gz_file.name}: {e}")
                            continue

                # Try Excel files (.xlsx, .xls) if MTX/CSV/TSV not available
                for pattern in ['*.xlsx', '*.xls']:
                    files = [f for f in supp_dir.glob(pattern)
                            if 'barcodes' not in f.name and 'features' not in f.name and 'genes' not in f.name]
                    if files:
                        xl_file = files[0]
                        logger.debug(f"    Loading {xl_file.name} (Excel format)...")
                        try:
                            df = pd.read_excel(xl_file, index_col=0)
                            # Convert to numeric, coerce errors to NaN
                            df = df.apply(pd.to_numeric, errors='coerce')
                            # Transpose if genes are rows
                            if df.shape[0] > df.shape[1]:
                                df = df.T
                            if df.shape[0] >= 2:
                                logger.debug(f"    Successfully loaded Excel: {df.shape[0]} samples × {df.shape[1]} genes")
                                return df
                        except Exception as e:
                            logger.debug(f"    Failed to load {xl_file.name}: {e}")
                            continue

                # Try HDF5 files (.h5, .hdf5, .loom) if above formats not available
                for pattern in ['*.h5', '*.hdf5', '*.loom']:
                    files = [f for f in supp_dir.glob(pattern)
                            if 'barcodes' not in f.name and 'features' not in f.name]
                    if files:
                        h5_file = files[0]
                        logger.debug(f"    Loading {h5_file.name} (HDF5 format)...")
                        try:
                            # Try pandas HDF5 reader first
                            df = pd.read_hdf(h5_file)
                            if isinstance(df, pd.DataFrame):
                                df = df.apply(pd.to_numeric, errors='coerce')
                                if df.shape[0] > df.shape[1]:
                                    df = df.T
                                if df.shape[0] >= 2:
                                    logger.debug(f"    Successfully loaded HDF5: {df.shape[0]} samples × {df.shape[1]} genes")
                                    return df
                        except Exception:
                            # Fallback: use h5py to find largest 2D dataset
                            try:
                                import h5py
                                with h5py.File(h5_file, 'r') as hf:
                                    # Find largest 2D numeric dataset (likely expression matrix)
                                    best_key, best_size = None, 0
                                    def find_matrix(name, obj):
                                        nonlocal best_key, best_size
                                        if isinstance(obj, h5py.Dataset) and len(obj.shape) == 2:
                                            size = obj.shape[0] * obj.shape[1]
                                            if size > best_size:
                                                best_key, best_size = name, size
                                    hf.visititems(find_matrix)
                                    if best_key:
                                        data = hf[best_key][()]
                                        df = pd.DataFrame(data)
                                        df = df.apply(pd.to_numeric, errors='coerce')
                                        if df.shape[0] > df.shape[1]:
                                            df = df.T
                                        if df.shape[0] >= 2:
                                            logger.debug(f"    Successfully loaded HDF5 (h5py): {df.shape[0]} samples × {df.shape[1]} genes")
                                            return df
                            except Exception as e2:
                                logger.debug(f"    Failed to load {h5_file.name} (h5py): {e2}")
                                continue

                # Try AnnData files (.h5ad) if above formats not available
                for pattern in ['*.h5ad']:
                    files = [f for f in supp_dir.glob(pattern)]
                    if files:
                        h5ad_file = files[0]
                        logger.debug(f"    Loading {h5ad_file.name} (AnnData format)...")
                        try:
                            import anndata as ad
                            adata = ad.read_h5ad(h5ad_file)
                            df = adata.to_df()  # samples × genes
                            df = df.apply(pd.to_numeric, errors='coerce')
                            if df.shape[0] >= 2:
                                logger.debug(f"    Successfully loaded h5ad: {df.shape[0]} samples × {df.shape[1]} genes")
                                return df
                        except ImportError:
                            logger.debug("    anndata not installed; skipping .h5ad files")
                        except Exception as e:
                            logger.debug(f"    Failed to load {h5ad_file.name}: {e}")
                            continue

                # Try Parquet files (.parquet) if above formats not available
                for pattern in ['*.parquet']:
                    files = [f for f in supp_dir.glob(pattern)
                            if 'barcodes' not in f.name and 'features' not in f.name]
                    if files:
                        pq_file = files[0]
                        logger.debug(f"    Loading {pq_file.name} (Parquet format)...")
                        try:
                            df = pd.read_parquet(pq_file)
                            df = df.apply(pd.to_numeric, errors='coerce')
                            if df.shape[0] > df.shape[1]:
                                df = df.T
                            if df.shape[0] >= 2:
                                logger.debug(f"    Successfully loaded Parquet: {df.shape[0]} samples × {df.shape[1]} genes")
                                return df
                        except Exception as e:
                            logger.debug(f"    Failed to load {pq_file.name}: {e}")
                            continue

                # Try Feather files (.feather) if above formats not available
                for pattern in ['*.feather']:
                    files = [f for f in supp_dir.glob(pattern)]
                    if files:
                        ft_file = files[0]
                        logger.debug(f"    Loading {ft_file.name} (Feather format)...")
                        try:
                            df = pd.read_feather(ft_file)
                            # Feather may have gene names as a column (not index)
                            if df.columns[0] in ('gene', 'Gene', 'gene_id', 'Gene_ID', 'feature', 'Feature'):
                                df = df.set_index(df.columns[0])
                            df = df.apply(pd.to_numeric, errors='coerce')
                            if df.shape[0] > df.shape[1]:
                                df = df.T
                            if df.shape[0] >= 2:
                                logger.debug(f"    Successfully loaded Feather: {df.shape[0]} samples × {df.shape[1]} genes")
                                return df
                        except Exception as e:
                            logger.debug(f"    Failed to load {ft_file.name}: {e}")
                            continue

            logger.debug(f"    No suitable supplementary files found")
            return None

        except Exception as e:
            logger.debug(f"    _try_supplementary failed: {e}")
            return None

    @classmethod
    def _translate_ensg_to_hgnc(cls, expression_df: pd.DataFrame) -> pd.DataFrame:
        """
        FIX D: Translate Ensembl gene IDs (ENSG/ENSMUSG) to HGNC symbols using mygene.info API.

        Handles both human (ENSG) and mouse (ENSMUSG) IDs.
        Unmapped genes are dropped to keep column space clean.

        Returns:
            DataFrame with HGNC symbols as column names (possibly fewer columns)
        """
        try:
            import mygene
        except ImportError:
            logger.warning("  mygene not installed; cannot translate ENSG IDs. Install with: pip install mygene")
            return expression_df

        mg = mygene.MyGeneInfo()
        ensg_ids = list(expression_df.columns)

        # Detect species
        is_mouse = any(str(c).startswith('ENSMUSG') for c in ensg_ids)
        species = 'mouse' if is_mouse else 'human'

        logger.debug(f"    Translating {len(ensg_ids)} {species} Ensembl IDs to HGNC symbols...")

        # Batch query: 1000 IDs per request
        batch_size = 1000
        translation_map = {}

        for i in range(0, len(ensg_ids), batch_size):
            batch = ensg_ids[i:i + batch_size]
            try:
                results = mg.querymany(
                    batch,
                    scopes='ensembl.gene',
                    fields='symbol',
                    species=species,
                    returnall=False,
                    verbose=False,
                )
                for result in results:
                    if 'symbol' in result and not result.get('notfound', False):
                        query_id = result.get('query', '')
                        symbol = str(result['symbol']).upper().strip()
                        if query_id and symbol:
                            translation_map[query_id] = symbol
            except Exception as e:
                logger.debug(f"    mygene batch {i//batch_size + 1} failed: {e}")
                continue

        if not translation_map:
            logger.warning("    mygene translation returned no results; keeping original ENSG column names")
            return expression_df

        # Rename columns using translation_map, drop unmapped columns
        mapped_cols = {old: new for old, new in translation_map.items() if old in expression_df.columns}
        expression_df = expression_df.rename(columns=mapped_cols)

        # Drop columns that still have ENSG IDs (unmapped)
        hgnc_cols = [c for c in expression_df.columns if not (str(c).startswith('ENSG') or str(c).startswith('ENSMUSG'))]
        expression_df = expression_df[hgnc_cols]

        # Collapse duplicates (multiple ENSG IDs can map to same HGNC symbol) by mean
        if expression_df.columns.duplicated().any():
            expression_df = expression_df.T.groupby(level=0).mean().T

        logger.debug(f"    Translated {len(mapped_cols)}/{len(ensg_ids)} ENSG IDs to {len(hgnc_cols)} unique HGNC symbols")
        return expression_df

    @classmethod
    def _translate_entrez_to_hgnc(cls, expression_df: pd.DataFrame) -> pd.DataFrame:
        """
        Translate Entrez gene IDs (numeric strings) to HGNC symbols using mygene.info API.
        Uses batch queries for efficiency and includes error recovery.

        Args:
            expression_df: Expression matrix with Entrez IDs as column names

        Returns:
            Expression matrix with HGNC symbols as column names (fewer columns if some IDs don't map)
        """
        try:
            import mygene
        except ImportError:
            logger.warning("  mygene not installed; cannot translate Entrez IDs. Install with: pip install mygene")
            return expression_df

        mg = mygene.MyGeneInfo()
        entrez_ids = list(expression_df.columns)

        logger.debug(f"    Translating {len(entrez_ids)} Entrez IDs to HGNC symbols...")

        # Batch query: 1000 IDs per request (mygene.info batch limit)
        batch_size = 1000
        translation_map = {}
        failed_ids = []

        for i in range(0, len(entrez_ids), batch_size):
            batch = entrez_ids[i:i + batch_size]
            try:
                results = mg.querymany(
                    batch,
                    scopes='entrezgene',
                    fields='symbol',
                    species='human',
                    returnall=False,
                    verbose=False,
                )
                for result in results:
                    if 'symbol' in result and not result.get('notfound', False):
                        query_id = result.get('query', '')
                        symbol = str(result['symbol']).upper().strip()
                        if query_id and symbol:
                            translation_map[query_id] = symbol
                    elif result.get('notfound', False):
                        failed_ids.append(result.get('query', ''))
            except Exception as e:
                logger.debug(f"    mygene Entrez batch {i//batch_size + 1}/{(len(entrez_ids)+batch_size-1)//batch_size} failed: {str(e)[:100]}")
                # Don't add to failed_ids; treat as temporary network error
                continue

        if not translation_map:
            logger.warning("    mygene Entrez translation returned no results; keeping original column names")
            return expression_df

        # Rename columns using translation_map, drop unmapped columns
        mapped_cols = {old: new for old, new in translation_map.items() if old in expression_df.columns}
        expression_df = expression_df.rename(columns=mapped_cols)

        # Drop columns that still have numeric Entrez IDs (unmapped)
        hgnc_cols = [c for c in expression_df.columns if not str(c).isdigit()]
        expression_df = expression_df[hgnc_cols]

        # Collapse duplicates (multiple Entrez IDs can map to same HGNC symbol) by mean
        if expression_df.columns.duplicated().any():
            expression_df = expression_df.T.groupby(level=0).mean().T

        logger.debug(f"    Translated {len(mapped_cols)}/{len(entrez_ids)} Entrez IDs to {len(hgnc_cols)} unique HGNC symbols")
        return expression_df

    @classmethod
    def _translate_refseq_to_hgnc(cls, expression_df: pd.DataFrame) -> pd.DataFrame:
        """
        Translate RefSeq gene IDs (NM_*, NR_*, etc.) to HGNC symbols using mygene.info API.

        Args:
            expression_df: Expression matrix with RefSeq IDs as column names

        Returns:
            Expression matrix with HGNC symbols as column names (fewer columns if some IDs don't map)
        """
        try:
            import mygene
        except ImportError:
            logger.warning("  mygene not installed; cannot translate RefSeq IDs. Install with: pip install mygene")
            return expression_df

        mg = mygene.MyGeneInfo()
        refseq_ids = list(expression_df.columns)

        logger.debug(f"    Translating {len(refseq_ids)} RefSeq IDs to HGNC symbols...")

        # Batch query: 1000 IDs per request
        batch_size = 1000
        translation_map = {}

        for i in range(0, len(refseq_ids), batch_size):
            batch = refseq_ids[i:i + batch_size]
            try:
                results = mg.querymany(
                    batch,
                    scopes='refseq',
                    fields='symbol',
                    species='human',
                    returnall=False,
                    verbose=False,
                )
                for result in results:
                    if 'symbol' in result and not result.get('notfound', False):
                        query_id = result.get('query', '')
                        symbol = str(result['symbol']).upper().strip()
                        if query_id and symbol:
                            translation_map[query_id] = symbol
            except Exception as e:
                logger.debug(f"    mygene RefSeq batch {i//batch_size + 1} failed: {e}")
                continue

        if not translation_map:
            logger.warning("    mygene RefSeq translation returned no results; keeping original column names")
            return expression_df

        # Rename columns using translation_map, drop unmapped columns
        REFSEQ_PREFIXES = ('NM_', 'NR_', 'XM_', 'XR_', 'NP_', 'NG_')
        mapped_cols = {old: new for old, new in translation_map.items() if old in expression_df.columns}
        expression_df = expression_df.rename(columns=mapped_cols)

        # Drop columns that still have RefSeq IDs (unmapped)
        hgnc_cols = [c for c in expression_df.columns if not str(c).startswith(REFSEQ_PREFIXES)]
        expression_df = expression_df[hgnc_cols]

        # Collapse duplicates (multiple RefSeq IDs can map to same HGNC symbol) by mean
        if expression_df.columns.duplicated().any():
            expression_df = expression_df.T.groupby(level=0).mean().T

        logger.debug(f"    Translated {len(mapped_cols)}/{len(refseq_ids)} RefSeq IDs to {len(hgnc_cols)} unique HGNC symbols")
        return expression_df

    @classmethod
    def _impute_nan(cls, expression_df: pd.DataFrame) -> pd.DataFrame:
        """
        Impute NaN values using column-wise median (fixes numpy.ndarray.median() bug).

        Args:
            expression_df: Expression matrix with possible NaN values

        Returns:
            Expression matrix with NaN values imputed
        """
        if not expression_df.isna().any().any():
            return expression_df  # No NaN values

        logger.debug(f"    Imputing NaN values...")

        expression_df = expression_df.copy()

        # Column-wise imputation
        for col in expression_df.columns:
            if expression_df[col].isna().any():
                col_median = expression_df[col].median()

                if pd.isna(col_median):
                    # If entire column is NaN, use global median (properly computed)
                    valid_vals = expression_df.values[~np.isnan(expression_df.values)]
                    if len(valid_vals) > 0:
                        col_median = np.nanmedian(valid_vals)
                    else:
                        col_median = 0.0

                expression_df[col].fillna(col_median, inplace=True)

        return expression_df

    @classmethod
    def _build_expression_matrix(cls, gse, accession_id: str, probe_to_symbol: dict) -> Optional[pd.DataFrame]:
        """
        Three-stage fallback for building expression matrix.

        1. Try SOFT table with probe-to-symbol mapping
        2. If <10 genes, try supplementary files
        3. If both fail, return None

        Returns:
            Expression DataFrame (samples × genes) or None
        """
        # Stage 1: SOFT table
        expression_df = cls._try_soft_table(gse, probe_to_symbol)

        # Check if SOFT table has sufficient genes (FIX 1: Make logic explicit)
        has_sufficient_genes = expression_df is not None and expression_df.shape[1] >= 10

        if has_sufficient_genes:
            logger.debug(f"    SOFT table: success ({expression_df.shape[1]} genes)")
        else:
            # Stage 2: Supplementary files (trigger on None OR insufficient genes)
            if expression_df is not None:
                logger.debug(f"    SOFT table: {expression_df.shape[1]} genes (insufficient); trying supplementary files...")
            else:
                logger.debug(f"    SOFT table: failed completely; trying supplementary files...")

            expression_df = cls._try_supplementary(gse, accession_id)

            if expression_df is not None:
                logger.debug(f"    Supplementary file: success ({expression_df.shape[1]} genes)")
            else:
                logger.debug(f"    Supplementary file: also failed; no fallback available")
                return None

        # FIX B: Normalize gene column names (uppercase, strip whitespace)
        expression_df.columns = expression_df.columns.str.upper().str.strip()

        # Collapse duplicate columns (same gene symbol after normalization) by mean
        if expression_df.columns.duplicated().any():
            logger.debug(f"    Collapsing {expression_df.columns.duplicated().sum()} duplicate gene columns by mean...")
            expression_df = expression_df.T.groupby(level=0).mean().T

        # FIX D: Detect and translate ENSG IDs to HGNC symbols
        ensg_cols = [c for c in expression_df.columns if str(c).startswith('ENSG') or str(c).startswith('ENSMUSG')]
        ensg_fraction = len(ensg_cols) / max(len(expression_df.columns), 1)

        if ensg_fraction > 0.5:
            logger.debug(f"    {len(ensg_cols)}/{len(expression_df.columns)} columns are Ensembl IDs; translating to HGNC...")
            expression_df = cls._translate_ensg_to_hgnc(expression_df)
            logger.debug(f"    After translation: {expression_df.shape[1]} HGNC gene columns")

        # Detect and translate Entrez IDs (numeric strings)
        entrez_cols = [c for c in expression_df.columns if str(c).isdigit()]
        entrez_fraction = len(entrez_cols) / max(len(expression_df.columns), 1)
        if entrez_fraction > 0.5:
            logger.debug(f"    {len(entrez_cols)}/{len(expression_df.columns)} columns are Entrez IDs; translating to HGNC...")
            expression_df = cls._translate_entrez_to_hgnc(expression_df)
            logger.debug(f"    After Entrez translation: {expression_df.shape[1]} HGNC gene columns")

        # Detect and translate RefSeq IDs (NM_*, NR_*, etc.)
        REFSEQ_PREFIXES = ('NM_', 'NR_', 'XM_', 'XR_', 'NP_', 'NG_')
        refseq_cols = [c for c in expression_df.columns if str(c).startswith(REFSEQ_PREFIXES)]
        refseq_fraction = len(refseq_cols) / max(len(expression_df.columns), 1)
        if refseq_fraction > 0.5:
            logger.debug(f"    {len(refseq_cols)}/{len(expression_df.columns)} columns are RefSeq IDs; translating to HGNC...")
            expression_df = cls._translate_refseq_to_hgnc(expression_df)
            logger.debug(f"    After RefSeq translation: {expression_df.shape[1]} HGNC gene columns")

        # Impute NaN values
        expression_df = cls._impute_nan(expression_df)

        # FIX 4: Verify no NaN remains after imputation
        if expression_df.isna().any().any():
            remaining_nans = expression_df.isna().sum().sum()
            logger.warning(f"    WARNING: {remaining_nans} NaN values remain after imputation; filling with 0")
            expression_df = expression_df.fillna(0.0)

        # FIX 5: Log2 transform if needed (only on non-NaN values)
        valid_vals = expression_df.values[~np.isnan(expression_df.values)]
        if len(valid_vals) > 0:
            max_val = np.nanmax(valid_vals)
            if max_val > 50:
                logger.debug(f"    Log2 transforming (max value: {max_val:.1f})...")
                expression_df = np.log2(expression_df + 1)

        return expression_df

    @classmethod
    def filter_expressed_genes(
        cls,
        expression_df: pd.DataFrame,
        min_expression_quantile: float = 0.25,
        remove_zero_variance: bool = True
    ) -> pd.DataFrame:
        """
        Filter genes by expression level and variance.

        Args:
            expression_df: Expression matrix (samples × genes)
            min_expression_quantile: Minimum expression threshold as quantile (0-1)
            remove_zero_variance: Whether to remove genes with zero variance

        Returns:
            Filtered expression matrix
        """
        n_genes_initial = expression_df.shape[1]

        # Filter by expression level: keep genes in top quantile by mean expression
        if min_expression_quantile > 0:
            mean_expr = expression_df.mean()
            threshold = mean_expr.quantile(min_expression_quantile)
            expressed_genes = mean_expr[mean_expr >= threshold].index
            expression_df = expression_df[expressed_genes]
            logger.debug(f"    Filtered by expression: {n_genes_initial} → {expression_df.shape[1]} genes")

        # Filter by variance: remove genes with zero or near-zero variance
        if remove_zero_variance:
            variance = expression_df.var()
            variable_genes = variance[variance > 1e-10].index
            expression_df = expression_df[variable_genes]
            logger.debug(f"    Filtered by variance: {expression_df.shape[1]} → {len(variable_genes)} genes")

        return expression_df


class BenchmarkDataLoaderFactory:
    """Factory for loading benchmark datasets based on manifest row."""

    @staticmethod
    def load(row: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load dataset based on manifest row (source + accession_id).

        Args:
            row: Pandas Series from benchmark manifest

        Returns:
            (gene_expression, phenotype_labels)
        """
        source = row['source']
        accession = row['accession_id']
        dataset_name = row['dataset_name']

        logger.info(f"Loading {dataset_name}...")

        try:
            if source == 'GDC':
                # Try local loader first (for manually downloaded TCGA data)
                # Check from repo root (since script is run from scripts/ directory)
                repo_root = Path(__file__).parent.parent.parent
                local_dir = repo_root / 'tcga_data' / accession
                logger.debug(f"  Checking for local data at: {local_dir}")
                if local_dir.exists() and list(local_dir.glob('*.tsv')):
                    logger.info(f"  Using local GDC loader for {accession}...")
                    try:
                        # Import with fallback for different import contexts
                        try:
                            from .local_gdc_loader import LocalGDCDataLoader
                        except ImportError:
                            import sys
                            sys.path.insert(0, str(Path(__file__).parent))
                            from local_gdc_loader import LocalGDCDataLoader
                        # Use absolute path to tcga_data directory
                        tcga_data_path = str(repo_root / 'tcga_data')
                        return LocalGDCDataLoader.load(accession, data_dir=tcga_data_path)
                    except Exception as e:
                        logger.warning(f"  Local GDC loader failed, falling back to API: {e}")
                        return GDCDataLoader.load(accession)
                else:
                    logger.info(f"  Using GDC API loader for {accession}...")
                    return GDCDataLoader.load(accession)

            elif source == 'GEO':
                return RobustGEOLoader.load(accession)

            elif source == 'CellxGene':
                return SingleCellDataLoader.load_cellxgene(accession)

            elif source == 'ArrayExpress':
                # For ArrayExpress, assume h5ad file is available locally
                local_path = Path(f"data/benchmarks/{accession}.h5ad")
                if not local_path.exists():
                    raise DataLoaderError(f"ArrayExpress file not found: {local_path}")
                return SingleCellDataLoader.load_h5ad(str(local_path))

            elif source == 'BioConductor':
                # BioConductor datasets would need rpy2 implementation
                logger.warning(f"  BioConductor loader not yet implemented; skipping {accession}")
                raise DataLoaderError(f"BioConductor loader not implemented")

            elif source == 'NCBI':
                # For microbial data, assume files are cached
                otu_path = Path(f"data/benchmarks/{accession}_otu.txt")
                tax_path = Path(f"data/benchmarks/{accession}_taxonomy.txt")
                if not otu_path.exists() or not tax_path.exists():
                    raise DataLoaderError(f"Microbial data files not found: {otu_path}, {tax_path}")
                return CustomDataLoader.load_microbial_16s(str(otu_path), str(tax_path))

            elif source == 'Figshare' or source == 'HCA':
                # Assume files are cached locally
                csv_path = Path(f"data/benchmarks/{accession}_expression.csv")
                if not csv_path.exists():
                    raise DataLoaderError(f"Data file not found: {csv_path}")
                phenotype_col = row.get('ground_truth_label_column')
                return CustomDataLoader.load_generic_csv(str(csv_path), phenotype_col)

            elif source == 'Custom':
                # Generic CSV loader
                csv_path = Path(f"data/benchmarks/{accession}_expression.csv")
                if not csv_path.exists():
                    raise DataLoaderError(f"Data file not found: {csv_path}")
                phenotype_col = row.get('ground_truth_label_column')
                return CustomDataLoader.load_generic_csv(str(csv_path), phenotype_col)

            else:
                raise DataLoaderError(f"Unknown data source: {source}")

        except DataLoaderError as e:
            logger.error(f"Data loading error: {e}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error loading {dataset_name}: {e}")
            raise DataLoaderError(f"Failed to load {dataset_name}: {e}")
