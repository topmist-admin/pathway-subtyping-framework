"""
Tests for the utils/performance module.

Tests cover:
- Chunked VCF reading (gzip support, multi-allelic expansion)
- Chunked gene burden computation
- Parallel pathway scoring
- Memory estimation
- Cohort downsampling
- Progress tracking
"""

import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.utils.performance import (
    ProcessingStats,
    ProgressTracker,
    chunked_vcf_reader,
    compute_gene_burdens_chunked,
    downsample_cohort,
    estimate_memory_usage,
    parallel_pathway_scores,
)


class TestChunkedVCFReader:
    """Tests for chunked VCF reading."""

    @pytest.fixture
    def simple_vcf(self, tmp_path):
        """Create a simple VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t100\tv1\tA\tG\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=missense;CADD=25\tGT\t0/1\t0/0\t1/1
chr1\t200\tv2\tC\tT\t99\tPASS\tGENE=CHD8;CONSEQUENCE=frameshift;CADD=35\tGT\t0/0\t0/1\t0/1
chr1\t300\tv3\tG\tA\t99\tPASS\tGENE=NRXN1;CONSEQUENCE=missense;CADD=20\tGT\t0/1\t0/1\t0/0
"""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    @pytest.fixture
    def gzipped_vcf(self, tmp_path):
        """Create a gzipped VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tv1\tA\tG\t99\tPASS\tGENE=SHANK3\tGT\t0/1\t0/0
chr1\t200\tv2\tC\tT\t99\tPASS\tGENE=CHD8\tGT\t0/0\t0/1
"""
        vcf_path = tmp_path / "test.vcf.gz"
        with gzip.open(vcf_path, "wt") as f:
            f.write(vcf_content)
        return vcf_path

    @pytest.fixture
    def multi_allelic_vcf(self, tmp_path):
        """Create a VCF with multi-allelic variants."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t100\tvar1\tA\tG,T\t99\tPASS\tGENE=TEST\tGT\t0/1\t0/2\t1/2
"""
        vcf_path = tmp_path / "multi.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    def test_read_simple_vcf(self, simple_vcf):
        """Test reading a simple VCF file."""
        chunks = list(chunked_vcf_reader(str(simple_vcf), chunk_size=10))

        assert len(chunks) == 1
        variants_df, genotypes_df, samples = chunks[0]

        assert len(variants_df) == 3
        assert len(samples) == 3
        assert samples == ["S1", "S2", "S3"]

    def test_read_gzipped_vcf(self, gzipped_vcf):
        """Test reading a gzipped VCF file."""
        chunks = list(chunked_vcf_reader(str(gzipped_vcf), chunk_size=10))

        assert len(chunks) == 1
        variants_df, genotypes_df, samples = chunks[0]

        assert len(variants_df) == 2
        assert len(samples) == 2

    def test_chunked_reading(self, simple_vcf):
        """Test that chunking works correctly."""
        chunks = list(chunked_vcf_reader(str(simple_vcf), chunk_size=2))

        assert len(chunks) == 2
        assert len(chunks[0][0]) == 2  # First chunk has 2 variants
        assert len(chunks[1][0]) == 1  # Second chunk has 1 variant

    def test_sample_subset(self, simple_vcf):
        """Test reading only a subset of samples."""
        chunks = list(
            chunked_vcf_reader(str(simple_vcf), chunk_size=10, sample_subset=["S1", "S3"])
        )

        variants_df, genotypes_df, samples = chunks[0]
        assert samples == ["S1", "S3"]
        assert len(genotypes_df.columns) == 2

    def test_multi_allelic_expansion(self, multi_allelic_vcf):
        """Test that multi-allelic variants are expanded."""
        chunks = list(chunked_vcf_reader(str(multi_allelic_vcf), chunk_size=10))

        variants_df, genotypes_df, samples = chunks[0]

        # One multi-allelic with 2 ALTs should become 2 records
        assert len(variants_df) == 2
        assert "var1_1" in variants_df["id"].values
        assert "var1_2" in variants_df["id"].values

    def test_multi_allelic_genotypes(self, multi_allelic_vcf):
        """Test genotype parsing for multi-allelic variants."""
        chunks = list(chunked_vcf_reader(str(multi_allelic_vcf), chunk_size=10))

        variants_df, genotypes_df, samples = chunks[0]

        # For S1 (0/1): allele 1 count=1, allele 2 count=0
        assert genotypes_df.loc[0, "S1"] == 1  # First alt
        assert genotypes_df.loc[1, "S1"] == 0  # Second alt

        # For S2 (0/2): allele 1 count=0, allele 2 count=1
        assert genotypes_df.loc[0, "S2"] == 0
        assert genotypes_df.loc[1, "S2"] == 1

        # For S3 (1/2): allele 1 count=1, allele 2 count=1
        assert genotypes_df.loc[0, "S3"] == 1
        assert genotypes_df.loc[1, "S3"] == 1

    def test_file_not_found(self, tmp_path):
        """Test error when VCF file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            list(chunked_vcf_reader(str(tmp_path / "nonexistent.vcf")))


class TestComputeGeneBurdensChunked:
    """Tests for chunked gene burden computation."""

    @pytest.fixture
    def annotated_vcf(self, tmp_path):
        """Create a VCF with gene annotations."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tv1\tA\tG\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=missense;CADD=25\tGT\t0/1\t0/0
chr1\t200\tv2\tC\tT\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=frameshift;CADD=35\tGT\t0/1\t0/0
chr1\t300\tv3\tG\tA\t99\tPASS\tGENE=CHD8;CONSEQUENCE=missense;CADD=20\tGT\t0/0\t0/1
"""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    @pytest.fixture
    def vcf_missing_cadd(self, tmp_path):
        """Create a VCF with missing CADD scores."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tv1\tA\tG\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=frameshift\tGT\t0/1\t0/0
chr1\t200\tv2\tC\tT\t99\tPASS\tGENE=CHD8;CONSEQUENCE=missense\tGT\t0/0\t0/1
"""
        vcf_path = tmp_path / "missing_cadd.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    def test_basic_burden_computation(self, annotated_vcf):
        """Test basic gene burden computation."""
        burdens = compute_gene_burdens_chunked(str(annotated_vcf))

        assert isinstance(burdens, pd.DataFrame)
        assert "SHANK3" in burdens.columns
        assert "CHD8" in burdens.columns
        assert len(burdens) == 2  # Two samples

    def test_burden_values_nonzero(self, annotated_vcf):
        """Test that burden values are computed correctly."""
        burdens = compute_gene_burdens_chunked(str(annotated_vcf))

        # S1 has variants in SHANK3, S2 has variant in CHD8
        assert burdens.loc["S1", "SHANK3"] > 0
        assert burdens.loc["S2", "CHD8"] > 0

    def test_missing_cadd_defaults(self, vcf_missing_cadd):
        """Test that missing CADD scores use defaults."""
        burdens = compute_gene_burdens_chunked(str(vcf_missing_cadd))

        # Should still compute burdens with default CADD values
        assert burdens.loc["S1", "SHANK3"] > 0
        assert burdens.loc["S2", "CHD8"] > 0

    def test_progress_callback(self, annotated_vcf):
        """Test progress callback is called."""
        progress_values = []

        def callback(current, total):
            progress_values.append(current)

        compute_gene_burdens_chunked(str(annotated_vcf), progress_callback=callback)

        assert len(progress_values) > 0


class TestParallelPathwayScores:
    """Tests for parallel pathway scoring."""

    @pytest.fixture
    def gene_burdens(self):
        """Create sample gene burden data with variance.

        Note: Gene values are chosen so that pathway means (averages of gene pairs)
        have variance across samples. Avoid complementary genes that sum to constant.
        """
        return pd.DataFrame(
            {
                "SHANK3": [1.0, 0.5, 0.0, 0.8],
                "CHD8": [0.2, 0.7, 0.4, 0.1],  # Not complementary to SHANK3
                "NRXN1": [0.1, 0.6, 0.3, 0.9],
                "SETD5": [0.2, 0.8, 0.3, 0.7],
                "ASH1L": [0.7, 0.2, 0.5, 0.4],
                "SCN2A": [0.3, 0.9, 0.1, 0.6],
            },
            index=["S1", "S2", "S3", "S4"],
        )

    @pytest.fixture
    def pathways(self):
        """Create sample pathways with distinct gene combinations."""
        return {
            "SYNAPTIC": ["SHANK3", "NRXN1"],
            "CHROMATIN": ["CHD8", "SETD5"],
            "ION_CHANNEL": ["ASH1L", "SCN2A"],  # Distinct genes
        }

    def test_basic_scoring(self, gene_burdens, pathways):
        """Test basic pathway scoring."""
        scores = parallel_pathway_scores(gene_burdens, pathways)

        assert isinstance(scores, pd.DataFrame)
        assert len(scores) == 4  # 4 samples
        assert len(scores.columns) >= 2  # At least 2 pathways after zero-variance filter

    def test_z_score_normalization(self, gene_burdens, pathways):
        """Test that scores are Z-normalized."""
        scores = parallel_pathway_scores(gene_burdens, pathways)

        # Z-scores should have mean ~0 and std ~1 (within tolerance)
        for col in scores.columns:
            assert abs(scores[col].mean()) < 0.1
            assert 0.9 < scores[col].std() < 1.1

    def test_zero_variance_handling(self):
        """Test handling of zero-variance pathways."""
        # Create burdens with one constant gene and varied genes
        # Note: GENE2 + GENE3 must NOT sum to a constant to avoid zero-variance pathway
        gene_burdens = pd.DataFrame(
            {
                "GENE1": [0.5, 0.5, 0.5, 0.5],  # Constant - will cause zero variance
                "GENE2": [1.0, 0.5, 0.0, 0.8],
                "GENE3": [0.2, 0.7, 0.4, 0.1],  # Not complementary to GENE2
                "GENE4": [0.2, 0.8, 0.3, 0.7],
            },
            index=["S1", "S2", "S3", "S4"],
        )

        pathways = {
            "VARIED1": ["GENE2", "GENE3"],  # Will have variance (genes not complementary)
            "VARIED2": ["GENE3", "GENE4"],  # Will have variance
            "CONSTANT": ["GENE1", "GENE1"],  # Zero variance - should be filtered
        }

        # Should not raise (we have 2 valid pathways after filtering)
        scores = parallel_pathway_scores(gene_burdens, pathways)
        assert len(scores.columns) >= 2
        assert "CONSTANT" not in scores.columns

    def test_pathway_with_few_genes(self, gene_burdens):
        """Test pathways with genes not in burden data."""
        pathways = {
            "VALID1": ["SHANK3", "CHD8"],  # Both genes in burden data
            "VALID2": ["NRXN1", "SETD5"],  # Both genes in burden data
            "VALID3": ["ASH1L", "SCN2A"],  # Both genes in burden data
            "INVALID": ["UNKNOWN1", "UNKNOWN2"],  # No overlap with burden data
        }

        scores = parallel_pathway_scores(gene_burdens, pathways)

        # Invalid pathway should be excluded
        assert "INVALID" not in scores.columns
        # At least 2 valid pathways should remain (needed for downstream clustering)
        assert len(scores.columns) >= 2

    def test_insufficient_pathways_error(self, gene_burdens):
        """Test error when no valid pathways."""
        pathways = {"INVALID": ["UNKNOWN1", "UNKNOWN2"]}

        with pytest.raises(ValueError, match="No valid pathways"):
            parallel_pathway_scores(gene_burdens, pathways)

    def test_parallel_execution(self, gene_burdens):
        """Test that parallel execution works."""
        # Create many pathways with different gene combinations
        all_genes = list(gene_burdens.columns)
        pathways = {}
        for i in range(10):
            # Use different gene pairs to avoid all-identical pathways
            g1 = all_genes[i % len(all_genes)]
            g2 = all_genes[(i + 1) % len(all_genes)]
            pathways[f"PATH_{i}"] = [g1, g2]

        scores = parallel_pathway_scores(gene_burdens, pathways, n_workers=2)

        # Should have at least some pathways after zero-variance filtering
        assert len(scores.columns) >= 2


class TestMemoryEstimation:
    """Tests for memory estimation."""

    def test_basic_estimation(self):
        """Test basic memory estimation."""
        estimate = estimate_memory_usage(n_samples=100, n_variants=1000, n_genes=500, n_pathways=50)

        assert "total_estimated_mb" in estimate
        assert estimate["total_estimated_mb"] > 0

    def test_scaling(self):
        """Test that estimates scale with data size."""
        small = estimate_memory_usage(100, 1000, 500, 50)
        large = estimate_memory_usage(1000, 10000, 5000, 500)

        assert large["total_estimated_mb"] > small["total_estimated_mb"]

    def test_chunk_size_recommendation(self):
        """Test chunk size recommendation."""
        estimate = estimate_memory_usage(10000, 100000, 20000, 100)

        assert "recommended_chunk_size" in estimate
        assert 100 <= estimate["recommended_chunk_size"] <= 1000


class TestDownsampleCohort:
    """Tests for cohort downsampling."""

    @pytest.fixture
    def phenotypes(self):
        """Create sample phenotypes."""
        return pd.DataFrame(
            {
                "age": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
                "group": ["A", "A", "A", "B", "B", "B", "C", "C", "C", "C"],
            },
            index=[f"S{i}" for i in range(10)],
        )

    def test_basic_downsampling(self, phenotypes):
        """Test basic random downsampling."""
        selected = downsample_cohort(phenotypes, target_size=5)

        assert len(selected) == 5
        assert all(s in phenotypes.index for s in selected)

    def test_stratified_downsampling(self, phenotypes):
        """Test stratified downsampling."""
        selected = downsample_cohort(phenotypes, target_size=6, stratify_by="group")

        assert len(selected) == 6
        # Should have samples from each group
        groups = phenotypes.loc[selected, "group"].value_counts()
        assert len(groups) == 3  # All 3 groups represented

    def test_reproducibility(self, phenotypes):
        """Test that downsampling is reproducible with seed."""
        selected1 = downsample_cohort(phenotypes, target_size=5, seed=42)
        selected2 = downsample_cohort(phenotypes, target_size=5, seed=42)

        assert selected1 == selected2

    def test_different_seeds(self, phenotypes):
        """Test that different seeds give different results."""
        selected1 = downsample_cohort(phenotypes, target_size=5, seed=42)
        selected2 = downsample_cohort(phenotypes, target_size=5, seed=123)

        # Very unlikely to be identical
        assert selected1 != selected2

    def test_target_larger_than_cohort(self, phenotypes):
        """Test when target is larger than cohort."""
        selected = downsample_cohort(phenotypes, target_size=100)

        assert len(selected) == len(phenotypes)


class TestProgressTracker:
    """Tests for progress tracking."""

    def test_basic_tracking(self, caplog):
        """Test basic progress tracking."""
        import logging

        caplog.set_level(logging.INFO)

        tracker = ProgressTracker(total=100, desc="Testing", log_interval=50)
        tracker.update(50)
        tracker.update(100)
        tracker.finish()

        assert "Testing" in caplog.text

    def test_unknown_total(self, caplog):
        """Test tracking with unknown total."""
        import logging

        caplog.set_level(logging.INFO)

        tracker = ProgressTracker(total=None, desc="Processing", log_interval=10)
        tracker.update(10)
        tracker.update(20)
        tracker.finish()

        assert "Processing" in caplog.text


class TestProcessingStats:
    """Tests for ProcessingStats dataclass."""

    def test_create_stats(self):
        """Test creating ProcessingStats."""
        stats = ProcessingStats(
            total_samples=100,
            total_variants=1000,
            chunks_processed=10,
            peak_memory_mb=256.0,
            processing_time_seconds=30.5,
        )

        assert stats.total_samples == 100
        assert stats.total_variants == 1000
        assert stats.chunks_processed == 10
