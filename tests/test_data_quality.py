"""
Tests for the data_quality module (Week 4 - Real-World Data Support).

Tests cover:
- Multi-allelic variant handling
- Missing annotation graceful handling
- Data quality reporting
- VCF validation
"""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from pathway_subtyping.data_quality import (
    DataQualityReport,
    VCFDataQualityError,
    expand_multi_allelic,
    load_vcf_with_quality_check,
    parse_genotype,
    parse_info_field,
    validate_vcf_for_pipeline,
)


class TestParseGenotype:
    """Tests for genotype parsing."""

    def test_standard_homozygous_ref(self):
        """Test parsing 0/0 genotype."""
        count, valid = parse_genotype("0/0")
        assert count == 0
        assert valid is True

    def test_standard_heterozygous(self):
        """Test parsing 0/1 genotype."""
        count, valid = parse_genotype("0/1")
        assert count == 1
        assert valid is True

    def test_standard_homozygous_alt(self):
        """Test parsing 1/1 genotype."""
        count, valid = parse_genotype("1/1")
        assert count == 2
        assert valid is True

    def test_phased_genotype(self):
        """Test parsing phased genotype (pipe separator)."""
        count, valid = parse_genotype("0|1")
        assert count == 1
        assert valid is True

    def test_multi_allelic_genotype_0_2_target_1(self):
        """Test parsing multi-allelic genotype 0/2 for allele 1."""
        count, valid = parse_genotype("0/2", target_allele=1)
        assert count == 0  # Zero copies of allele 1
        assert valid is True

    def test_multi_allelic_genotype_0_2_target_2(self):
        """Test parsing multi-allelic genotype 0/2 for allele 2."""
        count, valid = parse_genotype("0/2", target_allele=2)
        assert count == 1  # One copy of allele 2
        assert valid is True

    def test_multi_allelic_het_1_2_target_1(self):
        """Test parsing multi-allelic genotype 1/2 for allele 1."""
        count, valid = parse_genotype("1/2", target_allele=1)
        assert count == 1  # One copy of allele 1
        assert valid is True

    def test_multi_allelic_het_1_2_target_2(self):
        """Test parsing multi-allelic genotype 1/2 for allele 2."""
        count, valid = parse_genotype("1/2", target_allele=2)
        assert count == 1  # One copy of allele 2
        assert valid is True

    def test_homozygous_alt_2(self):
        """Test parsing multi-allelic homozygous 2/2 for allele 2."""
        count, valid = parse_genotype("2/2", target_allele=2)
        assert count == 2  # Two copies of allele 2
        assert valid is True

    def test_missing_genotype_dots(self):
        """Test parsing missing genotype ./."""
        count, valid = parse_genotype("./.")
        assert count == 0
        assert valid is False

    def test_missing_genotype_single_dot(self):
        """Test parsing single dot genotype."""
        count, valid = parse_genotype(".")
        assert count == 0
        assert valid is False

    def test_genotype_with_format_fields(self):
        """Test parsing genotype with additional format fields."""
        count, valid = parse_genotype("0/1:35:99")
        assert count == 1
        assert valid is True

    def test_empty_genotype(self):
        """Test parsing empty genotype."""
        count, valid = parse_genotype("")
        assert count == 0
        assert valid is False

    def test_malformed_genotype(self):
        """Test parsing malformed genotype."""
        count, valid = parse_genotype("X/Y")
        assert count == 0
        assert valid is False


class TestParseInfoField:
    """Tests for INFO field parsing."""

    def test_standard_info_field(self):
        """Test parsing standard INFO field."""
        info = parse_info_field("GENE=SHANK3;CONSEQUENCE=missense_variant;CADD=28.5")
        assert info["GENE"] == "SHANK3"
        assert info["CONSEQUENCE"] == "missense_variant"
        assert info["CADD"] == "28.5"

    def test_flag_info_field(self):
        """Test parsing INFO field with flags."""
        info = parse_info_field("GENE=CHD8;SOMATIC")
        assert info["GENE"] == "CHD8"
        assert info["SOMATIC"] is True

    def test_empty_info_field(self):
        """Test parsing empty INFO field."""
        info = parse_info_field(".")
        assert info == {}

    def test_multi_value_info_field(self):
        """Test parsing INFO field with multiple values."""
        info = parse_info_field("AC=10,5;AF=0.1,0.05")
        assert info["AC"] == ["10", "5"]
        assert info["AF"] == ["0.1", "0.05"]

    def test_missing_value(self):
        """Test parsing INFO field with missing value."""
        info = parse_info_field("GENE=.;CADD=25")
        assert info["GENE"] is None
        assert info["CADD"] == "25"


class TestExpandMultiAllelic:
    """Tests for multi-allelic variant expansion."""

    def test_expand_biallelic(self):
        """Test expanding a single alternate allele (no expansion needed)."""
        expanded = expand_multi_allelic(
            chrom="chr1",
            pos=100,
            vid="var1",
            ref="A",
            alts=["G"],
            info_dict={"GENE": "TEST"},
            sample_genotypes=["0/1", "0/0"],
            samples=["S1", "S2"],
        )
        assert len(expanded) == 1
        assert expanded[0]["alt"] == "G"
        assert expanded[0]["genotypes"]["S1"] == 1
        assert expanded[0]["genotypes"]["S2"] == 0

    def test_expand_multi_allelic(self):
        """Test expanding two alternate alleles."""
        expanded = expand_multi_allelic(
            chrom="chr1",
            pos=100,
            vid="var1",
            ref="A",
            alts=["G", "T"],
            info_dict={"GENE": "TEST"},
            sample_genotypes=["0/1", "0/2", "1/2"],
            samples=["S1", "S2", "S3"],
        )
        assert len(expanded) == 2

        # First allele (G)
        assert expanded[0]["alt"] == "G"
        assert expanded[0]["id"] == "var1_1"
        assert expanded[0]["genotypes"]["S1"] == 1  # 0/1 has allele 1
        assert expanded[0]["genotypes"]["S2"] == 0  # 0/2 has allele 2, not 1
        assert expanded[0]["genotypes"]["S3"] == 1  # 1/2 has allele 1

        # Second allele (T)
        assert expanded[1]["alt"] == "T"
        assert expanded[1]["id"] == "var1_2"
        assert expanded[1]["genotypes"]["S1"] == 0  # 0/1 has allele 1, not 2
        assert expanded[1]["genotypes"]["S2"] == 1  # 0/2 has allele 2
        assert expanded[1]["genotypes"]["S3"] == 1  # 1/2 has allele 2


class TestDataQualityReport:
    """Tests for DataQualityReport class."""

    def test_coverage_calculation(self):
        """Test annotation coverage percentage calculation."""
        report = DataQualityReport()
        report.parsed_variants = 100
        report.variants_with_gene = 80
        report.variants_with_consequence = 70
        report.variants_with_cadd = 50

        assert report.gene_coverage == 80.0
        assert report.consequence_coverage == 70.0
        assert report.cadd_coverage == 50.0

    def test_is_usable_pass(self):
        """Test is_usable returns True when coverage is sufficient."""
        report = DataQualityReport()
        report.parsed_variants = 100
        report.variants_with_gene = 60  # 60% coverage

        assert report.is_usable is True

    def test_is_usable_fail(self):
        """Test is_usable returns False when coverage is insufficient."""
        report = DataQualityReport()
        report.parsed_variants = 100
        report.variants_with_gene = 40  # 40% coverage

        assert report.is_usable is False

    def test_to_dict(self):
        """Test conversion to dictionary."""
        report = DataQualityReport()
        report.total_variants = 100
        report.parsed_variants = 98
        report.skipped_variants = 2
        report.variants_with_gene = 80

        result = report.to_dict()
        assert result["total_variants"] == 100
        assert result["parsed_variants"] == 98
        assert result["is_usable"] is True

    def test_add_warning(self):
        """Test adding warnings."""
        report = DataQualityReport()
        report.add_warning("Test warning 1")
        report.add_warning("Test warning 2")
        report.add_warning("Test warning 1")  # Duplicate

        assert len(report.warnings) == 2

    def test_summary(self):
        """Test summary generation."""
        report = DataQualityReport()
        report.total_variants = 100
        report.parsed_variants = 100
        report.variants_with_gene = 80

        summary = report.summary()
        assert "Total variants: 100" in summary
        assert "GENE:" in summary
        assert "PASS" in summary


class TestLoadVCFWithQualityCheck:
    """Tests for VCF loading with quality checks."""

    @pytest.fixture
    def vcf_with_annotations(self, tmp_path):
        """Create a well-annotated VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tvar1\tA\tG\t99\tPASS\tGENE=SHANK3;CONSEQUENCE=missense;CADD=28.5\tGT\t0/1\t0/0
chr1\t200\tvar2\tC\tT\t99\tPASS\tGENE=CHD8;CONSEQUENCE=frameshift;CADD=35\tGT\t0/0\t0/1
"""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    @pytest.fixture
    def vcf_without_annotations(self, tmp_path):
        """Create a VCF file without annotations."""
        vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tvar1\tA\tG\t99\tPASS\t.\tGT\t0/1\t0/0
chr1\t200\tvar2\tC\tT\t99\tPASS\t.\tGT\t0/0\t0/1
"""
        vcf_path = tmp_path / "no_anno.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    @pytest.fixture
    def vcf_with_multi_allelic(self, tmp_path):
        """Create a VCF file with multi-allelic variants."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t100\tvar1\tA\tG,T\t99\tPASS\tGENE=TEST;CONSEQUENCE=missense;CADD=25\tGT\t0/1\t0/2\t1/2
"""
        vcf_path = tmp_path / "multi.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    def test_load_annotated_vcf(self, vcf_with_annotations):
        """Test loading well-annotated VCF."""
        variants_df, genotypes_df, samples, report = load_vcf_with_quality_check(
            str(vcf_with_annotations)
        )

        assert len(variants_df) == 2
        assert len(samples) == 2
        assert report.gene_coverage == 100.0
        assert report.is_usable is True

    def test_load_unannotated_vcf_non_strict(self, vcf_without_annotations):
        """Test loading unannotated VCF in non-strict mode."""
        variants_df, genotypes_df, samples, report = load_vcf_with_quality_check(
            str(vcf_without_annotations), strict=False
        )

        assert len(variants_df) == 2
        assert report.gene_coverage == 0.0
        assert report.is_usable is False

    def test_load_unannotated_vcf_strict(self, vcf_without_annotations):
        """Test loading unannotated VCF in strict mode raises error."""
        with pytest.raises(VCFDataQualityError) as exc_info:
            load_vcf_with_quality_check(str(vcf_without_annotations), strict=True)

        assert "insufficient" in str(exc_info.value).lower()
        assert len(exc_info.value.fix_suggestions) > 0

    def test_load_multi_allelic(self, vcf_with_multi_allelic):
        """Test loading VCF with multi-allelic variants."""
        variants_df, genotypes_df, samples, report = load_vcf_with_quality_check(
            str(vcf_with_multi_allelic)
        )

        # Should expand 1 multi-allelic to 2 bi-allelic
        assert len(variants_df) == 2
        assert report.multi_allelic_variants == 1
        assert report.multi_allelic_expanded == 2

    def test_file_not_found(self, tmp_path):
        """Test appropriate error for missing file."""
        with pytest.raises(FileNotFoundError):
            load_vcf_with_quality_check(str(tmp_path / "nonexistent.vcf"))


class TestValidateVCFForPipeline:
    """Tests for the validation function."""

    @pytest.fixture
    def valid_vcf(self, tmp_path):
        """Create a valid VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Consequence">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10
chr1\t100\tv1\tA\tG\t99\tPASS\tGENE=A;CONSEQUENCE=missense;CADD=25\tGT\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0
chr1\t200\tv2\tC\tT\t99\tPASS\tGENE=B;CONSEQUENCE=stop;CADD=35\tGT\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1
"""
        vcf_path = tmp_path / "valid.vcf"
        vcf_path.write_text(vcf_content)
        return vcf_path

    def test_validate_valid_vcf(self, valid_vcf):
        """Test validation passes for valid VCF."""
        is_valid, report, suggestions = validate_vcf_for_pipeline(str(valid_vcf), verbose=False)

        assert is_valid is True
        assert report.is_usable is True
        assert len(suggestions) == 0 or report.cadd_coverage < 30

    def test_validate_missing_file(self, tmp_path):
        """Test validation handles missing file."""
        is_valid, report, suggestions = validate_vcf_for_pipeline(
            str(tmp_path / "missing.vcf"), verbose=False
        )

        assert is_valid is False
        assert len(report.errors) > 0


class TestVCFDataQualityError:
    """Tests for the VCFDataQualityError exception."""

    def test_error_message_format(self):
        """Test error message includes all components."""
        report = DataQualityReport()
        report.parsed_variants = 100
        report.variants_with_gene = 30  # Low coverage

        error = VCFDataQualityError(
            "Test error",
            report,
            ["Suggestion 1", "Suggestion 2"],
        )

        error_str = str(error)
        assert "Test error" in error_str
        assert "Suggestion 1" in error_str
        assert "Suggestion 2" in error_str

    def test_error_attributes(self):
        """Test error has correct attributes."""
        report = DataQualityReport()
        suggestions = ["Fix 1", "Fix 2"]

        error = VCFDataQualityError("Message", report, suggestions)

        assert error.message == "Message"
        assert error.report is report
        assert error.fix_suggestions == suggestions
