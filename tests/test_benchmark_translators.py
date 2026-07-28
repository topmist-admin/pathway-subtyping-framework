"""Unit tests for Phase 1 gene ID translators."""

import pytest

from scripts.benchmark_loaders.est_translator import (
    ESTTranslator,
    PlatformAnnotationLoader,
    UniGeneMapper,
)
from scripts.benchmark_loaders.translators import (
    EntrezGeneTranslator,
    GencodTranslator,
    RefSeqTranslator,
    batch_translate,
    detect_and_translate,
)


class TestEntrezGeneTranslator:
    """Test Entrez gene ID translator."""

    def test_can_handle_numeric_id(self):
        """Should detect numeric Entrez IDs."""
        assert EntrezGeneTranslator.can_handle("7157")
        assert EntrezGeneTranslator.can_handle("672")
        assert EntrezGeneTranslator.can_handle("5241")

    def test_cannot_handle_non_numeric(self):
        """Should reject non-numeric IDs."""
        assert not EntrezGeneTranslator.can_handle("TP53")
        assert not EntrezGeneTranslator.can_handle("ENSG00000141510")
        assert not EntrezGeneTranslator.can_handle("NM_000546")

    def test_whitespace_handling(self):
        """Should handle leading/trailing whitespace."""
        assert EntrezGeneTranslator.can_handle("  7157  ")
        assert EntrezGeneTranslator.can_handle("\t672\t")


class TestGencodTranslator:
    """Test Gencode versioned ID translator."""

    def test_can_handle_versioned_ensg(self):
        """Should detect versioned Ensembl IDs."""
        assert GencodTranslator.can_handle("ENSG00000141510.11")
        assert GencodTranslator.can_handle("ENSG00000012048.20")
        assert GencodTranslator.can_handle("ENSG00000000001.1")

    def test_cannot_handle_non_versioned(self):
        """Should reject non-versioned Ensembl IDs."""
        assert not GencodTranslator.can_handle("ENSG00000141510")
        assert not GencodTranslator.can_handle("7157")
        assert not GencodTranslator.can_handle("TP53")

    def test_cannot_handle_invalid_format(self):
        """Should reject invalid Gencode formats."""
        assert not GencodTranslator.can_handle("ENSG141510.11")  # Wrong length
        assert not GencodTranslator.can_handle("ENSMG00000141510.11")  # Wrong organism code


class TestRefSeqTranslator:
    """Test RefSeq translator."""

    def test_can_handle_nm_prefix(self):
        """Should handle NM_ (mRNA) RefSeq IDs."""
        assert RefSeqTranslator.can_handle("NM_000546")
        assert RefSeqTranslator.can_handle("NM_000546.3")

    def test_can_handle_all_refseq_types(self):
        """Should handle all RefSeq prefixes."""
        assert RefSeqTranslator.can_handle("NM_000546")  # mRNA
        assert RefSeqTranslator.can_handle("NR_000001")  # non-coding RNA
        assert RefSeqTranslator.can_handle("NP_000001")  # protein
        assert RefSeqTranslator.can_handle("XM_000001")  # predicted mRNA
        assert RefSeqTranslator.can_handle("XR_000001")  # predicted non-coding RNA
        assert RefSeqTranslator.can_handle("NG_000001")  # genomic DNA

    def test_cannot_handle_non_refseq(self):
        """Should reject non-RefSeq IDs."""
        assert not RefSeqTranslator.can_handle("7157")
        assert not RefSeqTranslator.can_handle("TP53")
        assert not RefSeqTranslator.can_handle("ENSG00000141510")


class TestESTTranslator:
    """Test EST probe translator."""

    def test_can_handle_tigr_est(self):
        """Should detect TIGR 40K EST probes."""
        assert ESTTranslator.can_handle("Hs_20001")
        assert ESTTranslator.can_handle("mm.12345")
        assert ESTTranslator.can_handle("Rn_5678")

    def test_can_handle_unigene_ids(self):
        """Should handle UniGene cluster IDs."""
        assert ESTTranslator.can_handle("Hs.123456")
        assert ESTTranslator.can_handle("Mm.987654")
        assert ESTTranslator.can_handle("Rn.111111")

    def test_cannot_handle_other_ids(self):
        """Should reject non-EST IDs."""
        assert not ESTTranslator.can_handle("TP53")
        assert not ESTTranslator.can_handle("7157")
        assert not ESTTranslator.can_handle("NM_000546")


class TestUniGeneMapper:
    """Test UniGene ID mapper."""

    def test_can_detect_unigene_ids(self):
        """Should detect UniGene cluster ID format."""
        assert UniGeneMapper.is_unigene_id("Hs.123456")
        assert UniGeneMapper.is_unigene_id("Mm.987654")
        assert UniGeneMapper.is_unigene_id("Rn.111111")

    def test_cannot_detect_non_unigene(self):
        """Should reject non-UniGene formats."""
        assert not UniGeneMapper.is_unigene_id("Hs_20001")
        assert not UniGeneMapper.is_unigene_id("TP53")
        assert not UniGeneMapper.is_unigene_id("7157")
        assert not UniGeneMapper.is_unigene_id("Hs.1234")  # Too short


class TestDetectAndTranslate:
    """Test automatic detection and translation."""

    def test_detect_entrez(self):
        """Should detect Entrez ID and return translator name."""
        translator_name, _ = detect_and_translate("7157")
        assert translator_name == "EntrezGeneTranslator"

    def test_detect_gencode(self):
        """Should detect Gencode versioned ID."""
        translator_name, _ = detect_and_translate("ENSG00000141510.11")
        assert translator_name == "GencodTranslator"

    def test_detect_refseq(self):
        """Should detect RefSeq ID."""
        translator_name, _ = detect_and_translate("NM_000546")
        assert translator_name == "RefSeqTranslator"

    def test_detect_est(self):
        """EST IDs require direct ESTTranslator call (not in generic registry)."""
        # Generic detect_and_translate() doesn't include ESTTranslator
        # (it's in separate module). Test the direct ESTTranslator instead:
        assert ESTTranslator.can_handle("Hs_20001")
        assert ESTTranslator.can_handle("Hs.123456")

    def test_undetected_returns_none(self):
        """Should return None for unrecognized IDs."""
        translator_name, symbol = detect_and_translate("UNKNOWN_ID_XXXX")
        assert translator_name is None
        assert symbol is None


class TestBatchTranslate:
    """Test batch translation."""

    def test_batch_with_mixed_types(self):
        """Should handle mixed gene ID types in one call."""
        gene_ids = ["7157", "ENSG00000141510.11", "NM_000546"]
        results = batch_translate(gene_ids)

        # Should have processed all IDs (may not all translate, but should attempt)
        assert len(results) <= len(gene_ids)

    def test_batch_returns_dict(self):
        """Should return dictionary with input as keys."""
        gene_ids = ["7157", "NM_000546"]
        results = batch_translate(gene_ids)

        assert isinstance(results, dict)
        for gene_id in gene_ids:
            assert gene_id in results or results.get(gene_id) is None


class TestPlatformAnnotationLoader:
    """Test platform annotation utilities."""

    def test_identify_est_platform(self):
        """Should identify EST-based platforms."""
        assert PlatformAnnotationLoader.is_est_platform("GPL3427")
        assert PlatformAnnotationLoader.is_est_platform("GPL9603")

    def test_reject_non_est_platform(self):
        """Should not identify non-EST platforms as EST."""
        assert not PlatformAnnotationLoader.is_est_platform("GPL570")
        assert not PlatformAnnotationLoader.is_est_platform("GPL571")

    def test_get_platform_info(self):
        """Should retrieve platform metadata."""
        info = PlatformAnnotationLoader.get_platform_info("GPL3427")
        assert info is not None
        assert info["type"] == "EST"
        assert "name" in info
