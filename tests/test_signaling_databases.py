"""Tests for pathway_subtyping.signaling_databases module."""

import json
from pathlib import Path
from unittest import mock

import pytest

from pathway_subtyping.signaling_databases import (
    SignalingDatabase,
    SignalingDatabaseResult,
    SignalingInteraction,
    _parse_cellphonedb_complexes,
    _parse_cellphonedb_genes,
    _resolve_partner_genes,
    convert_interactions_to_pathways,
    load_cellchatdb,
    load_cellphonedb,
    merge_signaling_databases,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_cellphonedb_interactions_csv(tmp_path):
    """Minimal CellPhoneDB interaction_input.csv."""
    content = (
        "id_cp_interaction,partner_a,partner_b,protein_name_a,protein_name_b,"
        "annotation_strategy,source,is_ppi,classification\n"
        "CPI-SS001,P12345,P67890,Protein_A,Protein_B,curated,literature,"
        "True,Adhesion by Cadherin\n"
        "CPI-SS002,P12345,complex_ABC,Protein_A,Complex_ABC,curated,"
        "literature,False,Adhesion by Cadherin\n"
        "CPI-SS003,P11111,P22222,Protein_C,Protein_D,curated,literature,"
        "True,WNT signaling\n"
        "CPI-SS004,P33333,P44444,Protein_E,Protein_F,curated,literature,"
        "True,Notch signaling\n"
        "CPI-SS005,P55555,P66666,Protein_G,Protein_H,curated,literature,"
        "True,\n"
    )
    path = tmp_path / "interaction_input.csv"
    path.write_text(content)
    return path


@pytest.fixture
def mock_cellphonedb_gene_csv(tmp_path):
    """Minimal CellPhoneDB gene_input.csv."""
    content = (
        "gene_name,uniprot,hgnc_symbol,ensembl\n"
        "CDH1,P12345,CDH1,ENSG00000039068\n"
        "CDH2,P67890,CDH2,ENSG00000170558\n"
        "WNT1,P11111,WNT1,ENSG00000125084\n"
        "FZD1,P22222,FZD1,ENSG00000157404\n"
        "DLL1,P33333,DLL1,ENSG00000198719\n"
        "NOTCH1,P44444,NOTCH1,ENSG00000148400\n"
        "ITGA1,P55555,ITGA1,ENSG00000213949\n"
        "ITGB1,P66666,ITGB1,ENSG00000150093\n"
        "SUBU1,P77777,SUBU1,ENSG00000000001\n"
        "SUBU2,P88888,SUBU2,ENSG00000000002\n"
        "SUBU3,P99999,SUBU3,ENSG00000000003\n"
    )
    path = tmp_path / "gene_input.csv"
    path.write_text(content)
    return path


@pytest.fixture
def mock_cellphonedb_complex_csv(tmp_path):
    """Minimal CellPhoneDB complex_input.csv."""
    content = (
        "complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n"
        "complex_ABC,P77777,P88888,P99999,,\n"
    )
    path = tmp_path / "complex_input.csv"
    path.write_text(content)
    return path


@pytest.fixture
def mock_cellchatdb_csv(tmp_path):
    """Minimal CellChatDB-exported CSV."""
    content = (
        "interaction_name,pathway_name,ligand,receptor,annotation,evidence\n"
        "COL1A1_CD44,COLLAGEN,COL1A1,CD44,Secreted Signaling,KEGG\n"
        "COL1A2_CD44,COLLAGEN,COL1A2,CD44,Secreted Signaling,KEGG\n"
        "WNT1_FZD1,WNT,WNT1,FZD1,Secreted Signaling,KEGG\n"
        "DLL1_NOTCH1,NOTCH,DLL1,NOTCH1,Cell-Cell Contact,KEGG\n"
        "JAG1_NOTCH1,NOTCH,JAG1,NOTCH1,Cell-Cell Contact,KEGG\n"
    )
    path = tmp_path / "cellchatdb_human.csv"
    path.write_text(content)
    return path


@pytest.fixture
def sample_interactions():
    """Pre-built list of SignalingInteraction objects."""
    return [
        SignalingInteraction(
            interaction_id="I1",
            partner_a="CDH1",
            partner_b="CDH2",
            genes_a=["CDH1"],
            genes_b=["CDH2"],
            pathway_name="Adhesion",
            source_database=SignalingDatabase.CELLPHONEDB,
            annotation="curated",
        ),
        SignalingInteraction(
            interaction_id="I2",
            partner_a="WNT1",
            partner_b="FZD1",
            genes_a=["WNT1"],
            genes_b=["FZD1"],
            pathway_name="WNT",
            source_database=SignalingDatabase.CELLPHONEDB,
            annotation="curated",
        ),
        SignalingInteraction(
            interaction_id="I3",
            partner_a="WNT2",
            partner_b="FZD2",
            genes_a=["WNT2"],
            genes_b=["FZD2"],
            pathway_name="WNT",
            source_database=SignalingDatabase.CELLPHONEDB,
            annotation="curated",
        ),
    ]


# ---------------------------------------------------------------------------
# Tests: SignalingDatabase enum
# ---------------------------------------------------------------------------


class TestSignalingDatabase:
    def test_cellphonedb_value(self):
        assert SignalingDatabase.CELLPHONEDB.value == "cellphonedb"

    def test_cellchatdb_value(self):
        assert SignalingDatabase.CELLCHATDB.value == "cellchatdb"

    def test_from_string(self):
        assert SignalingDatabase("cellphonedb") == SignalingDatabase.CELLPHONEDB
        assert SignalingDatabase("cellchatdb") == SignalingDatabase.CELLCHATDB

    def test_invalid_value(self):
        with pytest.raises(ValueError):
            SignalingDatabase("invalid_db")


# ---------------------------------------------------------------------------
# Tests: SignalingInteraction
# ---------------------------------------------------------------------------


class TestSignalingInteraction:
    def test_to_dict(self):
        interaction = SignalingInteraction(
            interaction_id="I1",
            partner_a="P12345",
            partner_b="P67890",
            genes_a=["CDH1"],
            genes_b=["CDH2"],
            pathway_name="Adhesion",
            source_database=SignalingDatabase.CELLPHONEDB,
            annotation="curated",
        )
        d = interaction.to_dict()
        assert d["interaction_id"] == "I1"
        assert d["source_database"] == "cellphonedb"
        assert d["genes_a"] == ["CDH1"]
        assert d["pathway_name"] == "Adhesion"

    def test_all_genes_property(self):
        interaction = SignalingInteraction(
            interaction_id="I1",
            partner_a="A",
            partner_b="B",
            genes_a=["CDH1", "CDH2"],
            genes_b=["CDH2", "CDH3"],
            pathway_name="Test",
            source_database=SignalingDatabase.CELLPHONEDB,
        )
        assert interaction.all_genes == ["CDH1", "CDH2", "CDH3"]

    def test_all_genes_sorted(self):
        interaction = SignalingInteraction(
            interaction_id="I1",
            partner_a="A",
            partner_b="B",
            genes_a=["ZZZ"],
            genes_b=["AAA"],
            pathway_name="Test",
            source_database=SignalingDatabase.CELLCHATDB,
        )
        assert interaction.all_genes == ["AAA", "ZZZ"]

    def test_json_serializable(self):
        interaction = SignalingInteraction(
            interaction_id="I1",
            partner_a="A",
            partner_b="B",
            genes_a=["CDH1"],
            genes_b=["CDH2"],
            pathway_name="Test",
            source_database=SignalingDatabase.CELLPHONEDB,
        )
        json_str = json.dumps(interaction.to_dict())
        assert isinstance(json_str, str)


# ---------------------------------------------------------------------------
# Tests: SignalingDatabaseResult
# ---------------------------------------------------------------------------


class TestSignalingDatabaseResult:
    def test_to_dict(self):
        result = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD1", "WNT1"]},
            n_interactions=5,
            n_pathways=1,
            n_unique_genes=2,
        )
        d = result.to_dict()
        assert d["database"] == "cellphonedb"
        assert d["n_interactions"] == 5
        assert d["n_pathways"] == 1
        assert "Signaling: WNT" in d["pathway_names"]

    def test_format_report(self):
        result = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD1", "WNT1"]},
            n_interactions=5,
            n_pathways=1,
            n_unique_genes=2,
        )
        report = result.format_report()
        assert "Signaling Database Report" in report
        assert "cellphonedb" in report
        assert "Signaling: WNT" in report

    def test_format_report_with_warnings(self):
        result = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={},
            warnings=["Some warning"],
        )
        report = result.format_report()
        assert "Warnings" in report
        assert "Some warning" in report

    def test_get_citations_cellphonedb(self):
        result = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={},
        )
        citations = result.get_citations()
        assert len(citations) == 1
        assert "Efremova" in citations[0]

    def test_get_citations_cellchatdb(self):
        result = SignalingDatabaseResult(
            database=SignalingDatabase.CELLCHATDB,
            pathway_gene_sets={},
        )
        citations = result.get_citations()
        assert len(citations) == 1
        assert "Jin" in citations[0]

    def test_json_serializable(self):
        result = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD1", "WNT1"]},
            n_interactions=5,
            n_pathways=1,
            n_unique_genes=2,
            warnings=["test"],
        )
        json_str = json.dumps(result.to_dict())
        assert isinstance(json_str, str)


# ---------------------------------------------------------------------------
# Tests: _parse_cellphonedb_genes
# ---------------------------------------------------------------------------


class TestParseCellPhoneDBGenes:
    def test_parse_gene_csv(self, mock_cellphonedb_gene_csv):
        mapping = _parse_cellphonedb_genes(mock_cellphonedb_gene_csv)
        assert len(mapping) == 11
        assert mapping["P12345"] == "CDH1"
        assert mapping["P77777"] == "SUBU1"

    def test_hgnc_preferred_over_gene_name(self, tmp_path):
        content = (
            "gene_name,uniprot,hgnc_symbol,ensembl\n" "OldName,P00001,NewName,ENSG00000000000\n"
        )
        path = tmp_path / "gene.csv"
        path.write_text(content)
        mapping = _parse_cellphonedb_genes(path)
        assert mapping["P00001"] == "NewName"

    def test_fallback_to_gene_name(self, tmp_path):
        content = "gene_name,uniprot,hgnc_symbol,ensembl\n" "FallbackGene,P00002,,ENSG00000000000\n"
        path = tmp_path / "gene.csv"
        path.write_text(content)
        mapping = _parse_cellphonedb_genes(path)
        assert mapping["P00002"] == "FallbackGene"

    def test_empty_file(self, tmp_path):
        path = tmp_path / "gene.csv"
        path.write_text("gene_name,uniprot,hgnc_symbol,ensembl\n")
        mapping = _parse_cellphonedb_genes(path)
        assert mapping == {}


# ---------------------------------------------------------------------------
# Tests: _parse_cellphonedb_complexes
# ---------------------------------------------------------------------------


class TestParseCellPhoneDBComplexes:
    def test_parse_complex_csv(self, mock_cellphonedb_complex_csv, mock_cellphonedb_gene_csv):
        gene_map = _parse_cellphonedb_genes(mock_cellphonedb_gene_csv)
        complexes = _parse_cellphonedb_complexes(mock_cellphonedb_complex_csv, gene_map)
        assert "complex_ABC" in complexes
        assert sorted(complexes["complex_ABC"]) == ["SUBU1", "SUBU2", "SUBU3"]

    def test_partial_subunits(self, tmp_path):
        content = (
            "complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n"
            "partial_complex,P12345,UNKNOWN_ID,,,\n"
        )
        path = tmp_path / "complex.csv"
        path.write_text(content)
        gene_map = {"P12345": "CDH1"}
        complexes = _parse_cellphonedb_complexes(path, gene_map)
        assert complexes["partial_complex"] == ["CDH1"]

    def test_empty_file(self, tmp_path):
        path = tmp_path / "complex.csv"
        path.write_text("complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n")
        complexes = _parse_cellphonedb_complexes(path, {})
        assert complexes == {}


# ---------------------------------------------------------------------------
# Tests: _resolve_partner_genes
# ---------------------------------------------------------------------------


class TestResolvePartnerGenes:
    def test_uniprot_lookup(self):
        genes = _resolve_partner_genes("P12345", {"P12345": "CDH1"}, {})
        assert genes == ["CDH1"]

    def test_complex_lookup(self):
        genes = _resolve_partner_genes("my_complex", {}, {"my_complex": ["GENE1", "GENE2"]})
        assert genes == ["GENE1", "GENE2"]

    def test_gene_name_fallback(self):
        genes = _resolve_partner_genes("MYC", {}, {})
        assert genes == ["MYC"]

    def test_empty_partner(self):
        genes = _resolve_partner_genes("", {}, {})
        assert genes == []

    def test_colon_prefix(self):
        genes = _resolve_partner_genes("simple:BRCA1", {}, {})
        assert genes == ["BRCA1"]

    def test_complex_preferred_over_uniprot(self):
        """Complex lookup takes priority over UniProt lookup."""
        genes = _resolve_partner_genes(
            "shared_id",
            {"shared_id": "SINGLE_GENE"},
            {"shared_id": ["GENE1", "GENE2"]},
        )
        assert genes == ["GENE1", "GENE2"]


# ---------------------------------------------------------------------------
# Tests: load_cellphonedb
# ---------------------------------------------------------------------------


class TestLoadCellPhoneDB:
    def test_load_from_files(
        self,
        mock_cellphonedb_interactions_csv,
        mock_cellphonedb_gene_csv,
        mock_cellphonedb_complex_csv,
    ):
        result = load_cellphonedb(
            interactions_path=mock_cellphonedb_interactions_csv,
            gene_path=mock_cellphonedb_gene_csv,
            complex_path=mock_cellphonedb_complex_csv,
        )
        assert isinstance(result, SignalingDatabaseResult)
        assert result.database == SignalingDatabase.CELLPHONEDB
        assert result.n_interactions > 0
        assert result.n_pathways > 0

    def test_pathway_gene_sets_format(
        self,
        mock_cellphonedb_interactions_csv,
        mock_cellphonedb_gene_csv,
        mock_cellphonedb_complex_csv,
    ):
        result = load_cellphonedb(
            interactions_path=mock_cellphonedb_interactions_csv,
            gene_path=mock_cellphonedb_gene_csv,
            complex_path=mock_cellphonedb_complex_csv,
        )
        for name, genes in result.pathway_gene_sets.items():
            assert isinstance(name, str)
            assert isinstance(genes, list)
            assert all(isinstance(g, str) for g in genes)

    def test_complex_resolution(
        self,
        mock_cellphonedb_interactions_csv,
        mock_cellphonedb_gene_csv,
        mock_cellphonedb_complex_csv,
    ):
        result = load_cellphonedb(
            interactions_path=mock_cellphonedb_interactions_csv,
            gene_path=mock_cellphonedb_gene_csv,
            complex_path=mock_cellphonedb_complex_csv,
        )
        # CPI-SS002 has complex_ABC as partner_b → should resolve to SUBU1-3
        cadherin_genes = result.pathway_gene_sets.get("Signaling: Adhesion by Cadherin", [])
        assert "SUBU1" in cadherin_genes or "SUBU2" in cadherin_genes

    def test_signaling_prefix(
        self,
        mock_cellphonedb_interactions_csv,
        mock_cellphonedb_gene_csv,
        mock_cellphonedb_complex_csv,
    ):
        result = load_cellphonedb(
            interactions_path=mock_cellphonedb_interactions_csv,
            gene_path=mock_cellphonedb_gene_csv,
            complex_path=mock_cellphonedb_complex_csv,
        )
        for name in result.pathway_gene_sets:
            assert name.startswith("Signaling: ")

    def test_empty_classification_skipped(
        self,
        mock_cellphonedb_interactions_csv,
        mock_cellphonedb_gene_csv,
        mock_cellphonedb_complex_csv,
    ):
        result = load_cellphonedb(
            interactions_path=mock_cellphonedb_interactions_csv,
            gene_path=mock_cellphonedb_gene_csv,
            complex_path=mock_cellphonedb_complex_csv,
        )
        # CPI-SS005 has empty classification → should be in warnings
        assert any("empty classification" in w for w in result.warnings)

    def test_min_genes_filter(
        self,
        mock_cellphonedb_interactions_csv,
        mock_cellphonedb_gene_csv,
        mock_cellphonedb_complex_csv,
    ):
        result = load_cellphonedb(
            interactions_path=mock_cellphonedb_interactions_csv,
            gene_path=mock_cellphonedb_gene_csv,
            complex_path=mock_cellphonedb_complex_csv,
            min_genes_per_pathway=10,
        )
        # With min_genes=10, most small pathways should be filtered out
        # Adhesion by Cadherin has CDH1, CDH2, SUBU1, SUBU2, SUBU3 = 5 genes
        assert result.n_pathways == 0 or all(
            len(genes) >= 10 for genes in result.pathway_gene_sets.values()
        )

    def test_download_with_mock(self, tmp_path):
        # Create mock files in the cache dir
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        interactions = (
            "id_cp_interaction,partner_a,partner_b,protein_name_a,"
            "protein_name_b,annotation_strategy,source,is_ppi,classification\n"
            "CPI-001,P11,P22,A,B,curated,lit,True,TestPathway\n"
        )
        genes = (
            "gene_name,uniprot,hgnc_symbol,ensembl\n"
            "GENE1,P11,GENE1,ENSG1\n"
            "GENE2,P22,GENE2,ENSG2\n"
        )
        complexes = "complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n"

        def mock_urlretrieve(url, path):
            if "interaction" in url:
                Path(path).write_text(interactions)
            elif "gene" in url:
                Path(path).write_text(genes)
            elif "complex" in url:
                Path(path).write_text(complexes)

        with mock.patch(
            "pathway_subtyping.signaling_databases.urllib.request.urlretrieve",
            side_effect=mock_urlretrieve,
        ):
            result = load_cellphonedb(cache_dir=cache_dir)
            assert result.n_interactions == 1
            assert "Signaling: TestPathway" in result.pathway_gene_sets

    def test_cache_hit_skips_download(self, tmp_path):
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        # Pre-populate cache
        interactions_content = (
            "id_cp_interaction,partner_a,partner_b,protein_name_a,"
            "protein_name_b,annotation_strategy,source,is_ppi,classification\n"
            "CPI-001,P11,P22,A,B,curated,lit,True,CachedPathway\n"
        )
        genes_content = "gene_name,uniprot,hgnc_symbol,ensembl\n" "G1,P11,G1,E1\nG2,P22,G2,E2\n"
        complexes_content = "complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n"
        (cache_dir / "cellphonedb_interaction_input.csv").write_text(interactions_content)
        (cache_dir / "cellphonedb_gene_input.csv").write_text(genes_content)
        (cache_dir / "cellphonedb_complex_input.csv").write_text(complexes_content)

        with mock.patch(
            "pathway_subtyping.signaling_databases.urllib.request.urlretrieve"
        ) as mock_dl:
            result = load_cellphonedb(cache_dir=cache_dir)
            mock_dl.assert_not_called()
            assert "Signaling: CachedPathway" in result.pathway_gene_sets

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_cellphonedb(
                interactions_path=tmp_path / "nonexistent.csv",
                gene_path=tmp_path / "gene.csv",
                complex_path=tmp_path / "complex.csv",
            )

    def test_missing_columns_raises(self, tmp_path):
        bad_csv = tmp_path / "bad.csv"
        bad_csv.write_text("col_a,col_b\nval1,val2\n")
        gene_csv = tmp_path / "gene.csv"
        gene_csv.write_text("gene_name,uniprot,hgnc_symbol,ensembl\n")
        complex_csv = tmp_path / "complex.csv"
        complex_csv.write_text("complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n")
        with pytest.raises(ValueError, match="Missing required columns"):
            load_cellphonedb(
                interactions_path=bad_csv,
                gene_path=gene_csv,
                complex_path=complex_csv,
            )

    def test_empty_interactions_raises(self, tmp_path):
        interactions_csv = tmp_path / "interactions.csv"
        interactions_csv.write_text(
            "id_cp_interaction,partner_a,partner_b,protein_name_a,"
            "protein_name_b,annotation_strategy,source,is_ppi,classification\n"
        )
        gene_csv = tmp_path / "gene.csv"
        gene_csv.write_text("gene_name,uniprot,hgnc_symbol,ensembl\n")
        complex_csv = tmp_path / "complex.csv"
        complex_csv.write_text("complex_name,uniprot_1,uniprot_2,uniprot_3,uniprot_4,uniprot_5\n")
        with pytest.raises(ValueError, match="No valid interactions"):
            load_cellphonedb(
                interactions_path=interactions_csv,
                gene_path=gene_csv,
                complex_path=complex_csv,
            )


# ---------------------------------------------------------------------------
# Tests: load_cellchatdb
# ---------------------------------------------------------------------------


class TestLoadCellChatDB:
    def test_load_from_csv(self, mock_cellchatdb_csv):
        result = load_cellchatdb(mock_cellchatdb_csv)
        assert isinstance(result, SignalingDatabaseResult)
        assert result.database == SignalingDatabase.CELLCHATDB
        assert result.n_interactions == 5

    def test_pathway_gene_sets_keys(self, mock_cellchatdb_csv):
        result = load_cellchatdb(mock_cellchatdb_csv)
        expected_keys = {
            "Signaling: COLLAGEN",
            "Signaling: WNT",
            "Signaling: NOTCH",
        }
        assert set(result.pathway_gene_sets.keys()) == expected_keys

    def test_gene_resolution(self, mock_cellchatdb_csv):
        result = load_cellchatdb(mock_cellchatdb_csv)
        collagen_genes = result.pathway_gene_sets["Signaling: COLLAGEN"]
        assert "COL1A1" in collagen_genes
        assert "COL1A2" in collagen_genes
        assert "CD44" in collagen_genes

    def test_signaling_prefix(self, mock_cellchatdb_csv):
        result = load_cellchatdb(mock_cellchatdb_csv)
        for name in result.pathway_gene_sets:
            assert name.startswith("Signaling: ")

    def test_min_genes_filter(self, mock_cellchatdb_csv):
        result = load_cellchatdb(mock_cellchatdb_csv, min_genes_per_pathway=5)
        # No pathway in mock has 5+ genes
        assert result.n_pathways == 0

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_cellchatdb(tmp_path / "nonexistent.csv")

    def test_missing_columns_raises(self, tmp_path):
        bad_csv = tmp_path / "bad.csv"
        bad_csv.write_text("col_a,col_b\nval1,val2\n")
        with pytest.raises(ValueError, match="Missing required columns"):
            load_cellchatdb(bad_csv)

    def test_empty_csv_raises(self, tmp_path):
        empty_csv = tmp_path / "empty.csv"
        empty_csv.write_text("interaction_name,pathway_name,ligand,receptor,annotation\n")
        with pytest.raises(ValueError, match="No valid interactions"):
            load_cellchatdb(empty_csv)

    def test_interactions_stored(self, mock_cellchatdb_csv):
        result = load_cellchatdb(mock_cellchatdb_csv)
        assert len(result.interactions) == 5
        assert all(i.source_database == SignalingDatabase.CELLCHATDB for i in result.interactions)


# ---------------------------------------------------------------------------
# Tests: convert_interactions_to_pathways
# ---------------------------------------------------------------------------


class TestConvertInteractionsToPathways:
    def test_group_by_pathway_name(self, sample_interactions):
        pathways = convert_interactions_to_pathways(sample_interactions)
        assert "Signaling: WNT" in pathways
        assert "Signaling: Adhesion" in pathways

    def test_group_by_annotation(self, sample_interactions):
        pathways = convert_interactions_to_pathways(
            sample_interactions, group_by="annotation", prefix=""
        )
        assert "curated" in pathways

    def test_min_genes_filter(self, sample_interactions):
        pathways = convert_interactions_to_pathways(sample_interactions, min_genes=3)
        # WNT has 4 genes (WNT1, FZD1, WNT2, FZD2), Adhesion has 2
        assert "Signaling: WNT" in pathways
        assert "Signaling: Adhesion" not in pathways

    def test_genes_sorted_and_unique(self, sample_interactions):
        pathways = convert_interactions_to_pathways(sample_interactions)
        for genes in pathways.values():
            assert genes == sorted(set(genes))

    def test_empty_list(self):
        pathways = convert_interactions_to_pathways([])
        assert pathways == {}

    def test_custom_prefix(self, sample_interactions):
        pathways = convert_interactions_to_pathways(sample_interactions, prefix="CellSignal")
        for name in pathways:
            assert name.startswith("CellSignal: ")

    def test_no_prefix(self, sample_interactions):
        pathways = convert_interactions_to_pathways(sample_interactions, prefix="")
        assert "WNT" in pathways
        assert "Adhesion" in pathways


# ---------------------------------------------------------------------------
# Tests: merge_signaling_databases
# ---------------------------------------------------------------------------


class TestMergeSignalingDatabases:
    def test_merge_disjoint(self):
        r1 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD1", "WNT1"]},
        )
        r2 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLCHATDB,
            pathway_gene_sets={"Signaling: NOTCH": ["DLL1", "NOTCH1"]},
        )
        merged = merge_signaling_databases(r1, r2)
        assert "Signaling: WNT" in merged
        assert "Signaling: NOTCH" in merged

    def test_merge_overlapping_unions_genes(self):
        r1 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD1", "WNT1"]},
        )
        r2 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLCHATDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD2", "WNT1"]},
        )
        merged = merge_signaling_databases(r1, r2)
        wnt_genes = merged["Signaling: WNT"]
        assert "FZD1" in wnt_genes
        assert "FZD2" in wnt_genes
        assert "WNT1" in wnt_genes
        assert len(wnt_genes) == 3

    def test_merge_single_result(self):
        r1 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: WNT": ["FZD1", "WNT1"]},
        )
        merged = merge_signaling_databases(r1)
        assert merged == {"Signaling: WNT": ["FZD1", "WNT1"]}

    def test_merge_empty(self):
        merged = merge_signaling_databases()
        assert merged == {}

    def test_merged_genes_sorted(self):
        r1 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLPHONEDB,
            pathway_gene_sets={"Signaling: P": ["ZZZ", "AAA"]},
        )
        r2 = SignalingDatabaseResult(
            database=SignalingDatabase.CELLCHATDB,
            pathway_gene_sets={"Signaling: P": ["MMM"]},
        )
        merged = merge_signaling_databases(r1, r2)
        assert merged["Signaling: P"] == ["AAA", "MMM", "ZZZ"]
