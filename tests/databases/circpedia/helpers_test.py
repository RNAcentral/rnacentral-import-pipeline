# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
Tests for CIRCpedia helper functions.
"""

import pytest
from rnacentral_pipeline.databases.circpedia import helpers


class TestParseGenomicLocation:
    """Test parsing of genomic location strings."""

    def test_parses_standard_location(self):
        """Test parsing standard chromosome location."""
        result = helpers.parse_genomic_location("chr1:12345-67890")
        assert result == {
            "chromosome": "1",
            "start": 12345,
            "end": 67890,
        }

    def test_parses_location_without_chr_prefix(self):
        """Test parsing location without 'chr' prefix."""
        result = helpers.parse_genomic_location("X:1000-2000")
        assert result == {
            "chromosome": "X",
            "start": 1000,
            "end": 2000,
        }

    def test_parses_mitochondrial_chromosome(self):
        """Test parsing mitochondrial chromosome."""
        result = helpers.parse_genomic_location("chrM:100-500")
        assert result == {
            "chromosome": "M",
            "start": 100,
            "end": 500,
        }

    def test_returns_none_for_invalid_format(self):
        """Test that invalid formats return None."""
        assert helpers.parse_genomic_location("invalid") is None
        assert helpers.parse_genomic_location("chr1:123") is None
        assert helpers.parse_genomic_location("") is None
        assert helpers.parse_genomic_location(None) is None


class TestParseExonPositions:
    """Test parsing of exon position strings."""

    def test_parses_single_exon(self):
        """Test parsing single exon."""
        result = helpers.parse_exon_positions("100-200")
        assert result == [(100, 200)]

    def test_parses_multiple_exons(self):
        """Test parsing multiple exons."""
        result = helpers.parse_exon_positions("100-200,300-400,500-600")
        assert result == [(100, 200), (300, 400), (500, 600)]

    def test_parses_with_spaces(self):
        """Test parsing with extra whitespace."""
        result = helpers.parse_exon_positions(" 100-200 , 300-400 ")
        assert result == [(100, 200), (300, 400)]

    def test_returns_empty_for_invalid(self):
        """Test that invalid input returns empty list."""
        assert helpers.parse_exon_positions("") == []
        assert helpers.parse_exon_positions(None) == []
        assert helpers.parse_exon_positions("invalid") == []


class TestParseStrand:
    """Test parsing of strand information."""

    def test_parses_forward_strand(self):
        """Test parsing forward strand."""
        assert helpers.parse_strand("+") == 1

    def test_parses_reverse_strand(self):
        """Test parsing reverse strand."""
        assert helpers.parse_strand("-") == -1

    def test_returns_zero_for_unknown(self):
        """Test that unknown strand returns 0."""
        assert helpers.parse_strand("") == 0
        assert helpers.parse_strand(None) == 0
        assert helpers.parse_strand(".") == 0


class TestPrimaryId:
    """Test primary ID generation."""

    def test_generates_primary_id(self):
        """Test that primary ID is generated correctly."""
        assert helpers.primary_id("circ123") == "circ123"


class TestAccession:
    """Test accession generation."""

    def test_generates_consistent_accession(self):
        """Test that accession is consistent for same input."""
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        acc1 = helpers.accession("circ123", location)
        acc2 = helpers.accession("circ123", location)
        assert acc1 == acc2
        assert acc1 == "CIRCPEDIA:circ123_1:1000-2000"

    def test_generates_different_accessions_for_different_locations(self):
        """Test that different locations produce different accessions."""
        loc1 = {"chromosome": "1", "start": 1000, "end": 2000}
        loc2 = {"chromosome": "1", "start": 3000, "end": 4000}
        acc1 = helpers.accession("circ123", loc1)
        acc2 = helpers.accession("circ123", loc2)
        assert acc1 != acc2
        assert acc1 == "CIRCPEDIA:circ123_1:1000-2000"
        assert acc2 == "CIRCPEDIA:circ123_1:3000-4000"

    def test_includes_database_prefix(self):
        """Test that accession is prefixed with database name."""
        location = {"chromosome": "X", "start": 5000, "end": 6000}
        acc = helpers.accession("hsa_circ_0001", location)
        assert acc.startswith("CIRCPEDIA:")
        assert "hsa_circ_0001" in acc


class TestUrl:
    """Test URL generation."""

    def test_generates_url(self):
        """Test that URL is generated correctly."""
        url = helpers.url("circ123")
        assert "circpediav3" in url.lower()
        assert "circ123" in url
        assert "/circrna/" in url.lower()

    def test_generates_direct_page_url(self):
        """Test that URL points to direct circRNA page, not search."""
        url = helpers.url("hsa_circ_0001")
        assert url == "https://bits.fudan.edu.cn/circpediav3/circrna/hsa_circ_0001"
        assert "search" not in url.lower()


class TestSeqVersion:
    """Test sequence version."""

    def test_returns_version_one(self):
        """Test that version 1 is returned."""
        row = {}
        assert helpers.seq_version(row) == "1"


class TestNoteData:
    """Test note data generation."""

    def test_includes_genomic_location(self):
        """Test that genomic location is included in notes."""
        row = {}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data(row, location)
        assert notes["genomic_location"] == "1:1000-2000"

    def test_includes_expression_data(self):
        """Test that expression data is included when available."""
        row = {"fpm": 10.5}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data(row, location)
        assert notes["expression_fpm"] == 10.5

    def test_includes_cell_line(self):
        """Test that cell line is included when available."""
        row = {"cell_line": "HeLa"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data(row, location)
        assert notes["cell_line"] == "HeLa"


class TestChromosome:
    """Test chromosome extraction."""

    def test_extracts_chromosome(self):
        """Test that chromosome is extracted correctly."""
        location = {"chromosome": "X", "start": 1000, "end": 2000}
        assert helpers.chromosome(location) == "X"


class TestGene:
    """Test gene extraction."""

    def test_extracts_gene_name(self):
        """Test that gene name is extracted when present."""
        row = {"gene": "TP53"}
        assert helpers.gene(row) == "TP53"

    def test_returns_none_when_missing(self):
        """Test that None is returned when gene is missing."""
        row = {}
        assert helpers.gene(row) is None

    def test_returns_none_for_empty_gene(self):
        """Test that None is returned for empty gene."""
        row = {"gene": None}
        assert helpers.gene(row) is None


class TestProduct:
    """Test product description generation."""

    def test_generates_product_with_gene(self):
        """Test product description with gene name."""
        row = {}
        product = helpers.product(row, "TP53")
        assert product == "TP53 circular RNA"

    def test_generates_product_without_gene(self):
        """Test product description without gene name."""
        row = {}
        product = helpers.product(row, None)
        assert product == "circular RNA"


class TestReferences:
    """Test reference generation."""

    def test_generates_references(self):
        """Test that references are generated."""
        refs = helpers.references()
        assert len(refs) > 0
        assert any("CIRCpedia" in ref.location for ref in refs)


class TestParseLocationField:
    """Test parsing of combined location field from CIRCpedia V3."""

    def test_parses_location_with_strand(self):
        """Test parsing location with chromosome, coordinates, and strand."""
        result = helpers.parse_location_field("V:15874634-15876408(-)")
        assert result == {
            "chromosome": "V",
            "start": 15874634,
            "end": 15876408,
            "strand": -1,
        }

    def test_parses_location_with_positive_strand(self):
        """Test parsing location with positive strand."""
        result = helpers.parse_location_field("chr1:1000-2000(+)")
        assert result == {
            "chromosome": "1",
            "start": 1000,
            "end": 2000,
            "strand": 1,
        }

    def test_parses_chromosome_without_chr_prefix(self):
        """Test parsing chromosome name without 'chr' prefix."""
        result = helpers.parse_location_field("III:5000-6000(+)")
        assert result == {
            "chromosome": "III",
            "start": 5000,
            "end": 6000,
            "strand": 1,
        }

    def test_parses_mitochondrial_chromosome(self):
        """Test parsing mitochondrial chromosome."""
        result = helpers.parse_location_field("chrM:100-500(-)")
        assert result == {
            "chromosome": "M",
            "start": 100,
            "end": 500,
            "strand": -1,
        }

    def test_returns_none_for_invalid_format(self):
        """Test that invalid format returns None."""
        assert helpers.parse_location_field("invalid") is None
        assert helpers.parse_location_field("chr1:100-200") is None  # Missing strand
        assert helpers.parse_location_field("chr1:100(+)") is None  # Missing end
        assert helpers.parse_location_field("") is None
        assert helpers.parse_location_field(None) is None


class TestNoteDataFromTsv:
    """Test note data extraction from TSV rows."""

    def test_includes_genomic_location(self):
        """Test that genomic location is included."""
        row = {}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["genomic_location"] == "1:1000-2000"

    def test_includes_circname(self):
        """Test that circname is included when present."""
        row = {"circname": "circ-TP53"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["circname"] == "circ-TP53"

    def test_includes_length(self):
        """Test that length is included when present."""
        row = {"length": "1500"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["length"] == 1500

    def test_includes_subcellular_location(self):
        """Test that subcellular location is included when present."""
        row = {"subcell_location": "cytoplasm"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["subcell_location"] == "cytoplasm"

    def test_filters_none_values(self):
        """Test that editing sites with 'none' are filtered."""
        row = {"editing_site": "none", "DIS3_signal": "none"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert "editing_site" not in notes
        assert "DIS3_signal" not in notes

    def test_includes_editing_sites(self):
        """Test that valid editing sites are included."""
        row = {"editing_site": "A>I;C>U"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["editing_site"] == "A>I;C>U"

    def test_includes_transcript_ids(self):
        """Test that transcript IDs are included when present."""
        row = {"transcript_Ensembl": "ENST123", "transcript_Refseq": "NM_123"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["transcript_ensembl"] == "ENST123"
        assert notes["transcript_refseq"] == "NM_123"

    def test_filters_na_transcript_ids(self):
        """Test that 'NA' transcript IDs are filtered."""
        row = {"transcript_Ensembl": "NA", "transcript_Refseq": "NA"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert "transcript_ensembl" not in notes
        assert "transcript_refseq" not in notes


class TestProductFromGene:
    """Test product description generation from gene name."""

    def test_generates_product_with_gene(self):
        """Test product description with gene name."""
        product = helpers.product_from_gene("TP53")
        assert product == "TP53 circular RNA"

    def test_generates_product_without_gene(self):
        """Test product description without gene name."""
        product = helpers.product_from_gene(None)
        assert product == "circular RNA"
