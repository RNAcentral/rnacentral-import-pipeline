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


class TestPrimaryId:
    """Test primary ID generation."""

    def test_generates_primary_id(self):
        """Test that primary ID matches the circRNA ID."""
        assert helpers.primary_id("hsa_circ_0001") == "hsa_circ_0001"


class TestAccession:
    """Test accession generation."""

    def test_generates_consistent_accession(self):
        """Test accession generation with database prefix."""
        location = {"chromosome": "1", "start": 100, "end": 200}
        acc = helpers.accession("hsa_circ_0001", location)
        assert acc == "CIRCPEDIA:hsa_circ_0001_1:100-200"

    def test_generates_different_accessions_for_different_locations(self):
        """Test that different locations produce different accessions."""
        location1 = {"chromosome": "1", "start": 100, "end": 200}
        location2 = {"chromosome": "1", "start": 300, "end": 400}
        acc1 = helpers.accession("hsa_circ_0001", location1)
        acc2 = helpers.accession("hsa_circ_0001", location2)
        assert acc1 != acc2

    def test_includes_database_prefix(self):
        """Test that accession includes CIRCPEDIA prefix."""
        location = {"chromosome": "1", "start": 100, "end": 200}
        acc = helpers.accession("hsa_circ_0001", location)
        assert acc.startswith("CIRCPEDIA:")


class TestUrl:
    """Test URL generation."""

    def test_generates_url(self):
        """Test URL generation for circRNA."""
        url = helpers.url("hsa_circ_0001")
        assert "https://bits.fudan.edu.cn/circpediav3/circrna/hsa_circ_0001" == url

    def test_generates_direct_page_url(self):
        """Test that URL points to direct circRNA page, not search."""
        url = helpers.url("hsa_circ_0001")
        assert "/circrna/" in url
        assert "search" not in url


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

    def test_includes_dis3_motif(self):
        """Test that DIS3 motif is included when present."""
        row = {"DIS3_motif": "AUUUA"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["DIS3_motif"] == "AUUUA"

    def test_includes_orthology(self):
        """Test that orthology information is included when present."""
        row = {"Orthology": "mouse;rat"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["orthology"] == "mouse;rat"

    def test_includes_tgs_support(self):
        """Test that TGS support is included when present."""
        row = {"TGS": "Nanopore"}
        location = {"chromosome": "1", "start": 1000, "end": 2000}
        notes = helpers.note_data_from_tsv(row, location)
        assert notes["TGS_support"] == "Nanopore"


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


class TestReferences:
    """Test reference generation."""

    def test_generates_references(self):
        """Test that references are generated."""
        refs = helpers.references()
        assert len(refs) > 0
        assert any("CIRCpedia" in ref.location for ref in refs)
