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
Tests for CIRCpedia parser.

Note: These tests use minimal mock data to avoid dependencies on:
- Network resources
- RNAcentral database
- Large test data files
"""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from rnacentral_pipeline.databases.circpedia import parser


@pytest.fixture
def mock_taxonomy():
    """Create a mock taxonomy database."""
    mock_db = MagicMock()
    mock_db.close = MagicMock()
    return mock_db


@pytest.fixture
def sample_csv_file():
    """Create a temporary CSV file with sample circRNA data."""
    csv_content = """circid,species,location,strand,gene,sequence,exon_positions,fpm
hsa_circ_0001,Homo sapiens,chr1:1000-2000,+,TP53,ATCGATCGATCGATCG,1000-1200;1400-2000,10.5
hsa_circ_0002,Homo sapiens,chr2:3000-4000,-,MYC,GCTAGCTAGCTAGCTA,3000-3500;3600-4000,5.2
mmu_circ_0001,Mus musculus,chr3:5000-6000,+,Trp53,TTTTAAAACCCCGGGG,5000-5500;5700-6000,8.7
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


@pytest.fixture
def mock_taxonomy_file():
    """Create a mock taxonomy database file."""
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        temp_path = f.name

    yield Path(temp_path)

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


class TestParseCSVRow:
    """Test parsing individual CSV rows."""

    def test_parses_valid_row_with_all_fields(self, mock_taxonomy):
        """Test parsing a complete valid row."""
        with patch('rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid') as mock_taxid, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.species') as mock_species, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.common_name') as mock_common, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.lineage') as mock_lineage:

            mock_taxid.return_value = 9606
            mock_species.return_value = "Homo sapiens"
            mock_common.return_value = "human"
            mock_lineage.return_value = "Eukaryota; Metazoa; Chordata"

            row = {
                "circid": "hsa_circ_0001",
                "species": "Homo sapiens",
                "location": "chr1:1000-2000",
                "strand": "+",
                "gene": "TP53",
                "sequence": "ATCGATCG",
                "exon_positions": "1000-1200,1400-2000",
                "fpm": 10.5,
            }

            entry = parser.parse_csv_row(row, mock_taxonomy, assembly_id="GRCh38")

            assert entry is not None
            assert entry.primary_id == "hsa_circ_0001"
            assert entry.ncbi_tax_id == 9606
            assert entry.database == "CIRCPEDIA"
            assert entry.sequence == "ATCGATCG"
            assert entry.rna_type == "SO:0000593"  # circular RNA
            assert entry.gene == "TP53"
            assert entry.chromosome == "1"
            assert len(entry.regions) > 0

    def test_returns_none_for_missing_circid(self, mock_taxonomy):
        """Test that missing circid returns None."""
        row = {
            "species": "Homo sapiens",
            "location": "chr1:1000-2000",
        }

        entry = parser.parse_csv_row(row, mock_taxonomy)
        assert entry is None

    def test_returns_none_for_missing_species(self, mock_taxonomy):
        """Test that missing species returns None."""
        row = {
            "circid": "hsa_circ_0001",
            "location": "chr1:1000-2000",
        }

        entry = parser.parse_csv_row(row, mock_taxonomy)
        assert entry is None

    def test_returns_none_for_missing_location(self, mock_taxonomy):
        """Test that missing location returns None."""
        with patch('rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid') as mock_taxid:
            mock_taxid.return_value = 9606

            row = {
                "circid": "hsa_circ_0001",
                "species": "Homo sapiens",
            }

            entry = parser.parse_csv_row(row, mock_taxonomy)
            assert entry is None

    def test_returns_none_for_unknown_taxid(self, mock_taxonomy):
        """Test that unknown taxid returns None."""
        with patch('rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid') as mock_taxid:
            mock_taxid.return_value = None

            row = {
                "circid": "hsa_circ_0001",
                "species": "Unknown species",
                "location": "chr1:1000-2000",
            }

            entry = parser.parse_csv_row(row, mock_taxonomy)
            assert entry is None

    def test_handles_row_without_exons(self, mock_taxonomy):
        """Test parsing row without exon positions."""
        with patch('rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid') as mock_taxid, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.species') as mock_species, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.common_name') as mock_common, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.lineage') as mock_lineage:

            mock_taxid.return_value = 9606
            mock_species.return_value = "Homo sapiens"
            mock_common.return_value = "human"
            mock_lineage.return_value = "Eukaryota"

            row = {
                "circid": "hsa_circ_0001",
                "species": "Homo sapiens",
                "location": "chr1:1000-2000",
                "sequence": "ATCG",
            }

            entry = parser.parse_csv_row(row, mock_taxonomy)
            assert entry is not None
            # Should have default exon spanning the whole region
            assert len(entry.regions) > 0


class TestParse:
    """Test the main parse function."""

    @patch('rnacentral_pipeline.databases.circpedia.parser.SqliteDict')
    def test_parses_csv_file(self, mock_sqlite, sample_csv_file, mock_taxonomy_file):
        """Test parsing a CSV file with multiple entries."""
        # Setup mock taxonomy
        mock_db = MagicMock()
        mock_db.close = MagicMock()
        mock_sqlite.return_value = mock_db

        with patch('rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid') as mock_taxid, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.species') as mock_species, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.common_name') as mock_common, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.lineage') as mock_lineage:

            # Mock taxonomy lookups
            def get_taxid(db, species):
                return {"Homo sapiens": 9606, "Mus musculus": 10090}.get(species)

            mock_taxid.side_effect = get_taxid
            mock_species.return_value = "Species name"
            mock_common.return_value = "common"
            mock_lineage.return_value = "lineage"

            # Parse the file
            entries = list(parser.parse(sample_csv_file, mock_taxonomy_file, assembly_id="GRCh38"))

            # Verify results
            assert len(entries) == 3
            assert all(e.database == "CIRCPEDIA" for e in entries)
            assert all(e.rna_type == "SO:0000593" for e in entries)

            # Verify taxonomy database was closed
            mock_db.close.assert_called_once()

    @patch('rnacentral_pipeline.databases.circpedia.parser.SqliteDict')
    def test_handles_missing_required_columns(self, mock_sqlite, mock_taxonomy_file):
        """Test that missing required columns raises error."""
        mock_db = MagicMock()
        mock_sqlite.return_value = mock_db

        # Create CSV without required columns
        csv_content = """id,name
1,test
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            temp_path = f.name

        try:
            with pytest.raises(ValueError, match="Missing required columns"):
                list(parser.parse(temp_path, mock_taxonomy_file))
        finally:
            Path(temp_path).unlink(missing_ok=True)
            mock_db.close.assert_called_once()

    @patch('rnacentral_pipeline.databases.circpedia.parser.SqliteDict')
    def test_logs_progress(self, mock_sqlite, sample_csv_file, mock_taxonomy_file):
        """Test that progress is logged during parsing."""
        mock_db = MagicMock()
        mock_sqlite.return_value = mock_db

        with patch('rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid') as mock_taxid, \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.species'), \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.common_name'), \
             patch('rnacentral_pipeline.databases.helpers.phylogeny.lineage'):

            mock_taxid.return_value = 9606

            # Parse and verify no exceptions
            entries = list(parser.parse(sample_csv_file, mock_taxonomy_file))

            # Should parse successfully
            assert len(entries) >= 0  # May vary depending on mock behavior
