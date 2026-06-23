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
def sample_annotation_file():
    """Create a temporary TSV annotation file with sample circRNA data."""
    tsv_content = """circID\tspecies\tLocation\tgene_Refseq\tgene_Ensembl\tcircname\tlength\tsubcell_location
hsa_circ_0001\tHomo sapiens\tchr1:1000-2000(+)\tTP53\tENSG00000141510\tcirc-TP53\t1000\tcytoplasm
hsa_circ_0002\tHomo sapiens\tchr2:3000-4000(-)\tMYC\tENSG00000136997\tcirc-MYC\t1000\tnucleus
mmu_circ_0001\tMus musculus\tchr3:5000-6000(+)\tTrp53\tENSMUSG00000059552\tcirc-Trp53\t1000\tcytoplasm
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(tsv_content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


@pytest.fixture
def sample_fasta_file():
    """Create a temporary FASTA file with sample sequences."""
    fasta_content = """>hsa_circ_0001
ATCGATCGATCGATCG
>hsa_circ_0002
GCTAGCTAGCTAGCTA
>mmu_circ_0001
TTTTAAAACCCCGGGG
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        f.write(fasta_content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


def _mock_seq_record(sequence):
    """Create a mock SeqRecord with a .seq attribute."""
    record = MagicMock()
    record.seq = sequence
    return record


class TestParseTSVRow:
    """Test parsing individual TSV rows."""

    def test_parses_valid_row_with_all_fields(self):
        """Test parsing a complete valid row."""
        sequences = {"hsa_circ_0001": _mock_seq_record("ATCGATCG")}

        with patch(
            "rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid"
        ) as mock_taxid, patch(
            "rnacentral_pipeline.databases.helpers.phylogeny.species"
        ) as mock_species, patch(
            "rnacentral_pipeline.databases.helpers.phylogeny.common_name"
        ) as mock_common, patch(
            "rnacentral_pipeline.databases.helpers.phylogeny.lineage"
        ) as mock_lineage:

            mock_taxid.return_value = 9606
            mock_species.return_value = "Homo sapiens"
            mock_common.return_value = "human"
            mock_lineage.return_value = "Eukaryota; Metazoa; Chordata"

            row = {
                "circID": "hsa_circ_0001",
                "species": "Homo sapiens",
                "Location": "chr1:1000-2000(+)",
                "gene_Refseq": "TP53",
                "circname": "circ-TP53",
                "length": "8",
            }

            entry = parser.parse_tsv_row(row, sequences, assembly_id="GRCh38")

            assert entry is not None
            assert entry.primary_id == "hsa_circ_0001"
            assert entry.ncbi_tax_id == 9606
            assert entry.database == "CIRCPEDIA"
            assert entry.sequence == "ATCGATCG"
            assert entry.rna_type == "SO:0002291"  # circular RNA
            assert entry.gene == "TP53"
            assert entry.chromosome == "1"
            assert len(entry.regions) > 0

    def test_returns_none_for_missing_circid(self):
        """Test that missing circID returns None."""
        sequences = {}
        row = {
            "species": "Homo sapiens",
            "Location": "chr1:1000-2000(+)",
        }

        entry = parser.parse_tsv_row(row, sequences)
        assert entry is None

    def test_returns_none_for_missing_species(self):
        """Test that missing species returns None."""
        sequences = {"hsa_circ_0001": _mock_seq_record("ATCG")}
        row = {
            "circID": "hsa_circ_0001",
            "Location": "chr1:1000-2000(+)",
        }

        entry = parser.parse_tsv_row(row, sequences)
        assert entry is None

    def test_returns_none_for_missing_location(self):
        """Test that missing location returns None."""
        sequences = {"hsa_circ_0001": _mock_seq_record("ATCG")}

        with patch(
            "rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid"
        ) as mock_taxid:
            mock_taxid.return_value = 9606

            row = {
                "circID": "hsa_circ_0001",
                "species": "Homo sapiens",
            }

            entry = parser.parse_tsv_row(row, sequences)
            assert entry is None

    def test_returns_none_for_missing_sequence(self):
        """Test that missing sequence returns None."""
        sequences = {}  # No sequence for this circID

        with patch(
            "rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid"
        ) as mock_taxid:
            mock_taxid.return_value = 9606

            row = {
                "circID": "hsa_circ_0001",
                "species": "Homo sapiens",
                "Location": "chr1:1000-2000(+)",
            }

            entry = parser.parse_tsv_row(row, sequences)
            assert entry is None

    def test_returns_none_for_unknown_taxid(self):
        """Test that unknown taxid returns None."""
        sequences = {"hsa_circ_0001": _mock_seq_record("ATCG")}

        with patch(
            "rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid"
        ) as mock_taxid:
            mock_taxid.return_value = None

            row = {
                "circID": "hsa_circ_0001",
                "species": "Unknown species",
                "Location": "chr1:1000-2000(+)",
            }

            entry = parser.parse_tsv_row(row, sequences)
            assert entry is None


class TestParse:
    """Test the main parse function."""

    def test_parses_tsv_and_fasta_files(
        self, sample_annotation_file, sample_fasta_file
    ):
        """Test parsing TSV annotation and FASTA sequence files."""
        with patch(
            "rnacentral_pipeline.databases.circpedia.helpers.get_ncbi_taxid"
        ) as mock_taxid, patch(
            "rnacentral_pipeline.databases.helpers.phylogeny.species"
        ) as mock_species, patch(
            "rnacentral_pipeline.databases.helpers.phylogeny.common_name"
        ) as mock_common, patch(
            "rnacentral_pipeline.databases.helpers.phylogeny.lineage"
        ) as mock_lineage:

            # Mock taxonomy lookups
            def get_taxid(species):
                return {"Homo sapiens": 9606, "Mus musculus": 10090}.get(species)

            mock_taxid.side_effect = get_taxid
            mock_species.return_value = "Species name"
            mock_common.return_value = "common"
            mock_lineage.return_value = "lineage"

            # Parse the files
            entries = list(
                parser.parse(
                    sample_annotation_file,
                    sample_fasta_file,
                    assembly_id="GRCh38",
                )
            )

            # Verify results
            assert len(entries) == 3
            assert all(e.database == "CIRCPEDIA" for e in entries)
            assert all(e.rna_type == "SO:0002291" for e in entries)

    def test_handles_missing_required_columns(self, sample_fasta_file):
        """Test that missing required columns raises error."""
        # Create TSV without required columns
        tsv_content = """id\tname
1\ttest
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(tsv_content)
            temp_path = f.name

        try:
            with pytest.raises(ValueError, match="Missing required column"):
                list(parser.parse(temp_path, sample_fasta_file))
        finally:
            Path(temp_path).unlink(missing_ok=True)
