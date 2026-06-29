# -*- coding: utf-8 -*-

import io
from pathlib import Path

import pytest

from rnacentral_pipeline.databases.data import RelatedCoordinate, RelatedSequence
from rnacentral_pipeline.databases.mirtrondb.parser import find_coords, parse

DATA = Path(__file__).parent / "data.tsv"


@pytest.fixture
def entries():
    with open(DATA) as f:
        result = list(parse(f))
    return {e.accession: e for e in result}


class TestParse:
    def test_yields_correct_count(self, entries):
        assert len(entries) == 4

    def test_precursor_entry_fields(self, entries):
        e = entries["hsa-mirtron-1230"]
        assert e.primary_id == "MIRTRONDB:218"
        assert e.ncbi_tax_id == 9606
        assert e.database == "MIRTRONDB"
        assert e.rna_type == "SO:0001244"
        assert e.seq_version == "1"
        assert e.gene == "ANKRD24(14)"
        assert e.sequence == "GACACACGGGGGCTGAGAGCAGAACCCCGAGCTGGACTCTGCTCCCTCCCCCCAG"
        assert "precursor" in e.description
        assert "ANKRD24(14)" in e.description
        assert (
            e.url
            == "http://mirtrondb.cp.utfpr.edu.br/fetch_details.php?mrt_details=hsa-mirtron-1230"
        )

    def test_mature_entry_fields(self, entries):
        e = entries["hsa-mirtron-1230-3p"]
        assert e.primary_id == "MIRTRONDB:1501"
        assert e.ncbi_tax_id == 9606
        assert e.rna_type == "SO:0000276"
        assert e.sequence == "ACTCTGCTCCCTCCCCCCAG"

    def test_newly_added_species(self, entries):
        e = entries["rno-mirtron-1224"]
        assert e.ncbi_tax_id == 10116
        assert e.database == "MIRTRONDB"

    def test_empty_host_gene_is_none(self, entries):
        e = entries["rno-mirtron-1224"]
        assert e.gene is None
        assert "(" not in e.description

    def test_precursor_has_mature_relationships(self, entries):
        pre = entries["hsa-mirtron-1230"]
        assert len(pre.related_sequences) == 2
        rel_ids = {r.sequence_id for r in pre.related_sequences}
        assert rel_ids == {"hsa-mirtron-1230-3p", "hsa-mirtron-1230-5p"}
        for rel in pre.related_sequences:
            assert rel.relationship == "mature_product"

    def test_mature_has_precursor_relationship(self, entries):
        mat = entries["hsa-mirtron-1230-3p"]
        assert len(mat.related_sequences) == 1
        assert mat.related_sequences[0].sequence_id == "hsa-mirtron-1230"
        assert mat.related_sequences[0].relationship == "precursor"

    def test_mature_coordinates_in_precursor(self, entries):
        pre = entries["hsa-mirtron-1230"]
        rel_3p = [
            r for r in pre.related_sequences if r.sequence_id == "hsa-mirtron-1230-3p"
        ][0]
        assert len(rel_3p.coordinates) == 1
        coord = rel_3p.coordinates[0]
        assert pre.sequence[coord.start : coord.stop] == "ACTCTGCTCCCTCCCCCCAG"

    def test_precursor_without_matures_has_no_relationships(self, entries):
        pre = entries["rno-mirtron-1224"]
        assert len(pre.related_sequences) == 0


class TestFindCoords:
    def test_found(self):
        coords = find_coords("test", "ABCDEFGH", "CDE")
        assert len(coords) == 1
        assert coords[0] == RelatedCoordinate(start=2, stop=5)

    def test_not_found(self):
        coords = find_coords("test", "ABCDEFGH", "XYZ")
        assert coords == []


class TestParseErrors:
    def test_invalid_header(self):
        data = "not a valid header\n"
        with pytest.raises(AssertionError, match="Invalid header line"):
            list(parse(io.StringIO(data)))

    def test_unknown_rna_type(self):
        data = (
            "##mirtronDB tabular format\n"
            "specie\tid\tname\tpaper\ttype\thairpin arm\thost gene\tchromosome\tstart\tend\tstrand\tsequence\n"
            "H. sapiens \t1 \ttest-1 \t123 \tunknown \t- \tGENE1 \t1 \t100 \t200 \t+ \tACGUACGUACGUACGUACGU\n"
        )
        with pytest.raises(KeyError):
            list(parse(io.StringIO(data)))

    def test_unknown_species(self):
        data = (
            "##mirtronDB tabular format\n"
            "specie\tid\tname\tpaper\ttype\thairpin arm\thost gene\tchromosome\tstart\tend\tstrand\tsequence\n"
            "Z. unknown \t1 \ttest-1 \t123 \tmature \t3p \tGENE1 \t1 \t100 \t200 \t+ \tACGUACGUACGUACGUACGU\n"
        )
        with pytest.raises(KeyError):
            list(parse(io.StringIO(data)))

    def test_duplicate_accession(self):
        data = (
            "##mirtronDB tabular format\n"
            "specie\tid\tname\tpaper\ttype\thairpin arm\thost gene\tchromosome\tstart\tend\tstrand\tsequence\n"
            "H. sapiens \t1 \tdup-name \t123 \tmature \t3p \tGENE1 \t1 \t100 \t200 \t+ \tACGUACGUACGUACGUACGU\n"
            "H. sapiens \t2 \tdup-name \t123 \tmature \t5p \tGENE1 \t1 \t100 \t200 \t+ \tACGUACGUACGUACGUACGU\n"
        )
        with pytest.raises(ValueError, match="Duplicate accession"):
            list(parse(io.StringIO(data)))
