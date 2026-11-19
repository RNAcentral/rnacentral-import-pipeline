# -*- coding: utf-8 -*-

"""
Expected values are what RNAcentral production holds for these LncBook
transcripts (rnc_accessions / rnc_sequence_regions / rnc_sequence_exons / rna),
so these tests pin the parser to the data that actually gets loaded.
"""

import pytest

from rnacentral_pipeline.databases.data import Strand
from rnacentral_pipeline.databases.lncbook import parser

SAMPLE = "data/lncbook/grch38-sample.json"


@pytest.fixture(scope="module")
def entries():
    with open(SAMPLE) as raw:
        return {e.accession: e for e in parser.parse(raw)}


def test_only_entries_with_a_grch38_location_are_kept(entries):
    # HSALNT0000002 only has a GRCh37 location
    assert sorted(entries) == ["LncBook:HSALNT0276659", "LncBook:HSALNT0279498"]


def test_grch37_locations_are_dropped(entries):
    # The fixture gives HSALNT0276659 an extra GRCh37 location as well
    regions = entries["LncBook:HSALNT0276659"].regions
    assert [r.assembly_id for r in regions] == ["GRCh38"]


def test_sequences_are_upper_cased(entries):
    seq = entries["LncBook:HSALNT0276659"].sequence
    assert seq == seq.upper()
    assert len(seq) == 815  # rna.len for URS0000E9A9DD
    assert seq.startswith("TAGGAAAGGTTCAACTTGGCAAATGTAGTG")


@pytest.mark.parametrize(
    "accession,chromosome,exons",
    [
        (
            "LncBook:HSALNT0276659",
            "22",
            [(28500011, 28500795), (28514757, 28514786)],
        ),
        ("LncBook:HSALNT0279498", "MT", [(2, 262), (1699, 1767)]),
    ],
)
def test_regions_match_production(entries, accession, chromosome, exons):
    (region,) = entries[accession].regions
    assert region.chromosome == chromosome
    assert region.strand == Strand.reverse  # rnc_sequence_regions.strand = -1
    assert [(e.start, e.stop) for e in region.exons] == exons


def test_entries_match_production_metadata(entries):
    entry = entries["LncBook:HSALNT0279498"]
    assert entry.primary_id == "HSALNT0279498"  # rnc_accessions.external_id
    assert entry.ncbi_tax_id == 9606
    assert entry.rna_type == "SO:0001904"
    assert entry.gene == "HSALNG0136336"
    assert len(entry.sequence) == 330  # rna.len for URS0001BCE0A3


def test_a_file_with_no_grch38_data_is_an_error():
    with open("data/lncbook/grch37-only.json") as raw:
        with pytest.raises(ValueError):
            list(parser.parse(raw))
