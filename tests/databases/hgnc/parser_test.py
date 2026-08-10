# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import csv
import os
from pathlib import Path

import attr
import pytest

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.hgnc import helpers, parser
from rnacentral_pipeline.databases.hgnc.data import Context, HgncEntry

DATA_PATH = Path("data/hgnc/current-data.json")

# data/hgnc/current-mapping.tsv is a regression baseline regenerated against
# the database. It previously carried seven xfails for entries that "used the
# wrong Ensembl gene" - those were the genomic sequence lookup failing to match
# spliced transcripts, and they resolve correctly now.


@pytest.fixture(scope="module")
def raw_entries():
    return helpers.load(DATA_PATH)


@pytest.fixture(scope="module")
def current_data(raw_entries):
    return {entry.hgnc_id: entry for entry in raw_entries}


@pytest.fixture(scope="module")
def context(raw_entries):
    return Context.build(os.environ["PGDATABASE"], raw_entries)


@pytest.fixture(scope="module")
def parsed(context, raw_entries):
    return {e.primary_id: e for e in parser.as_entries(context, raw_entries)}


def known_mappings():
    with open("data/hgnc/current-mapping.tsv", "r") as raw:
        return list(csv.reader(raw, delimiter="\t"))


def hgnc_entry(**kwargs) -> HgncEntry:
    defaults = {
        "symbol": "EXAMPLE1",
        "name": "example gene",
        "hgnc_id": "HGNC:1",
        "ucsc_id": None,
        "hgnc_rna_type": "RNA, long non-coding",
        "agr_id": None,
        "ensembl_gene_id": None,
        "lncipedia_id": None,
        "rnacentral_id": None,
        "previous_names": [],
        "previous_symbols": [],
        "refseq_id": None,
        "ena_ids": [],
        "gene_groups": [],
    }
    defaults.update(kwargs)
    return HgncEntry(**defaults)


@pytest.mark.hgnc
def test_so_term_handles_known_types():
    assert helpers.so_term(hgnc_entry(hgnc_rna_type="RNA, micro")) == "SO:0000276"
    assert helpers.so_term(hgnc_entry(hgnc_rna_type="RNA, transfer")) == "SO:0000253"


@pytest.mark.hgnc
def test_so_term_uses_gene_groups_for_ribosomal():
    ribosomal = {"hgnc_rna_type": "RNA, ribosomal"}
    assert helpers.so_term(hgnc_entry(**ribosomal)) == "SO:0000252"
    assert (
        helpers.so_term(hgnc_entry(gene_groups=["5S ribosomal RNAs"], **ribosomal))
        == "SO:0000652"
    )
    assert (
        helpers.so_term(hgnc_entry(gene_groups=["12S RNA"], **ribosomal))
        == "SO:0002344"
    )


@pytest.mark.hgnc
def test_so_term_skips_unknown_types_instead_of_raising():
    """A new HGNC locus type must cost one entry, not the whole run."""
    assert helpers.so_term(hgnc_entry(hgnc_rna_type="RNA, newly invented")) is None


@pytest.mark.hgnc
def test_longest_sequence_wins():
    rows = [
        ("key", "URS0000000001", 100),
        ("key", "URS0000000002", 300),
        ("key", "URS0000000003", 200),
        ("other", "URS0000000004", 10),
    ]
    assert helpers._longest(rows) == {
        "key": "URS0000000002",
        "other": "URS0000000004",
    }


@pytest.mark.db
@pytest.mark.hgnc
@pytest.mark.parametrize("urs,hgnc_id", known_mappings())
def test_maps_sequences_correctly(current_data, context, urs, hgnc_id):
    if hgnc_id not in current_data:
        return
    entry = current_data[hgnc_id]
    assert parser.rnacentral_id(context, entry) == urs


@pytest.mark.db
@pytest.mark.hgnc
def test_maps_the_expected_number_of_entries(parsed):
    """
    Guards against silent degradation. A failing Ensembl batch, or a lost
    partition filter, drops entries without raising anything.
    """
    assert len(parsed) > 7900


@pytest.mark.db
@pytest.mark.hgnc
def test_produces_one_entry_per_gene(context, raw_entries):
    ids = [e.primary_id for e in parser.as_entries(context, raw_entries)]
    assert len(ids) == len(set(ids))


@pytest.mark.db
@pytest.mark.hgnc
def test_produces_expected_data(parsed):
    assert attr.asdict(parsed["HGNC:34365"]) == attr.asdict(
        Entry(
            primary_id="HGNC:34365",
            accession="HGNC:34365",
            ncbi_tax_id=9606,
            database="HGNC",
            sequence="GTCTACGGCCATACCACCCTGAACGCGCCCGATCTCGTCTGATCTCGGAAGCTAAGCAGGGTCGGGCCTGGTTAGTACTTGGATGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTTT",
            regions=[],
            rna_type="SO:0000652",
            url="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:34365",
            seq_version="1",
            description="Homo sapiens (human) RNA, 5S ribosomal 4",
            locus_tag="RNA5S4",
            gene="RNA5S4",
            gene_synonyms=["RN5S4"],
        )
    )
