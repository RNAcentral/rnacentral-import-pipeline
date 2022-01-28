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
from rnacentral_pipeline.databases.hgnc import parser
from rnacentral_pipeline.databases.hgnc import helpers
from rnacentral_pipeline.databases.hgnc.data import Context


@pytest.fixture(scope="module")
def current_data():
    path = Path("data/hgnc/current-data.json")
    data = {}
    for entry in helpers.load(path):
        data[entry.hgnc_id] = entry
    return data


@pytest.fixture(scope="module")
def context():
    return Context.build(os.environ["PGDATABASE"])


@pytest.fixture(scope="module")
def parsed():
    data_path = Path("data/hgnc/current-data.json")
    data = {}
    for entry in parser.parse(data_path, os.environ["PGDATABASE"]):
        data[entry.primary_id] = entry
    return data


# These appear to have something wrong with them. It looks like either they are
# using the wrong Ensembl gene (though everything matches in the db) or
# annotatikons changed in a way I can't track.
ISSUES = {
    "HGNC:19380",
    "HGNC:51945",
    "HGNC:54172",
    "HGNC:54171",
    "HGNC:52553",
    "HGNC:51654",
    "HGNC:52379",
}


def known_mappings():
    params = []
    with open("data/hgnc/current-mapping.tsv", "r") as raw:
        data = list(csv.reader(raw, delimiter="\t"))
        for entry in data:
            if entry[1] in ISSUES:
                entry = pytest.param(
                    *entry, marks=pytest.mark.xfail(reason="Unknown issue")
                )
            params.append(entry)
    return params


@pytest.mark.db
@pytest.mark.parametrize("urs,hgnc_id", known_mappings())
def test_maps_sequences_correctly(current_data, context, urs, hgnc_id):
    if hgnc_id not in current_data:
        return None
    entry = current_data[hgnc_id]
    assert parser.rnacentral_id(context, entry) == urs


@pytest.mark.skip(reason="Not sure of correct count")
@pytest.mark.db
def test_can_parse_all_data(parsed):
    data_path = Path("data/hgnc/current-data.json")
    with data_path.open("r") as raw:
        expected = len(raw.readlines())
    assert len(parsed) == expected


@pytest.mark.db
def test_produces_one_entry_pre_sequence():
    data_path = Path("data/hgnc/current-data.json")
    seen = set()
    for entry in parser.parse(data_path, os.environ["PGDATABASE"]):
        assert entry.primary_id not in seen
        seen.add(entry.primary_id)


@pytest.mark.skip
@pytest.mark.db
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
