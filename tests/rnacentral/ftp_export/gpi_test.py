# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

import os

import psycopg2
import pytest
from pypika import Table

from rnacentral_pipeline.rnacentral.ftp_export import gpi


@pytest.fixture(scope="module")
def conn():
    db_url = os.getenv("PGDATABASE")
    with psycopg2.connect(db_url) as conn:
        yield conn


def entry(conn, urs_taxid: str) -> gpi.GpiEntry:
    pre = Table("rnc_rna_precomputed")
    related = Table("rnc_related_sequences")
    generic_query = gpi.generic_query().where(pre.id == urs_taxid)
    mirbase_query = gpi.mirbase_info_query().where(
        related.source_urs_taxid == urs_taxid
    )
    entries = gpi.load(conn, generic_query=generic_query, mirbase_query=mirbase_query)
    entries = list(entries)
    assert len(entries) == 1
    return entries[0]


@pytest.mark.db
@pytest.mark.parametrize(
    "urs_taxid,expected",
    [
        (
            "URS000012F9EC_9606",
            gpi.GpiEntry(
                urs_taxid="URS000012F9EC_9606",
                description="Homo sapiens (human) hsa-miR-4691-3p",
                rna_type="miRNA",
                symbol="hsa-miR-4691-3p",
                precursors={"URS000075C981_9606"},
            ),
        ),
        (
            "URS00001AE6A4_10090",
            gpi.GpiEntry(
                urs_taxid="URS00001AE6A4_10090",
                description="Mus musculus (house mouse) mmu-miR-700-5p",
                rna_type="miRNA",
                symbol="mmu-miR-700-5p",
                precursors={"URS000075A134_10090"},
            ),
        ),
    ],
)
def test_builds_correct_data(conn, urs_taxid, expected):
    assert entry(conn, urs_taxid) == expected
