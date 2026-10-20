# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import attr
import pytest

from rnacentral_pipeline.databases import data as dat
from rnacentral_pipeline.databases.ensembl import protists
from rnacentral_pipeline.databases.helpers import publications as pubs

from . import helpers

# Trimmed to 6 representative ncRNA genes (of the original chromosome's 71)
# covering distinct rna types, rather than shipping the full ~5MB chromosome.

FIXTURE = "test-data/ensembl_protists/Leishmania_major.ASM272v2.chromosome.36.dat"


@pytest.fixture(scope="module")  # pylint: disable=no-member
def leis_36():
    return helpers.parse(protists.parse, FIXTURE)


@pytest.mark.parametrize(
    "filename,count",
    [
        (FIXTURE, 6),
    ],
)
def test_can_parse_all_entries(filename, count):
    assert len(helpers.parse(protists.parse, filename)) == count


def test_can_parse_expected_accounts(leis_36):
    assert set(d.primary_id for d in leis_36) == {
        "EPrT00049906997",
        "ENSRNAT00049766004",
        "EPrT00049906982",
        "EPrT00049907002",
        "EPrT00049906963",
        "EPrT00049906996",
    }


@pytest.mark.parametrize(
    "accession,rna_type",
    [
        ("ENSEMBL_PROTISTS:EPrT00049906997", "SO:0000275"),
        ("ENSEMBL_PROTISTS:ENSRNAT00049766004", "SO:0000673"),
        ("ENSEMBL_PROTISTS:EPrT00049907002", "SO:0000673"),
        ("ENSEMBL_PROTISTS:EPrT00049906963", "SO:0000274"),
    ],
)
def test_can_get_expected_rna_types(leis_36, accession, rna_type):
    val = helpers.entry_for(leis_36, accession)
    assert val.rna_type == rna_type


def test_only_has_expected_ncrna_as_rna_type(leis_36):
    for entry in leis_36:
        assert entry.rna_type != "SO:0000655"


def test_can_parse_expected_data(leis_36):
    assert attr.asdict(
        helpers.entry_for(leis_36, "ENSEMBL_PROTISTS:EPrT00049906982")
    ) == attr.asdict(
        dat.Entry(
            primary_id="EPrT00049906982",
            accession="ENSEMBL_PROTISTS:EPrT00049906982",
            ncbi_tax_id=347515,
            database="ENSEMBL_PROTISTS",
            sequence="GCGCAGTTGGTCTAGTGGTAGAATTCTCGCCTGCCACGCGGGAGGCCCGGGTTCGATTCCCGGACTGCGCA",
            regions=[
                dat.SequenceRegion(
                    chromosome="36",
                    strand=1,
                    exons=[dat.Exon(start=164, stop=234)],
                    assembly_id="ASM272v2",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                ),
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            note_data={},
            xref_data={},
            species="Leishmania major",
            common_name=None,
            lineage=(
                "Eukaryota; Discoba; Euglenozoa; Kinetoplastea; Metakinetoplastina; "
                "Trypanosomatida; Trypanosomatidae; Leishmaniinae; Leishmania; "
                "Leishmania major strain Friedlin"
            ),
            gene="LMJF_36_TRNAGLY_01",
            description="Leishmania major tRNA-Gly",
            references=[pubs.reference("doi:10.1093/nar/gkx1011")],
        )
    )
