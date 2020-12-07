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
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases.ensembl.genomes import protists

from . import helpers


@pytest.fixture(scope='module')  # pylint: disable=no-member
def leis_36():
    return helpers.parse(protists.parse, 'data/ensembl_protists/Leishmania_major.ASM272v2.chromosome.36.dat')


@pytest.mark.parametrize('filename,count', [
    ('data/ensembl_protists/Leishmania_major.ASM272v2.chromosome.36.dat', 71),
])
def test_can_parse_all_entries(filename, count):
    assert len(helpers.parse(protists.parse, filename)) == count


def test_can_parse_expected_accounts(leis_36):
    assert set(d.primary_id for d in leis_36) == {
        "ENSRNAT00049766034",
        "EPrT00049906992",
        "EPrT00049906958",
        "ENSRNAT00049766030",
        "EPrT00049906968",
        "ENSRNAT00049766010",
        "ENSRNAT00049766044",
        "EPrT00049906987",
        "EPrT00049906997",
        "ENSRNAT00049766004",
        "EPrT00049906979",
        "ENSRNAT00049765998",
        "ENSRNAT00049765996",
        "EPrT00049906985",
        "EPrT00049906983",
        "EPrT00049906993",
        "ENSRNAT00049766042",
        "EPrT00049906969",
        "EPrT00049906989",
        "EPrT00049906962",
        "ENSRNAT00049766018",
        "EPrT00049906991",
        "ENSRNAT00049766028",
        "EPrT00049906970",
        "ENSRNAT00049766022",
        "EPrT00049906971",
        "EPrT00049906967",
        "EPrT00049906982",
        "ENSRNAT00049766038",
        "ENSRNAT00049766020",
        "EPrT00049907001",
        "EPrT00049906974",
        "ENSRNAT00049766012",
        "EPrT00049907002",
        "EPrT00049906963",
        "EPrT00049906994",
        "EPrT00049906986",
        "ENSRNAT00049766006",
        "EPrT00049906973",
        "EPrT00049906981",
        "EPrT00049906965",
        "EPrT00049906988",
        "EPrT00049906959",
        "ENSRNAT00049766026",
        "EPrT00049906978",
        "EPrT00049906977",
        "EPrT00049906990",
        "EPrT00049906960",
        "EPrT00049906996",
        "ENSRNAT00049766008",
        "EPrT00049906964",
        "EPrT00049906998",
        "EPrT00049906966",
        "ENSRNAT00049766040",
        "EPrT00049906984",
        "EPrT00049906972",
        "EPrT00049907000",
        "ENSRNAT00049766032",
        "EPrT00049906976",
        "ENSRNAT00049766000",
        "ENSRNAT00049766036",
        "ENSRNAT00049766014",
        "EPrT00049906975",
        "ENSRNAT00049765994",
        "EPrT00049906980",
        "EPrT00049906961",
        "ENSRNAT00049766016",
        "ENSRNAT00049766024",
        "EPrT00049906995",
        "ENSRNAT00049766002",
        "EPrT00049906999",
    }


@pytest.mark.parametrize('accession,rna_type', [
    ('ENSEMBL_PROTISTS:EPrT00049906997', 'SO:0000275'),
    ('ENSEMBL_PROTISTS:ENSRNAT00049766004', 'SO:0000274'),
    ('ENSEMBL_PROTISTS:EPrT00049907002', 'SO:0000253'),
    ('ENSEMBL_PROTISTS:EPrT00049906963', 'SO:0000274'),
])
def test_can_get_expected_rna_types(leis_36, accession, rna_type):
    val = helpers.entry_for(leis_36, accession)
    assert val.rna_type == rna_type


def test_only_has_expected_ncrna_as_rna_type(leis_36):
    for entry in leis_36:
        if entry.primary_id == 'EPrT00049906996':
            continue
        assert entry.rna_type != 'SO:0000655'


def test_can_parse_expected_data(leis_36):
    assert attr.asdict(helpers.entry_for(leis_36, 'ENSEMBL_PROTISTS:EPrT00049906982')) == attr.asdict(dat.Entry(
        primary_id='EPrT00049906982',
        accession='ENSEMBL_PROTISTS:EPrT00049906982',
        ncbi_tax_id=347515,
        database='ENSEMBL_PROTISTS',
        sequence='GCGCAGTTGGTCTAGTGGTAGAATTCTCGCCTGCCACGCGGGAGGCCCGGGTTCGATTCCCGGACTGCGCA',
        regions=[
            dat.SequenceRegion(
                chromosome='36',
                strand=1,
                exons=[dat.Exon(start=1032512, stop=1032582)],
                assembly_id='ASM272v2',
                coordinate_system=dat.CoordinateSystem.one_based(),
            ),
        ],
        rna_type='SO:0000253',
        url='',
        seq_version='1',
        note_data={},
        xref_data={},
        species='Leishmania major',
        common_name=None,
        lineage=(
            'Eukaryota; Discoba; Euglenozoa; Kinetoplastea; Metakinetoplastina; '
            'Trypanosomatida; Trypanosomatidae; Leishmaniinae; Leishmania; '
            'Leishmania major strain Friedlin'
        ),
        gene='LMJF_36_TRNAGLY_01',
        locus_tag='LMJF_36_TRNAGLY_01',
        description='Leishmania major tRNA-Gly',
        references=[pubs.reference('doi:10.1093/nar/gkx1011')],
    ))
