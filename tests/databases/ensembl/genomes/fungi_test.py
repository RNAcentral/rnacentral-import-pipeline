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
from rnacentral_pipeline.databases.ensembl.genomes import fungi

from . import helpers


@pytest.fixture(scope='module')  # pylint: disable=no-member
def asp_1():
    return helpers.parse(fungi.parse, 'data/ensembl_fungi/Aspergillus_nidulans.chromosome.I.dat')


@pytest.mark.parametrize('filename,count', [
    ('data/ensembl_fungi/Aspergillus_nidulans.chromosome.I.dat', 24),
])
def test_can_parse_all_entries(filename, count):
    assert len(helpers.parse(fungi.parse, filename)) == count


def test_can_parse_expected_ids(asp_1):
    assert set(d.primary_id for d in asp_1) == {
        'EBT00005235647',
        'EBT00005235679',
        'EBT00005235661',
        'EBT00005235639',
        'EBT00005235669',
        'EBT00005235667',
        'EBT00005235653',
        'EBT00005235637',
        'EBT00005235663',
        'EBT00005235683',
        'EBT00005235643',
        'EBT00005235641',
        'EBT00005235671',
        'EBT00005235649',
        'EBT00005235681',
        'EBT00005235645',
        'EBT00005235673',
        'EBT00005235677',
        'EBT00005235665',
        'EBT00005235659',
        'EBT00005235651',
        'EBT00005235655',
        'EBT00005235675',
        'EBT00005235657',
    }


def test_can_produce_correct_data(asp_1):
    assert attr.asdict(helpers.entry_for(asp_1, "ENSEMBL_FUNGI:EBT00005235663")) == attr.asdict(dat.Entry(
        primary_id="EBT00005235663",
        accession='ENSEMBL_FUNGI:EBT00005235663',
        ncbi_tax_id=227321,
        database='ENSEMBL_FUNGI',
        sequence='AAGGGTTCGCGTCGATTGGGGGTAGTGTAAGCCCAATGACACTATGTGGCAAGGCATGCGCAAGCACAGAGGCCCTGCCACGAGATCCAAAAATAAGATTAATCAATGTTGTGTTCCTTAGAACACTTTCAGAGATTAAACTGAAAAGAACAGATACTACACTTGATCTAAGCCAAAAGGCAAAGAGAGCTGG',
        regions=[
            dat.SequenceRegion(
                chromosome='I',
                strand=-1,
                exons=[dat.Exon(start=1978646, stop=1978838)],
                assembly_id='ASM1142v1',
                coordinate_system=dat.CoordinateSystem.one_based(),
            ),
        ],
        rna_type='snRNA',
        url='',
        seq_version='1',
        note_data={},
        xref_data={'RFAM': ['RF00004']},
        species='Aspergillus nidulans FGSC A4',
        common_name=None,
        lineage=(
            'Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; '
            'Eurotiomycetes; Eurotiomycetidae; Eurotiales; Aspergillaceae; '
            'Aspergillus; Aspergillus subgen. Nidulantes; '
            'Aspergillus nidulans FGSC A4'
        ),
        gene="EBG00005235662",
        locus_tag='U2',
        description='Aspergillus nidulans FGSC A4 snRNA U2',
        references=[pubs.reference('doi:10.1093/nar/gkx1011')],
    ))


def test_does_not_create_ncRNA_rna_type(asp_1):
    for entry in asp_1:
        assert entry.rna_type != 'ncRNA'
