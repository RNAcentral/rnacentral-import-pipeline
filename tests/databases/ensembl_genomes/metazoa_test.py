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
from rnacentral_pipeline.databases.ensembl_genomes import metazoa

from . import helpers


@pytest.fixture(scope='module')  # pylint: disable=no-member
def apis_1():
    return helpers.parse(metazoa.parse, 'data/ensembl_metazoa/Apis_mellifera.Amel_4.5.chromosome_group.1.dat')


@pytest.mark.parametrize('filename,count', [
    ('data/ensembl_metazoa/Apis_mellifera.Amel_4.5.chromosome_group.1.dat', 36),
])
def test_can_parse_all_entries(filename, count):
    assert len(helpers.parse(metazoa.parse, filename)) == count


def test_gets_all_ids(apis_1):
    assert set(d.primary_id for d in apis_1) == {
        'ENSRNA022717521-T1',
        "ENSRNA022717522-T1",
        "ENSRNA022712925-T1",
        "ENSRNA022712939-T1",
        "ENSRNA022717520-T1",
        "ENSRNA022717518-T1",
        "ENSRNA022717519-T1",
        "ENSRNA022717517-T1",
        "ENSRNA022717515-T1",
        "ENSRNA022717516-T1",
        "ENSRNA022717513-T1",
        "ENSRNA022717514-T1",
        "ENSRNA022712978-T1",
        "ENSRNA022712985-T1",
        "ENSRNA022717880-T1",
        "ENSRNA022717863-T1",
        "ENSRNA022717861-T1",
        "ENSRNA022717512-T1",
        "ENSRNA022717510-T1",
        "ENSRNA022717511-T1",
        "ENSRNA022717829-T1",
        "ENSRNA022717825-T1",
        "ENSRNA022717509-T1",
        "ENSRNA022712955-T1",
        "ENSRNA022717507-T1",
        "ENSRNA022717508-T1",
        "ENSRNA022717506-T1",
        "ENSRNA022717504-T1",
        "ENSRNA022717505-T1",
        "ENSRNA022712934-T1",
        "ENSRNA022712958-T1",
        "ENSRNA022717503-T1",
        "ENSRNA022717501-T1",
        "ENSRNA022717502-T1",
        "ENSRNA022717500-T1",
        "ENSRNA022712948-T1",
        "ENSRNA022717521-T1",
    }


@pytest.mark.parametrize('accession,rna_type', [
    ('ENSEMBL_METAZOA:ENSRNA022717521-T1', 'precursor_RNA'),  # TODO: pre_miRNA
    ("ENSEMBL_METAZOA:ENSRNA022712925-T1", 'tRNA'),
    ("ENSEMBL_METAZOA:ENSRNA022717829-T1", 'snoRNA'),
])
def test_can_assign_expected_rna_types(apis_1, accession, rna_type):
    assert helpers.entry_for(apis_1, accession).rna_type == rna_type


@pytest.mark.parametrize('accession,description', [
    ('ENSEMBL_METAZOA:ENSRNA022717825-T1', 'Apis mellifera (European honey bee) U4 spliceosomal RNA'),
])
def test_can_get_expected_descriptions(apis_1, accession, description):
    assert helpers.entry_for(apis_1, accession).description == description


def test_can_get_expected_data(apis_1):
    assert attr.asdict(helpers.entry_for(apis_1, 'ENSEMBL_METAZOA:ENSRNA022717521-T1')) == attr.asdict(dat.Entry(
        primary_id='ENSRNA022717521-T1',
        accession='ENSEMBL_METAZOA:ENSRNA022717521-T1',
        ncbi_tax_id=7460,
        database='ENSEMBL_METAZOA',
        sequence='AGGGTCGAAGAGTGAGTAAATGGCCGAGGGTGATTTGGGCCTTAGTGGTCCTGTGGTGGCTGCGTACGAATCCTACTGGCCTGCTAAGTCCCAAGTGATTCTCGGCTCGCGCTGCGATA',
        regions=[
            dat.SequenceRegion(
                chromosome='1',
                strand=1,
                exons=[dat.Exon(start=204902, stop=205020)],
                assembly_id='Amel_4.5',
                coordinate_system=dat.CoordinateSystem.one_based(),
            ),
        ],
        rna_type='precursor_RNA',
        url='',
        seq_version='1',
        note_data={},
        xref_data={},
        species='Apis mellifera',
        common_name='European honey bee',
        lineage=(
            'Eukaryota; Metazoa; Ecdysozoa; Arthropoda; '
            'Hexapoda; Insecta; Pterygota; Neoptera; Holometabola; '
            'Hymenoptera; Apocrita; Aculeata; Apoidea; Apidae; Apis; '
            'Apis mellifera'
        ),
        gene='ENSRNA022717521',
        locus_tag='ame-mir-193',
        description='Apis mellifera (European honey bee) precursor RNA ame-mir-193',
        references=[pubs.reference('doi:10.1093/nar/gkx1011')],
    ))
