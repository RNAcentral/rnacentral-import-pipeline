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

from rnacentral_pipeline.databases.data import regions
from rnacentral_pipeline.rnacentral.ftp_export.coordinates import bed

from .helpers import fetch_coord


def fetch_data(rna_id, assembly):
    data = [bed.BedEntry.from_coordinate(c) for c in fetch_coord(rna_id, assembly)]
    assert len(data) == 1
    return data[0]


def test_can_build_bed_from_region():
    data = fetch_data('URS000082BE64_9606', "GRCh38")
    assert attr.asdict(data) == attr.asdict(bed.BedEntry(
        rna_id='URS000082BE64_9606',
        rna_type='snoRNA',
        databases='snOPY',
        region=region.SequenceRegion(
            assembly_id='GRCh38',
            chromosome='3',
            strand=1,
            exons=[regions.Exon(start=32676136, stop=32676226)],
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
    ))


def test_can_build_entry_with_several_databases():
    data = fetch_data('URS0000368518_9606', "GRCh38")
    assert attr.asdict(data) == attr.asdict(bed.BedEntry(
        rna_id='URS0000368518_9606',
        rna_type='lncRNA',
        databases='Ensembl,GENCODE',
        region=region.SequenceRegion(
            assembly_id='GRCh38',
            chromosome='2',
            strand=1,
            exons=[regions.Exon(start=37208875, stop=37212677)],
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
    ))


def test_can_build_entry_from_pig():
    data = fetch_data('URS000099C6E5_9598', 'Pan_tro_3.0')
    assert attr.asdict(data) == attr.asdict(bed.BedEntry(
        rna_id='URS000099C6E5_9598',
        rna_type='SRP_RNA',
        databases='Rfam',
        region=region.SequenceRegion(
            assembly_id='Pan_tro_3.0',
            chromosome='X',
            strand=-1,
            exons=[regions.Exon(start=73819430, stop=73819577)],
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
    ))


def test_can_build_entry_with_several_exons():
    data = fetch_data('URS0000000055_9606', 'GRCh38')
    assert attr.asdict(data) == attr.asdict(bed.BedEntry(
        rna_id='URS0000000055_9606',
        rna_type='lncRNA',
        databases='Ensembl,GENCODE,LNCipedia,NONCODE',
        region=region.SequenceRegion(
            assembly_id='GRCh38',
            chromosome='6',
            strand=-1,
            exons=[
                regions.Exon(start=57171005, stop=57171098),
                regions.Exon(start=57173376, stop=57173480),
                regions.Exon(start=57173737, stop=57174236),
            ],
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
    ))


@pytest.mark.parametrize('upi,assembly,expected', [
    ('URS0000368518_9606', "GRCh38", 'chr2'),
    ('URS000001B2EC_9606', "GRCh38", 'chrM'),
    ('URS000071014B_112509', 'IBSC_v2', 'chr3H'),
])
def test_gets_correct_chromosome(upi, assembly, expected):
    assert fetch_data(upi, assembly).bed_chromosome == expected


@pytest.mark.parametrize('upi,assembly,expected', [
    ('URS0000368518_9606', "GRCh38", [1]),
    ('URS000001B2EC_9606', "GRCh38", [1]),
    ('URS000071014B_112509', 'IBSC_v2', [1]),
    ('URS0000000055_9606', 'GRCh38', [1, 1, 1]),
])
def test_gets_correct_bed_sizes(upi, assembly, expected):
    assert fetch_data(upi, assembly).sizes() == expected


@pytest.mark.parametrize('upi,assembly,expected', [
    ('URS0000368518_9606', "GRCh38", [0]),
    ('URS000001B2EC_9606', "GRCh38", [0]),
    ('URS000071014B_112509', 'IBSC_v2', [0]),
    ('URS0000000055_9606', 'GRCh38', [1, 1, 1]),
])
def test_gets_correct_bed_starts(upi, assembly, expected):
    assert fetch_data(upi, assembly).starts() == expected


@pytest.mark.skip
@pytest.mark.parametrize('upi,assembly,expected', [
    ('URS0000368518_9606', "GRCh38", []),
    ('URS000001B2EC_9606', "GRCh38", []),
    ('URS000071014B_112509', 'IBSC_v2', []),
    ('URS000099C6E5_9598', 'Pan_tro_3.0', []),
    ('URS0000000055_9606', 'GRCh38', [1, 1, 1]),
])
def test_gets_generates_expected_writeable(upi, assembly, expected):
    assert fetch_data(upi, assembly).writeable() == expected
