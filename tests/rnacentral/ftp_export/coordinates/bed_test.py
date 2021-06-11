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
    data = fetch_data("URS000082BE64_9606", "GRCh38")
    assert attr.asdict(data) == attr.asdict(
        bed.BedEntry(
            rna_id="URS000082BE64_9606",
            rna_type="snoRNA",
            databases="snOPY",
            region=regions.SequenceRegion(
                assembly_id="GRCh38",
                chromosome="3",
                strand=1,
                exons=[regions.Exon(start=32676135, stop=32676226)],
                coordinate_system=regions.CoordinateSystem.zero_based(),
            ),
        )
    )


def test_can_build_entry_with_several_databases():
    data = fetch_data("URS0000368518_9606", "GRCh38")
    assert attr.asdict(data) == attr.asdict(
        bed.BedEntry(
            rna_id="URS0000368518_9606",
            rna_type="lncRNA",
            databases="Ensembl,GENCODE",
            region=regions.SequenceRegion(
                assembly_id="GRCh38",
                chromosome="2",
                strand=1,
                exons=[regions.Exon(start=37208874, stop=37212677)],
                coordinate_system=regions.CoordinateSystem.zero_based(),
            ),
        )
    )


def test_can_build_entry_from_pig():
    data = fetch_data("URS000099C6E5_9598", "Pan_tro_3.0")
    assert attr.asdict(data) == attr.asdict(
        bed.BedEntry(
            rna_id="URS000099C6E5_9598",
            rna_type="SRP_RNA",
            databases="Rfam",
            region=regions.SequenceRegion(
                assembly_id="Pan_tro_3.0",
                chromosome="X",
                strand=-1,
                exons=[regions.Exon(start=73819429, stop=73819577)],
                coordinate_system=regions.CoordinateSystem.zero_based(),
            ),
        )
    )


def test_can_build_entry_with_several_exons():
    data = fetch_data("URS0000000055_9606", "GRCh38")
    assert attr.asdict(data) == attr.asdict(
        bed.BedEntry(
            rna_id="URS0000000055_9606",
            rna_type="lncRNA",
            databases="Ensembl,GENCODE,LNCipedia,NONCODE",
            region=regions.SequenceRegion(
                assembly_id="GRCh38",
                chromosome="6",
                strand=-1,
                exons=[
                    regions.Exon(start=57171004, stop=57171098),
                    regions.Exon(start=57173375, stop=57173480),
                    regions.Exon(start=57173736, stop=57174236),
                ],
                coordinate_system=regions.CoordinateSystem.zero_based(),
            ),
        )
    )


@pytest.mark.parametrize(
    "upi,assembly,expected",
    [
        ("URS0000368518_9606", "GRCh38", "chr2"),
        ("URS000001B2EC_9606", "GRCh38", "chrM"),
        pytest.param(
            "URS000071014B_112509", "IBSC_v2", "chr3H", marks=pytest.mark.xfail
        ),
    ],
)
def test_gets_correct_chromosome(upi, assembly, expected):
    assert fetch_data(upi, assembly).bed_chromosome == expected


@pytest.mark.parametrize(
    "upi,assembly,expected",
    [
        ("URS0000368518_9606", "GRCh38", [3803]),
        ("URS000001B2EC_9606", "GRCh38", [68]),
        pytest.param("URS000071014B_112509", "IBSC_v2", [1], marks=pytest.mark.xfail),
        ("URS0000000055_9606", "GRCh38", [94, 105, 500]),
    ],
)
def test_gets_correct_bed_sizes(upi, assembly, expected):
    assert fetch_data(upi, assembly).sizes() == expected


@pytest.mark.parametrize(
    "upi,assembly,expected",
    [
        ("URS0000368518_9606", "GRCh38", [0]),
        ("URS000001B2EC_9606", "GRCh38", [0]),
        pytest.param("URS000071014B_112509", "IBSC_v2", [0], marks=pytest.mark.xfail),
        ("URS0000000055_9606", "GRCh38", [0, 2371, 2732]),
    ],
)
def test_gets_correct_bed_starts(upi, assembly, expected):
    assert fetch_data(upi, assembly).starts() == expected


@pytest.mark.parametrize(
    "upi,assembly,expected",
    [
        (
            "URS0000368518_9606",
            "GRCh38",
            [
                "chr2",
                37208874,
                37212677,
                "URS0000368518_9606",
                0,
                "+",
                37208874,
                37212677,
                "63,125,151",
                1,
                "3803",
                "0",
                ".",
                "lncRNA",
                "Ensembl,GENCODE",
            ],
        ),
        (
            "URS000001B2EC_9606",
            "GRCh38",
            [
                "chrM",
                9990,
                10058,
                "URS000001B2EC_9606",
                0,
                "+",
                9990,
                10058,
                "63,125,151",
                1,
                "68",
                "0",
                ".",
                "tRNA",
                "ENA",
            ],
        ),
        # ('URS000071014B_112509', 'IBSC_v2', []),
        (
            "URS000099C6E5_9598",
            "Pan_tro_3.0",
            [
                "chrX",
                73819429,
                73819577,
                "URS000099C6E5_9598",
                0,
                "-",
                73819429,
                73819577,
                "63,125,151",
                1,
                "148",
                "0",
                ".",
                "SRP_RNA",
                "Rfam",
            ],
        ),
        (
            "URS0000000055_9606",
            "GRCh38",
            [
                "chr6",
                57171004,
                57174236,
                "URS0000000055_9606",
                0,
                "-",
                57171004,
                57174236,
                "63,125,151",
                3,
                "94,105,500",
                "0,2371,2732",
                ".",
                "lncRNA",
                "Ensembl,GENCODE,LNCipedia,NONCODE",
            ],
        ),
    ],
)
def test_gets_generates_expected_writeable(upi, assembly, expected):
    assert fetch_data(upi, assembly).writeable() == expected
