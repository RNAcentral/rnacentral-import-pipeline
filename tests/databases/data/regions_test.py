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

import six.moves.cPickle as pickle

import pytest

from rnacentral_pipeline.databases.data.regions import *


@pytest.mark.parametrize(
    "raw,expected",
    [
        (1, Strand.forward),
        (1.0, Strand.forward),
        ("1", Strand.forward),
        ("+", Strand.forward),
        (-1, Strand.reverse),
        ("-1", Strand.reverse),
        ("-", Strand.reverse),
        (0, Strand.unknown),
        (0.0, Strand.unknown),
        (".", Strand.unknown),
    ],
)
def test_can_convert_a_strand(raw, expected):
    assert Strand.build(raw) is expected


@pytest.mark.parametrize(
    "raw",
    [
        1.1,
        -0.2,
        "bob",
    ],
)
def test_fails_with_bad_strands(raw):
    with pytest.raises(UnknownStrand):
        Strand.build(raw)


@pytest.mark.parametrize(
    "name,expected",
    [
        (
            "0-start, half-open",
            CoordinateSystem(basis=CoordinateStart.zero, close_status=CloseStatus.open),
        ),
        (
            "1-start, fully-closed",
            CoordinateSystem(
                basis=CoordinateStart.one, close_status=CloseStatus.closed
            ),
        ),
    ],
)
def test_can_build_coordinate_system_from_name(name, expected):
    assert CoordinateSystem.from_name(name) == expected


@pytest.mark.parametrize(
    "name",
    [
        "0-start, half-open",
        "1-start, fully-closed",
    ],
)
def test_can_generate_coordinate_name(name):
    assert CoordinateSystem.from_name(name).name() == name


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=12), Exon(start=11, stop=12)),
        ("1-start, fully-closed", Exon(start=10, stop=12), Exon(start=10, stop=12)),
    ],
)
def test_can_correctly_normalize_an_exon(name, exon, expected):
    assert CoordinateSystem.from_name(name).normalize(exon) == expected


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=10), Exon(start=10, stop=10)),
        ("0-start, half-open", Exon(start=10, stop=11), Exon(start=10, stop=11)),
        ("0-start, half-open", Exon(start=10, stop=12), Exon(start=10, stop=12)),
        ("1-start, fully-closed", Exon(start=10, stop=10), Exon(start=9, stop=10)),
        ("1-start, fully-closed", Exon(start=10, stop=11), Exon(start=9, stop=11)),
        ("1-start, fully-closed", Exon(start=10, stop=12), Exon(start=9, stop=12)),
    ],
)
def test_can_correctly_switch_to_zero_based(name, exon, expected):
    assert CoordinateSystem.from_name(name).as_zero_based(exon) == expected


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=11), Exon(start=11, stop=11)),
        ("0-start, half-open", Exon(start=10, stop=12), Exon(start=11, stop=12)),
        ("1-start, fully-closed", Exon(start=10, stop=10), Exon(start=10, stop=10)),
        ("1-start, fully-closed", Exon(start=10, stop=11), Exon(start=10, stop=11)),
        ("1-start, fully-closed", Exon(start=10, stop=12), Exon(start=10, stop=12)),
    ],
)
def test_can_correctly_to_one_based(name, exon, expected):
    assert CoordinateSystem.from_name(name).normalize(exon) == expected


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=10), 0),
        ("0-start, half-open", Exon(start=10, stop=11), 1),
        ("0-start, half-open", Exon(start=10, stop=12), 2),
        ("1-start, fully-closed", Exon(start=10, stop=10), 1),
        ("1-start, fully-closed", Exon(start=10, stop=11), 2),
        ("1-start, fully-closed", Exon(start=10, stop=12), 3),
    ],
)
def test_can_correctly_compute_lengths(name, exon, expected):
    assert CoordinateSystem.from_name(name).size(exon) == expected


@pytest.mark.parametrize(
    "upi,strand,coordinate_system,expected",
    [
        (
            "",
            Strand.forward,
            CoordinateSystem.from_name("0-start, half-open"),
            "@4/11-22,31-40:+",
        ),
        (
            "U1",
            Strand.forward,
            CoordinateSystem.from_name("0-start, half-open"),
            "U1@4/11-22,31-40:+",
        ),
        (
            "",
            Strand.reverse,
            CoordinateSystem.from_name("0-start, half-open"),
            "@4/11-22,31-40:-",
        ),
        (
            "U2",
            Strand.reverse,
            CoordinateSystem.from_name("0-start, half-open"),
            "U2@4/11-22,31-40:-",
        ),
        (
            "",
            Strand.forward,
            CoordinateSystem.from_name("1-start, fully-closed"),
            "@4/10-22,30-40:+",
        ),
        (
            "U3",
            Strand.forward,
            CoordinateSystem.from_name("1-start, fully-closed"),
            "U3@4/10-22,30-40:+",
        ),
        (
            "",
            Strand.reverse,
            CoordinateSystem.from_name("1-start, fully-closed"),
            "@4/10-22,30-40:-",
        ),
        (
            "U2",
            Strand.reverse,
            CoordinateSystem.from_name("1-start, fully-closed"),
            "U2@4/10-22,30-40:-",
        ),
    ],
)
def test_can_generate_correct_region_names(upi, strand, coordinate_system, expected):
    region = SequenceRegion(
        assembly_id="GRCh38",
        chromosome="4",
        strand=strand,
        exons=[
            Exon(start=10, stop=22),
            Exon(start=30, stop=40),
        ],
        coordinate_system=coordinate_system,
    )
    assert region.name(upi=upi) == expected


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=10), 0),
        ("0-start, half-open", Exon(start=10, stop=11), 1),
        ("0-start, half-open", Exon(start=10, stop=12), 2),
        ("1-start, fully-closed", Exon(start=10, stop=10), 1),
        ("1-start, fully-closed", Exon(start=10, stop=11), 2),
        ("1-start, fully-closed", Exon(start=10, stop=12), 3),
    ],
)
def test_region_can_correctly_compute_lengths(name, exon, expected):
    region = SequenceRegion(
        assembly_id="GRCh38",
        chromosome="4",
        strand=-1,
        exons=[exon],
        coordinate_system=CoordinateSystem.from_name(name),
    )
    assert region.sizes() == [expected]


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=12), Exon(start=10, stop=12)),
        ("1-start, fully-closed", Exon(start=10, stop=12), Exon(start=9, stop=12)),
    ],
)
def test_region_can_become_zero_based(name, exon, expected):
    region = SequenceRegion(
        assembly_id="GRCh38",
        chromosome="4",
        strand=-1,
        exons=[exon],
        coordinate_system=CoordinateSystem.from_name(name),
    )
    assert region.as_zero_based() == attr.evolve(
        region,
        exons=[expected],
        coordinate_system=CoordinateSystem.from_name("0-start, half-open"),
    )


@pytest.mark.parametrize(
    "name,exon,expected",
    [
        ("0-start, half-open", Exon(start=10, stop=12), Exon(start=11, stop=12)),
        ("1-start, fully-closed", Exon(start=10, stop=12), Exon(start=10, stop=12)),
    ],
)
def test_region_can_become_one_based(name, exon, expected):
    region = SequenceRegion(
        assembly_id="GRCh38",
        chromosome="4",
        strand=1,
        exons=[exon],
        coordinate_system=CoordinateSystem.from_name(name),
    )
    assert region.as_one_based() == attr.evolve(
        region,
        exons=[expected],
        coordinate_system=CoordinateSystem.from_name("1-start, fully-closed"),
    )


@pytest.mark.parametrize(
    "accession,strand,coordinate_system,expected",
    [
        ("a1", Strand.unknown, CoordinateSystem.from_name("0-start, half-open"), []),
        (
            "a1",
            Strand.forward,
            CoordinateSystem.from_name("0-start, half-open"),
            [["a1", "@4/31-40:+", "4", 1, "GRCh38", 1, 31, 40]],
        ),
        (
            "U1",
            Strand.forward,
            CoordinateSystem.from_name("0-start, half-open"),
            [["U1", "U1@4/31-40:+", "4", 1, "GRCh38", 1, 31, 40]],
        ),
        ("U1", Strand.unknown, CoordinateSystem.from_name("0-start, half-open"), []),
        (
            "a3",
            Strand.reverse,
            CoordinateSystem.from_name("0-start, half-open"),
            [["a3", "@4/31-40:-", "4", -1, "GRCh38", 1, 31, 40]],
        ),
        (
            "U2",
            Strand.reverse,
            CoordinateSystem.from_name("0-start, half-open"),
            [["U2", "U2@4/31-40:-", "4", -1, "GRCh38", 1, 31, 40]],
        ),
        ("U2", Strand.unknown, CoordinateSystem.from_name("0-start, half-open"), []),
        (
            "a2",
            Strand.forward,
            CoordinateSystem.from_name("1-start, fully-closed"),
            [["a2", "@4/30-40:+", "4", 1, "GRCh38", 1, 30, 40]],
        ),
        ("a2", Strand.unknown, CoordinateSystem.from_name("1-start, fully-closed"), []),
        (
            "U3",
            Strand.forward,
            CoordinateSystem.from_name("1-start, fully-closed"),
            [["U3", "U3@4/30-40:+", "4", 1, "GRCh38", 1, 30, 40]],
        ),
        ("U3", Strand.unknown, CoordinateSystem.from_name("1-start, fully-closed"), []),
        (
            "A1",
            Strand.reverse,
            CoordinateSystem.from_name("1-start, fully-closed"),
            [["A1", "@4/30-40:-", "4", -1, "GRCh38", 1, 30, 40]],
        ),
        ("A1", Strand.unknown, CoordinateSystem.from_name("1-start, fully-closed"), []),
        (
            "U1",
            Strand.reverse,
            CoordinateSystem.from_name("1-start, fully-closed"),
            [["U1", "U1@4/30-40:-", "4", -1, "GRCh38", 1, 30, 40]],
        ),
        ("U1", Strand.unknown, CoordinateSystem.from_name("1-start, fully-closed"), []),
    ],
)
def test_can_generate_correct_writeable(accession, strand, coordinate_system, expected):
    region = SequenceRegion(
        assembly_id="GRCh38",
        chromosome="4",
        strand=strand,
        exons=[Exon(start=30, stop=40)],
        coordinate_system=coordinate_system,
    )
    upi = accession.startswith("U")
    assert list(region.writeable(accession, is_upi=upi)) == expected


@pytest.mark.parametrize(
    "data",
    [
        (Strand.reverse),
        (Strand.unknown),
        (Strand.forward),
        (CoordinateStart.zero),
        (CoordinateStart.one),
        (CloseStatus.closed),
        (CloseStatus.open),
    ],
)
def test_can_serialize_enums_to_and_from_pickle(data):
    assert pickle.loads(pickle.dumps(data)) is data


@pytest.mark.parametrize(
    "region",
    [
        (
            SequenceRegion(
                assembly_id="GRCh38",
                chromosome="1",
                strand=Strand.unknown,
                exons=[Exon(start=100000, stop=2000000), Exon(start=-1, stop=10000)],
                coordinate_system=CoordinateSystem.from_name("0-start, half-open"),
            )
        ),
        (
            SequenceRegion(
                assembly_id="GRCh38",
                chromosome="MT",
                strand=Strand.forward,
                exons=[Exon(start=10, stop=200), Exon(start=-1, stop=10000)],
                coordinate_system=CoordinateSystem.from_name("1-start, fully-closed"),
            )
        ),
    ],
)
def test_can_serialize_regions_to_and_from_pickle(region):
    data = attr.asdict(region)
    loaded = pickle.loads(pickle.dumps(data))
    assert loaded == data
    assert SequenceRegion(**loaded) == region
    assert pickle.loads(pickle.dumps(region)) == region
