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

import shutil
from pathlib import Path

import pytest

from rnacentral_pipeline.databases.helpers.hashes import md5
from rnacentral_pipeline.rnacentral.r2dt import data, parser

BASE = Path("data/r2dt/output")
MODEL_INFO = BASE / "model-info.json"


def parse(base=BASE, allow_missing=False):
    with MODEL_INFO.open("r") as handle:
        return list(parser.parse(handle, Path(base), allow_missing=allow_missing))


@pytest.fixture(scope="module")
def parsed():
    return {r.urs: r for r in parse()}


@pytest.mark.r2dt
def test_load_model_info_indexes_by_name():
    with MODEL_INFO.open("r") as handle:
        info = parser.load_model_info(handle)
    assert info["SAM"].alias == "RF00162"
    assert info["SAM"].source is data.Source.rfam
    # gtrnadb names are registered with both dashes and underscores
    assert info["E-Gln"] is info["E_Gln"]


@pytest.mark.r2dt
def test_can_process_a_directory(parsed):
    assert set(parsed) == {
        "URS00000F9D45_9606",
        "URS0000C5FF65",
        "URS0000A7635A",
        "URS0000A0BF23",
    }


@pytest.mark.r2dt
def test_parse_rejects_a_missing_directory():
    with pytest.raises(ValueError):
        parse(base="data/r2dt/no-such-directory")


@pytest.mark.r2dt
def test_parse_is_empty_without_metadata(tmp_path):
    assert parse(base=tmp_path) == []


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "urs,source,model_id",
    [
        ("URS00000F9D45_9606", data.Source.crw, 1),
        ("URS0000C5FF65", data.Source.ribovision, 2),
        ("URS0000A7635A", data.Source.rfam, 3),
        ("URS0000A0BF23", data.Source.gtrnadb, 4),
    ],
)
def test_assigns_the_expected_model(parsed, urs, source, model_id):
    assert parsed[urs].source is source
    assert parsed[urs].model_id == model_id


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "urs,svg,fasta",
    [
        (
            "URS00000F9D45_9606",
            "svg/URS00000F9D45_9606-d.16.b.C.perfringens.colored.svg",
            "fasta/URS00000F9D45_9606-d.16.b.C.perfringens.fasta",
        ),
        (
            "URS0000C5FF65",
            "svg/URS0000C5FF65-EC_LSU_3D.colored.svg",
            "fasta/URS0000C5FF65-EC_LSU_3D.fasta",
        ),
        (
            "URS0000A7635A",
            "svg/URS0000A7635A.colored.svg",
            "fasta/URS0000A7635A-RF00162.fasta",
        ),
        (
            "URS0000A0BF23",
            "svg/URS0000A0BF23-E-Gln.colored.svg",
            "fasta/URS0000A0BF23-E_Gln.fasta",
        ),
    ],
)
def test_can_produce_correct_paths(parsed, urs, svg, fasta):
    info = parsed[urs].info
    assert info.svg == BASE / "results" / svg
    assert info.fasta == BASE / "results" / fasta


@pytest.mark.r2dt
def test_attaches_ribovore_hits_only_where_they_exist(parsed):
    assert parsed["URS00000F9D45_9606"].hit_info.model == "d.16.b.C.perfringens"
    assert parsed["URS0000C5FF65"].hit_info.model == "EC_LSU_3D"
    # gtrnadb has no ribotyper stage, and there is no hit for the Rfam URS
    assert parsed["URS0000A0BF23"].hit_info is None
    assert parsed["URS0000A7635A"].hit_info is None


@pytest.mark.r2dt
def test_gets_correct_count(parsed):
    val = parsed["URS0000A7635A"]
    assert val.overlap_count() == 0
    assert val.basepair_count() == 24


@pytest.mark.r2dt
def test_produces_valid_data_for_rfam(parsed):
    val = parsed["URS0000A7635A"]
    writeable = val.writeable()

    # The SVG is too long to include here, so compare it to the file it reads.
    with (BASE / "results/svg/URS0000A7635A.colored.svg").open("r") as raw:
        assert raw.read().replace("\n", "") == val.svg()

    assert writeable == [
        "URS0000A7635A",
        3,
        "((((((((......(((...(((.....)))......))).(((.(((......))))))........((((......))))...)))))))).",
        0,
        24,
        None,
        None,
        None,
        None,
        None,
        True,
    ]


@pytest.mark.r2dt
def test_produces_valid_data_for_ribovision(parsed):
    val = parsed["URS0000C5FF65"]
    assert val.writeable()[3:] == [3, 854, 2, 2902, 2, 2890, 0.999, True]


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "urs,md5_hash",
    [
        ("URS00000F9D45_9606", "2204b2f0ac616b8366a3b5f37aa123b8"),
        ("URS0000A7635A", "9504c4b9a1cea77fa2c4ef8082d7b996"),
    ],
)
def test_can_extract_expected_svg_data(parsed, urs, md5_hash):
    svg = parsed[urs].svg()
    assert "\n" not in svg
    assert svg.startswith("<svg")
    assert md5(svg.encode()) == md5_hash


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "urs,secondary,bp_count",
    [
        (
            "URS00000F9D45_9606",
            "(((((((((....((((((((.....((((((............))))..))....)))))).)).(((((......((.((.(((....))))).)).....))))).)))))))))...",
            35,
        ),
        (
            "URS0000A7635A",
            "((((((((......(((...(((.....)))......))).(((.(((......))))))........((((......))))...)))))))).",
            24,
        ),
        (
            "URS0000A0BF23",
            "(((((((..((((........)))).(((((.......))))).....(((((.......)))))))))))).",
            21,
        ),
    ],
)
def test_can_extract_expected_dot_bracket_data(parsed, urs, secondary, bp_count):
    val = parsed[urs]
    assert val.dot_bracket() == secondary
    assert val.basepair_count() == bp_count
    assert val.modeled_length() == len(secondary)


@pytest.mark.r2dt
def test_skips_models_that_are_not_in_the_mapping(tmp_path):
    shutil.copytree(BASE, tmp_path / "output")
    metadata = tmp_path / "output/results/tsv/metadata.tsv"
    metadata.write_text(metadata.read_text() + "URS0000A7635B\tno-such-model\trfam\n")
    assert "URS0000A7635B" not in {r.urs for r in parse(base=tmp_path / "output")}


@pytest.mark.r2dt
def test_prefers_a_non_rfam_hit_over_an_rfam_one(tmp_path):
    shutil.copytree(BASE, tmp_path / "output")
    metadata = tmp_path / "output/results/tsv/metadata.tsv"
    # the same URS drawn against both an Rfam and a CRW model
    metadata.write_text(
        "URS00000F9D45_9606\tSAM\trfam\n"
        "URS00000F9D45_9606\td.16.b.C.perfringens\tcrw\n"
    )
    results = parse(base=tmp_path / "output")
    assert [(r.urs, r.source) for r in results] == [
        ("URS00000F9D45_9606", data.Source.crw)
    ]


@pytest.mark.r2dt
def test_missing_ribotyper_data_is_fatal_unless_allowed(tmp_path):
    shutil.copytree(BASE, tmp_path / "output")
    (tmp_path / "output/rfam/.ribotyper.long.out").unlink()

    with pytest.raises(ValueError):
        parse(base=tmp_path / "output")

    assert len(parse(base=tmp_path / "output", allow_missing=True)) == 4


@pytest.mark.r2dt
def test_missing_result_files_are_fatal_unless_allowed(tmp_path):
    shutil.copytree(BASE, tmp_path / "output")
    (tmp_path / "output/results/svg/URS0000A7635A.colored.svg").unlink()

    with pytest.raises(ValueError):
        parse(base=tmp_path / "output")

    results = parse(base=tmp_path / "output", allow_missing=True)
    assert "URS0000A7635A" not in {r.urs for r in results}
