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

from pathlib import Path

import attr
import pytest

from rnacentral_pipeline import ribovore
from rnacentral_pipeline.databases.data import RibovoreResult
from rnacentral_pipeline.rnacentral.r2dt import ribovore as r2dt_ribovore


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/r2dt/crw/crw.ribotyper.long.out", 1),
        ("data/r2dt/failed-data.ribotyper.long.out", 16),
        ("data/r2dt/ribovision/.ribotyper.long.out", 300),
    ],
)
def test_can_parse_whole_file(filename, count):
    assert len(list(ribovore.parse_file(Path(filename)))) == count


@pytest.mark.r2dt
def test_can_parse_a_simple_result():
    path = Path("data/r2dt/crw/crw.ribotyper.long.out")
    results = list(ribovore.parse_file(path))
    assert len(results) == 1
    assert attr.asdict(results[0]) == attr.asdict(
        RibovoreResult(
            target="URS00000F9D45_9606",
            status="PASS",
            length=1588,
            fm=1,
            fam="SSU",
            domain="Bacteria",
            model="d.16.b.C.perfringens",
            strand=1,
            ht=1,
            tscore=1093.0,
            bscore=1093.0,
            bevalue=0.0,
            tcov=0.999,
            bcov=0.999,
            bfrom=3,
            bto=1588,
            mfrom=3,
            mto=1512,
        )
    )


@pytest.mark.r2dt
def test_parse_directory_finds_prefixed_and_hidden_files():
    # ribotyper writes either <dirname>.ribotyper.long.out or .ribotyper.long.out
    assert len(list(ribovore.parse_directory(Path("data/r2dt/crw")))) == 1
    assert len(list(ribovore.parse_directory(Path("data/r2dt/ribovision")))) == 300


@pytest.mark.r2dt
def test_parse_directory_raises_when_there_is_no_result_file(tmp_path):
    with pytest.raises(ribovore.MissingRibotyperDataException):
        list(ribovore.parse_directory(tmp_path))


@pytest.mark.r2dt
def test_can_produce_dict_of_results():
    assert r2dt_ribovore.as_dict(Path("data/r2dt/crw")) == {
        "URS00000F9D45_9606": RibovoreResult(
            target="URS00000F9D45_9606",
            status="PASS",
            length=1588,
            fm=1,
            fam="SSU",
            domain="Bacteria",
            model="d.16.b.C.perfringens",
            strand=1,
            ht=1,
            tscore=1093.0,
            bscore=1093.0,
            bevalue=0.0,
            tcov=0.999,
            bcov=0.999,
            bfrom=3,
            bto=1588,
            mfrom=3,
            mto=1512,
        )
    }


@pytest.mark.r2dt
def test_as_dict_is_keyed_by_target_and_drops_failures():
    results = r2dt_ribovore.as_dict(Path("data/r2dt/ribovision"))
    assert len(results) == 248
    assert all(target == hit.target for target, hit in results.items())
    assert all(hit.status != "FAIL" for hit in results.values())
    assert "URS0000C5FF65" in results
    assert "URS000083FA52" not in results  # FAIL in the fixture


@pytest.mark.r2dt
def test_as_dict_respects_allow_missing(tmp_path):
    with pytest.raises(ValueError):
        r2dt_ribovore.as_dict(tmp_path)
    assert r2dt_ribovore.as_dict(tmp_path, allow_missing=True) == {}
