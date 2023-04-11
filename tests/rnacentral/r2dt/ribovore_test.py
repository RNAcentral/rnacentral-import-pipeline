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

from rnacentral_pipeline.rnacentral.r2dt import ribovore


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/traveler/crw/crw.ribotyper.long.out", 1),
        ("data/traveler/failed-data.ribotyper.long.out", 12),
    ],
)
def test_can_parse_whole_file(filename, count):
    assert len(list(ribovore.parse(Path(filename)))) == count


@pytest.mark.r2dt
def test_can_parse_a_simple_result():
    path = Path("data/traveler/crw/crw.ribotyper.long.out")
    data = list(ribovore.parse(path))
    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(
        ribovore.RibovoreResult(
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
def test_can_produce_dict_of_results():
    data = ribovore.as_dict(Path("data/traveler/crw"))
    assert data == {
        "URS00000F9D45_9606": ribovore.RibovoreResult(
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
