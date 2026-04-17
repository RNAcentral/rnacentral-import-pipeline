# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import pytest

from rnacentral_pipeline.databases.generic.parser import parse


@pytest.mark.skip()
def test_runs_validation_on_data():
    pass


@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/json-schema/v020/circatlas.json", 3),
        ("data/json-schema/v020/flybase.json", 5),
        ("data/json-schema/v020/flybase-scaRNA.json", 1),
        ("data/json-schema/v020/japonicusdb.json", 3),
        ("data/json-schema/v020/lincipedia.json", 1),
        ("data/json-schema/v020/lncbook.json", 3),
        ("data/json-schema/v020/lncipedia-5.0.json", 1),
        ("data/json-schema/v020/lncipedia-with-isoforms.json", 5),
        ("data/json-schema/v020/missing-mirbase.json", 2),
        ("data/json-schema/v020/modomics-modifications.json", 1),
        ("data/json-schema/v020/pombase.json", 1),
        ("data/json-schema/v020/shift-mirbase.json", 1),
        ("data/json-schema/v020/shift-mirbase-2.json", 1),
        ("data/json-schema/v020/tarbase.json", 1),
    ],
)
def test_can_parse_v0_2_0_data(filename, count):
    with open(filename, "rb") as raw:
        assert len(list(parse(raw))) == count
