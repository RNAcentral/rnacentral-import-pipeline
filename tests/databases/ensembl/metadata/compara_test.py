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

import pytest

from rnacentral_pipeline.databases.ensembl.metadata import compara


@pytest.fixture
def simple():
    with open("data/ensembl/compara/simple.fasta", "r") as raw:
        yield raw


@pytest.mark.ensembl
def test_can_produce_expected_data(simple):
    assert list(compara.data(simple)) == [
        (
            "707b0ba6ac8febddca12a3126ff8bcfd7fd2838411b4316235162c43d0854907",
            [
                "ENSMUST00000082483",
                "ENSRNOT00000078435",
                "ENSRNOT00000082313",
                "ENSRNOT00000082656",
                "ENSRNOT00000082946",
                "ENSRNOT00000083118",
                "ENSRNOT00000083605",
                "ENSRNOT00000086169",
                "ENSRNOT00000086419",
                "ENSRNOT00000090157",
            ],
        ),
        (
            "2be126e36b38bb326dcde5358cfd1d9b76b8b8e46ea2b9876dacbca8ef894943",
            [
                "ENSMUST00000179813",
                "ENSMUST00000180265",
                "ENSRNOT00000078860",
                "ENSRNOT00000079197",
                "ENSRNOT00000087954",
                "ENSRNOT00000088084",
                "ENSRNOT00000091107",
            ],
        ),
    ]
