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

from rnacentral_pipeline.databases.rgd.helpers import indexed
from rnacentral_pipeline.databases.rgd import parser as rgd


@pytest.fixture
def sequences():
    with indexed("data/rgd/sequences.fa.gz") as sequences:
        yield sequences


def test_can_find_version():
    with open("data/rgd/rat_genes.txt", "r") as raw:
        assert rgd.get_version(raw) == "genes-version-2.2.5"


@pytest.mark.xfail()
def test_can_parse_data(sequences):
    with open("data/rgd/rat_genes.txt", "r") as raw:
        entries = list(rgd.parse(raw, sequences))
    assert len(entries) == 14
