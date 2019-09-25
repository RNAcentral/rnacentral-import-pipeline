# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import tempfile

import pytest

from rnacentral_pipeline.databases.ncbi.gene import fetch
from rnacentral_pipeline.databases.ncbi.gene import helpers
from rnacentral_pipeline.databases.ncbi.gene import parser


def parse(filename):
    with tempfile.NamedTemporaryFile() as tmp:
        with open(filename, 'r') as raw:
            ncrnas = list(helpers.ncrnas(raw))
            fetch.fetch_and_write(ncrnas, tmp)
        tmp.flush()
        tmp.seek(0)
        return list(parser.parse(tmp))


@pytest.mark.parametrize('filename,count', [
    ('data/ncbi_gene/simple.txt', 11),
])
def test_parses_all_data(filename, count):
    assert len(parse(filename)) == count
