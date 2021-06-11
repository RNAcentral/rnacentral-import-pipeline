# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
from furl import furl

from rnacentral_pipeline.databases.gtrnadb import urls as fetch

REMOTE = furl("http://trna.ucsc.edu/download/RNAcentral/export2019/")


def test_gets_correct_urls():
    assert fetch.urls_for(REMOTE) == [
        "http://trna.ucsc.edu/download/RNAcentral/export2019/archaea_tRNAs.json.gz",
        "http://trna.ucsc.edu/download/RNAcentral/export2019/bacteria_tRNAs.tar.gz",
        "http://trna.ucsc.edu/download/RNAcentral/export2019/eukaryotes_tRNAs.tar.gz",
    ]


def test_can_extract_ursl_from_text():
    with open("data/gtrnadb/files.html", "r") as raw:
        urls = list(fetch.extract_urls(REMOTE, raw.read()))
        assert urls == [
            furl(
                "http://trna.ucsc.edu/download/RNAcentral/export2019/archaea_tRNAs.json.gz"
            ),
            furl(
                "http://trna.ucsc.edu/download/RNAcentral/export2019/bacteria_tRNAs.tar.gz"
            ),
            furl(
                "http://trna.ucsc.edu/download/RNAcentral/export2019/eukaryotes_tRNAs.tar.gz"
            ),
        ]
