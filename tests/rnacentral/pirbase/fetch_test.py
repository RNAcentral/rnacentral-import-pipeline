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

from rnacentral_pipeline.rnacentral.pirbase import fetch


def test_can_extract_base_url():
    assert fetch.base_url(
        furl("http://www.regulatoryrna.org/database/piRNA/forJSON.html")
    ) == furl("http://www.regulatoryrna.org/database/piRNA")


def test_can_extract_urls():
    base = furl("http://www.regulatoryrna.org/database/piRNA/")

    with open("data/pirbase/links.html") as raw:
        urls = fetch.extract_urls(base, raw.read())
    assert urls == [
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_hsa.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_mmu.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_rno.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dme.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_cel.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dre.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_gga.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_xtr.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_bmo.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_nve.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_bta.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_mfa.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_mml.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_cja.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_aca.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_tbe.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_ssc.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_der.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dya.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dvi.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_ocu.json.gz"),
    ]


def test_can_preform_full_fetch():
    urls = fetch.find_urls(furl("http://www.regulatoryrna.org/database/piRNA/forJSON.html"))
    assert urls == [
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_hsa.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_mmu.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_rno.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dme.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_cel.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dre.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_gga.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_xtr.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_bmo.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_nve.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_bta.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_mfa.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_mml.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_cja.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_aca.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_tbe.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_ssc.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_der.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dya.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_dvi.json.gz"),
        furl("http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/json/piRBase_ocu.json.gz"),
    ]
