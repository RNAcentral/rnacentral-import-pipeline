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

import json
from io import StringIO

import pytest

from rnacentral_pipeline.databases.lncbook import parser

RAW_DATA = "data/json-schema/v020/lncbook.json"


def load_raw():
    with open(RAW_DATA, "r") as raw:
        return json.load(raw)


def test_rejects_data_with_no_grch38_entries():
    with open(RAW_DATA, "r") as raw:
        with pytest.raises(ValueError):
            list(parser.parse(raw))


def test_keeps_only_grch38_regions_and_remaps_mitochondrial_chromosome():
    data = load_raw()
    entry = data["data"][0]
    entry["genomeLocations"][0]["assembly"] = "GRCh38"
    entry["genomeLocations"][0]["exons"][0]["chromosome"] = "chrM"
    entry["sequence"] = entry["sequence"].lower()

    handle = StringIO(json.dumps(data))
    entries = list(parser.parse(handle))

    assert len(entries) == 1
    (parsed,) = entries
    assert parsed.sequence == entry["sequence"].upper()
    assert len(parsed.regions) == 1
    assert parsed.regions[0].chromosome == "MT"


def test_drops_entries_with_no_grch38_regions_but_keeps_others():
    data = load_raw()
    data["data"][0]["genomeLocations"][0]["assembly"] = "GRCh38"

    handle = StringIO(json.dumps(data))
    entries = list(parser.parse(handle))

    assert len(entries) == 1
