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

import importlib.util
import io
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "bin" / "json2fasta.py"
SPEC = importlib.util.spec_from_file_location("json2fasta", MODULE_PATH)
json2fasta = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(json2fasta)


@pytest.mark.utils
def test_parse_can_read_multiline_json_entry():
    raw = io.StringIO('{"id":"URS0001",\n"description":"sample","sequence":"ACGUN"}\n')

    parsed = list(json2fasta.parse(raw))

    assert parsed == [
        {"id": "URS0001", "description": "sample", "sequence": "ACGUN"},
    ]


@pytest.mark.utils
def test_sequences_emits_fasta_records_from_multiline_json():
    raw = io.StringIO(
        '{"id":"URS0001",\n"description":"sample description","sequence":"ACGUN"}\n'
    )

    records = list(json2fasta.sequences(raw))

    assert len(records) == 1
    assert records[0].id == "URS0001"
    assert str(records[0].seq) == "ACGUN"
    assert records[0].description == "sample description"


@pytest.mark.utils
def test_parse_raises_for_incomplete_json_at_eof():
    raw = io.StringIO('{"id":"URS0001","sequence":"ACGUN"')

    with pytest.raises(json2fasta.json.JSONDecodeError):
        list(json2fasta.parse(raw))
