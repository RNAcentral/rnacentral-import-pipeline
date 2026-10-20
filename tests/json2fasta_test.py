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
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "bin" / "json2fasta.py"
SPEC = importlib.util.spec_from_file_location("json2fasta", MODULE_PATH)
json2fasta = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(json2fasta)


@pytest.mark.utils
def test_parse_reads_one_entry_per_line():
    raw = io.StringIO(
        '{"id": "URS0001", "description": "sample description", "sequence": "ACGUN"}\n'
    )

    parsed = list(json2fasta.parse(raw))

    assert parsed == [
        {"id": "URS0001", "description": "sample description", "sequence": "ACGUN"},
    ]


@pytest.mark.utils
def test_sequences_emits_fasta_records():
    raw = io.StringIO(
        '{"id": "URS0001", "description": "sample description", "sequence": "ACGUN"}\n'
    )

    records = list(json2fasta.sequences(raw))

    assert len(records) == 1
    assert records[0].id == "URS0001"
    assert str(records[0].seq) == "ACGUN"
    assert records[0].description == "sample description"


@pytest.mark.utils
def test_parse_fails_fast_on_an_entry_split_across_lines(caplog):
    # A dump must emit exactly one compact JSON object per line - buffering
    # across lines to tolerate a split entry is what caused a 256GB OOM in
    # production (a single bad line would otherwise accumulate the rest of
    # the multi-hundred-GB dump into one string). parse() must fail on the
    # first line instead, not wait for the object to complete.
    raw = io.StringIO('{"id":"URS0001",\n"description":"sample","sequence":"ACGUN"}\n')

    with caplog.at_level(logging.ERROR):
        with pytest.raises(json.JSONDecodeError):
            list(json2fasta.parse(raw))

    assert "Undecodable JSON on line 1" in caplog.text


@pytest.mark.utils
def test_parse_raises_for_incomplete_json_at_eof():
    raw = io.StringIO('{"id":"URS0001","sequence":"ACGUN"')

    with pytest.raises(json.JSONDecodeError):
        list(json2fasta.parse(raw))


@pytest.mark.utils
def test_script_runs_without_pythonpath(tmp_path):
    # Nextflow runs the script from a work dir in a container where the package
    # is not installed, so it must find rnacentral_pipeline on its own.
    raw = tmp_path / "in.json"
    raw.write_text('{"id": "URS0001", "description": "d", "sequence": "ACGU"}\n')
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    subprocess.run(
        [sys.executable, str(MODULE_PATH), str(raw), "out.fasta"],
        cwd=tmp_path,
        env=env,
        check=True,
    )
    assert (tmp_path / "out.fasta").read_text().startswith(">URS0001")
