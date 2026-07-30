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

import csv
import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import r2dt
from rnacentral_pipeline.rnacentral.r2dt import s3

BASE = Path("data/r2dt/output").absolute()
MODEL_INFO = str(BASE / "model-info.json")


@pytest.mark.cli
@pytest.mark.r2dt
def test_can_process_svgs():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            r2dt.cli, ["process-svgs", MODEL_INFO, str(BASE), "r2dt.csv"]
        )
        assert result.exit_code == 0, result.output

        with open("r2dt.csv", "r") as raw:
            rows = list(csv.reader(raw))
        assert {r[0] for r in rows} == {
            "URS00000F9D45_9606",
            "URS0000C5FF65",
            "URS0000A7635A",
            "URS0000A0BF23",
        }


@pytest.mark.cli
@pytest.mark.r2dt
def test_can_prepare_s3_files():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            r2dt.cli,
            ["prepare-s3", MODEL_INFO, str(BASE), "svgs", "files.txt"],
        )
        assert result.exit_code == 0, result.output

        listed = Path("files.txt").read_text().split()
        assert len(listed) == 4
        for path in listed:
            with gzip.open(path, "rt") as raw:
                assert raw.read().startswith("<svg")


@pytest.mark.cli
@pytest.mark.r2dt
def test_upload_s3_fails_loudly(monkeypatch, tmp_path):
    class _Broken:
        def put_object(self, **kwargs):
            raise RuntimeError("simulated S3 500")

    monkeypatch.setattr(s3, "client", lambda *a, **k: _Broken())
    listing = tmp_path / "files.txt"
    svg = tmp_path / "URS0000F7F700.svg.gz"
    svg.write_bytes(b"body")
    listing.write_text(f"{svg}\n")

    result = CliRunner().invoke(r2dt.cli, ["upload-s3", str(listing)])
    assert result.exit_code != 0
    assert isinstance(result.exception, RuntimeError)


@pytest.mark.cli
@pytest.mark.r2dt
def test_verify_s3_exits_non_zero_when_objects_are_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(s3, "verify", lambda *a, **k: 3)
    listing = tmp_path / "files.txt"
    listing.write_text("")

    result = CliRunner().invoke(r2dt.cli, ["verify-s3", str(listing)])
    assert result.exit_code != 0
    assert "3 files missing or mismatched" in result.output


@pytest.mark.cli
@pytest.mark.r2dt
def test_list_s3_writes_urs_to_output(monkeypatch, tmp_path):
    def fake_list(env, out, **kwargs):
        out.write("URS0000F7F700\nURS0000F7F701\n")
        return 2

    monkeypatch.setattr(s3, "list_svgs", fake_list)

    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(r2dt.cli, ["list-s3", "listing.txt"])
        assert result.exit_code == 0, result.output
        assert Path("listing.txt").read_text().split() == [
            "URS0000F7F700",
            "URS0000F7F701",
        ]
