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

import gzip

import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.genome_mapping import compare


def write_gz(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as out:
        out.write(text)


GFF_LINE = "1\tRNAcentral\tncRNA\t1\t9\t.\t+\t.\tID=URS0000000001_9606\n"
BED_LINE = "1\t0\t9\tURS0000000001_9606\t0\t+\n"


@pytest.mark.parametrize(
    "kind,name,line",
    [("gff", "a.gff3.gz", GFF_LINE), ("bed", "a.bed.gz", BED_LINE)],
)
def test_written_ids_reads_matching_files_only(tmp_path, kind, name, line):
    write_gz(tmp_path / name, "# header URS0000000009_1\n" + line)
    write_gz(tmp_path / "other.txt.gz", "URS0000000002_1\n")
    found = compare.written_ids(tmp_path, kind)
    assert found["urs_taxid"].to_list() == ["URS0000000001_9606"]


def test_missing_from_files(tmp_path, monkeypatch):
    write_gz(tmp_path / "a.bed.gz", BED_LINE)
    ids = ["URS0000000001_9606", "URS0000000002_9606"]
    monkeypatch.setattr(
        compare, "database_ids", lambda _: pl.DataFrame({"urs_taxid": ids})
    )
    missing = compare.missing_from_files(tmp_path, "ignored", "bed")
    assert missing["urs_taxid"].to_list() == ["URS0000000002_9606"]


def test_written_ids_empty_directory(tmp_path):
    found = compare.written_ids(tmp_path, "gff")
    assert found.is_empty()
    assert found.columns == ["urs_taxid"]


def test_written_ids_deduplicates_across_files(tmp_path):
    write_gz(tmp_path / "a.bed.gz", BED_LINE + BED_LINE)
    write_gz(tmp_path / "b.bed.gz", BED_LINE)
    assert compare.written_ids(tmp_path, "bed").height == 1


def test_missing_from_files_nothing_missing(tmp_path, monkeypatch):
    write_gz(tmp_path / "a.bed.gz", BED_LINE)
    monkeypatch.setattr(
        compare,
        "database_ids",
        lambda _: pl.DataFrame({"urs_taxid": ["URS0000000001_9606"]}),
    )
    assert compare.missing_from_files(tmp_path, "ignored", "bed").is_empty()


def test_database_ids_selects_active_with_coordinates(monkeypatch):
    class Cursor:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def execute(self, query):
            self.query = query

        def fetchall(self):
            return [("URS0000000001_9606",), ("URS0000000001_9606",)]

    class Conn:
        closed = False
        cursor_ = Cursor()

        def cursor(self):
            return self.cursor_

        def close(self):
            self.closed = True

    conn = Conn()
    monkeypatch.setattr(compare.pg, "connect", lambda url: conn)
    ids = compare.database_ids("dsn")
    assert ids["urs_taxid"].to_list() == ["URS0000000001_9606"]
    assert "is_active AND has_coordinates" in conn.cursor_.query
    assert conn.closed
