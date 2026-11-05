# -*- coding: utf-8 -*-

import pytest

pl = pytest.importorskip("polars")

from rnacentral_pipeline.rnacentral.ftp_export import parquet  # noqa: E402


def dump(n):
    return "".join(f'{{"id": "u{i}", "sequence": "ACGT"}}\n' for i in range(n))


def test_ndjson_is_written_across_several_batches(tmp_path):
    src = tmp_path / "dump.json"
    out = tmp_path / "out.parquet"
    src.write_text(dump(5))
    parquet.ndjson_2_parquet(src, out, batch_size=2)
    assert pl.read_parquet(out)["id"].to_list() == [f"u{i}" for i in range(5)]


def test_copy_sink_reassembles_lines_split_across_chunks(tmp_path):
    out = tmp_path / "sink.parquet"
    sink = parquet._ParquetLineSink(out, batch_size=2)
    blob = dump(5)
    for i in range(0, len(blob), 7):
        sink.write(blob[i : i + 7].encode())
    sink.close()
    assert pl.read_parquet(out).height == 5


def test_copy_sink_keeps_a_final_line_with_no_newline(tmp_path):
    out = tmp_path / "sink.parquet"
    sink = parquet._ParquetLineSink(out, batch_size=10)
    sink.write(dump(2).rstrip("\n"))
    sink.close()
    assert pl.read_parquet(out).height == 2


def test_copy_sink_writes_no_file_when_nothing_arrives(tmp_path):
    out = tmp_path / "sink.parquet"
    parquet._ParquetLineSink(out, batch_size=2).close()
    assert not out.exists()


def test_copy_query_streams_the_query_through_the_sink(tmp_path):
    class Cursor:
        def __enter__(self):
            return self

        def __exit__(self, *_):
            return False

        def copy_expert(self, sql, sink):
            self.sql = sql
            sink.write(dump(3))

    class Conn:
        def cursor(self):
            return Cursor()

    query = tmp_path / "q.sql"
    query.write_text("COPY (select 1) TO STDOUT")
    out = tmp_path / "out.parquet"
    parquet.copy_query_2_parquet(query, Conn(), out, batch_size=2)
    assert pl.read_parquet(out).height == 3
