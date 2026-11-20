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
