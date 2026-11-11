# -*- coding: utf-8 -*-

import io
import json

import pytest

from rnacentral_pipeline.rnacentral.ftp_export import sequences_json as sj


def test_parse_reads_one_entry_per_line():
    lines = "".join('{"id": "u%d"}\n' % i for i in range(3))
    assert [e["id"] for e in sj.parse(io.StringIO(lines))] == ["u0", "u1", "u2"]


def test_parse_skips_blank_lines():
    assert list(sj.parse(io.StringIO('\n{"id": "a"}\n  \n'))) == [{"id": "a"}]


def test_decode_entry_undoes_the_doubled_backslashes_from_copy():
    # COPY text format doubles backslashes, which makes `\\n` an invalid JSON escape pair here.
    assert sj.decode_entry('{"d": "x\\\\ny"}') == {"d": "x\\ny"}
    assert sj.decode_entry('{"d": "a\\\\"b"}'.replace('\\\\"', '\\\\\\"')) == {
        "d": 'a\\"b'
    }


def test_decode_entry_reports_the_original_error_when_both_attempts_fail():
    with pytest.raises(json.JSONDecodeError):
        sj.decode_entry("{bad")


def test_parse_fails_on_a_bad_line_without_reading_on(caplog):
    handle = io.StringIO('{"id": "u0"}\n{bad\n{"id": "u2"}\n')
    entries = sj.parse(handle)
    assert next(entries) == {"id": "u0"}
    with pytest.raises(json.JSONDecodeError):
        next(entries)
    assert "line 2" in caplog.text
