# -*- coding: utf-8 -*-

import io
import json

from rnacentral_pipeline.rnacentral.ftp_export.sequences_json import filter_dump

# json_build_object output spaces the colon, and COPY doubles backslashes.
DUMP = (
    '{"id" : "URS1", "description" : "x", "sequence" : "AC"}\n'
    '{"id" : "URS2", "description" : "a\\\\nb", "sequence" : "GG"}\n'
    "\n"
    '{"id":"URS3","sequence":"UU"}\n'
)


def run(ids):
    out = io.StringIO()
    filter_dump(io.StringIO(DUMP), ids, out)
    return out.getvalue()


def test_keeps_matching_lines_unchanged():
    assert run({"URS2"}) == DUMP.splitlines(keepends=True)[1]


def test_matches_compact_json_too():
    assert run({"URS3"}) == '{"id":"URS3","sequence":"UU"}\n'


def test_preserves_dump_order_and_skips_blank_lines():
    kept = run({"URS1", "URS3"}).splitlines()
    assert [json.loads(line)["id"] for line in kept] == ["URS1", "URS3"]


def test_no_ids_keeps_nothing():
    assert run(set()) == ""
