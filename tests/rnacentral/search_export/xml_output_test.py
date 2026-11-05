# -*- coding: utf-8 -*-

import io
import json
from datetime import date

import pytest
from lxml import etree

from rnacentral_pipeline.rnacentral.search_export import exporter


def entries(n):
    return [etree.Element("entry", {"id": f"URS{i}_9606"}) for i in range(n)]


def test_write_wraps_the_entries_in_a_database_document():
    out, count = io.StringIO(), io.StringIO()

    exporter.write(iter(entries(2)), out, count)

    root = etree.fromstring(out.getvalue())
    assert root.tag == "database"
    assert root.findtext("name") == "RNAcentral"
    assert root.findtext("release_date") == date.today().strftime("%d/%m/%Y")
    assert [e.get("id") for e in root.find("entries")] == ["URS0_9606", "URS1_9606"]
    assert root.findtext("entry_count") == "2"
    assert count.getvalue() == "2"


def test_write_refuses_to_produce_an_empty_export():
    with pytest.raises(ValueError):
        exporter.write(iter([]), io.StringIO(), io.StringIO())


def test_as_xml_builds_one_entry_per_json_line(monkeypatch):
    monkeypatch.setattr(
        exporter, "raw_builder", lambda e: etree.Element("entry", {"id": e["id"]})
    )
    source = io.StringIO("".join(json.dumps({"id": f"U{i}"}) + "\n" for i in range(3)))
    out, count = io.StringIO(), io.StringIO()

    exporter.as_xml(source, out, count)

    assert count.getvalue() == "3"
    assert etree.fromstring(out.getvalue()).find("entries").__len__() == 3


def test_release_note_sums_the_counts_from_every_chunk():
    out = io.StringIO()

    exporter.release_note(out, "28", [io.StringIO("10"), io.StringIO("5")])

    lines = out.getvalue().splitlines()
    assert lines[0] == "release=28"
    assert lines[1] == "release_date=" + date.today().strftime("%d-%B-%Y")
    assert lines[2] == "entries=15"
