# -*- coding: utf-8 -*-

import io
import json
import sqlite3

import pytest

from rnacentral_pipeline.rnacentral.search_export import metadata


def lines(*rows):
    return io.StringIO("".join(json.dumps(r) + "\n" for r in rows))


def test_merge_turns_each_row_into_one_entry_per_metadata_type():
    handle = lines({"rna_id": "URS1_9606", "go": {"a": 1}, "r2dt": {"b": 2}})

    assert list(metadata.merge(handle)) == [
        ("go", "URS1_9606", {"a": 1}),
        ("r2dt", "URS1_9606", {"b": 2}),
    ]


def test_merge_of_nothing_is_nothing():
    assert list(metadata.merge(io.StringIO(""))) == []


def test_load_so_terms_reads_type_and_tree():
    handle = lines({"so_rna_type": "SO:1", "so_term_tree": [["SO:1", "x"]]})

    assert list(metadata.load_so_terms(handle)) == [("SO:1", [["SO:1", "x"]])]


def test_the_cache_returns_what_was_merged(tmp_path):
    db = str(tmp_path / "cache.db")
    rows = lines(
        {"rna_id": "URS1_9606", "go": {"a": 1}},
        {"rna_id": "URS2_9606", "r2dt": {"b": 2}},
    )
    tree = lines({"so_rna_type": "SO:1", "so_term_tree": [["SO:1", "x"]]})

    metadata.write_merge(rows, tree, db)

    with metadata.open(db) as cache:
        assert sorted(cache.known_metadata_types) == ["go", "r2dt"]
        assert cache.lookup("URS1_9606") == {
            "go": {"a": 1},
            "r2dt": {},
            "secondary": {},
        }
        assert cache.lookup("URS9_1", missing=None)["go"] == {}
        assert cache.lookup_so("SO:1") == [["SO:1", "x"]]


def test_the_cache_connection_is_closed_when_the_block_exits(tmp_path):
    with metadata.open(str(tmp_path / "cache.db")) as cache:
        pass
    with pytest.raises(sqlite3.ProgrammingError):
        cache.conn.execute("select 1")


def test_a_missing_so_term_falls_back_to_ncrna_and_an_unknown_one_fails(tmp_path):
    with metadata.open(str(tmp_path / "cache.db")) as cache:
        assert cache.lookup_so(None) == [("SO:0000655", "ncRNA")]
        with pytest.raises(AssertionError):
            cache.lookup_so("SO:404")


def test_the_first_value_wins_for_a_repeated_metadata_entry(tmp_path):
    db = str(tmp_path / "cache.db")
    rows = lines(
        {"rna_id": "URS1_9606", "go": {"a": 1}},
        {"rna_id": "URS1_9606", "go": {"a": 2}},
    )

    metadata.write_merge(rows, io.StringIO(""), db)

    with metadata.open(db) as cache:
        assert cache.lookup("URS1_9606")["go"] == {"a": 1}


def test_so_term_tree_writes_one_json_line_per_type(monkeypatch):
    monkeypatch.setattr(metadata.so, "load_ontology", lambda url: "ontology")
    monkeypatch.setattr(metadata.so, "rna_type_tree", lambda ont, term: [[term, ont]])
    out = io.StringIO()

    metadata.write_so_term_tree(io.StringIO("SO:1\nSO:2\n"), "url", out)

    assert [json.loads(line) for line in out.getvalue().splitlines()] == [
        {"so_rna_type": "SO:1", "so_term_tree": [["SO:1", "ontology"]]},
        {"so_rna_type": "SO:2", "so_term_tree": [["SO:2", "ontology"]]},
    ]
