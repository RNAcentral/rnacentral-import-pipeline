import json

import pytest

from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.hgnc import parser


def _write(tmp_path, docs):
    path = tmp_path / "raw.json"
    path.write_text(json.dumps({"response": {"docs": docs}}))
    return path


DOCS = [
    {"hgnc_id": "HGNC:1", "symbol": "AAA", "name": "a", "locus_type": "RNA, micro"},
    {"hgnc_id": "HGNC:2", "symbol": "BBB", "name": "b", "locus_type": "RNA, micro"},
    {"hgnc_id": "HGNC:3", "symbol": "CCC", "name": "c", "locus_type": "RNA, micro"},
]


def _mapped(monkeypatch):
    """Record which accessions actually reach the (expensive) mapping step."""
    seen = []

    def fake_map(ctx, entry):
        seen.append(entry.hgnc_id)
        return "URS_" + entry.hgnc_id

    monkeypatch.setattr(parser, "rnacentral_id", fake_map)
    monkeypatch.setattr(parser, "as_entry", lambda ctx, entry, urs: (entry.hgnc_id, urs))
    # Context.build must never be called when there is nothing to parse. It takes
    # the entries to look up in bulk, so the stub has to accept them too.
    monkeypatch.setattr(parser.Context, "build", lambda db_url, entries: object())
    return seen


def test_only_changed_and_new_are_parsed(tmp_path, monkeypatch):
    seen = _mapped(monkeypatch)
    path = _write(tmp_path, DOCS)

    sigs = {d["hgnc_id"]: manifest.record_signature(d) for d in DOCS}
    # HGNC:2 changed, HGNC:3 gone, HGNC:4 will be new below.
    previous = dict(sigs)
    previous["HGNC:2"] = "stale"
    previous["HGNC:9"] = "was-here"          # dropped
    del previous["HGNC:1"]                    # so HGNC:1 looks new too

    result = parser.parse(path, "postgres://unused", previous)
    entries = list(result.entries)

    assert sorted(seen) == ["HGNC:1", "HGNC:2"]          # new + changed only
    assert sorted(e[0] for e in entries) == ["HGNC:1", "HGNC:2"]
    assert result.deletions == ["HGNC:9"]                # dropped accession
    assert set(result.signatures) == {"HGNC:1", "HGNC:2", "HGNC:3"}


def test_no_change_parses_nothing_and_never_connects(tmp_path, monkeypatch):
    seen = _mapped(monkeypatch)

    def explode(db_url, entries):
        raise AssertionError("Context.build must not run when nothing changed")

    monkeypatch.setattr(parser.Context, "build", explode)
    path = _write(tmp_path, DOCS)
    previous = {d["hgnc_id"]: manifest.record_signature(d) for d in DOCS}

    result = parser.parse(path, "postgres://unused", previous)
    assert list(result.entries) == []
    assert seen == []
    assert result.deletions == []


def test_bootstrap_parses_everything(tmp_path, monkeypatch):
    seen = _mapped(monkeypatch)
    path = _write(tmp_path, DOCS)

    result = parser.parse(path, "postgres://unused", {})
    list(result.entries)
    assert sorted(seen) == ["HGNC:1", "HGNC:2", "HGNC:3"]
    assert result.deletions == []
