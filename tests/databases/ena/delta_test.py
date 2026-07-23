"""
Unit tests for ENA delta parsing (docs/incremental-parsing-ena.md).

These cover the two pure, database-free passes -- signatures and filter -- and the
signature/diff behaviour. The database-side diff (manifest.diff_via_db) needs a live
Postgres and is exercised by tests/xref-incremental-parity, not here.
"""

from pathlib import Path

from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.ena import delta


def _record(accession, version, seq="acgtacgtac", description="test"):
    return (
        f"ID   {accession}; SV {version}; linear; RNA; STD; PRO; {len(seq)} BP.\n"
        "XX\n"
        f"AC   {accession};\n"
        "XX\n"
        f"DE   {description}\n"
        "XX\n"
        f"FT   source          1..{len(seq)}\n"
        f"FT   ncRNA           1..{len(seq)}\n"
        "XX\n"
        f"SQ   Sequence {len(seq)} BP;\n"
        f"     {seq}        {len(seq)}\n"
        "//\n"
    )


def _chunk(tmp_path, *records):
    path = tmp_path / "chunk.ncr"
    path.write_text("".join(records))
    return path


def test_signature_key_is_accession_dot_version(tmp_path):
    path = _chunk(tmp_path, _record("AB111111", 1), _record("AB222222", 3))
    sigs = dict(delta.iter_signatures(path))
    assert set(sigs) == {"AB111111.1", "AB222222.3"}


def test_signatures_are_stable_and_content_sensitive(tmp_path):
    path_a = _chunk(tmp_path, _record("AB111111", 1, description="one"))
    first = dict(delta.iter_signatures(path_a))
    second = dict(delta.iter_signatures(path_a))
    assert first == second

    path_b = _chunk(tmp_path, _record("AB111111", 1, description="changed"))
    changed = dict(delta.iter_signatures(path_b))
    assert changed["AB111111.1"] != first["AB111111.1"]


def test_filter_keeps_only_listed_accessions(tmp_path):
    path = _chunk(
        tmp_path,
        _record("AB111111", 1),
        _record("AB222222", 3),
        _record("AB333333", 2),
    )
    only = tmp_path / "to_parse.txt"
    only.write_text("AB222222.3\nAB333333.2\n")
    out = tmp_path / "filtered.ncr"

    kept = delta.filter_records(path, only, out)
    assert kept == 2
    assert set(dict(delta.iter_signatures(out))) == {"AB222222.3", "AB333333.2"}


def test_filter_keep_all_sentinel_copies_everything(tmp_path):
    path = _chunk(tmp_path, _record("AB111111", 1), _record("AB222222", 3))
    only = tmp_path / "to_parse.txt"
    only.write_text(delta.KEEP_ALL + "\n")
    out = tmp_path / "filtered.ncr"

    assert delta.filter_records(path, only, out) == 2


def test_filter_empty_to_parse_keeps_nothing(tmp_path):
    path = _chunk(tmp_path, _record("AB111111", 1))
    only = tmp_path / "to_parse.txt"
    only.write_text("")
    out = tmp_path / "filtered.ncr"

    assert delta.filter_records(path, only, out) == 0
    assert out.read_text() == ""


def test_version_bump_is_new_plus_dropped(tmp_path):
    """A sequence-version bump changes the accession key: .N is new, .N-1 dropped."""
    old = _chunk(tmp_path, _record("AB111111", 1))
    old_sigs = dict(delta.iter_signatures(old))

    new = _chunk(tmp_path, _record("AB111111", 2))
    new_sigs = dict(delta.iter_signatures(new))

    diff = manifest.compute_diff(new_sigs, old_sigs)
    assert diff.new == frozenset({"AB111111.2"})
    assert diff.dropped == frozenset({"AB111111.1"})
