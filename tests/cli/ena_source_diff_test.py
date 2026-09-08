# -*- coding: utf-8 -*-

"""
Source-level pruning for ENA: a source whose archives have not changed is never
fetched, and — the part that can lose data if it is wrong — its records must not
then look dropped. See docs/incremental-parsing-ena.md.
"""

import csv

from click.testing import CliRunner

from rnacentral_pipeline.cli import ena as ena_cli
from rnacentral_pipeline.databases.ena import delta


def _boom(*args, **kwargs):
    raise AssertionError("the stored signatures must not be read under --force-full")


STORED = {
    "/ena/wgs/aaa": "sig-aaa",
    "/ena/wgs/aab": "sig-aab",
    "/ena/wgs/gone": "sig-gone",
}


def _source_diff(tmp_path, current, monkeypatch):
    listing = tmp_path / "sources.tsv"
    listing.write_text("".join(f"{p}\t{s}\n" for p, s in current.items()))

    monkeypatch.setattr(
        ena_cli.manifest, "load_file_signatures_for", lambda url, db: dict(STORED)
    )

    out = {name: tmp_path / name for name in ("to_fetch.txt", "scanned.txt", "s.csv")}
    result = CliRunner().invoke(
        ena_cli.cli,
        ["source-diff", str(listing)] + [str(out[n]) for n in out],
    )
    assert result.exit_code == 0, result.output
    return out


def test_only_changed_and_new_sources_are_fetched(tmp_path, monkeypatch):
    out = _source_diff(
        tmp_path,
        {
            "/ena/wgs/aaa": "sig-aaa",  # unchanged
            "/ena/wgs/aab": "sig-moved",  # changed
            "/ena/wgs/aac": "sig-aac",  # new
        },
        monkeypatch,
    )

    assert out["to_fetch.txt"].read_text().split() == [
        "/ena/wgs/aab",
        "/ena/wgs/aac",
    ]


def test_scanned_covers_changed_and_vanished_but_never_unchanged(tmp_path, monkeypatch):
    """
    scanned is what the record diff may retire. An unchanged source must stay out of
    it, or skipping the fetch would retire every record it holds.
    """
    out = _source_diff(
        tmp_path,
        {"/ena/wgs/aaa": "sig-aaa", "/ena/wgs/aab": "sig-moved"},
        monkeypatch,
    )

    scanned = out["scanned.txt"].read_text().split()
    assert "/ena/wgs/aaa" not in scanned
    assert sorted(scanned) == ["/ena/wgs/aab", "/ena/wgs/gone"]


def test_sources_artifact_lists_every_current_source(tmp_path, monkeypatch):
    """What gets promoted after the load: the current set, so gone sources go too."""
    out = _source_diff(
        tmp_path, {"/ena/wgs/aaa": "sig-aaa", "/ena/wgs/aab": "sig-moved"}, monkeypatch
    )

    with out["s.csv"].open() as handle:
        assert list(csv.reader(handle)) == [
            ["ENA", "/ena/wgs/aaa", "sig-aaa"],
            ["ENA", "/ena/wgs/aab", "sig-moved"],
        ]


def test_delta_diff_resolves_each_label_back_to_its_source(tmp_path, monkeypatch):
    """manifest.csv has to carry the real path, not the chunk-name label."""
    source = "/ena/wgs/aab"
    signatures = tmp_path / "signatures.csv"
    signatures.write_text(f"{delta.source_label(source)},AB111111.1,sig-one\n")

    scanned = tmp_path / "scanned.txt"
    scanned.write_text(source + "\n")

    to_parse = tmp_path / "to_parse.txt"
    deletions = tmp_path / "deletions.csv"
    manifest_csv = tmp_path / "manifest.csv"

    result = CliRunner().invoke(
        ena_cli.cli,
        [
            "delta-diff",
            "--force-full",
            str(signatures),
            str(scanned),
            str(to_parse),
            str(deletions),
            str(manifest_csv),
        ],
    )

    assert result.exit_code == 0, result.output
    with manifest_csv.open() as handle:
        assert list(csv.reader(handle)) == [["ENA", "AB111111.1", "sig-one", source]]


def test_force_full_fetches_every_source(tmp_path, monkeypatch):
    """
    A forced full run releases with FULL, which retires whatever the load does not
    carry. Skipping an unchanged source there would retire every record in it.
    """
    listing = tmp_path / "sources.tsv"
    listing.write_text("/ena/wgs/aaa\tsig-aaa\n/ena/wgs/aab\tsig-aab\n")

    monkeypatch.setattr(ena_cli.manifest, "load_file_signatures_for", _boom)

    to_fetch = tmp_path / "to_fetch.txt"
    scanned = tmp_path / "scanned.txt"
    sources = tmp_path / "sources.csv"

    result = CliRunner().invoke(
        ena_cli.cli,
        [
            "source-diff",
            "--force-full",
            str(listing),
            str(to_fetch),
            str(scanned),
            str(sources),
        ],
    )

    assert result.exit_code == 0, result.output
    assert to_fetch.read_text().split() == ["/ena/wgs/aaa", "/ena/wgs/aab"]
    assert scanned.read_text().split() == ["/ena/wgs/aaa", "/ena/wgs/aab"]
