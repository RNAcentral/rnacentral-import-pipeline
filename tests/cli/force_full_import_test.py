# -*- coding: utf-8 -*-

"""
--force-full is the escape hatch for the delta import: it must parse everything
without consulting the stored manifest, so a drifted manifest cannot hold back a
reimport. See docs/incremental-parsing.md.
"""

import json
from pathlib import Path

from click.testing import CliRunner

from rnacentral_pipeline.cli import ena as ena_cli
from rnacentral_pipeline.cli import hgnc as hgnc_cli
from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.ena import delta


def _boom(*args, **kwargs):
    raise AssertionError("the stored manifest must not be read under --force-full")


def test_hgnc_force_full_ignores_the_stored_manifest(tmp_path, monkeypatch):
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({"response": {"docs": []}}))

    monkeypatch.setattr(manifest, "load_signatures_for", _boom)
    seen = {}

    def fake_parse(path, db_url, previous):
        seen["previous"] = previous
        return hgnc_cli.parser.ParseResult(
            entries=iter([]), signatures={}, deletions=[]
        )

    monkeypatch.setattr(hgnc_cli.parser, "parse", fake_parse)

    result = CliRunner().invoke(
        hgnc_cli.cli, ["map", "--force-full", str(raw), str(tmp_path)]
    )

    assert result.exit_code == 0, result.output
    assert seen["previous"] == {}


def test_ena_force_full_keeps_everything_without_a_diff(tmp_path, monkeypatch):
    label = delta.source_label("/ena/wgs/aaa")
    signatures = tmp_path / "signatures.csv"
    signatures.write_text(f"{label},AB111111.1,sig-one\n{label},AB222222.3,sig-two\n")

    scanned = tmp_path / "scanned.txt"
    scanned.write_text("/ena/wgs/aaa\n")

    monkeypatch.setattr(ena_cli.manifest, "dump_signatures", _boom)
    monkeypatch.setattr(ena_cli.manifest, "diff_via_polars", _boom)

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
    assert to_parse.read_text().strip() == delta.KEEP_ALL
    assert deletions.read_text() == ""
    assert manifest_csv.read_text().splitlines() == [
        "ENA,AB111111.1,sig-one,/ena/wgs/aaa",
        "ENA,AB222222.3,sig-two,/ena/wgs/aaa",
    ]
