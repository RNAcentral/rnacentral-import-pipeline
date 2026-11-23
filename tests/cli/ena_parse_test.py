# -*- coding: utf-8 -*-

"""
`rnac ena parse`, run the way workflows/databases/ena.nf runs it:

    rnac ena parse --counts <raw>-counts.txt filtered.ncr <tpa> ribotyper-results <model_lengths> .

on the ribotyper cases in data/ena/to-exclude.
"""

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import ena as ena_cli
from rnacentral_pipeline.databases.ena import parser

CASES = Path("data/ena/to-exclude").absolute()
TPA = Path("data/ena/tpa/combined/mapping.tsv").absolute()


def run(case, tmp_path, ribotyper=None):
    ribotyper = ribotyper or CASES / case / "ribotyper-results"
    out = tmp_path / "out"
    counts = tmp_path / "counts.txt"
    result = CliRunner().invoke(
        ena_cli.cli,
        [
            "parse",
            "--counts",
            str(counts),
            str(CASES / case / "sequences.dat"),
            str(TPA),
            str(ribotyper),
            str(CASES / "model-lengths.csv"),
            str(out),
        ],
    )
    return result, out, counts


@pytest.fixture(autouse=True)
def csv_output(monkeypatch):
    monkeypatch.setenv("RNAC_OUTPUT_FORMAT", "csv")


@pytest.fixture
def notifications(monkeypatch):
    sent = []
    monkeypatch.setattr(ena_cli, "send_notification", lambda *a: sent.append(a))
    return sent


def test_an_included_rrna_is_written_with_counts(tmp_path, notifications):
    result, out, counts = run("lazr01007255", tmp_path)
    assert result.exit_code == 0, result.output
    accessions = (out / "accessions.csv").read_text()
    assert "LAZR01007255.1:2..3402:rRNA" in accessions
    assert counts.read_text()
    assert notifications == []


def test_a_file_with_only_excluded_entries_is_not_an_error(tmp_path, notifications):
    # ribotyper rejects the only rRNA in kbtv01008406
    result, out, counts = run("kbtv01008406", tmp_path)
    assert result.exit_code == 0, result.output
    assert "KBTV01008406" not in "".join(p.read_text() for p in out.glob("*.csv"))
    assert counts.read_text()
    assert notifications == []


@pytest.mark.xfail(
    reason="a ValueError while parsing is caught as 'no entries could be written': "
    "it is posted to Slack and the command exits 0 with partial output"
)
def test_a_parse_error_fails_the_command(tmp_path, monkeypatch, notifications):
    real = parser.parse_with_context

    def broken(ctx, path):
        yield from real(ctx, path)
        raise ValueError("bad record")

    monkeypatch.setattr(ena_cli.parser, "parse_with_context", broken)
    # A real ribotyper run leaves a .log, which the error branch attaches
    ribotyper = tmp_path / "ribotyper-results"
    shutil.copytree(CASES / "lazr01007255" / "ribotyper-results", ribotyper)
    (ribotyper / "ribotyper-results.ribotyper.log").write_text("log\n")
    result, _, _ = run("lazr01007255", tmp_path, ribotyper)
    assert notifications == []
    assert result.exit_code != 0
