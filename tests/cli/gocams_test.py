import csv

import pyarrow.parquet as pq
from click.testing import CliRunner

from rnacentral_pipeline.cli import gocams as gocams_cli


def _fake_find(urs_taxids):
    # Consumes the input lazily, like the real generator, so a file closed too
    # early would surface here.
    for urs_taxid in urs_taxids:
        if urs_taxid == "URS_HIT_9606":
            yield urs_taxid, "abc123"
            yield urs_taxid, "def456"  # one URS in two CAMs gives two rows


def _run(tmp_path, monkeypatch, fmt, lines="URS_HIT_9606\nURS_MISS_9606\n\n"):
    # format_option writes RNAC_OUTPUT_FORMAT into os.environ; setenv first so
    # monkeypatch restores it afterwards.
    monkeypatch.setenv("RNAC_OUTPUT_FORMAT", "csv")
    monkeypatch.setattr(gocams_cli, "find_gocam", _fake_find)
    infile = tmp_path / "urs.txt"
    infile.write_text(lines)
    out = tmp_path / "out"
    result = CliRunner().invoke(
        gocams_cli.cli, ["search", "--format", fmt, str(infile), str(out)]
    )
    return result, out


def test_search_writes_csv_rows_for_hits_only(tmp_path, monkeypatch):
    result, out = _run(tmp_path, monkeypatch, "csv")

    assert result.exit_code == 0, result.output
    with (out / "go_cam_memberships.csv").open() as handle:
        rows = list(csv.reader(handle))
    assert rows == [["URS_HIT_9606", "abc123"], ["URS_HIT_9606", "def456"]]


def test_search_writes_parquet(tmp_path, monkeypatch):
    result, out = _run(tmp_path, monkeypatch, "parquet")

    assert result.exit_code == 0, result.output
    table = pq.read_table(out / "go_cam_memberships.parquet")
    assert table.column_names == ["urs_taxid", "go_cam_id"]
    assert table.to_pylist() == [
        {"urs_taxid": "URS_HIT_9606", "go_cam_id": "abc123"},
        {"urs_taxid": "URS_HIT_9606", "go_cam_id": "def456"},
    ]


def test_search_empty_input_writes_empty_csv(tmp_path, monkeypatch):
    result, out = _run(tmp_path, monkeypatch, "csv", lines="\n")

    assert result.exit_code == 0, result.output
    assert (out / "go_cam_memberships.csv").read_text() == ""


def test_search_api_error_fails(tmp_path, monkeypatch):
    def boom(urs_taxids):
        raise RuntimeError("API down")
        yield  # pragma: no cover - makes this a generator

    monkeypatch.setenv("RNAC_OUTPUT_FORMAT", "csv")
    monkeypatch.setattr(gocams_cli, "find_gocam", boom)
    infile = tmp_path / "urs.txt"
    infile.write_text("URS_X_1\n")

    result = CliRunner().invoke(
        gocams_cli.cli,
        ["search", "--format", "csv", str(infile), str(tmp_path / "out")],
    )

    assert result.exit_code != 0
