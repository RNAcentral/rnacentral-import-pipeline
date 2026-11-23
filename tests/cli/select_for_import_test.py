# -*- coding: utf-8 -*-

"""
`rnac scan-imports select-for-import`, as workflows/databases/select.nf runs it:
latest_md5s.csv holds "<db_name>,<md5>" lines, compared against
rnc_import_tracker. TRACKER is production's table: it has one row.
"""

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import scan_imports

TRACKER = [
    {
        "id": 1,
        "db_name": "intact",
        "db_id": 42,
        "last_import_date": "2026-01-20 13:58:04.07019",
        "file_md5": "44cc7e39dadc3101343b20efed8535ef",
    }
]


class Cursor:
    def __init__(self, rows):
        self.rows = rows

    def execute(self, sql):
        assert sql.strip().rstrip(";") == "SELECT * FROM rnc_import_tracker"

    def fetchall(self):
        return self.rows

    def close(self):
        pass


class Connection:
    def __init__(self, rows):
        self.rows = rows

    def cursor(self, cursor_factory=None):
        return Cursor(self.rows)

    def close(self):
        pass


def select(tmp_path, monkeypatch, md5s, tracker=TRACKER):
    monkeypatch.setattr(
        scan_imports.psycopg2, "connect", lambda url: Connection(tracker)
    )
    latest = tmp_path / "latest_md5s.csv"
    latest.write_text("".join(f"{name},{md5}\n" for name, md5 in md5s))
    output = tmp_path / "db_selection.config"
    result = CliRunner().invoke(
        scan_imports.cli,
        [
            "select-for-import",
            str(latest),
            str(output),
            "--db-url",
            "postgres://unused",
        ],
    )
    assert result.exit_code == 0, result.output
    return {
        line.strip() for line in output.read_text().splitlines() if ".run = " in line
    }


def test_an_unchanged_file_is_switched_off(tmp_path, monkeypatch):
    lines = select(tmp_path, monkeypatch, [("intact", TRACKER[0]["file_md5"])])
    assert lines == {"intact.run = false"}


def test_a_changed_file_is_switched_on(tmp_path, monkeypatch):
    lines = select(tmp_path, monkeypatch, [("intact", "0" * 32)])
    assert lines == {"intact.run = true"}


def test_a_database_never_imported_is_switched_on(tmp_path, monkeypatch):
    lines = select(
        tmp_path,
        monkeypatch,
        [("intact", TRACKER[0]["file_md5"]), ("mirbase", "f" * 32)],
    )
    assert lines == {"intact.run = false", "mirbase.run = true"}


def test_the_output_is_a_nextflow_params_block(tmp_path, monkeypatch):
    select(tmp_path, monkeypatch, [("intact", "0" * 32)])
    text = (tmp_path / "db_selection.config").read_text()
    assert text.split()[:3] == ["params", "{", "databases"]


@pytest.mark.xfail(
    reason="an empty rnc_import_tracker gives a DataFrame with no db_name column, "
    "so set_index('db_name') raises KeyError instead of selecting everything"
)
def test_an_empty_tracker_selects_everything(tmp_path, monkeypatch):
    lines = select(tmp_path, monkeypatch, [("intact", "0" * 32)], tracker=[])
    assert lines == {"intact.run = true"}
