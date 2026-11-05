# -*- coding: utf-8 -*-

import io
from contextlib import contextmanager
from datetime import date

from rnacentral_pipeline.rnacentral.ftp_export import release_note


class FakeCursor:
    def __init__(self, counts):
        self.counts = counts
        self.sql = None

    def execute(self, sql):
        self.sql = sql

    def fetchone(self):
        return (self.counts[self.sql],)


def patch_cursor(monkeypatch, counts):
    cur = FakeCursor(counts)

    @contextmanager
    def fake(_db_url):
        yield cur

    monkeypatch.setattr(release_note, "cursor", fake)


def test_counts_are_formatted_with_thousands_separators(monkeypatch):
    patch_cursor(
        monkeypatch,
        {
            "select count(*) from rna": 1234567,
            "select count(*) from xref where deleted = 'N'": 89,
            "select count(*) from rnc_database": 1000,
        },
    )
    assert release_note.unique_sequences("db") == "1,234,567"
    assert release_note.active_xrefs("db") == "89"
    assert release_note.active_databases("db") == "1,000"


def test_current_date_is_day_month_year():
    assert release_note.current_date() == date.today().strftime("%d/%m/%Y")


def test_write_fills_the_template(monkeypatch):
    patch_cursor(
        monkeypatch,
        {
            "select count(*) from rna": 10,
            "select count(*) from xref where deleted = 'N'": 20,
            "select count(*) from rnc_database": 30,
        },
    )
    out = io.StringIO()
    template = io.StringIO(
        "{release} {unique_sequences} {active_xrefs} {databases} {date}"
    )
    release_note.write(template, 28, out, "db")
    assert out.getvalue() == f"28 10 20 30 {release_note.current_date()}"
