# -*- coding: utf-8 -*-

import io

from click.testing import CliRunner

from rnacentral_pipeline.cli import ftp_export
from rnacentral_pipeline.rnacentral.ftp_export import ensembl


def test_dense_regions_are_cut_into_smaller_ranges(monkeypatch):
    # Equal id windows of four would give 1..4 a thousand xrefs and 5..8 almost none.
    counts = [(1, 400), (2, 400), (3, 100), (4, 100), (5, 1), (6, 1), (7, 1), (8, 1)]
    monkeypatch.setattr(ensembl, "xref_counts", lambda db_url: iter(counts))
    out = io.StringIO()

    ensembl.write_ranges(out, 4, 500)

    assert out.getvalue().splitlines() == ["rna,1,1", "rna,2,3", "rna,4,7", "rna,8,8"]


def test_the_command_passes_its_limits_through(monkeypatch):
    seen = {}

    def fake(output, max_rna, max_xrefs, db_url=None):
        seen.update(max_rna=max_rna, max_xrefs=max_xrefs, db_url=db_url)
        output.write("rna,1,2\n")

    monkeypatch.setattr(ensembl, "write_ranges", fake)

    result = CliRunner().invoke(
        ftp_export.cli, ["ensembl-ranges", "--db-url", "x", "100000", "400000", "-"]
    )

    assert result.exit_code == 0, result.output
    assert seen == {"max_rna": 100000, "max_xrefs": 400000, "db_url": "x"}
    assert result.output == "rna,1,2\n"


class FakeCursor:
    def __init__(self, rows):
        self.rows = rows
        self.executed = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def execute(self, sql):
        self.executed = sql

    def __iter__(self):
        return iter(self.rows)


class FakeConnection:
    def __init__(self, cursor):
        self.fake = cursor
        self.cursor_name = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def cursor(self, name=None):
        self.cursor_name = name
        return self.fake


def test_xref_counts_streams_the_grouped_counts_through_a_named_cursor(monkeypatch):
    cursor = FakeCursor([(1, 3), (2, 5)])
    conn = FakeConnection(cursor)
    monkeypatch.setattr(ensembl.psycopg2, "connect", lambda url: conn)

    assert list(ensembl.xref_counts("db")) == [(1, 3), (2, 5)]
    assert conn.cursor_name == "ensembl_xref_counts"
    assert "xref.deleted = 'N'" in cursor.executed
    assert "GROUP BY rna.id ORDER BY rna.id" in cursor.executed
